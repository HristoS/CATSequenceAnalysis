using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class SearchLogTests
    {
        private static readonly SearchRequest similar = new(SearchKind.Similar, new CatalogueFilter(Kind: "virus", Provider: "log"), 20, LengthRatio: 0.5);

        // Random DNA with a replication-like GC skew and its origin at the middle, so the query is rotated (OriginRotationTests).
        private static string SkewedFromMiddle(int length, int seed)
        {
            var chars = Generator.GetRandomDNA(length, seed).ToCharArray();
            var rng = new Random(seed + 1);
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < 0.1)
                    chars[i] = i < length / 2 ? (chars[i] == 'C' ? 'G' : chars[i]) : (chars[i] == 'G' ? 'C' : chars[i]);
            var dna = new string(chars);
            return dna[(length / 2)..] + dna[..(length / 2)];
        }

        private static void AssertSameLadder(ZCurveProfileLadder? actual, ZCurveProfileLadder? expected)
        {
            Assert.That(actual is null, Is.EqualTo(expected is null));
            if (expected is null) return;
            foreach (var factor in expected.Factors)
                Assert.That(actual!.Level(factor).Windows, Is.EqualTo(expected.Level(factor).Windows), $"level x{factor}");
        }

        [Test]
        public async Task JobRebuildsTheQueryExactlyWithoutTheDna()
        {
            var query = QueryProfile.FromSequence(SkewedFromMiddle(600_123, seed: 1));
            Assert.That(query.Rotated, Is.Not.Null, "the test needs a rotated query");

            var id = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, similar, query);
            var restored = (await SearchLog.QueryAsync(ThrowawayDatabase.DataSource, id))!;
            var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!;

            Assert.That((restored.LengthBp, restored.Sha256, restored.ReverseSha256), Is.EqualTo((600_123L, query.Sha256, query.ReverseSha256)));
            AssertSameLadder(restored.Forward, query.Forward);
            AssertSameLadder(restored.Reverse, query.Reverse);
            AssertSameLadder(restored.Rotated, query.Rotated);
            AssertSameLadder(restored.RotatedReverse, query.RotatedReverse);
            Assert.That((run.Status, run.Request, run.QueryLengthBp), Is.EqualTo((SearchStatus.Queued, similar, 600_123L)));
            Assert.That((run.ExpiresAt - run.CreatedAt).TotalDays, Is.EqualTo(30).Within(1 / 24.0), "D7; a DST change shifts a calendar day by an hour");
            Assert.That((run.StartedAt, run.CompletedAt, run.WindowSize), Is.EqualTo(((DateTimeOffset?)null, (DateTimeOffset?)null, (int?)null)));
        }

        [Test]
        public async Task ProfileUploadJobHasNoHashAndOneStrand()
        {
            var forward = new ZCurveWindowedProfile(Generator.GetRandomDNA(20_000, seed: 2));
            var request = new SearchRequest(SearchKind.Localize, new CatalogueFilter(), 5, WithinSequenceId: 42);

            var id = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, request, QueryProfile.FromProfiles(forward, null));
            var restored = (await SearchLog.QueryAsync(ThrowawayDatabase.DataSource, id))!;

            Assert.That((restored.LengthBp, restored.Sha256, restored.ReverseSha256, restored.Reverse), Is.EqualTo((20_000L, (byte[]?)null, (byte[]?)null, (ZCurveProfileLadder?)null)));
            Assert.That(restored.Forward.Fine.Windows, Is.EqualTo(forward.ToHalfPrecision().Windows), "rounded as stored levels are");
            Assert.That((await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!.Request, Is.EqualTo(request));
        }

        [Test]
        public async Task SimilarityResultReadsBackListByListInPages()
        {
            var ids = new List<long>();
            for (int i = 0; i < 5; i++) ids.Add(await ThrowawayDatabase.InsertRecordAsync("log-similar", $"s{i}", 30_000));
            var ranked = new[]
            {
                new SimilarityHit(ids[0], 1, Strand.Reverse, new PValue(-4.5, PValueKind.Bound), -0.2, SharedNeighbours: 17, Rotated: true),
                new SimilarityHit(ids[1], 0.91, Strand.Forward, new PValue(-2.25, PValueKind.Empirical), 1.5, SharedNeighbours: 3),
                new SimilarityHit(ids[2], 0.88, Strand.Forward),
            };
            var result = new SimilarityResult(2_000, new LengthRange(15_000, 60_001), ranked, 214,
                new[] { (ids[0], Strand.Reverse) }, new[] { ids[3], ids[4] }, FewWindows: true, NeighbourRanked: true);
            var id = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, similar, QueryProfile.FromSequence(Generator.GetRandomDNA(30_000, seed: 3)));

            Assert.That(await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, id, result), Is.True);

            var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!;
            Assert.That((run.Status, run.WindowSize, run.Lengths, run.RankedTotal, run.Candidates, run.WeakSignal, run.NeighbourRanked),
                Is.EqualTo((SearchStatus.Done, (int?)2_000, (LengthRange?)new LengthRange(15_000, 60_001), (int?)214, (int?)null, (bool?)true, (bool?)true)));
            Assert.That(run.CompletedAt, Is.Not.Null);

            var all = await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, id, HitList.Ranked);
            Assert.That(all.Select(h => new SimilarityHit(h.SequenceId, h.Score!.Value, h.Strand!.Value, h.P, h.Log10E, h.SharedNeighbours, h.Rotated)),
                Is.EqualTo(ranked));
            Assert.That(all.Select(h => h.Position), Is.EqualTo(new[] { 1, 2, 3 }));
            var page = await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, id, HitList.Ranked, afterPosition: 1, limit: 1);
            Assert.That(page.Select(h => h.SequenceId), Is.EqualTo(new[] { ids[1] }));
            Assert.That(await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, id, HitList.Identical),
                Is.EqualTo(new[] { new StoredHit(HitList.Identical, 1, ids[0], Strand.Reverse) }));
            Assert.That(await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, id, HitList.NotComparable),
                Is.EqualTo(new[] { new StoredHit(HitList.NotComparable, 1, ids[3]), new StoredHit(HitList.NotComparable, 2, ids[4]) }));
            Assert.That(await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, id, result), Is.False, "a done run isn't written twice");
        }

        [Test]
        public async Task LocalizationResultKeepsOffsetsAndCoarseScores()
        {
            var source = await ThrowawayDatabase.InsertRecordAsync("log-localize", "l0", 400_000, "bacterium");
            var other = await ThrowawayDatabase.InsertRecordAsync("log-localize", "l1", 300_000, "bacterium");
            var hits = new[]
            {
                new LocalizationHit(source, 123_500, Strand.Reverse, 0.97, 0.93, new PValue(-3, PValueKind.Empirical), -1),
                new LocalizationHit(other, 0, Strand.Forward, 0.71, 0.75),
            };
            var request = new SearchRequest(SearchKind.Localize, new CatalogueFilter(Kind: "bacterium"), 20);
            var id = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, request, QueryProfile.FromSequence(Generator.GetRandomDNA(10_000, seed: 4)));

            await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, id, new LocalizationResult(2_000, 812, hits, Array.Empty<long>(), WeakSignal: false));

            var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!;
            Assert.That((run.WindowSize, run.Candidates, run.Lengths, run.RankedTotal, run.WeakSignal, run.NeighbourRanked),
                Is.EqualTo(((int?)2_000, (int?)812, (LengthRange?)null, (int?)null, (bool?)false, (bool?)null)));
            var stored = await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, id, HitList.Ranked);
            Assert.That(stored.Select(h => new LocalizationHit(h.SequenceId, h.OffsetBp!.Value, h.Strand!.Value, h.Score!.Value, h.CoarseScore!.Value, h.P, h.Log10E)),
                Is.EqualTo(hits));
            Assert.That(await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, id, HitList.NotComparable), Is.Empty);
        }

        [Test]
        public async Task ExpiryKeepsOnlyHashLengthAndSummary()
        {
            var record = await ThrowawayDatabase.InsertRecordAsync("log-expiry", "e0", 30_000);
            var query = QueryProfile.FromSequence(Generator.GetRandomDNA(30_000, seed: 5));
            var result = new SimilarityResult(500, new LengthRange(15_000, 60_001), new[] { new SimilarityHit(record, 0.9, Strand.Forward) }, 1,
                Array.Empty<(long, Strand)>(), Array.Empty<long>(), FewWindows: false);
            var due = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, similar, query);
            var deleted = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, similar, query);
            var kept = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, similar, query);
            foreach (var id in new[] { due, deleted, kept }) await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, id, result);
            await using (var backdate = ThrowawayDatabase.DataSource.CreateCommand("update search_run set expires_at = now() - interval '1 second' where public_id = $1"))
            {
                backdate.Parameters.AddWithValue(due);
                await backdate.ExecuteNonQueryAsync();
            }

            Assert.That(await SearchLog.ExpireAsync(ThrowawayDatabase.DataSource, deleted), Is.True);
            Assert.That(await SearchLog.ExpireAsync(ThrowawayDatabase.DataSource, deleted), Is.False, "already expired");
            Assert.That(await SearchLog.ExpireDueAsync(ThrowawayDatabase.DataSource), Is.EqualTo(1));

            foreach (var id in new[] { due, deleted })
            {
                var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!;
                Assert.That((run.Status, run.QueryLengthBp, run.WindowSize, run.RankedTotal), Is.EqualTo((SearchStatus.Expired, 30_000L, (int?)500, (int?)1)));
                Assert.That(await SearchLog.QueryAsync(ThrowawayDatabase.DataSource, id), Is.Null);
                Assert.That(await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, id, HitList.Ranked), Is.Empty);
            }
            await using (var hash = ThrowawayDatabase.DataSource.CreateCommand(
                "select query_sha256 = $2 and query_reverse_sha256 = $3 from search_run where public_id = $1"))
            {
                hash.Parameters.AddWithValue(due);
                hash.Parameters.AddWithValue(query.Sha256!);
                hash.Parameters.AddWithValue(query.ReverseSha256!);
                Assert.That(await hash.ExecuteScalarAsync(), Is.True, "D7 keeps the hashes");
            }
            Assert.That((await SearchLog.GetAsync(ThrowawayDatabase.DataSource, kept))!.Status, Is.EqualTo(SearchStatus.Done));
            Assert.That(await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, kept, HitList.Ranked), Has.Count.EqualTo(1));
        }

        [Test]
        public async Task ExpiredOrFailedJobTakesNoResult()
        {
            var query = QueryProfile.FromSequence(Generator.GetRandomDNA(30_000, seed: 6));
            var result = new SimilarityResult(500, new LengthRange(15_000, 60_001), Array.Empty<SimilarityHit>(), 0,
                Array.Empty<(long, Strand)>(), Array.Empty<long>(), FewWindows: false);
            var expired = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, similar, query);
            var failed = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, similar, query);

            await SearchLog.ExpireAsync(ThrowawayDatabase.DataSource, expired);
            Assert.That(await SearchLog.FailAsync(ThrowawayDatabase.DataSource, failed, "no stored level"), Is.True);

            Assert.That(await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, expired, result), Is.False, "a DELETE during the run wins");
            Assert.That(await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, failed, result), Is.False);
            var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, failed))!;
            Assert.That((run.Status, run.Error), Is.EqualTo((SearchStatus.Failed, "no stored level")));
            Assert.That(await SearchLog.GetAsync(ThrowawayDatabase.DataSource, Guid.NewGuid()), Is.Null);
        }
    }
}
