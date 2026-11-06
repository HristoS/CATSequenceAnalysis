using System.Diagnostics;
using System.Net.Http.Json;
using System.Text;
using System.Text.Json;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// The search worker in the API process (S4.3). Each test starts its own host, so only one worker runs at a time; jobs other
    /// fixtures left queued may run too, and the assertions look only at this fixture's own.
    /// </summary>
    public class WorkerTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "worker";

        private static StringContent Fasta(string dna) => new($">q\n{dna}\n", Encoding.ASCII, "text/x-fasta");

        private static async Task<long> StoreAsync(string provider, string accession, string dna, long? ownerId = null)
        {
            var id = ownerId is long owner
                ? await ThrowawayDatabase.InsertUploadAsync(owner, accession, dna.Length, "bacterium")
                : await ThrowawayDatabase.InsertRecordAsync(provider, accession, dna.Length, "bacterium");
            var ladder = RecordStore.Ladder(new ZCurveWindowedProfile(dna));
            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var factor in ladder.Factors)
                if (ladder.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, id, method, methodVersion, factor, ladder.Level(factor));
            await transaction.CommitAsync();
            return id;
        }

        private static async Task<SearchRun> WaitForAsync(Guid id, params SearchStatus[] statuses)
        {
            var watch = Stopwatch.StartNew();
            while (watch.Elapsed < TimeSpan.FromSeconds(60))
            {
                var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!;
                if (statuses.Contains(run.Status)) return run;
                await Task.Delay(50);
            }
            Assert.Fail($"Search {id} didn't reach {string.Join(" or ", statuses)} within 60 s.");
            return null!;
        }

        private static async Task<Guid> PostAsync(HttpClient client, string route, string dna)
        {
            var response = await client.PostAsync(route, Fasta(dna));
            Assert.That((int)response.StatusCode, Is.EqualTo(202), await response.Content.ReadAsStringAsync());
            return (await response.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("id").GetGuid();
        }

        private static async Task ExecuteAsync(string sql, Guid id)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(sql);
            command.Parameters.AddWithValue(id);
            await command.ExecuteNonQueryAsync();
        }

        [Test]
        public async Task QueuedSearchesRunOnTheirOwnAndSeeRecordsStoredSinceTheCacheLoaded()
        {
            const string provider = "worker-run";
            var genomes = Enumerable.Range(0, 3).Select(i => Generator.GetRandomDNA(120_000, seed: 90 + i)).ToArray();
            var ids = new List<long>();
            for (int i = 0; i < genomes.Length; i++) ids.Add(await StoreAsync(provider, $"w{i}", genomes[i]));
            using var factory = ApiHost.Create(methodVersion);
            var (_, client) = await factory.NewClientAsync("alice");

            var watch = Stopwatch.StartNew();
            var similar = await PostAsync(client, $"/api/searches/similar?provider={provider}", genomes[0]);
            var first = await WaitForAsync(similar, SearchStatus.Done, SearchStatus.Failed);
            var similarMs = watch.ElapsedMilliseconds;
            watch.Restart();
            var localize = await PostAsync(client, $"/api/searches/localize?provider={provider}", genomes[2].Substring(40_000, 12_000));
            await WaitForAsync(localize, SearchStatus.Done, SearchStatus.Failed);
            var localizeMs = watch.ElapsedMilliseconds;
            TestContext.Out.WriteLine($"POST to Done: find-similar {similarMs} ms (first, cache cold), localize {localizeMs} ms (streamed)");

            Assert.That((first.Status, first.RankedTotal), Is.EqualTo((SearchStatus.Done, 3)));
            var ranked = await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, similar, HitList.Ranked);
            Assert.That((ranked[0].SequenceId, ranked[0].Score), Is.EqualTo((ids[0], (double?)1.0)));
            var placed = await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, localize, HitList.Ranked);
            Assert.That(placed[0].SequenceId, Is.EqualTo(ids[2]));
            Assert.That(placed[0].OffsetBp, Is.EqualTo(40_000).Within(500));

            // Stored after the worker's cache loaded this band: ranked with its profile, not "not comparable".
            var late = await StoreAsync(provider, "late", genomes[1]);
            var again = await PostAsync(client, $"/api/searches/similar?provider={provider}", genomes[1]);
            Assert.That((await WaitForAsync(again, SearchStatus.Done, SearchStatus.Failed)).Status, Is.EqualTo(SearchStatus.Done));
            Assert.That((await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, again, HitList.Ranked)).Select(h => h.SequenceId), Does.Contain(late));
            Assert.That(await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, again, HitList.NotComparable), Is.Empty);
        }

        [Test]
        public async Task AFailingSearchKeepsItsReason()
        {
            var (alice, _) = await Users.CreateAsync(ThrowawayDatabase.DataSource, $"alice-{Guid.NewGuid():N}");
            var levelless = await ThrowawayDatabase.InsertUploadAsync(alice.Id, "no-levels", 50_000);
            var id = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource,
                new SearchRequest(SearchKind.Localize, new CatalogueFilter(ViewerId: alice.Id), 20, WithinSequenceId: levelless),
                QueryProfile.FromSequence(Generator.GetRandomDNA(3_000, seed: 95)));

            using var factory = ApiHost.Create(methodVersion);
            _ = factory.Server;

            var run = await WaitForAsync(id, SearchStatus.Done, SearchStatus.Failed);
            Assert.That((run.Status, run.Error), Is.EqualTo((SearchStatus.Failed, $"Record {levelless} has no stored 500 bp level.")));
        }

        [Test]
        public async Task OnStartTheWorkerRequeuesInterruptedSearchesFailsRepeatOffendersAndExpiresOldResults()
        {
            var query = QueryProfile.FromSequence(Generator.GetRandomDNA(20_000, seed: 96));
            var request = new SearchRequest(SearchKind.Similar, new CatalogueFilter(Provider: "worker-nobody"), 10, LengthRatio: 0.5);
            var interrupted = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, request, query);
            var offender = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, request, query);
            var old = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, request, query);
            await ExecuteAsync("update search_run set status = 'running', started_at = now() - interval '1 minute', attempts = 1 where public_id = $1", interrupted);
            await ExecuteAsync("update search_run set status = 'running', started_at = now() - interval '1 minute', attempts = 3 where public_id = $1", offender);
            await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, old, new SimilarityResult(500, new LengthRange(10_000, 40_001), [], 0, [], [], FewWindows: false));
            await ExecuteAsync("update search_run set expires_at = now() - interval '1 second' where public_id = $1", old);

            using var factory = ApiHost.Create(methodVersion);
            _ = factory.Server;

            var rerun = await WaitForAsync(interrupted, SearchStatus.Done, SearchStatus.Failed);
            var gaveUp = await WaitForAsync(offender, SearchStatus.Done, SearchStatus.Failed);
            var expired = await WaitForAsync(old, SearchStatus.Expired);
            Assert.That(rerun.Status, Is.EqualTo(SearchStatus.Done));
            Assert.That((gaveUp.Status, gaveUp.Error), Is.EqualTo((SearchStatus.Failed, "The search stopped its worker 3 times; it was not retried.")));
            Assert.That(expired.RankedTotal, Is.Zero, "the summary outlives the results");
        }
    }
}
