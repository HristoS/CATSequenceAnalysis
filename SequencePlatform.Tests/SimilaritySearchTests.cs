using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class SimilaritySearchTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "search";

        // Stores the ladder levels and the sequence_sha256 a real ingest would, under this class's method version.
        private static async Task<long> StoreAsync(string provider, string accession, string dna, string kind = "bacterium")
        {
            var sequenceId = await ThrowawayDatabase.InsertRecordAsync(provider, accession, dna.Length, kind);
            var query = QueryProfile.FromSequence(dna);

            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var update = new Npgsql.NpgsqlCommand("update sequence_record set sequence_sha256 = $1 where id = $2", connection, transaction))
            {
                update.Parameters.AddWithValue(query.Sha256!);
                update.Parameters.AddWithValue(sequenceId);
                await update.ExecuteNonQueryAsync();
            }
            foreach (var factor in query.Forward.Factors)
                if (query.Forward.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, query.Forward.Level(factor));
            await transaction.CommitAsync();
            return sequenceId;
        }

        private static SimilaritySearch NewSearch(bool useCache = true) => new(ThrowawayDatabase.DataSource, useCache, method: method, methodVersion: methodVersion);

        [Test]
        public async Task RanksOnlyComparableLengthsAtTheBandLevelOnTheBetterStrand()
        {
            const string provider = "search-band";
            var lengths = new[] { 900_000, 600_000, 300_000, 200_000, 40_000 };
            var stored = new Dictionary<long, (int Length, ZCurveProfileLadder Ladder)>();
            for (int i = 0; i < lengths.Length; i++)
            {
                var dna = Generator.GetRandomDNA(lengths[i], seed: 100 + i);
                stored[await StoreAsync(provider, $"r{i}", dna)] = (lengths[i], RecordStore.Ladder(new ZCurveWindowedProfile(dna)));
            }
            var query = QueryProfile.FromSequence(Generator.GetRandomDNA(500_000, seed: 1));

            var cached = await NewSearch().SearchAsync(query, new CatalogueFilter(Provider: provider));
            var streamed = await NewSearch(useCache: false).SearchAsync(query, new CatalogueFilter(Provider: provider));

            Assert.That(cached.QueryWindowSize, Is.EqualTo(8_000), "15 windows at 32 kbp leave a half-length candidate only 7");
            Assert.That(cached.FewWindows, Is.False);
            Assert.That(cached.Lengths, Is.EqualTo(new LengthRange(250_000, 1_000_001)));
            Assert.That(cached.Ranked.Select(h => stored[h.SequenceId].Length), Is.EquivalentTo(new[] { 900_000, 600_000, 300_000 }));
            foreach (var hit in cached.Ranked)
            {
                var candidate = stored[hit.SequenceId].Ladder.Level(16);
                var forward = ZCurveWindowedProfile.Compare(query.Forward.Level(16), candidate);
                var reverse = ZCurveWindowedProfile.Compare(query.Reverse!.Level(16), candidate);
                Assert.That(hit.Score, Is.EqualTo(Math.Max(forward, reverse)));
                Assert.That(hit.Strand, Is.EqualTo(reverse > forward ? Strand.Reverse : Strand.Forward));
            }
            Assert.That(cached.Ranked.Select(h => h.Score), Is.Ordered.Descending);
            Assert.That((streamed.QueryWindowSize, streamed.Lengths, streamed.RankedTotal), Is.EqualTo((cached.QueryWindowSize, cached.Lengths, cached.RankedTotal)));
            Assert.That(streamed.Ranked, Is.EqualTo(cached.Ranked), "the cache and the database give the same result");
            Assert.That(streamed.NotComparable, Is.EqualTo(cached.NotComparable).And.Empty);
        }

        [Test]
        public async Task RecordStoredOnTheOtherStrandIsFoundOnTheReverseStrandAndAsIdentical()
        {
            const string provider = "search-strand";
            var dna = Generator.GetRandomDNA(400_123, seed: 7);   // not a multiple of 500: a profile flip wouldn't be exact
            var id = await StoreAsync(provider, "plus", dna);
            await StoreAsync(provider, "other", Generator.GetRandomDNA(400_000, seed: 8));

            var reversed = await NewSearch().SearchAsync(QueryProfile.FromSequence(QueryProfile.ReverseComplement(dna)), new CatalogueFilter(Provider: provider));
            var same = await NewSearch().SearchAsync(QueryProfile.FromSequence(dna.ToLowerInvariant()), new CatalogueFilter(Provider: provider));

            Assert.That(reversed.Ranked[0], Is.EqualTo(new SimilarityHit(id, 1, Strand.Reverse)));
            Assert.That(reversed.Identical, Is.EqualTo(new[] { (id, Strand.Reverse) }));
            Assert.That(same.Ranked[0], Is.EqualTo(new SimilarityHit(id, 1, Strand.Forward)));
            Assert.That(same.Identical, Is.EqualTo(new[] { (id, Strand.Forward) }), "letter case doesn't change the hash");
        }

        [Test]
        public async Task FilterAndRatioNarrowTheCandidates()
        {
            const string provider = "search-filter";
            var virus = await StoreAsync(provider, "virus", Generator.GetRandomDNA(30_000, seed: 1), kind: "virus");
            var longerVirus = await StoreAsync(provider, "longer", Generator.GetRandomDNA(50_000, seed: 2), kind: "virus");
            await StoreAsync(provider, "bacterium", Generator.GetRandomDNA(31_000, seed: 3));
            var query = QueryProfile.FromSequence(Generator.GetRandomDNA(30_000, seed: 4));

            var half = await NewSearch().SearchAsync(query, new CatalogueFilter(Kind: "virus", Provider: provider));
            var close = await NewSearch().SearchAsync(query, new CatalogueFilter(Kind: "virus", Provider: provider), lengthRatio: 0.75);

            Assert.That(half.Ranked.Select(h => h.SequenceId), Is.EquivalentTo(new[] { virus, longerVirus }));
            Assert.That(close.Ranked.Select(h => h.SequenceId), Is.EqualTo(new[] { virus }), "50 kbp is outside [22.5, 40] kbp");
            Assert.That(() => NewSearch().SearchAsync(query, lengthRatio: 0.2), Throws.InstanceOf<ArgumentOutOfRangeException>());
        }

        [Test]
        public async Task QueryShorterThanOneWindowComparesWithNothing()
        {
            var result = await NewSearch().SearchAsync(QueryProfile.FromSequence(Generator.GetRandomDNA(499, seed: 1)),
                new CatalogueFilter(Provider: "search-nothing"));

            Assert.That(result.QueryWindowSize, Is.Null);
            Assert.That(result.Ranked, Is.Empty);
            Assert.That(result.FewWindows, Is.True);
        }

        [TestCase(7_906, 0.5, true)]     // HPV16's length: 15 windows, a half-length candidate keeps 7
        [TestCase(8_000, 0.5, false)]
        [TestCase(8_000, 0.25, true)]
        public async Task ShortQueriesAreFlaggedAsMostlyComposition(int length, double ratio, bool expected)
        {
            var result = await NewSearch().SearchAsync(QueryProfile.FromSequence(Generator.GetRandomDNA(length, seed: 1)),
                new CatalogueFilter(Provider: "search-nothing"), lengthRatio: ratio);

            Assert.That(result.FewWindows, Is.EqualTo(expected));
        }

        // The cache holds each level only for the record lengths CachedLengths gives; a candidate outside it
        // would be silently missing. Every query length and allowed ratio must keep its whole band inside.
        [Test]
        public void EveryBandFitsTheCachedLengthsOfItsLevel()
        {
            var search = NewSearch();
            var factors = ZCurveProfileLadder.DefaultFactors.ToList();
            foreach (var ratio in new[] { 0.25, 0.3, 0.5, 0.75, 1.0 })
                for (long length = 500; length < 20_000_000; length = length * 21 / 20 + 1)
                {
                    int windows = (int)(length / 500);
                    var levelWindows = factors.Select(f => windows / f).ToList();
                    int level = Enumerable.Range(0, factors.Count).LastOrDefault(i => ratio * levelWindows[i] >= ZCurveProfileLadder.DefaultMinWindows);
                    var range = search.CachedLengths(level);
                    long low = (long)Math.Ceiling(ratio * length), high = (long)Math.Floor(length / ratio);
                    Assert.That(range.Contains(low) && range.Contains(high), Is.True, $"length {length}, ratio {ratio}, level {level}: [{low}, {high}] vs {range}");
                }
        }

        [Test]
        public void ReverseComplementKeepsOnlyLettersAndComplementsIupacCodes()
        {
            Assert.That(QueryProfile.ReverseComplement("acgT\nRYKM-BVDH 1SWN"), Is.EqualTo("NWSDHBVKMRYACGT"));
        }
    }
}
