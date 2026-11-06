using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class LocalizationTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "localize";

        // Stores the ladder levels a real ingest would, under this class's method version.
        private static async Task<long> StoreAsync(string provider, string accession, string dna)
        {
            var sequenceId = await ThrowawayDatabase.InsertRecordAsync(provider, accession, dna.Length, "bacterium");
            var ladder = RecordStore.Ladder(new ZCurveWindowedProfile(dna));

            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var factor in ladder.Factors)
                if (ladder.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, ladder.Level(factor));
            await transaction.CommitAsync();
            return sequenceId;
        }

        private static Localization NewLocalization(bool useCache = true) =>
            new(ThrowawayDatabase.DataSource, useCache, method, methodVersion);

        [TestCase(5_000, 500)]
        [TestCase(20_000, 2_000)]
        [TestCase(50_000, 8_000)]
        [TestCase(150_000, 32_000)]
        public async Task FragmentIsPlacedInItsRecordOnEitherStrand(int length, int coarseWindowSize)
        {
            string provider = $"localize-place-{length}";
            var genomes = Enumerable.Range(0, 3).Select(i => Generator.GetRandomDNA(300_000, seed: 10 + i)).ToArray();
            var ids = new List<long>();
            for (int i = 0; i < genomes.Length; i++) ids.Add(await StoreAsync(provider, $"g{i}", genomes[i]));
            const int start = 123_456;   // off the 500 bp grid
            var piece = genomes[1].Substring(start, length);
            var filter = new CatalogueFilter(Provider: provider);

            var forward = await NewLocalization().LocalizeAsync(QueryProfile.FromSequence(piece), filter);
            var reverse = await NewLocalization().LocalizeAsync(QueryProfile.FromSequence(QueryProfile.ReverseComplement(piece)), filter);
            var streamed = await NewLocalization(useCache: false).LocalizeAsync(QueryProfile.FromSequence(piece), filter);

            Assert.That((forward.CoarseWindowSize, forward.Candidates, forward.WeakSignal), Is.EqualTo((coarseWindowSize, 3, length < 8_000)));
            Assert.That(forward.Hits.Select(h => h.SequenceId), Is.EquivalentTo(ids), "each record is placed once");
            Assert.That(forward.Hits[0].SequenceId, Is.EqualTo(ids[1]));
            Assert.That(forward.Hits[0].Strand, Is.EqualTo(Strand.Forward));
            Assert.That(forward.Hits[0].OffsetBp, Is.EqualTo(start).Within(500));
            Assert.That(forward.Hits.Select(h => h.Score), Is.Ordered.Descending);
            Assert.That(reverse.Hits[0].SequenceId, Is.EqualTo(ids[1]));
            Assert.That(reverse.Hits[0].Strand, Is.EqualTo(Strand.Reverse));
            Assert.That(reverse.Hits[0].OffsetBp, Is.EqualTo(start).Within(500));
            Assert.That(streamed.Hits, Is.EqualTo(forward.Hits), "the cache and the database give the same result");
        }

        [Test]
        public async Task RefinedScoreIsTheFineAlignmentNearTheCoarsePlacement()
        {
            const string provider = "localize-refine";
            var dna = Generator.GetRandomDNA(400_000, seed: 3);
            var id = await StoreAsync(provider, "g", dna);
            var piece = dna.Substring(200_321, 40_000);

            var result = await NewLocalization().LocalizeAsync(QueryProfile.FromSequence(piece), new CatalogueFilter(Provider: provider));
            var (score, offset) = ZCurveWindowedProfile.BestAlignment(new ZCurveWindowedProfile(piece).ToHalfPrecision(), new ZCurveWindowedProfile(dna).ToHalfPrecision());

            Assert.That(result.Hits.Single(), Is.EqualTo(result.Hits[0] with { SequenceId = id, OffsetBp = offset * 500L, Score = score }),
                "the refine finds the full 500 bp scan's best");
        }

        [Test]
        public async Task OnlyRecordsAtLeastAsLongAsTheFragmentAreCandidates()
        {
            const string provider = "localize-length";
            await StoreAsync(provider, "short", Generator.GetRandomDNA(30_000, seed: 1));
            var longer = await StoreAsync(provider, "long", Generator.GetRandomDNA(100_000, seed: 2));

            var result = await NewLocalization().LocalizeAsync(QueryProfile.FromSequence(Generator.GetRandomDNA(50_000, seed: 3)),
                new CatalogueFilter(Provider: provider));

            Assert.That(result.Candidates, Is.EqualTo(1));
            Assert.That(result.Hits.Select(h => h.SequenceId), Is.EqualTo(new[] { longer }));
        }

        [Test]
        public async Task ShortFragmentsArePlacedOnlyInsideANamedRecord()
        {
            const string provider = "localize-short";
            var dna = Generator.GetRandomDNA(100_000, seed: 4);
            var id = await StoreAsync(provider, "g", dna);
            // On the 500 bp grid: off it, 4 windows of random DNA are too few to place reliably, which is the weak signal
            // the result reports (here 2 kbp from 61,234 landed at 45,500).
            var piece = QueryProfile.FromSequence(QueryProfile.ReverseComplement(dna.Substring(61_000, 2_000)));

            Assert.That(() => NewLocalization().LocalizeAsync(piece, new CatalogueFilter(Provider: provider)), Throws.ArgumentException);
            var inside = await NewLocalization().LocalizeInRecordAsync(piece, id);
            var tooShort = await NewLocalization().LocalizeInRecordAsync(QueryProfile.FromSequence(dna.Substring(0, 499)), id);

            Assert.That(inside.WeakSignal, Is.True);
            Assert.That(inside.Hits.Single().Strand, Is.EqualTo(Strand.Reverse));
            Assert.That((inside.Hits.Single().OffsetBp, inside.Hits.Single().Score), Is.EqualTo((61_000L, 1.0)));
            Assert.That(tooShort.Hits, Is.Empty);
        }

        // The cache holds each level only for the lengths CachedLengths gives; a fragment scanned at a level below its range
        // would silently miss records. Every fragment length must have its candidates (records at least as long) inside it.
        [Test]
        public void EveryFragmentLevelHoldsItsCandidates()
        {
            var factors = ZCurveProfileLadder.DefaultFactors.ToList();
            for (long length = Localization.MinCatalogueLengthBp; length < 20_000_000; length = length * 21 / 20 + 1)
            {
                long windows = length / 500;
                int level = Enumerable.Range(0, factors.Count).Last(i => windows / factors[i] >= ZCurveProfileLadder.DefaultMinFragmentWindows);
                Assert.That(Localization.CachedLengths(level).Contains(length), Is.True, $"length {length}, level {level}");
            }
        }
    }
}
