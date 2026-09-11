using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// All other real-genome results in this document use ZCurveWindowedProfile's default 500 bp
    /// window, inherited unchanged from the unrelated k-mer-based WindowedProfile method - never
    /// swept against real data on its own merits. This reruns RealGenomeDistanceMatrixBenchmark's
    /// exact relationship-ordering check at several window sizes on the same 30-genome corpus, to see
    /// whether the ordering (and the default) holds up.
    ///
    /// Each genome's sequence is read from disk once and reused across every window size (profiles at
    /// different window sizes are cheap, linear re-scans of the same in-memory string - re-reading and
    /// re-decompressing 30 gzipped FASTA files per window size would be pure waste). Window sizes are
    /// chosen to bound total runtime: BestAlignment's O(m*(n-m+1)) cost (§6, §8.7) scales with window
    /// *count*, so halving window size roughly doubles both m and n for every mismatched-length pair
    /// in this corpus (genomes range 580 kbp-8.67 Mbp) - a 100 bp window was measured to cost an
    /// order of magnitude more than the 250 bp point below and was dropped for exactly that reason.
    /// </summary>
    public class WindowSizeSensitivityBenchmark
    {
        private static readonly GenomeMeta[] genomes = GenomeTaxonomy.Genomes;
        private static readonly int[] windowSizes = { 250, 500, 1000, 2000 };

        private static double Mean(List<double> v) => v.Average();

        private static double StdDev(List<double> v)
        {
            if (v.Count < 2) return 0;
            var mean = v.Average();
            return Math.Sqrt(v.Sum(x => (x - mean) * (x - mean)) / (v.Count - 1));
        }

        [Test]
        public void RelationshipOrderingHoldsAcrossWindowSizes()
        {
            var sequences = genomes.Select(g =>
                (Meta: g, Sequence: FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence))
                .ToList();

            var relationshipLevelsInOrder = new[]
            {
                RelationshipLevel.SameSpecies, RelationshipLevel.SameGenus,
                RelationshipLevel.SamePhylum, RelationshipLevel.DifferentPhylum,
            };

            using var report = File.CreateText("WindowSizeSensitivity.csv");
            report.WriteLine("WindowSize, RelationshipLevel, PairCount, Mean, StdDev");

            foreach (var windowSize in windowSizes)
            {
                var profiles = sequences.Select(s => (s.Meta, Profile: new ZCurveWindowedProfile(s.Sequence, windowSize))).ToList();

                var rows = new List<(RelationshipLevel Level, double Score)>();
                for (int a = 0; a < profiles.Count; a++)
                    for (int b = a + 1; b < profiles.Count; b++)
                    {
                        var score = ZCurveWindowedProfile.Compare(profiles[a].Profile, profiles[b].Profile);
                        rows.Add((GenomeTaxonomy.Classify(profiles[a].Meta, profiles[b].Meta), score));
                    }

                var byLevel = rows.GroupBy(r => r.Level).ToDictionary(g => g.Key, g => g.Select(r => r.Score).ToList());

                foreach (var level in relationshipLevelsInOrder)
                {
                    if (!byLevel.TryGetValue(level, out var scores)) continue;
                    report.WriteLine($"{windowSize}, {level}, {scores.Count}, {Mean(scores)}, {StdDev(scores)}");
                }

                var meanSameSpecies = Mean(byLevel[RelationshipLevel.SameSpecies]);
                var meanSameGenus = Mean(byLevel[RelationshipLevel.SameGenus]);
                var meanSamePhylum = Mean(byLevel[RelationshipLevel.SamePhylum]);
                var meanDifferentPhylum = Mean(byLevel[RelationshipLevel.DifferentPhylum]);

                TestContext.Out.WriteLine(
                    $"w={windowSize,5}  SameSpecies={meanSameSpecies:F4}  SameGenus={meanSameGenus:F4}  " +
                    $"SamePhylum={meanSamePhylum:F4}  DifferentPhylum={meanDifferentPhylum:F4}");

                // The central claim §8.6 makes at w=500 - mean similarity strictly decreases with
                // taxonomic distance - should hold at every window size tested, not just the default.
                Assert.That(meanSameSpecies, Is.GreaterThan(meanSameGenus),
                    $"ordering should hold at window size {windowSize}");
                Assert.That(meanSameGenus, Is.GreaterThan(meanSamePhylum),
                    $"ordering should hold at window size {windowSize}");
                Assert.That(meanSamePhylum, Is.GreaterThan(meanDifferentPhylum),
                    $"ordering should hold at window size {windowSize}");
            }
        }
    }
}
