using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// The systematic version of RealGenomeProfileTests.cs's six hand-picked pairs: every one of
    /// the C(30,2) = 435 unique pairs among data/genomes/'s 30 real bacterial chromosomes, each
    /// pair labeled by a taxonomically conservative relationship level (see
    /// <see cref="RelationshipLevel"/>) and scored with <see cref="ZCurveWindowedProfile.Compare"/>.
    /// Produces the evidence for the manuscript's real-data section: whether mean similarity is
    /// monotonically ordered by genuine phylogenetic distance across a real, unselected sample,
    /// not just a handful of pairs picked to make the point.
    ///
    /// Genus/phylum assignments (in <see cref="GenomeTaxonomy"/>) are deliberately conservative -
    /// only genus and phylum are asserted (both essentially undisputed for every organism here), not
    /// family or order, since
    /// several family-level placements in this set are unsettled or have shifted recently (e.g.
    /// Yersinia's family, Mycoplasmoides' phylum). Under-claiming relatedness (grouping two
    /// close-but-not-identical-genus organisms as merely "different phylum" would be the wrong
    /// direction of error) is avoided by checking genus first, phylum second - never a
    /// family-level guess.
    /// </summary>
    public class RealGenomeDistanceMatrixBenchmark
    {
        // Taxonomy labels (Genomes, RelationshipLevel, Classify) live in GenomeTaxonomy.cs, shared
        // with ExternalBaselineComparisonBenchmark.cs so both benchmarks compare against the exact
        // same relationship labels rather than two copies that could silently drift apart.
        private static readonly GenomeMeta[] Genomes = GenomeTaxonomy.Genomes;
        private static RelationshipLevel Classify(GenomeMeta a, GenomeMeta b) => GenomeTaxonomy.Classify(a, b);

        private static double Mean(List<double> v) => v.Average();

        private static double StdDev(List<double> v)
        {
            if (v.Count < 2) return 0;
            var mean = v.Average();
            return Math.Sqrt(v.Sum(x => (x - mean) * (x - mean)) / (v.Count - 1));
        }

        [Test]
        public void PairwiseDistanceMatrixAcrossAllRealGenomes()
        {
            var profiles = Genomes.Select(g =>
            {
                var primary = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName));
                return (Meta: g, Profile: new ZCurveWindowedProfile(primary.Sequence));
            }).ToList();

            var rows = new List<(GenomeMeta A, GenomeMeta B, RelationshipLevel Level, double Score)>();
            for (int i = 0; i < profiles.Count; i++)
                for (int j = i + 1; j < profiles.Count; j++)
                {
                    var score = ZCurveWindowedProfile.Compare(profiles[i].Profile, profiles[j].Profile);
                    var level = Classify(profiles[i].Meta, profiles[j].Meta);
                    rows.Add((profiles[i].Meta, profiles[j].Meta, level, score));
                }

            Assert.That(rows.Count, Is.EqualTo(Genomes.Length * (Genomes.Length - 1) / 2));

            using (var raw = File.CreateText("RealGenomeDistanceMatrix_Raw.csv"))
            {
                raw.WriteLine("OrganismA, OrganismB, RelationshipLevel, ZCurveCompare");
                foreach (var row in rows)
                    raw.WriteLine($"{row.A.OrganismName}, {row.B.OrganismName}, {row.Level}, {row.Score}");
            }

            var byLevel = rows.GroupBy(r => r.Level)
                .ToDictionary(g => g.Key, g => g.Select(r => r.Score).ToList());

            using (var summary = File.CreateText("RealGenomeDistanceMatrix_Summary.csv"))
            {
                summary.WriteLine("RelationshipLevel, PairCount, Mean, StdDev, Min, Max");
                foreach (var level in new[] { RelationshipLevel.SameSpecies, RelationshipLevel.SameGenus, RelationshipLevel.SamePhylum, RelationshipLevel.DifferentPhylum })
                {
                    if (!byLevel.TryGetValue(level, out var scores)) continue;
                    summary.WriteLine($"{level}, {scores.Count}, {Mean(scores)}, {StdDev(scores)}, {scores.Min()}, {scores.Max()}");
                    TestContext.Out.WriteLine($"{level,-16} n={scores.Count,4}  mean={Mean(scores):F4}  stdev={StdDev(scores):F4}  min={scores.Min():F4}  max={scores.Max():F4}");
                }
            }

            // The central claim: mean similarity strictly decreases as genuine taxonomic distance
            // increases, across the full unselected set of pairs - not just the hand-picked ones in
            // RealGenomeProfileTests.cs.
            var meanSameSpecies = Mean(byLevel[RelationshipLevel.SameSpecies]);
            var meanSameGenus = Mean(byLevel[RelationshipLevel.SameGenus]);
            var meanSamePhylum = Mean(byLevel[RelationshipLevel.SamePhylum]);
            var meanDifferentPhylum = Mean(byLevel[RelationshipLevel.DifferentPhylum]);

            Assert.That(meanSameSpecies, Is.GreaterThan(meanSameGenus));
            Assert.That(meanSameGenus, Is.GreaterThan(meanSamePhylum));
            Assert.That(meanSamePhylum, Is.GreaterThan(meanDifferentPhylum));
        }
    }
}
