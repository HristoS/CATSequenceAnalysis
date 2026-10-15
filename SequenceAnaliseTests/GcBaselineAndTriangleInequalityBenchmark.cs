using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Two cheap checks on the 30-genome corpus behind §8.6/§8.8, run before the platform's
    /// organism-organizing work (P5) builds on <see cref="ZCurveWindowedProfile.Compare"/>:
    ///
    ///  - C6, GC baseline: how much of the taxonomic signal is just GC content? The z axis of the
    ///    Z-curve is strong/weak H-bond, i.e. local GC, so a GC-only score is the obvious null. It is
    ///    scored as -|GC_a - GC_b| (higher = more similar) with the same AUC helper and the same two
    ///    thresholds as <see cref="RetrievalClassificationBenchmark"/>. If Z-curve's AUC is within
    ///    0.03 of GC-only, P5 has to show what the x and y axes add.
    ///  - C7, triangle inequality: is d = 1 - Compare usable as a tree distance? Compare is a
    ///    best-offset containment score, so nothing guarantees it is a metric. Every triple is checked
    ///    on all three inequalities. Mash distance is checked the same way as a reference point.
    /// </summary>
    public class GcBaselineAndTriangleInequalityBenchmark
    {
        private const double violationTolerance = 1e-12;

        private static readonly GenomeMeta[] genomes = GenomeTaxonomy.Genomes;

        [Test]
        public void GcBaselineAndTriangleInequality()
        {
            var primaries = genomes.Select(g => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence).ToList();
            var gc = primaries.Select(GcFraction).ToArray();
            var profiles = primaries.Select(s => new ZCurveWindowedProfile(s)).ToArray();
            var mash = ExternalBaselineData.ParseMashDistances(ExternalBaselineData.RawResultsPath("mash_dist_raw.tsv"));

            int n = genomes.Length;
            var zDistance = new double[n, n];
            var mashDistance = new double[n, n];
            double maxSelfDistance = 0;
            for (int a = 0; a < n; a++)
            {
                maxSelfDistance = Math.Max(maxSelfDistance, 1 - ZCurveWindowedProfile.Compare(profiles[a], profiles[a]));
                for (int b = a + 1; b < n; b++)
                {
                    zDistance[a, b] = zDistance[b, a] = 1 - ZCurveWindowedProfile.Compare(profiles[a], profiles[b]);
                    mashDistance[a, b] = mashDistance[b, a] = mash[ExternalBaselineData.PairKey(genomes[a].FileName, genomes[b].FileName)];
                }
            }

            WriteGcBaseline(gc, zDistance);
            WriteTriangleInequality(zDistance, mashDistance, maxSelfDistance);
        }

        // Over A, C, G, T only, so IUPAC ambiguity codes don't count, the same as both profile types.
        private static double GcFraction(string sequence)
        {
            long strong = 0, total = 0;
            foreach (var c in sequence)
                switch (char.ToUpperInvariant(c))
                {
                    case 'G': case 'C': strong++; total++; break;
                    case 'A': case 'T': total++; break;
                }
            return strong / (double)total;
        }

        private static void WriteGcBaseline(double[] gc, double[,] zDistance)
        {
            var rows = new List<(RelationshipLevel Level, double ZCurve, double Gc)>();
            for (int a = 0; a < genomes.Length; a++)
                for (int b = a + 1; b < genomes.Length; b++)
                    rows.Add((GenomeTaxonomy.Classify(genomes[a], genomes[b]), 1 - zDistance[a, b], -Math.Abs(gc[a] - gc[b])));

            var thresholds = new (string Name, Func<RelationshipLevel, bool> IsPositive)[]
            {
                ("SameGenusOrCloser_vs_Rest", l => l is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus),
                ("SamePhylumOrCloser_vs_DifferentPhylum", l => l is not RelationshipLevel.DifferentPhylum),
            };

            using var report = File.CreateText("GcBaseline_Retrieval.csv");
            report.WriteLine("Threshold, Method, PositiveCount, NegativeCount, AUC, AucMinusZCurve");
            foreach (var (name, isPositive) in thresholds)
            {
                var positives = rows.Where(r => isPositive(r.Level)).ToList();
                var negatives = rows.Where(r => !isPositive(r.Level)).ToList();
                var zAuc = RetrievalClassificationBenchmark.ComputeAuc(
                    positives.Select(r => r.ZCurve).ToList(), negatives.Select(r => r.ZCurve).ToList());
                var gcAuc = RetrievalClassificationBenchmark.ComputeAuc(
                    positives.Select(r => r.Gc).ToList(), negatives.Select(r => r.Gc).ToList());

                report.WriteLine($"{name}, ZCurveWindowedProfile, {positives.Count}, {negatives.Count}, {zAuc}, 0");
                report.WriteLine($"{name}, AbsGcDifference, {positives.Count}, {negatives.Count}, {gcAuc}, {gcAuc - zAuc}");
                TestContext.Out.WriteLine($"C6 {name} (n={positives.Count}/{negatives.Count}): AUC ZCurve={zAuc:F4} GC-only={gcAuc:F4} diff={gcAuc - zAuc:+0.0000;-0.0000}");
            }
        }

        private static void WriteTriangleInequality(double[,] zDistance, double[,] mashDistance, double maxSelfDistance)
        {
            using var report = File.CreateText("DistanceTriangleInequality.csv");
            report.WriteLine("Distance, Genomes, Triples, TriplesViolating, InequalitiesChecked, InequalitiesViolated, " +
                "ViolatedShare, MaxViolation, MeanViolation, MaxSelfDistance");

            foreach (var (name, distance, selfDistance) in new[]
            {
                ("OneMinusZCurveCompare", zDistance, maxSelfDistance),
                ("MashDistance", mashDistance, double.NaN),
            })
            {
                var t = CountViolations(distance);
                double share = t.InequalitiesViolated / (double)t.InequalitiesChecked;
                report.WriteLine($"{name}, {genomes.Length}, {t.Triples}, {t.TriplesViolating}, {t.InequalitiesChecked}, " +
                    $"{t.InequalitiesViolated}, {share}, {t.MaxViolation}, {t.MeanViolation}, {selfDistance}");
                TestContext.Out.WriteLine($"C7 {name}: {t.TriplesViolating}/{t.Triples} triples violate, " +
                    $"{t.InequalitiesViolated}/{t.InequalitiesChecked} inequalities ({share:P2}), max {t.MaxViolation:F4}, mean {t.MeanViolation:F4}, self {selfDistance:F6}");
            }
        }

        // A violation of d(a,b) <= d(a,c) + d(c,b) is measured as d(a,b) - d(a,c) - d(c,b), in distance units.
        private static (int Triples, int TriplesViolating, int InequalitiesChecked, int InequalitiesViolated, double MaxViolation, double MeanViolation)
            CountViolations(double[,] d)
        {
            int n = d.GetLength(0);
            int triples = 0, triplesViolating = 0, violated = 0;
            double maxViolation = 0, violationSum = 0;
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    for (int c = b + 1; c < n; c++)
                    {
                        triples++;
                        bool any = false;
                        foreach (var excess in new[] { d[a, b] - d[a, c] - d[c, b], d[a, c] - d[a, b] - d[b, c], d[b, c] - d[b, a] - d[a, c] })
                        {
                            if (excess <= violationTolerance) continue;
                            any = true;
                            violated++;
                            violationSum += excess;
                            maxViolation = Math.Max(maxViolation, excess);
                        }
                        if (any) triplesViolating++;
                    }
            return (triples, triplesViolating, triples * 3, violated, maxViolation, violated == 0 ? 0 : violationSum / violated);
        }
    }
}
