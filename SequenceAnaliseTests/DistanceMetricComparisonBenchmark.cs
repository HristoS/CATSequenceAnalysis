using System.Diagnostics;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Diagnostic, not production: focused on SquaredEuclidean specifically now - the one alternative
    /// (of Euclidean/SquaredEuclidean/Chebyshev/Manhattan) that measured better than the real baseline
    /// on every axis tested: strong correlation (Pearson 0.984), preserved taxonomic ordering, and a
    /// real ~1.3x speedup. Chebyshev and Manhattan were also tested against the same corpus and are
    /// not repeated here - Chebyshev in particular measured *slower* than the sqrt-based baseline
    /// despite having fewer nominal operations (likely branching in Math.Max's IEEE-754 NaN/signed-zero
    /// handling), so neither outperformed SquaredEuclidean on any axis.
    ///
    /// Compares three ways of turning per-window squared distances into a score:
    /// - Euclidean: today's real production formula (sqrt per window, then mean) - the baseline.
    /// - SquaredEuclidean: no sqrt at all (mean of squared distances).
    /// - SquaredEuclideanSqrtOfSum: sqrt applied once to the *summed* squared distance, then divided
    ///   by window count (`1 - sqrt(sum)/m`). NOT the same as RMSD (`sqrt(mean)` = `sqrt(sum/m)`) -
    ///   this has `1/m`, RMSD has `1/sqrt(m)` - so this formula's score carries an extra length-
    ///   dependent factor RMSD wouldn't: two pairs with identical average per-window divergence but
    ///   different query window counts m would NOT get the same score. Real query window counts in
    ///   this corpus range from ~1,160 to ~17,340 (15x), so this is a real, not theoretical, concern -
    ///   flagged before running, then checked empirically rather than assumed either way.
    ///
    /// ZCurveWindowedProfile itself is never modified: the baseline is the actual, unmodified
    /// Compare/BestAlignment, and the alternatives are independently reimplemented here against the
    /// same real 30-genome corpus and taxonomy labelsalready established in
    /// RealGenomeDistanceMatrixBenchmark.cs / ExternalBaselineComparisonBenchmark.cs.
    /// </summary>
    public class DistanceMetricComparisonBenchmark
    {
        // Three independent, formula-inlined reimplementations of BestAlignment's sliding search - no
        // Func<> delegate indirection, so a timing comparison reflects each formula's real cost.
        // CompareEuclidean is arithmetically and order-of-operations identical to
        // ZCurveWindowedProfile.Compare/BestAlignment, so it's a fair, same-shape stand-in for the
        // real baseline in a timing comparison, not just an approximation.

        private static double CompareEuclidean((double X, double Y, double Z)[] x, (double X, double Y, double Z)[] y)
        {
            var (q, r) = x.Length <= y.Length ? (x, y) : (y, x);
            if (q.Length == 0 || r.Length < q.Length) return 0;

            double best = double.NegativeInfinity;
            for (int offset = 0; offset + q.Length <= r.Length; offset++)
            {
                double sum = 0;
                for (int j = 0; j < q.Length; j++)
                {
                    var a = q[j];
                    var b = r[offset + j];
                    var dx = a.X - b.X;
                    var dy = a.Y - b.Y;
                    var dz = a.Z - b.Z;
                    sum += Math.Sqrt(dx * dx + dy * dy + dz * dz);
                }
                double score = 1 - sum / q.Length;
                if (score > best) best = score;
            }
            return best;
        }

        private static double CompareSquaredEuclidean((double X, double Y, double Z)[] x, (double X, double Y, double Z)[] y)
        {
            var (q, r) = x.Length <= y.Length ? (x, y) : (y, x);
            if (q.Length == 0 || r.Length < q.Length) return 0;

            double best = double.NegativeInfinity;
            for (int offset = 0; offset + q.Length <= r.Length; offset++)
            {
                double sum = 0;
                for (int j = 0; j < q.Length; j++)
                {
                    var a = q[j];
                    var b = r[offset + j];
                    var dx = a.X - b.X;
                    var dy = a.Y - b.Y;
                    var dz = a.Z - b.Z;
                    sum += dx * dx + dy * dy + dz * dz;
                }
                double score = 1 - sum / q.Length;
                if (score > best) best = score;
            }
            return best;
        }

        // The project owner's suggested modification: same squared-distance accumulation as
        // CompareSquaredEuclidean above, but sqrt applied once to the summed value rather than never
        // at all. See the class-level remarks for why this is not the same formula as RMSD despite
        // the surface resemblance (1/m here, not RMSD's 1/sqrt(m)).
        private static double CompareSquaredEuclideanSqrtOfSum((double X, double Y, double Z)[] x, (double X, double Y, double Z)[] y)
        {
            var (q, r) = x.Length <= y.Length ? (x, y) : (y, x);
            if (q.Length == 0 || r.Length < q.Length) return 0;

            double best = double.NegativeInfinity;
            for (int offset = 0; offset + q.Length <= r.Length; offset++)
            {
                double sum = 0;
                for (int j = 0; j < q.Length; j++)
                {
                    var a = q[j];
                    var b = r[offset + j];
                    var dx = a.X - b.X;
                    var dy = a.Y - b.Y;
                    var dz = a.Z - b.Z;
                    sum += dx * dx + dy * dy + dz * dz;
                }
                double score = 1 - Math.Sqrt(sum) / q.Length;
                if (score > best) best = score;
            }
            return best;
        }

        private static double Mean(IReadOnlyList<double> v) => v.Average();

        private static double StdDev(IReadOnlyList<double> v)
        {
            if (v.Count < 2) return 0;
            var mean = v.Average();
            return Math.Sqrt(v.Sum(x => (x - mean) * (x - mean)) / (v.Count - 1));
        }

        private static double Pearson(IReadOnlyList<double> x, IReadOnlyList<double> y)
        {
            var mx = x.Average();
            var my = y.Average();
            double sxy = 0, sxx = 0, syy = 0;
            for (int i = 0; i < x.Count; i++)
            {
                var dx = x[i] - mx;
                var dy = y[i] - my;
                sxy += dx * dy;
                sxx += dx * dx;
                syy += dy * dy;
            }
            return sxy / Math.Sqrt(sxx * syy);
        }

        // Ordinary least squares, y ~ slope*x + intercept. x is the cheap/predictor metric, y the
        // real baseline being approximated - so calling as LinearRegression(cheapScores,
        // euclideanScores) gives the direction an actual calibration would use. R-squared for simple
        // (one-predictor) linear regression is exactly Pearson's r squared - no separate derivation.
        private static (double Slope, double Intercept) LinearRegression(IReadOnlyList<double> x, IReadOnlyList<double> y)
        {
            var mx = x.Average();
            var my = y.Average();
            double sxy = 0, sxx = 0;
            for (int i = 0; i < x.Count; i++)
            {
                var dx = x[i] - mx;
                sxy += dx * (y[i] - my);
                sxx += dx * dx;
            }
            double slope = sxy / sxx;
            double intercept = my - slope * mx;
            return (slope, intercept);
        }

        // 1-based ranks, ties broken by averaging - the standard Spearman convention. Ties are
        // vanishingly unlikely for continuous real-genome scores but handled correctly regardless.
        private static double[] Rank(IReadOnlyList<double> values)
        {
            var indexed = values.Select((v, i) => (Value: v, Index: i)).OrderBy(p => p.Value).ToArray();
            var ranks = new double[values.Count];
            int i = 0;
            while (i < indexed.Length)
            {
                int j = i;
                while (j + 1 < indexed.Length && indexed[j + 1].Value == indexed[i].Value) j++;
                double avgRank = (i + 1 + j + 1) / 2.0;
                for (int k = i; k <= j; k++)
                    ranks[indexed[k].Index] = avgRank;
                i = j + 1;
            }
            return ranks;
        }

        private static double Spearman(IReadOnlyList<double> x, IReadOnlyList<double> y) =>
            Pearson(Rank(x), Rank(y));

        [Test]
        public void CompareSquaredEuclideanVariantsAgainstEuclideanBaseline()
        {
            var genomes = GenomeTaxonomy.Genomes;
            int expectedPairs = genomes.Length * (genomes.Length - 1) / 2;

            var profiles = new (GenomeMeta Meta, ZCurveWindowedProfile Profile, (double X, double Y, double Z)[] Windows)[genomes.Length];
            for (int g = 0; g < genomes.Length; g++)
            {
                var primary = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genomes[g].FileName));
                var profile = new ZCurveWindowedProfile(primary.Sequence);
                profiles[g] = (genomes[g], profile, profile.Windows.ToArray());
            }

            var rows = new (GenomeMeta A, GenomeMeta B, RelationshipLevel Level,
                double Euclidean, double SquaredEuclidean, double SquaredEuclideanSqrtOfSum)[expectedPairs];
            var euclideanStopwatch = new Stopwatch();
            var squaredStopwatch = new Stopwatch();
            var squaredSqrtOfSumStopwatch = new Stopwatch();
            int rowIndex = 0;
            for (int i = 0; i < profiles.Length; i++)
                for (int j = i + 1; j < profiles.Length; j++)
                {
                    var a = profiles[i].Meta;
                    var b = profiles[j].Meta;
                    var level = GenomeTaxonomy.Classify(a, b);
                    var x = profiles[i].Windows;
                    var y = profiles[j].Windows;

                    euclideanStopwatch.Start();
                    var euclidean = CompareEuclidean(x, y);
                    euclideanStopwatch.Stop();

                    squaredStopwatch.Start();
                    var squared = CompareSquaredEuclidean(x, y);
                    squaredStopwatch.Stop();

                    squaredSqrtOfSumStopwatch.Start();
                    var squaredSqrtOfSum = CompareSquaredEuclideanSqrtOfSum(x, y);
                    squaredSqrtOfSumStopwatch.Stop();

                    rows[rowIndex++] = (a, b, level, euclidean, squared, squaredSqrtOfSum);
                }

            Assert.That(rowIndex, Is.EqualTo(expectedPairs));
            Assert.That(rows.All(r =>
                    !double.IsNaN(r.Euclidean) && !double.IsInfinity(r.Euclidean) &&
                    !double.IsNaN(r.SquaredEuclidean) && !double.IsInfinity(r.SquaredEuclidean) &&
                    !double.IsNaN(r.SquaredEuclideanSqrtOfSum) && !double.IsInfinity(r.SquaredEuclideanSqrtOfSum)),
                "no metric should produce NaN/Infinity for any real pair");

            using (var raw = File.CreateText("DistanceMetricComparison_Raw.csv"))
            {
                raw.WriteLine("OrganismA, OrganismB, RelationshipLevel, Euclidean, SquaredEuclidean, SquaredEuclideanSqrtOfSum");
                for (int i = 0; i < rows.Length; i++)
                {
                    var row = rows[i];
                    raw.WriteLine($"{row.A.OrganismName}, {row.B.OrganismName}, {row.Level}, {row.Euclidean}, {row.SquaredEuclidean}, {row.SquaredEuclideanSqrtOfSum}");
                }
            }

            var byLevel = rows.GroupBy(r => r.Level).ToDictionary(g => g.Key, g => g.ToList());

            using (var summary = File.CreateText("DistanceMetricComparison_Summary.csv"))
            {
                summary.WriteLine("RelationshipLevel, PairCount, EuclideanMean, EuclideanStdDev, SquaredEuclideanMean, SquaredEuclideanStdDev, SquaredEuclideanSqrtOfSumMean, SquaredEuclideanSqrtOfSumStdDev");
                foreach (var level in new[] { RelationshipLevel.SameSpecies, RelationshipLevel.SameGenus, RelationshipLevel.SamePhylum, RelationshipLevel.DifferentPhylum })
                {
                    if (!byLevel.TryGetValue(level, out var levelRows)) continue;
                    var euclideanScores = levelRows.Select(r => r.Euclidean).ToList();
                    var squaredScores = levelRows.Select(r => r.SquaredEuclidean).ToList();
                    var squaredSqrtOfSumScores = levelRows.Select(r => r.SquaredEuclideanSqrtOfSum).ToList();

                    summary.WriteLine($"{level}, {levelRows.Count}, {Mean(euclideanScores)}, {StdDev(euclideanScores)}, {Mean(squaredScores)}, {StdDev(squaredScores)}, {Mean(squaredSqrtOfSumScores)}, {StdDev(squaredSqrtOfSumScores)}");
                    TestContext.Out.WriteLine($"{level,-16} n={levelRows.Count,4}  Euclidean={Mean(euclideanScores):F4}  SquaredEuclidean={Mean(squaredScores):F4}  SquaredEuclideanSqrtOfSum={Mean(squaredSqrtOfSumScores):F4}");
                }
            }

            var euclideanAll = rows.Select(r => r.Euclidean).ToList();
            var squaredAll = rows.Select(r => r.SquaredEuclidean).ToList();
            var squaredSqrtOfSumAll = rows.Select(r => r.SquaredEuclideanSqrtOfSum).ToList();

            var squaredPearson = Pearson(euclideanAll, squaredAll);
            var squaredSpearman = Spearman(euclideanAll, squaredAll);
            var (squaredSlope, squaredIntercept) = LinearRegression(squaredAll, euclideanAll);
            var squaredRSquared = squaredPearson * squaredPearson;

            var sqrtOfSumPearson = Pearson(euclideanAll, squaredSqrtOfSumAll);
            var sqrtOfSumSpearman = Spearman(euclideanAll, squaredSqrtOfSumAll);
            var (sqrtOfSumSlope, sqrtOfSumIntercept) = LinearRegression(squaredSqrtOfSumAll, euclideanAll);
            var sqrtOfSumRSquared = sqrtOfSumPearson * sqrtOfSumPearson;

            using (var corr = File.CreateText("DistanceMetricComparison_Correlation.csv"))
            {
                corr.WriteLine("Comparison, PairCount, PearsonR, SpearmanR, Slope, Intercept, RSquared");
                corr.WriteLine($"SquaredEuclidean_vs_Euclidean, {rows.Length}, {squaredPearson}, {squaredSpearman}, {squaredSlope}, {squaredIntercept}, {squaredRSquared}");
                corr.WriteLine($"SquaredEuclideanSqrtOfSum_vs_Euclidean, {rows.Length}, {sqrtOfSumPearson}, {sqrtOfSumSpearman}, {sqrtOfSumSlope}, {sqrtOfSumIntercept}, {sqrtOfSumRSquared}");
            }

            TestContext.Out.WriteLine($"SquaredEuclidean vs Euclidean:           Pearson={squaredPearson:F4}  Spearman={squaredSpearman:F4}");
            TestContext.Out.WriteLine($"  Calibration: Euclidean ≈ {squaredSlope:F4} × SquaredEuclidean + {squaredIntercept:F4}  (R²={squaredRSquared:F4})");
            TestContext.Out.WriteLine($"SquaredEuclideanSqrtOfSum vs Euclidean:  Pearson={sqrtOfSumPearson:F4}  Spearman={sqrtOfSumSpearman:F4}");
            TestContext.Out.WriteLine($"  Calibration: Euclidean ≈ {sqrtOfSumSlope:F4} × SquaredEuclideanSqrtOfSum + {sqrtOfSumIntercept:F4}  (R²={sqrtOfSumRSquared:F4})");

            // Real per-formula timing: each Stopwatch covers only that metric's own 435 CompareX
            // calls, interleaved pair-by-pair with the others (not sequential passes), so all three
            // see identical machine load/JIT-warmup conditions.
            double euclideanMs = euclideanStopwatch.Elapsed.TotalMilliseconds;
            double squaredMs = squaredStopwatch.Elapsed.TotalMilliseconds;
            double squaredSqrtOfSumMs = squaredSqrtOfSumStopwatch.Elapsed.TotalMilliseconds;

            using (var runtime = File.CreateText("DistanceMetricComparison_Runtime.csv"))
            {
                runtime.WriteLine("Metric, PairCount, TotalMs, SpeedupVsEuclidean");
                runtime.WriteLine($"Euclidean, {rows.Length}, {euclideanMs}, 1.0");
                runtime.WriteLine($"SquaredEuclidean, {rows.Length}, {squaredMs}, {euclideanMs / squaredMs}");
                runtime.WriteLine($"SquaredEuclideanSqrtOfSum, {rows.Length}, {squaredSqrtOfSumMs}, {euclideanMs / squaredSqrtOfSumMs}");
            }

            TestContext.Out.WriteLine($"Runtime (435 pairs): Euclidean={euclideanMs:F0}ms  SquaredEuclidean={squaredMs:F0}ms ({euclideanMs / squaredMs:F2}x)  SquaredEuclideanSqrtOfSum={squaredSqrtOfSumMs:F0}ms ({euclideanMs / squaredSqrtOfSumMs:F2}x)");

            // Sanity check that the harness itself is wired correctly: the real baseline should still
            // reproduce the established taxonomic ordering. Deliberately not asserting anything about the
            // two alternatives' ordering or correlation - that is the unknown this diagnostic exists
            // to discover, including whether the length-confound flagged above actually bites.
            var meanSameSpecies = Mean(byLevel[RelationshipLevel.SameSpecies].Select(r => r.Euclidean).ToList());
            var meanSameGenus = Mean(byLevel[RelationshipLevel.SameGenus].Select(r => r.Euclidean).ToList());
            var meanSamePhylum = Mean(byLevel[RelationshipLevel.SamePhylum].Select(r => r.Euclidean).ToList());
            var meanDifferentPhylum = Mean(byLevel[RelationshipLevel.DifferentPhylum].Select(r => r.Euclidean).ToList());
            Assert.That(meanSameSpecies, Is.GreaterThan(meanSameGenus));
            Assert.That(meanSameGenus, Is.GreaterThan(meanSamePhylum));
            Assert.That(meanSamePhylum, Is.GreaterThan(meanDifferentPhylum));
        }
    }
}
