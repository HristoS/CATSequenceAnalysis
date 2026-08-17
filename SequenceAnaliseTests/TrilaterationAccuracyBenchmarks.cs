using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Small-scale accuracy evidence: compares CATProfile (old trilateration), ZCurveWindowedProfile
    /// (new trilateration), and Needleman-Wunsch (an independent alignment-based ground truth)
    /// against each other across a range of known mutation rates. NW is O(n*m), so this file first
    /// calibrates NW's actual per-call cost on this machine (NWCalibrationTiming) and only commits
    /// to a full 5-seed grid at lengths where the resulting (8 mutation rates + 1 unrelated) x
    /// 5-seed = 45-call grid is affordable - see the Verdict column it writes.
    ///
    /// NWGroundTruthLengthProvider/NWSpotCheckLengthProvider are populated from those verdicts
    /// after NWCalibrationTiming has actually been run once (a small manual edit - NUnit has no
    /// clean way for one test's runtime output to feed another test's TestCaseSource within a
    /// single run).
    ///
    /// CAT column note: this uses CATProfile.CompareTrilateration, not the default
    /// CATProfile.Compare() (CompareDistances). Compare() currently returns NaN unconditionally
    /// under the working tree's uncommitted Benchmark.Match rewrite - CanContribute() (the
    /// triangle-inequality precondition CompareDistances' sqrt depends on) was measured False at
    /// every length from 10 to 20,000, a regression from the committed behavior CLAUDE.md
    /// documents (held up to length 10,000,000). CompareTrilateration has no such sqrt and remains
    /// numerically sound, so it's used here as the working CAT baseline pending a fix upstream.
    /// </summary>
    public class TrilaterationAccuracyBenchmarks
    {
        private const int GridConditionCount = 9; // 8 mutation rates + 1 "unrelated"
        private const double GoThresholdSeconds = 60;
        private const double SpotCheckThresholdSeconds = 600;

        public static IEnumerable<int> NWCalibrationLengthProvider()
        {
            yield return 500;
            yield return 1_000;
            yield return 2_000;
            yield return 5_000;
            yield return 10_000;
            yield return 20_000;
        }

        // Lengths cheap enough to always run the full 5-seed grid, regardless of calibration.
        public static IEnumerable<int> NWGroundTruthLengthProvider()
        {
            yield return 500;
            yield return 1_000;
            yield return 2_000;
            yield return 5_000;
        }

        // Populated from NWCalibrationTiming's actual measured verdicts (ReportNWCalibration.csv):
        // 10,000 -> 101.3s projected, 20,000 -> 454.5s projected, both "SpotCheck" (<=600s) while
        // 500-5,000 projected <=39.3s each ("Go", handled by NWGroundTruthLengthProvider instead).
        public static IEnumerable<int> NWSpotCheckLengthProvider()
        {
            yield return 10_000;
            yield return 20_000;
        }

        [Test, TestCaseSource(nameof(NWCalibrationLengthProvider))]
        public void NWCalibrationTiming(int length)
        {
            var refSeq = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "nwcal-ref"));
            var alignSeq = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "nwcal-align"));

            var ms = BenchmarkTiming.Measure(BenchmarkTiming.DefaultRepetitions, () => NeedlemanWunsch.Calculate(refSeq, alignSeq));
            var avg = ms.Average();
            var projectedFullGridSeconds = avg * GridConditionCount * AccuracyBenchmarkShared.SeedCount / 1000d;

            string verdict = projectedFullGridSeconds <= GoThresholdSeconds ? "Go"
                : projectedFullGridSeconds <= SpotCheckThresholdSeconds ? "SpotCheck"
                : "Skip";

            TestContext.Out.WriteLine($"length={length,7}  NW avg={avg,10:F3} ms  projectedFullGridSeconds={projectedFullGridSeconds,10:F1}  verdict={verdict}");

            using var report = File.CreateText($"NWCalibration_{length}.csv");
            report.WriteLine("Length, Rep1 ms, Rep2 ms, Rep3 ms, Rep4 ms, Rep5 ms, AvgMs, ProjectedFullGridSeconds, Verdict");
            report.WriteLine($"{length}, {string.Join(", ", ms)}, {avg}, {projectedFullGridSeconds}, {verdict}");
            report.Flush();

            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(NWGroundTruthLengthProvider))]
        public void AccuracySweep(int length) => RunAccuracyGrid(length, tier: "Full", seedCount: AccuracyBenchmarkShared.SeedCount);

        [Test, TestCaseSource(nameof(NWSpotCheckLengthProvider))]
        public void AccuracySpotCheck(int length) => RunAccuracyGrid(length, tier: "SpotCheck", seedCount: 1);

        private static void RunAccuracyGrid(int length, string tier, int seedCount)
        {
            var rawRows = new List<string>();
            var byCondition = new Dictionary<string, (List<double> Cat, List<double> ZCurve, List<double> NW)>();

            void Record(string condition, double? mutationRate, int seed, double cat, double zcurve, double nw)
            {
                rawRows.Add($"{length}, {tier}, {seed}, {(mutationRate.HasValue ? mutationRate.Value.ToString() : "")}, {condition}, {cat}, {zcurve}, {nw}");

                if (!byCondition.TryGetValue(condition, out var lists))
                {
                    lists = (new List<double>(), new List<double>(), new List<double>());
                    byCondition[condition] = lists;
                }
                lists.Cat.Add(cat);
                lists.ZCurve.Add(zcurve);
                lists.NW.Add(nw);
            }

            for (int seed = 0; seed < seedCount; seed++)
            {
                var baseDna = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "base", seed));
                var baseCat = new CATProfile(baseDna);
                var baseZc = new ZCurveWindowedProfile(baseDna);

                foreach (var rate in AccuracyBenchmarkShared.MutationRates)
                {
                    var mutated = AccuracyBenchmarkShared.Mutate(baseDna, rate, seed: BenchmarkTiming.Seed(length, "mut", seed, rate));
                    // CATProfile.Compare() (CompareDistances) returns NaN under the current
                    // uncommitted Benchmark.Match rewrite - CanContribute() is false at every
                    // length tested (10-20,000), violating the triangle-inequality precondition
                    // CompareDistances' sqrt relies on. CompareTrilateration avoids that sqrt
                    // entirely (pure linear algebra) and is the numerically-sound CAT baseline
                    // until that regression is fixed - see the article's evidence section.
                    var cat = CATProfile.CompareTrilateration(baseCat, new CATProfile(mutated));
                    var zc = ZCurveWindowedProfile.Compare(baseZc, new ZCurveWindowedProfile(mutated));
                    var nw = NeedlemanWunsch.Calculate(baseDna, mutated);
                    Record($"Mutation_{rate:F2}", rate, seed, cat, zc, nw);
                }

                var unrelated = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "unrelated", seed));
                var catU = CATProfile.CompareTrilateration(baseCat, new CATProfile(unrelated));
                var zcU = ZCurveWindowedProfile.Compare(baseZc, new ZCurveWindowedProfile(unrelated));
                var nwU = NeedlemanWunsch.Calculate(baseDna, unrelated);
                Record("Unrelated", null, seed, catU, zcU, nwU);
            }

            using (var raw = File.CreateText($"AccuracyGrid_Raw_{tier}_{length}.csv"))
            {
                raw.WriteLine("Length, Tier, Seed, MutationRate, Condition, CATTrilaterationCompare, ZCurveCompare, NWIdentity");
                foreach (var row in rawRows) raw.WriteLine(row);
                raw.Flush();
            }

            using (var summary = File.CreateText($"AccuracyGrid_Summary_{tier}_{length}.csv"))
            {
                summary.WriteLine("Length, Tier, Condition, SeedsUsed, CAT(Trilateration) Mean, CAT(Trilateration) StdDev, ZCurve Mean, ZCurve StdDev, NW Mean, NW StdDev");
                foreach (var kv in byCondition)
                {
                    var (cat, zc, nw) = kv.Value;
                    summary.WriteLine($"{length}, {tier}, {kv.Key}, {cat.Count}, {Mean(cat)}, {StdDev(cat)}, {Mean(zc)}, {StdDev(zc)}, {Mean(nw)}, {StdDev(nw)}");
                }
                summary.Flush();
            }

            TestContext.Out.WriteLine($"length={length,7} tier={tier} seeds={seedCount}: wrote AccuracyGrid_Raw_{tier}_{length}.csv and AccuracyGrid_Summary_{tier}_{length}.csv");

            Assert.Pass();
        }

        private static double Mean(List<double> values) => values.Average();

        private static double StdDev(List<double> values)
        {
            if (values.Count < 2) return 0;
            var mean = values.Average();
            var sumSquares = values.Sum(v => (v - mean) * (v - mean));
            return Math.Sqrt(sumSquares / (values.Count - 1));
        }

        [Test]
        public void CombineTrilaterationAccuracy()
        {
            CombineReports("AccuracyGrid_Raw_*.csv", "ReportTrilaterationAccuracy_Raw.csv",
                "Length, Tier, Seed, MutationRate, Condition, CATTrilaterationCompare, ZCurveCompare, NWIdentity");
            CombineReports("AccuracyGrid_Summary_*.csv", "ReportTrilaterationAccuracy_Summary.csv",
                "Length, Tier, Condition, SeedsUsed, CAT(Trilateration) Mean, CAT(Trilateration) StdDev, ZCurve Mean, ZCurve StdDev, NW Mean, NW StdDev");
            CombineReports("NWCalibration_*.csv", "ReportNWCalibration.csv",
                "Length, Rep1 ms, Rep2 ms, Rep3 ms, Rep4 ms, Rep5 ms, AvgMs, ProjectedFullGridSeconds, Verdict");

            Assert.Pass();
        }

        private static void CombineReports(string sourceGlob, string outputName, string header)
        {
            using var report = File.CreateText(outputName);
            report.WriteLine(header);
            foreach (var file in Directory.GetFiles(Directory.GetCurrentDirectory(), sourceGlob).OrderBy(f => f))
                foreach (var line in File.ReadAllLines(file).Skip(1))
                    report.WriteLine(line);

            report.Flush();
        }
    }
}
