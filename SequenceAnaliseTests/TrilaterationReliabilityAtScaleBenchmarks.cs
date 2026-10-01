using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Reliability evidence: the known-mutation-rate grid (same procedure as
    /// TrilaterationAccuracyBenchmarks) at lengths where Needleman-Wunsch is no longer affordable
    /// as a ground truth. The mutation rate itself is the ground truth here, exactly known by
    /// construction, which is what makes it possible to test reliability at the scale where the
    /// original CAT method's large-N collapse was diagnosed, not just at the small scale NW can
    /// validate.
    /// </summary>
    public class TrilaterationReliabilityAtScaleBenchmarks
    {
        [Test, TestCaseSource(typeof(TrilaterationPerformanceBenchmarks), nameof(TrilaterationPerformanceBenchmarks.PerformanceLengthProvider))]
        public void KnownMutationRateSweep(int length)
        {
            var rawRows = new List<string>();
            var byCondition = new Dictionary<string, List<double>>();

            void Record(string condition, double? mutationRate, int seed, double zcurve)
            {
                rawRows.Add($"{length}, {seed}, {(mutationRate.HasValue ? mutationRate.Value.ToString() : "")}, {condition}, {zcurve}");

                if (!byCondition.TryGetValue(condition, out var scores))
                    byCondition[condition] = scores = new List<double>();
                scores.Add(zcurve);
            }

            for (int seed = 0; seed < AccuracyBenchmarkShared.SeedCount; seed++)
            {
                var baseDna = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "base", seed));
                var baseZc = new ZCurveWindowedProfile(baseDna);

                foreach (var rate in AccuracyBenchmarkShared.MutationRates)
                {
                    var mutated = AccuracyBenchmarkShared.Mutate(baseDna, rate, seed: BenchmarkTiming.Seed(length, "mut", seed, rate));
                    var zc = ZCurveWindowedProfile.Compare(baseZc, new ZCurveWindowedProfile(mutated));
                    Record($"Mutation_{rate:F2}", rate, seed, zc);
                }

                var unrelated = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "unrelated", seed));
                var zcU = ZCurveWindowedProfile.Compare(baseZc, new ZCurveWindowedProfile(unrelated));
                Record("Unrelated", null, seed, zcU);
            }

            using (var raw = File.CreateText($"ReliabilityAtScale_Raw_{length}.csv"))
            {
                raw.WriteLine("Length, Seed, MutationRate, Condition, ZCurveCompare");
                foreach (var row in rawRows) raw.WriteLine(row);
                raw.Flush();
            }

            using (var summary = File.CreateText($"ReliabilityAtScale_Summary_{length}.csv"))
            {
                summary.WriteLine("Length, Condition, SeedsUsed, ZCurve Mean, ZCurve StdDev");
                foreach (var kv in byCondition)
                    summary.WriteLine($"{length}, {kv.Key}, {kv.Value.Count}, {Mean(kv.Value)}, {StdDev(kv.Value)}");
                summary.Flush();
            }

            TestContext.Out.WriteLine($"length={length,9}: wrote ReliabilityAtScale_Raw_{length}.csv and ReliabilityAtScale_Summary_{length}.csv");

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
        public void CombineTrilaterationReliabilityAtScale()
        {
            CombineReports("ReliabilityAtScale_Raw_*.csv", "ReportReliabilityAtScale_Raw.csv",
                "Length, Seed, MutationRate, Condition, ZCurveCompare");
            CombineReports("ReliabilityAtScale_Summary_*.csv", "ReportReliabilityAtScale_Summary.csv",
                "Length, Condition, SeedsUsed, ZCurve Mean, ZCurve StdDev");

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
