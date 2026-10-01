using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Measures ZCurveWindowedProfile construction time separately from comparison time - proving
    /// the two are genuinely separable, cacheable operations, not just conceptually so. Includes
    /// the fragment-localization scenario (sliding-alignment of a short query against a long
    /// reference).
    /// </summary>
    public class TrilaterationPerformanceBenchmarks
    {
        private const int repetitions = BenchmarkTiming.DefaultRepetitions;
        private const int localizationQueryLength = 5_000;
        private const string buildTimeHeader = "Length, ZCurve Rep1 ms, ZCurve Rep2 ms, ZCurve Rep3 ms, ZCurve Rep4 ms, ZCurve Rep5 ms, ZCurve Avg ms, ZCurve Bases per sec";
        private const string compareTimeHeader = "Length, WindowCount, ZCurve Rep1 ms, ZCurve Rep2 ms, ZCurve Rep3 ms, ZCurve Rep4 ms, ZCurve Rep5 ms, ZCurve Avg ms";

        public static IEnumerable<int> PerformanceLengthProvider()
        {
            yield return 20_000;
            yield return 40_000;
            yield return 80_000;
            yield return 160_000;
            yield return 1_280_000;
            yield return 2_560_000;
            yield return 5_120_000;
            yield return 10_240_000;
        }

        [Test, TestCaseSource(nameof(PerformanceLengthProvider))]
        public void BuildTime(int length)
        {
            var dna = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "build"));

            var zcMs = BenchmarkTiming.Measure(repetitions, () => new ZCurveWindowedProfile(dna));
            var zcAvg = zcMs.Average();

            TestContext.Out.WriteLine(
                $"length={length,9}  ZCurve build avg={zcAvg,10:F3} ms ({length / (zcAvg / 1000d),14:N0} bases/s)");

            using var report = File.CreateText($"TrilaterationBuildTime_{length}.csv");
            report.WriteLine(buildTimeHeader);
            report.WriteLine($"{length}, {string.Join(", ", zcMs)}, {zcAvg}, {length / (zcAvg / 1000d)}");
            report.Flush();

            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(PerformanceLengthProvider))]
        public void CompareTimeWholeSequence(int length)
        {
            var dnaX = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "compareX"));
            var dnaY = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "compareY"));

            var zcX = new ZCurveWindowedProfile(dnaX);
            var zcY = new ZCurveWindowedProfile(dnaY);

            var zcMs = BenchmarkTiming.Measure(repetitions, () => ZCurveWindowedProfile.Compare(zcX, zcY));
            var zcAvg = zcMs.Average();

            TestContext.Out.WriteLine(
                $"length={length,9}  ZCurve compare avg={zcAvg,10:F5} ms  windowCount={zcX.WindowCount}");

            using var report = File.CreateText($"TrilaterationCompareTime_WholeSequence_{length}.csv");
            report.WriteLine(compareTimeHeader);
            report.WriteLine($"{length}, {zcX.WindowCount}, {string.Join(", ", zcMs)}, {zcAvg}");
            report.Flush();

            Assert.Pass();
        }

        // A short query trajectory searched against a much longer reference trajectory via sliding
        // alignment - a single whole-sequence point has nothing analogous to localize with.
        [Test, TestCaseSource(nameof(PerformanceLengthProvider))]
        public void NewTrilaterationLocalizationTime(int referenceLength)
        {
            var query = Generator.GetRandomDNA(localizationQueryLength, seed: BenchmarkTiming.Seed(referenceLength, "query"));
            var reference = Generator.GetRandomDNA(referenceLength, seed: BenchmarkTiming.Seed(referenceLength, "reference"));

            var queryProfile = new ZCurveWindowedProfile(query);
            var referenceProfile = new ZCurveWindowedProfile(reference);

            var ms = BenchmarkTiming.Measure(repetitions, () => ZCurveWindowedProfile.BestAlignment(queryProfile, referenceProfile));
            var (score, offset) = ZCurveWindowedProfile.BestAlignment(queryProfile, referenceProfile);

            var avg = ms.Average();
            int offsetsSearched = referenceProfile.WindowCount - queryProfile.WindowCount + 1;
            double windowComparisonsPerSecond = (double)offsetsSearched * queryProfile.WindowCount / (avg / 1000d);

            TestContext.Out.WriteLine(
                $"referenceLength={referenceLength,9}  queryWindows={queryProfile.WindowCount}  referenceWindows={referenceProfile.WindowCount}  " +
                $"offsetsSearched={offsetsSearched}  avg={avg,10:F3} ms  ({windowComparisonsPerSecond,14:N0} window-comparisons/s)  offset={offset}  score={score}");

            using var report = File.CreateText($"TrilaterationLocalizationTime_{referenceLength}.csv");
            report.WriteLine("ReferenceLength, QueryLength, QueryWindowCount, ReferenceWindowCount, OffsetsSearched, Rep1 ms, Rep2 ms, Rep3 ms, Rep4 ms, Rep5 ms, Avg ms, WindowComparisonsPerSec, Offset");
            report.WriteLine($"{referenceLength}, {localizationQueryLength}, {queryProfile.WindowCount}, {referenceProfile.WindowCount}, {offsetsSearched}, {string.Join(", ", ms)}, {avg}, {windowComparisonsPerSecond}, {offset}");
            report.Flush();

            Assert.Pass();
        }

        // One-off spot check, not part of the systematic length sweep: how close to genome scale
        // does the current string-based representation actually reach? A literal 3,100,000,000 bp
        // (~human genome scale) request doesn't even compile against Generator.GetRandomDNA's
        // `int length` parameter (int.MaxValue = 2,147,483,647), and empirically, allocating a
        // char[3_100_000_000] itself throws OverflowException on this runtime - before the
        // documented .NET string cap would even apply. That cap was independently confirmed here:
        // string construction succeeds through 1,070,000,000 chars but throws OutOfMemoryException
        // starting at 1,073,741,823 (2^30 - 1), matching CLAUDE.md's documented ~1.07 billion char
        // ceiling almost exactly.
        [Test]
        public void GenomeScaleCeilingProbe()
        {
            const int length = 1_070_000_000; // just under the empirically-confirmed ~1.073B string cap

            var dna = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "genome-scale"));

            var buildMs = BenchmarkTiming.Measure(1, () => new ZCurveWindowedProfile(dna));
            var profile = new ZCurveWindowedProfile(dna);
            var other = new ZCurveWindowedProfile(Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "genome-scale-other")));
            var compareMs = BenchmarkTiming.Measure(1, () => ZCurveWindowedProfile.Compare(profile, other));

            TestContext.Out.WriteLine(
                $"length={length:N0}  build={buildMs[0]:F1} ms ({length / (buildMs[0] / 1000d):N0} bases/s)  " +
                $"windows={profile.WindowCount:N0}  wholeVsWholeCompare={compareMs[0]:F3} ms");

            using var report = File.CreateText("GenomeScaleCeilingProbe.csv");
            report.WriteLine("Length, BuildMs, WindowCount, CompareMs");
            report.WriteLine($"{length}, {buildMs[0]}, {profile.WindowCount}, {compareMs[0]}");
            report.Flush();

            Assert.Pass();
        }

        [Test]
        public void CombineTrilaterationPerformanceReports()
        {
            CombineReports("TrilaterationBuildTime_*.csv", "ReportTrilaterationBuildTime.csv", buildTimeHeader);
            CombineReports("TrilaterationCompareTime_WholeSequence_*.csv", "ReportTrilaterationCompareTime_WholeSequence.csv", compareTimeHeader);
            CombineReports("TrilaterationLocalizationTime_*.csv", "ReportTrilaterationLocalizationTime.csv",
                "ReferenceLength, QueryLength, QueryWindowCount, ReferenceWindowCount, OffsetsSearched, Rep1 ms, Rep2 ms, Rep3 ms, Rep4 ms, Rep5 ms, Avg ms, WindowComparisonsPerSec, Offset");

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
