using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Measures profile-construction time separately from comparison time - proving the two are
    /// genuinely separable, cacheable operations, not just conceptually so - for both the old
    /// single-point trilateration
    /// (CATProfile) and the new windowed Z-curve trilateration (ZCurveWindowedProfile). Includes
    /// the fragment-localization scenario that has no CAT equivalent - sliding-alignment
    /// localization is itself the new capability.
    /// </summary>
    public class TrilaterationPerformanceBenchmarks
    {
        private const int repetitions = BenchmarkTiming.DefaultRepetitions;
        private const int localizationQueryLength = 5_000;

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

            var catMs = BenchmarkTiming.Measure(repetitions, () => new CATProfile(dna));
            var zcMs = BenchmarkTiming.Measure(repetitions, () => new ZCurveWindowedProfile(dna));

            var catAvg = catMs.Average();
            var zcAvg = zcMs.Average();

            TestContext.Out.WriteLine(
                $"length={length,9}  CAT build avg={catAvg,10:F3} ms ({length / (catAvg / 1000d),14:N0} bases/s)  " +
                $"ZCurve build avg={zcAvg,10:F3} ms ({length / (zcAvg / 1000d),14:N0} bases/s)  ratio(CAT/ZCurve)={catAvg / zcAvg,6:F1}x");

            using var report = File.CreateText($"TrilaterationBuildTime_{length}.csv");
            report.WriteLine("Length, CAT Rep1 ms, CAT Rep2 ms, CAT Rep3 ms, CAT Rep4 ms, CAT Rep5 ms, CAT Avg ms, CAT Bases per sec, ZCurve Rep1 ms, ZCurve Rep2 ms, ZCurve Rep3 ms, ZCurve Rep4 ms, ZCurve Rep5 ms, ZCurve Avg ms, ZCurve Bases per sec");
            report.WriteLine($"{length}, {string.Join(", ", catMs)}, {catAvg}, {length / (catAvg / 1000d)}, {string.Join(", ", zcMs)}, {zcAvg}, {length / (zcAvg / 1000d)}");
            report.Flush();

            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(PerformanceLengthProvider))]
        public void CompareTimeWholeSequence(int length)
        {
            var dnaX = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "compareX"));
            var dnaY = Generator.GetRandomDNA(length, seed: BenchmarkTiming.Seed(length, "compareY"));

            var catX = new CATProfile(dnaX);
            var catY = new CATProfile(dnaY);
            var zcX = new ZCurveWindowedProfile(dnaX);
            var zcY = new ZCurveWindowedProfile(dnaY);

            var catMs = BenchmarkTiming.Measure(repetitions, () => CATProfile.Compare(catX, catY));
            var zcMs = BenchmarkTiming.Measure(repetitions, () => ZCurveWindowedProfile.Compare(zcX, zcY));

            var catAvg = catMs.Average();
            var zcAvg = zcMs.Average();

            TestContext.Out.WriteLine(
                $"length={length,9}  CAT compare avg={catAvg,10:F5} ms  ZCurve compare avg={zcAvg,10:F5} ms  windowCount={zcX.WindowCount}");

            using var report = File.CreateText($"TrilaterationCompareTime_WholeSequence_{length}.csv");
            report.WriteLine("Length, WindowCount, CAT Rep1 ms, CAT Rep2 ms, CAT Rep3 ms, CAT Rep4 ms, CAT Rep5 ms, CAT Avg ms, ZCurve Rep1 ms, ZCurve Rep2 ms, ZCurve Rep3 ms, ZCurve Rep4 ms, ZCurve Rep5 ms, ZCurve Avg ms");
            report.WriteLine($"{length}, {zcX.WindowCount}, {string.Join(", ", catMs)}, {catAvg}, {string.Join(", ", zcMs)}, {zcAvg}");
            report.Flush();

            Assert.Pass();
        }

        // The scenario with no CAT equivalent: a short query trajectory searched against a much
        // longer reference trajectory via sliding alignment. This IS "the new trilateration" -
        // CAT's single point has nothing analogous to localize with, only whole-vs-whole distance.
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
        // ceiling almost exactly. ZCurve-only (no CATProfile - see class remarks) since CAT would
        // add several minutes for data outside current scope with no bearing on this ceiling.
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
            CombineReports("TrilaterationBuildTime_*.csv", "ReportTrilaterationBuildTime.csv",
                "Length, CAT Rep1 ms, CAT Rep2 ms, CAT Rep3 ms, CAT Rep4 ms, CAT Rep5 ms, CAT Avg ms, CAT Bases per sec, ZCurve Rep1 ms, ZCurve Rep2 ms, ZCurve Rep3 ms, ZCurve Rep4 ms, ZCurve Rep5 ms, ZCurve Avg ms, ZCurve Bases per sec");
            CombineReports("TrilaterationCompareTime_WholeSequence_*.csv", "ReportTrilaterationCompareTime_WholeSequence.csv",
                "Length, WindowCount, CAT Rep1 ms, CAT Rep2 ms, CAT Rep3 ms, CAT Rep4 ms, CAT Rep5 ms, CAT Avg ms, ZCurve Rep1 ms, ZCurve Rep2 ms, ZCurve Rep3 ms, ZCurve Rep4 ms, ZCurve Rep5 ms, ZCurve Avg ms");
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
