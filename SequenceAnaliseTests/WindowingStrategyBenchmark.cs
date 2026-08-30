using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Reproduces the manuscript's §5.3 comparison table: the same unrelated-pair drift experiment
    /// (a short length vs. a much longer length) run against three different ways of compressing a
    /// per-base score into a profile - the raw endpoint (ZCurveProfile.Compare, §4.2), a fixed
    /// checkpoint *count* (ZCurveProfile.CompareTrajectory, §5.1), and a fixed window *size*
    /// (ZCurveWindowedProfile.Compare, §5.2). Previously these three numbers only existed as
    /// TestContext.Out.WriteLine output and code comments spread across ZCurveProfileTests.cs and
    /// ZCurveWindowedProfileTests.cs (docs/zcurve-windowed-profile-factcheck.md, discrepancy 4);
    /// this consolidates them into one checked-in CSV, mirroring the fix already applied to §8.5 in
    /// GenomeScaleStreamingTests.cs. Lengths and seeds match the pre-existing tests exactly:
    /// UnrelatedRandomSequencesOfSameLengthScoreNearOne (endpoint),
    /// TrajectoryWithFixedCheckpointCountStillCollapsesAtLargeEnoughLength (checkpoint count), and
    /// FixedWindowSizeAvoidsTheLargeNCollapse (window size).
    ///
    /// Also covers §5.1's separate "~0.66" shuffled-sequence figure (same discrepancy 4: "same
    /// applies to §5.1's k = 20 shuffled-sequence score"), matching
    /// ZCurveProfileTests.TrajectoryDetectsRegionalStructureThatTheEndpointCannotSee exactly.
    /// </summary>
    public class WindowingStrategyBenchmark
    {
        private static string BiasedDNA(char[] weightedBases, int length, int seed)
        {
            var rng = new Random(seed);
            return new string(Enumerable.Range(0, length).Select(_ => weightedBases[rng.Next(weightedBases.Length)]).ToArray());
        }

        // Fisher-Yates shuffle: same multiset of bases, different order/positions.
        private static string Shuffle(string s, int seed)
        {
            var rng = new Random(seed);
            var chars = s.ToCharArray();
            for (int i = chars.Length - 1; i > 0; i--)
            {
                int j = rng.Next(i + 1);
                (chars[i], chars[j]) = (chars[j], chars[i]);
            }
            return new string(chars);
        }

        [Test]
        public void UnrelatedPairDriftAcrossThreeWindowingStrategies()
        {
            var rows = new List<(string Construction, int N, double Score)>();

            // §4.2: raw endpoint, no windowing at all - a 10x length range (this row's own range in
            // the manuscript; unlike the other two rows below, no wider synthetic range was claimed
            // for it).
            var endpointShort = ZCurveProfile.Compare(
                new ZCurveProfile(Generator.GetRandomDNA(5_000, seed: 1)),
                new ZCurveProfile(Generator.GetRandomDNA(5_000, seed: 2)));
            var endpointLong = ZCurveProfile.Compare(
                new ZCurveProfile(Generator.GetRandomDNA(50_000, seed: 1)),
                new ZCurveProfile(Generator.GetRandomDNA(50_000, seed: 2)));
            rows.Add(("Endpoint", 5_000, endpointShort));
            rows.Add(("Endpoint", 50_000, endpointLong));

            // §5.1: fixed checkpoint *count* = 20, a 100x length range.
            const int checkpointCount = 20;
            var checkpointShort = ZCurveProfile.CompareTrajectory(
                new ZCurveProfile(Generator.GetRandomDNA(2_000, seed: 1), checkpointCount),
                new ZCurveProfile(Generator.GetRandomDNA(2_000, seed: 2), checkpointCount));
            var checkpointLong = ZCurveProfile.CompareTrajectory(
                new ZCurveProfile(Generator.GetRandomDNA(200_000, seed: 1), checkpointCount),
                new ZCurveProfile(Generator.GetRandomDNA(200_000, seed: 2), checkpointCount));
            rows.Add(("FixedCheckpointCount20", 2_000, checkpointShort));
            rows.Add(("FixedCheckpointCount20", 200_000, checkpointLong));

            // §5.2: fixed window *size* = 500 bp, the same 100x length range as the row above.
            var windowShort = ZCurveWindowedProfile.Compare(
                new ZCurveWindowedProfile(Generator.GetRandomDNA(2_000, seed: 1)),
                new ZCurveWindowedProfile(Generator.GetRandomDNA(2_000, seed: 2)));
            var windowLong = ZCurveWindowedProfile.Compare(
                new ZCurveWindowedProfile(Generator.GetRandomDNA(200_000, seed: 1)),
                new ZCurveWindowedProfile(Generator.GetRandomDNA(200_000, seed: 2)));
            rows.Add(("FixedWindowSize500", 2_000, windowShort));
            rows.Add(("FixedWindowSize500", 200_000, windowLong));

            using (var report = File.CreateText("ReportWindowingStrategy.csv"))
            {
                report.WriteLine("Construction, N, UnrelatedPairScore");
                foreach (var row in rows)
                    report.WriteLine($"{row.Construction}, {row.N}, {row.Score}");
            }

            foreach (var row in rows)
                TestContext.Out.WriteLine($"{row.Construction,-24} n={row.N,7}  Compare={row.Score:F4}");

            var endpointDrift = endpointLong - endpointShort;
            var checkpointDrift = checkpointLong - checkpointShort;
            var windowDrift = windowLong - windowShort;
            TestContext.Out.WriteLine(
                $"drift: endpoint={endpointDrift:+0.000;-0.000}  " +
                $"checkpointCount20={checkpointDrift:+0.000;-0.000}  " +
                $"windowSize500={windowDrift:+0.000;-0.000}");

            // The central §5.3 claim: decoupling positional resolution from total sequence length
            // (fixed window *size*) keeps the drift an order of magnitude smaller than fixing
            // checkpoint *count* does, over the same 100x length range.
            Assert.That(Math.Abs(windowDrift), Is.LessThan(Math.Abs(checkpointDrift)),
                "fixed window size should drift far less than fixed checkpoint count over the same length range");
        }

        [Test]
        public void StructuredSequenceShuffleDetectionAtFixedCheckpointCount()
        {
            // §5.1: an AT-rich-then-GC-rich two-block sequence, shuffled. The endpoint (§4.2)
            // cannot see the shuffle at all (same composition survives any permutation); the
            // k=20 trajectory (§5.1) does, partially - the "~0.66" figure this reproduces.
            char[] atRich = { 'A', 'A', 'T', 'T', 'A', 'G' }; // ~67% AT
            char[] gcRich = { 'G', 'C', 'G', 'C', 'A', 'T' }; // ~67% GC
            var structured = BiasedDNA(atRich, 2500, seed: 10) + BiasedDNA(gcRich, 2500, seed: 11);
            var shuffled = Shuffle(structured, seed: 20);

            const int checkpointCount = 20;
            var original = new ZCurveProfile(structured, checkpointCount);
            var reordered = new ZCurveProfile(shuffled, checkpointCount);

            var endpointScore = ZCurveProfile.Compare(original, reordered);
            var trajectoryScore = ZCurveProfile.CompareTrajectory(original, reordered);

            using (var report = File.CreateText("ReportRegionalStructureDetection.csv"))
            {
                report.WriteLine("Metric, CheckpointCount, Score");
                report.WriteLine($"EndpointCompare, 1, {endpointScore}");
                report.WriteLine($"TrajectoryCompare, {checkpointCount}, {trajectoryScore}");
            }

            TestContext.Out.WriteLine($"endpoint Compare (blind to order): {endpointScore}");
            TestContext.Out.WriteLine($"trajectory CompareTrajectory ({checkpointCount} checkpoints): {trajectoryScore}");

            Assert.That(endpointScore, Is.EqualTo(1.0).Within(1e-9),
                "endpoint still can't see the shuffle - same ceiling as before");
            Assert.That(trajectoryScore, Is.LessThan(endpointScore),
                "the trajectory should notice the regional AT->GC trend was destroyed, unlike the endpoint");
        }
    }
}
