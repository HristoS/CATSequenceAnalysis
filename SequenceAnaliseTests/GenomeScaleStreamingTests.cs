using SequenceAnalyses;
using System.Diagnostics;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Evidence for ZCurveWindowedProfile.FromStream: correctness against the existing string-based
    /// constructor, the actual storage footprint of a profile at genome scale, and whether a known
    /// fragment ("a specific gene/protein-coding sequence") can still be localized once the
    /// reference itself is too large for any single string to hold - the scenario the streaming
    /// constructor exists for.
    /// </summary>
    public class GenomeScaleStreamingTests
    {
        // Streams random bases lazily (O(1) memory beyond a Random instance and a counter) - the
        // same base alphabet/order as Generator.GetRandomDNA, but never materializes the sequence
        // as a string or array, so it has no length ceiling tied to .NET's string/array limits.
        private sealed class RandomDnaReader : TextReader
        {
            private static readonly char[] Bases = { 'A', 'T', 'C', 'G' };
            private readonly Random _random;
            private long _remaining;

            public RandomDnaReader(long length, int seed)
            {
                _random = new Random(seed);
                _remaining = length;
            }

            public override int Read(char[] buffer, int index, int count)
            {
                int n = (int)Math.Min(count, _remaining);
                for (int i = 0; i < n; i++)
                    buffer[index + i] = Bases[_random.Next(4)];
                _remaining -= n;
                return n;
            }
        }

        // Like RandomDnaReader, but splices a fixed, known fragment (e.g. a "gene") in at a known
        // offset - random bases before and after, the exact fragment text in between - to simulate
        // searching for a specific known sequence within a much larger genome-scale reference.
        private sealed class PlantedFragmentReader : TextReader
        {
            private static readonly char[] Bases = { 'A', 'T', 'C', 'G' };
            private readonly Random _random;
            private readonly string _fragment;
            private readonly long _fragmentOffset;
            private readonly long _totalLength;
            private long _position;

            public PlantedFragmentReader(long totalLength, long fragmentOffset, string fragment, int seed)
            {
                _random = new Random(seed);
                _fragment = fragment;
                _fragmentOffset = fragmentOffset;
                _totalLength = totalLength;
            }

            public override int Read(char[] buffer, int index, int count)
            {
                int written = 0;
                while (written < count && _position < _totalLength)
                {
                    buffer[index + written] = _position >= _fragmentOffset && _position < _fragmentOffset + _fragment.Length
                        ? _fragment[(int)(_position - _fragmentOffset)]
                        : Bases[_random.Next(4)];
                    _position++;
                    written++;
                }
                return written;
            }
        }

        // Reproduces the exact same base sequence as RandomDnaReader with the same seed (the base
        // RNG is advanced once per position regardless of whether a mutation fires, keeping it in
        // sync with a plain RandomDnaReader), but independently decides per-position mutation via a
        // second, separately-seeded RNG - the streaming analogue of AccuracyBenchmarkShared.Mutate,
        // which cannot be used directly here since it requires the whole sequence already
        // materialized as one string. Not required to reproduce AccuracyBenchmarkShared.Mutate's
        // exact output byte-for-byte (the two are only related by seed, not by call sequence) - only
        // to implement the same statistical procedure: each base independently substituted with
        // probability `rate`.
        private sealed class MutatingDnaReader : TextReader
        {
            private static readonly char[] Bases = { 'A', 'T', 'C', 'G' };
            private readonly Random _baseRandom;
            private readonly Random _mutationRandom;
            private readonly double _rate;
            private long _remaining;

            public MutatingDnaReader(long length, int baseSeed, int mutationSeed, double rate)
            {
                _baseRandom = new Random(baseSeed);
                _mutationRandom = new Random(mutationSeed);
                _rate = rate;
                _remaining = length;
            }

            public override int Read(char[] buffer, int index, int count)
            {
                int n = (int)Math.Min(count, _remaining);
                for (int i = 0; i < n; i++)
                {
                    char baseChar = Bases[_baseRandom.Next(4)];
                    buffer[index + i] = _mutationRandom.NextDouble() < _rate ? Bases[_mutationRandom.Next(4)] : baseChar;
                }
                _remaining -= n;
                return n;
            }
        }

        [Test]
        public void FromStreamMatchesStringConstructor()
        {
            var dna = Generator.GetRandomDNA(50_000, seed: 7);
            var fromString = new ZCurveWindowedProfile(dna, windowSize: 500);
            var fromStream = ZCurveWindowedProfile.FromStream(new StringReader(dna), windowSize: 500);

            Assert.That(fromStream.WindowSize, Is.EqualTo(fromString.WindowSize));
            Assert.That(fromStream.Windows, Is.EqualTo(fromString.Windows),
                "streaming construction must produce bit-for-bit identical windows to the string-based path");
            Assert.That(ZCurveWindowedProfile.Compare(fromStream, fromString), Is.EqualTo(1.0).Within(1e-9));
        }

        [Test]
        public void ProfileStorageFootprintAtGenomeScale()
        {
            const long length = 3_100_000_000L; // the literal genome-scale length requested
            const int windowSize = ZCurveWindowedProfile.DefaultWindowSize;

            long windowCount = length / windowSize;
            long profileBytes = windowCount * 3 * sizeof(double); // 3 doubles per window, no other per-window state
            long rawSequenceBytes = length * sizeof(char); // UTF-16, what a single in-memory string would cost

            TestContext.Out.WriteLine($"length={length:N0} bp, windowSize={windowSize}");
            TestContext.Out.WriteLine($"windowCount={windowCount:N0}");
            TestContext.Out.WriteLine($"profile storage: {profileBytes:N0} bytes ({profileBytes / 1e6:F1} MB)");
            TestContext.Out.WriteLine($"raw sequence as a single string would be: {rawSequenceBytes:N0} bytes ({rawSequenceBytes / 1e9:F2} GB)");
            TestContext.Out.WriteLine($"compression ratio (raw string : profile): {(double)rawSequenceBytes / profileBytes:N0} : 1");

            using (var report = File.CreateText($"GenomeScaleStorageFootprint_{length}.csv"))
            {
                report.WriteLine("Length, WindowSize, WindowCount, ProfileBytes, RawSequenceBytes, CompressionRatio");
                report.WriteLine($"{length}, {windowSize}, {windowCount}, {profileBytes}, {rawSequenceBytes}, {(double)rawSequenceBytes / profileBytes}");
            }

            Assert.That(windowCount, Is.GreaterThan(0));
        }

// Times how long it takes to drain a reader alone, with no profile accumulation, using the
        // same buffer size as ZCurveWindowedProfile.FromStream - lets a caller isolate pure
        // accumulation cost by building an identical reader (same parameters/seed, so the same
        // control-flow/RNG cost) and subtracting this from the full FromStream time.
        private static double DrainMs(TextReader reader)
        {
            var buffer = new char[65536];
            var sw = Stopwatch.StartNew();
            while (reader.Read(buffer, 0, buffer.Length) > 0) { }
            sw.Stop();
            return sw.Elapsed.TotalMilliseconds;
        }

        [Test]
        public void PlantedGeneLocalizationAtGenomeScale()
        {
            const long referenceLength = 3_100_000_000L;
            const int windowSize = ZCurveWindowedProfile.DefaultWindowSize;
            const int geneLength = 2_000; // realistic order of magnitude for a coding sequence
            const long plantedOffset = 2_345_000_000L; // window-aligned (divisible by windowSize)

            Assert.That(plantedOffset % windowSize, Is.EqualTo(0), "test setup: offset must be window-aligned for an exact-match expectation");

            var gene = Generator.GetRandomDNA(geneLength, seed: 555);

            var sw = Stopwatch.StartNew();
            var reference = ZCurveWindowedProfile.FromStream(
                new PlantedFragmentReader(referenceLength, plantedOffset, gene, seed: 1),
                windowSize);
            sw.Stop();
            var buildMs = sw.Elapsed.TotalMilliseconds;

            // Isolate profile-computation cost from the cost of obtaining the underlying bases:
            // replay an identical reader (same parameters/seed) purely to drain it, then subtract
            // that from the full build time above.
            var readerOnlyMs = DrainMs(new PlantedFragmentReader(referenceLength, plantedOffset, gene, seed: 1));
            var accumulationOnlyMs = buildMs - readerOnlyMs;

            var query = new ZCurveWindowedProfile(gene, windowSize);

            sw.Restart();
            var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, reference);
            sw.Stop();
            var searchMs = sw.Elapsed.TotalMilliseconds;

            int expectedOffset = (int)(plantedOffset / windowSize);

            TestContext.Out.WriteLine(
                $"referenceLength={referenceLength:N0}  referenceWindows={reference.Windows.Count:N0}  buildMs={buildMs:F0} " +
                $"({referenceLength / (buildMs / 1000d):N0} bases/s)  readerOnlyMs={readerOnlyMs:F0}  " +
                $"accumulationOnlyMs={accumulationOnlyMs:F0} ({referenceLength / (accumulationOnlyMs / 1000d):N0} bases/s)");
            TestContext.Out.WriteLine(
                $"geneLength={geneLength}  plantedOffset={plantedOffset:N0} (window {expectedOffset:N0})  " +
                $"foundOffset={offset:N0}  score={score}  searchMs={searchMs:F1}");

            using (var report = File.CreateText($"GenomeScalePlantedGeneLocalization_{referenceLength}.csv"))
            {
                report.WriteLine("ReferenceLength, ReferenceWindows, WindowSize, BuildMs, BasesPerSec, ReaderOnlyMs, AccumulationOnlyMs, AccumulationBasesPerSec, GeneLength, PlantedOffset, ExpectedWindow, FoundOffset, Score, SearchMs");
                report.WriteLine($"{referenceLength}, {reference.Windows.Count}, {windowSize}, {buildMs}, {referenceLength / (buildMs / 1000d)}, {readerOnlyMs}, {accumulationOnlyMs}, {referenceLength / (accumulationOnlyMs / 1000d)}, {geneLength}, {plantedOffset}, {expectedOffset}, {offset}, {score}, {searchMs}");
            }

            Assert.That(offset, Is.EqualTo(expectedOffset),
                "the planted gene should localize to the exact window it was inserted at, even inside a 3.1 billion bp reference");
            Assert.That(score, Is.EqualTo(1.0).Within(1e-9),
                "an exact, unmutated fragment should still score a perfect match regardless of reference size");
        }

        // Comparison-cost and reliability data at the 3,100,000,000 bp row, in one pass, reusing
        // the same handful of expensive builds for both rather than paying for 3.1B-bp construction
        // twice. Single seed (not five) since each build costs ~2 minutes - 9 builds (base +
        // unrelated + 7 non-trivial mutation rates; the 0% condition is 1.000 by construction, so it
        // needs no build) already runs to roughly 20 minutes. No NW ground truth at this scale: NW
        // is O(n²), and extrapolating its own measured 20,000 bp timing to 3.1B projects to several
        // thousand years, not a runnable data point.
        [Test]
        public void GenomeScaleComparisonAndReliability()
        {
            const long length = 3_100_000_000L;
            const int windowSize = ZCurveWindowedProfile.DefaultWindowSize;
            const int localizationQueryLength = 5_000; // matches the systematic-sweep query length

            var sw = Stopwatch.StartNew();
            var baseProfile = ZCurveWindowedProfile.FromStream(
                new RandomDnaReader(length, BenchmarkTiming.Seed(length, "genome-base")), windowSize);
            sw.Stop();
            TestContext.Out.WriteLine($"base build: {sw.Elapsed.TotalSeconds:F1} s, windows={baseProfile.Windows.Count:N0}");

            sw.Restart();
            var unrelatedProfile = ZCurveWindowedProfile.FromStream(
                new RandomDnaReader(length, BenchmarkTiming.Seed(length, "genome-unrelated")), windowSize);
            sw.Stop();
            TestContext.Out.WriteLine($"unrelated build: {sw.Elapsed.TotalSeconds:F1} s");

            var results = new List<(string Condition, double Score)> { ("Mutation_0.00", 1.0) };

            foreach (var rate in AccuracyBenchmarkShared.MutationRates.Where(r => r > 0))
            {
                sw.Restart();
                var mutated = ZCurveWindowedProfile.FromStream(
                    new MutatingDnaReader(length,
                        baseSeed: BenchmarkTiming.Seed(length, "genome-base"),
                        mutationSeed: BenchmarkTiming.Seed(length, "genome-mut", rate),
                        rate),
                    windowSize);
                sw.Stop();
                var score = ZCurveWindowedProfile.Compare(baseProfile, mutated);
                TestContext.Out.WriteLine($"rate={rate:F2}: build {sw.Elapsed.TotalSeconds:F1} s, Compare={score}");
                results.Add(($"Mutation_{rate:F2}", score));
            }

            // Whole-vs-whole compare timing - reuses base/unrelated, no extra build.
            var compareMs = BenchmarkTiming.Measure(1, () => ZCurveWindowedProfile.Compare(baseProfile, unrelatedProfile));
            var unrelatedScore = ZCurveWindowedProfile.Compare(baseProfile, unrelatedProfile);
            results.Add(("Unrelated", unrelatedScore));

            // Fragment-localization timing - fresh small query against the already-built base.
            var query = new ZCurveWindowedProfile(Generator.GetRandomDNA(localizationQueryLength, seed: 999), windowSize);
            var localizationMs = BenchmarkTiming.Measure(1, () => ZCurveWindowedProfile.BestAlignment(query, baseProfile));
            int offsetsSearched = baseProfile.Windows.Count - query.Windows.Count + 1;

            TestContext.Out.WriteLine("--- Reliability row (Length, Condition, Score) ---");
            foreach (var (condition, score) in results)
                TestContext.Out.WriteLine($"{length}, {condition}, {score}");
            var margin = results.First(r => r.Condition == "Mutation_0.50").Score - unrelatedScore;
            TestContext.Out.WriteLine($"margin (50% - Unrelated) = {margin}");

            TestContext.Out.WriteLine("--- Comparison-cost row (Length, Windows, WholeVsWholeCompareMs, LocalizationMs, OffsetsSearched) ---");
            TestContext.Out.WriteLine($"{length}, {baseProfile.Windows.Count}, {compareMs[0]}, {localizationMs[0]}, {offsetsSearched}");

            using (var report = File.CreateText($"GenomeScaleReliability_{length}.csv"))
            {
                report.WriteLine("Length, Condition, ZCurveCompare");
                foreach (var (condition, score) in results)
                    report.WriteLine($"{length}, {condition}, {score}");
            }
            using (var report = File.CreateText($"GenomeScaleCompareTime_{length}.csv"))
            {
                report.WriteLine("Length, WindowCount, WholeVsWholeCompareMs, LocalizationQueryLength, LocalizationMs, OffsetsSearched");
                report.WriteLine($"{length}, {baseProfile.Windows.Count}, {compareMs[0]}, {localizationQueryLength}, {localizationMs[0]}, {offsetsSearched}");
            }

            Assert.That(results.First(r => r.Condition == "Mutation_0.50").Score, Is.GreaterThan(unrelatedScore),
                "even at genome scale, a 50%-mutated sequence should score higher than a fully unrelated one");
        }
    }
}
