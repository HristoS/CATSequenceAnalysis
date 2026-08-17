using System.Diagnostics;
using System.Globalization;

namespace SequenceAnaliseTests
{
    internal static class BenchmarkTiming
    {
        public const int DefaultRepetitions = 5;

        // System.HashCode.Combine (and string.GetHashCode()) are deliberately randomized per
        // process - a DoS-resistance measure, not a bug - so using either to derive a seed for
        // `new Random(...)` gives a DIFFERENT effective seed on every test run despite looking
        // "seeded". This is a plain deterministic combinator (FNV-1a over each part's invariant
        // string form, with an explicit delimiter between parts) used only for seed derivation,
        // so the same (length, tag, ...) always produces the same Random sequence across runs -
        // required for the reproducibility this benchmark suite's numbers are reported under.
        public static int Seed(params object[] parts)
        {
            unchecked
            {
                const uint fnvPrime = 16777619;
                uint hash = 2166136261;
                foreach (var part in parts)
                {
                    var text = part is double d ? d.ToString("R", CultureInfo.InvariantCulture) : part?.ToString() ?? "";
                    foreach (var c in text)
                    {
                        hash ^= c;
                        hash *= fnvPrime;
                    }
                    hash ^= '|';
                    hash *= fnvPrime;
                }
                return (int)hash;
            }
        }

        // One untimed warm-up call (removes first-call JIT-tiering skew from the measurement),
        // then `repetitions` individually-timed calls. Returns every rep's elapsed ms so callers
        // can report both the average and the raw per-rep spread - same shape as
        // GeneratorTests.GeneratorSpeed's convention, plus this warm-up call.
        public static double[] Measure(int repetitions, Func<object?> action)
        {
            action();

            var elapsedMs = new double[repetitions];
            for (int i = 0; i < repetitions; i++)
            {
                var sw = Stopwatch.StartNew();
                action();
                sw.Stop();
                elapsedMs[i] = sw.Elapsed.TotalMilliseconds;
            }
            return elapsedMs;
        }
    }
}
