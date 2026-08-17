using SequenceAnalyses;
using System.Diagnostics;

namespace SequenceAnaliseTests
{
    public class GeneratorTests
    {
        public static IEnumerable<int> GeneratorLengthProvider()
        {
            yield return 1_000;
            yield return 10_000;
            yield return 100_000;
            yield return 1_000_000;
        }

        [Test, TestCaseSource(nameof(GeneratorLengthProvider))]
        public void GeneratorSpeed(int length)
        {
            const int repetitions = 5;
            var elapsedMs = new List<double>();

            for (int i = 0; i < repetitions; i++)
            {
                var sw = Stopwatch.StartNew();
                var dna = Generator.GetRandomDNA(length);
                sw.Stop();

                Assert.That(dna.Length, Is.EqualTo(length));
                Assert.That(dna.All(c => c is 'A' or 'C' or 'G' or 'T'), Is.True);

                elapsedMs.Add(sw.Elapsed.TotalMilliseconds);
            }

            var averageMs = elapsedMs.Average();
            var basesPerSecond = length / (averageMs / 1000d);

            TestContext.Out.WriteLine($"length={length,9}  avg={averageMs,10:F3} ms  throughput={basesPerSecond,14:N0} bases/s  (per run: {string.Join(", ", elapsedMs.Select(x => x.ToString("F2")))})");
        }
    }
}
