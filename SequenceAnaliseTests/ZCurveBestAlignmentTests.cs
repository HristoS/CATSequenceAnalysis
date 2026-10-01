using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    // Characterization tests for ZCurveWindowedProfile.BestAlignment/Compare: pin down the current
    // behavior of both SIMD shapes (wide query = across windows, narrow query = across offsets) and
    // their scalar tails before the kernel is changed. Vector256<float>.Count is 8, so query window
    // counts below/at/above 8 and offset counts that are/aren't multiples of 8 are what matter.
    public class ZCurveBestAlignmentTests
    {
        private const double tolerance = 1e-6;

        // Plain double-precision reference: no SIMD, no float accumulation, first best wins on ties.
        private static (double Score, int Offset) OracleBestAlignment(ZCurveWindowedProfile query, ZCurveWindowedProfile reference)
        {
            var scores = OracleScores(query, reference);
            double best = double.NegativeInfinity;
            int bestOffset = -1;
            for (int offset = 0; offset < scores.Length; offset++)
            {
                if (scores[offset] > best)
                {
                    best = scores[offset];
                    bestOffset = offset;
                }
            }
            return bestOffset < 0 ? (double.NaN, -1) : (best, bestOffset);
        }

        private static double[] OracleScores(ZCurveWindowedProfile query, ZCurveWindowedProfile reference)
        {
            var q = query.Windows;
            var r = reference.Windows;
            int m = q.Count;
            int n = r.Count;
            if (m == 0 || n < m) return Array.Empty<double>();

            var scores = new double[n - m + 1];
            for (int offset = 0; offset < scores.Length; offset++)
            {
                double sum = 0;
                for (int j = 0; j < m; j++)
                {
                    double dx = q[j].X - r[offset + j].X;
                    double dy = q[j].Y - r[offset + j].Y;
                    double dz = q[j].Z - r[offset + j].Z;
                    sum += Math.Sqrt(dx * dx + dy * dy + dz * dz);
                }
                scores[offset] = 1 - sum / m;
            }
            return scores;
        }

        private static ZCurveWindowedProfile RandomProfile(int windowCount, int windowSize, int seed) =>
            new ZCurveWindowedProfile(Generator.GetRandomDNA(windowCount * windowSize, seed), windowSize);

        private static string Mutate(string s, double rate, int seed)
        {
            var rng = new Random(seed);
            char[] bases = { 'A', 'C', 'G', 'T' };
            var chars = s.ToCharArray();
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < rate) chars[i] = bases[rng.Next(bases.Length)];
            return new string(chars);
        }

        // Covers narrow (m < 8) and wide (m >= 8) shapes, m % 8 != 0 window tails, and offset counts
        // (n - m + 1) that are 1, below 8, exactly 8/16, and not a multiple of 8.
        [TestCase(1, 1)] [TestCase(1, 2)] [TestCase(1, 8)] [TestCase(1, 9)] [TestCase(1, 40)]
        [TestCase(3, 3)] [TestCase(3, 10)] [TestCase(3, 11)] [TestCase(3, 18)] [TestCase(3, 50)]
        [TestCase(7, 7)] [TestCase(7, 14)] [TestCase(7, 15)] [TestCase(7, 26)]
        [TestCase(8, 8)] [TestCase(8, 9)] [TestCase(8, 15)] [TestCase(8, 16)] [TestCase(8, 29)]
        [TestCase(9, 9)] [TestCase(9, 17)] [TestCase(9, 32)]
        [TestCase(15, 15)] [TestCase(15, 50)]
        [TestCase(16, 16)] [TestCase(16, 23)] [TestCase(16, 53)]
        [TestCase(17, 17)] [TestCase(17, 56)]
        [TestCase(64, 64)] [TestCase(64, 197)]
        [TestCase(100, 100)] [TestCase(100, 305)]
        public void MatchesScalarOracle(int queryWindows, int referenceWindows)
        {
            const int windowSize = 20;
            int seed = BenchmarkTiming.Seed(nameof(MatchesScalarOracle), queryWindows, referenceWindows);
            var query = RandomProfile(queryWindows, windowSize, seed);
            var reference = RandomProfile(referenceWindows, windowSize, seed + 1);

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, reference);
            var oracleScores = OracleScores(query, reference);
            var (oracleScore, _) = OracleBestAlignment(query, reference);

            Assert.That(score, Is.EqualTo(oracleScore).Within(tolerance));
            // The kernel accumulates in float, so an offset within tolerance of the oracle's best is
            // an acceptable winner; anything clearly worse is a bug.
            Assert.That(offset, Is.InRange(0, oracleScores.Length - 1));
            Assert.That(oracleScores[offset], Is.GreaterThanOrEqualTo(oracleScore - tolerance));
        }

        [TestCase(8, 37)]
        [TestCase(13, 5)]
        [TestCase(16, 0)]
        [TestCase(40, 111)]
        [TestCase(3, 197)]
        public void PlantedExactFragmentIsFoundWithPerfectScore(int fragmentWindows, int plantedAt)
        {
            const int windowSize = 100;
            var reference = RandomProfile(200, windowSize, seed: 314 + fragmentWindows);
            var fragment = new ZCurveWindowedProfile(
                reference.DnaString!.Substring(plantedAt * windowSize, fragmentWindows * windowSize), windowSize);

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(fragment, reference);

            Assert.That(offset, Is.EqualTo(plantedAt));
            Assert.That(score, Is.EqualTo(1.0));
        }

        [TestCase(4)]
        [TestCase(12)]
        public void PlantedMutatedFragmentLocalizesToItsOrigin(int fragmentWindows)
        {
            const int windowSize = 500;
            const int plantedAt = 57;
            var reference = RandomProfile(300, windowSize, seed: 2718);
            var original = reference.DnaString!.Substring(plantedAt * windowSize, fragmentWindows * windowSize);
            var fragment = new ZCurveWindowedProfile(Mutate(original, 0.05, seed: 99), windowSize);

            var (_, offset) = ZCurveWindowedProfile.BestAlignment(fragment, reference);

            Assert.That(offset, Is.EqualTo(plantedAt));
        }

        [Test]
        public void EmptyQueryReturnsNaNAndMinusOne()
        {
            var empty = new ZCurveWindowedProfile("ACGT", windowSize: 500);
            var reference = RandomProfile(10, 500, seed: 1);

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(empty, reference);

            Assert.That(empty.WindowCount, Is.Zero);
            Assert.That(score, Is.NaN);
            Assert.That(offset, Is.EqualTo(-1));
        }

        [Test]
        public void QueryLongerThanReferenceReturnsNaNAndMinusOne()
        {
            var query = RandomProfile(11, 500, seed: 2);
            var reference = RandomProfile(10, 500, seed: 3);

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, reference);

            Assert.That(score, Is.NaN);
            Assert.That(offset, Is.EqualTo(-1));
        }

        [Test]
        public void MismatchedWindowSizesThrow()
        {
            var dna = Generator.GetRandomDNA(5000, seed: 4);
            var a = new ZCurveWindowedProfile(dna, windowSize: 500);
            var b = new ZCurveWindowedProfile(dna, windowSize: 250);

            Assert.Throws<ArgumentException>(() => ZCurveWindowedProfile.BestAlignment(a, b));
        }

        [TestCase(4)]
        [TestCase(20)]
        public void EqualLengthHasSingleOffsetZero(int windows)
        {
            var a = RandomProfile(windows, 500, seed: 5);
            var b = RandomProfile(windows, 500, seed: 6);

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(a, b);

            Assert.That(offset, Is.EqualTo(0));
            Assert.That(score, Is.EqualTo(OracleBestAlignment(a, b).Score).Within(tolerance));
        }

        [TestCase(3, 40)]
        [TestCase(20, 57)]
        public void CompareIsSymmetric(int aWindows, int bWindows)
        {
            var a = RandomProfile(aWindows, 500, seed: 7);
            var b = RandomProfile(bWindows, 500, seed: 8);

            Assert.That(ZCurveWindowedProfile.Compare(a, b), Is.EqualTo(ZCurveWindowedProfile.Compare(b, a)));
        }

        [Test]
        public void CompareWithEmptyProfileIsZero()
        {
            var empty = new ZCurveWindowedProfile("", windowSize: 500);
            var other = RandomProfile(10, 500, seed: 9);

            Assert.That(ZCurveWindowedProfile.Compare(empty, other), Is.EqualTo(0));
            Assert.That(ZCurveWindowedProfile.Compare(other, empty), Is.EqualTo(0));
        }

        // Reference = [noise block of `period` windows] + the same `period`-window pattern repeated 3
        // times, and the query is a prefix of that pattern, so offsets period, 2*period, 3*period all
        // score exactly 1.0. The lowest (period) must win. period 3 puts two perfect offsets in the same
        // 8-lane narrow batch; period 9 splits them across batches; queries >= 8 take the wide path.
        [TestCase(3, 3)]
        [TestCase(2, 5)]
        [TestCase(3, 9)]
        [TestCase(8, 8)]
        [TestCase(10, 12)]
        [TestCase(9, 9)]
        public void TiesResolveToLowestOffset(int queryWindows, int period)
        {
            var rng = new Random(BenchmarkTiming.Seed(nameof(TiesResolveToLowestOffset), queryWindows, period));
            (double, double, double) RandomWindow() =>
                (rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1);

            var noise = Enumerable.Range(0, period).Select(_ => RandomWindow()).ToList();
            var pattern = Enumerable.Range(0, period).Select(_ => RandomWindow()).ToList();
            var referenceWindows = noise.Concat(pattern).Concat(pattern).Concat(pattern).ToList();
            var reference = ZCurveWindowedProfile.FromPrecomputedWindows(500, referenceWindows);
            var query = ZCurveWindowedProfile.FromPrecomputedWindows(500, pattern.Take(queryWindows).ToList());

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, reference);

            Assert.That(score, Is.EqualTo(1.0));
            Assert.That(offset, Is.EqualTo(period));
        }

        // Bit-exact snapshot of the current kernel's output. Any numeric change to BestAlignment -
        // summation order, sqrt approximation, precision - fails here, and the failure message shows
        // how far the score drifted. Regenerate deliberately (and say so in the commit) if a change is
        // accepted.
        public static IEnumerable<TestCaseData> SnapshotCases()
        {
            // query bases, reference bases, window size, query seed, reference seed, planted mutated?, expected score bits, expected offset
            yield return new TestCaseData(2_000, 200_000, 500, 11, 12, false, 4606801935499001856L, 264).SetName("Snapshot_Narrow_Unrelated");
            yield return new TestCaseData(2_000, 200_000, 500, 13, 14, true, 4606962744258723840L, 133).SetName("Snapshot_Narrow_PlantedMutated");
            yield return new TestCaseData(6_500, 100_000, 500, 15, 16, false, 4606590295707121191L, 187).SetName("Snapshot_WideTail_Unrelated");
            yield return new TestCaseData(20_000, 1_000_000, 500, 17, 18, true, 4606898128015903949L, 666).SetName("Snapshot_Wide_PlantedMutated");
            yield return new TestCaseData(250_000, 750_000, 500, 19, 20, false, 4606330636836794270L, 86).SetName("Snapshot_WideLarge_Unrelated");
            yield return new TestCaseData(50_000, 50_000, 500, 21, 22, false, 4606325446876562719L, 0).SetName("Snapshot_EqualLength");
            yield return new TestCaseData(3_300, 99_900, 100, 23, 24, false, 4605536778452134074L, 666).SetName("Snapshot_Window100_Unrelated");
        }

        private static (ZCurveWindowedProfile Query, ZCurveWindowedProfile Reference) SnapshotProfiles(
            int queryBases, int referenceBases, int windowSize, int querySeed, int referenceSeed, bool plantedMutated)
        {
            var referenceDna = Generator.GetRandomDNA(referenceBases, referenceSeed);
            // Planted fragments start on a window boundary about a third of the way in.
            var queryDna = plantedMutated
                ? Mutate(referenceDna.Substring(referenceBases / 3 / windowSize * windowSize, queryBases), 0.1, querySeed)
                : Generator.GetRandomDNA(queryBases, querySeed);
            return (new ZCurveWindowedProfile(queryDna, windowSize), new ZCurveWindowedProfile(referenceDna, windowSize));
        }

        [TestCaseSource(nameof(SnapshotCases))]
        public void MatchesGoldenSnapshot(int queryBases, int referenceBases, int windowSize, int querySeed, int referenceSeed,
            bool plantedMutated, long expectedScoreBits, int expectedOffset)
        {
            var (query, reference) = SnapshotProfiles(queryBases, referenceBases, windowSize, querySeed, referenceSeed, plantedMutated);

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, reference);

            double expectedScore = BitConverter.Int64BitsToDouble(expectedScoreBits);
            Assert.Multiple(() =>
            {
                Assert.That(offset, Is.EqualTo(expectedOffset));
                Assert.That(BitConverter.DoubleToInt64Bits(score), Is.EqualTo(expectedScoreBits),
                    $"score {score:R} vs snapshot {expectedScore:R} (drift {score - expectedScore:E3})");
            });
        }
    }
}
