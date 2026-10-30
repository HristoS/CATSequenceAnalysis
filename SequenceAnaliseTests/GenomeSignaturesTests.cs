using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class GenomeSignaturesTests
    {
        // 200 kbp in four blocks of different composition (GC-rich, AT-rich, purine-rich, amino-rich), so the cumulative curve
        // bends and its principal directions are well defined; a multiple of 500, so the profile drops no partial window.
        private static string Structured(int seed)
        {
            var rng = new Random(seed);
            var weights = new[] { new[] { 1, 3, 3, 1 }, new[] { 3, 1, 1, 3 }, new[] { 3, 1, 3, 1 }, new[] { 3, 3, 1, 1 } };   // A, C, G, T
            var chars = new char[200_000];
            for (int i = 0; i < chars.Length; i++)
            {
                var w = weights[i / 50_000];
                int r = rng.Next(w.Sum()), k = 0;
                while (r >= w[k]) r -= w[k++];
                chars[i] = "ACGT"[k];
            }
            return new string(chars);
        }

        private static string ReverseComplement(string sequence) =>
            new string(sequence.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', var o => o }).ToArray());

        [Test]
        public void NaturalVectorMatchesTheDefinitionOnATinySequence()
        {
            // Each base once at positions 1..4: counts 1, means 1..4, no spread. Then "AAT": A at 1 and 2 (mean 1.5,
            // D = (0.25 + 0.25) / (2 * 3)), T at 3.
            Assert.That(GenomeSignatures.NaturalVector("ACGT"), Is.EqualTo(new double[] { 1, 1, 1, 1, 1, 2, 3, 4, 0, 0, 0, 0 }));
            var aat = GenomeSignatures.NaturalVector("AAT");
            Assert.That(aat[0], Is.EqualTo(2));
            Assert.That(aat[4], Is.EqualTo(1.5));
            Assert.That(aat[8], Is.EqualTo(0.5 / 6).Within(1e-12));
            Assert.That(aat[7], Is.EqualTo(3));
        }

        [Test]
        public void KarlinRhoIsTheSameOnBothStrands()
        {
            var sequence = Structured(1);
            var forward = GenomeSignatures.KarlinRho(GenomeSignatures.DinucleotideCounts(sequence));
            var reverse = GenomeSignatures.KarlinRho(GenomeSignatures.DinucleotideCounts(ReverseComplement(sequence)));
            Assert.That(forward, Is.EqualTo(reverse).Within(1e-12));
            Assert.That(GenomeSignatures.KarlinDelta(forward, reverse), Is.LessThan(1e-12));
            // Independent draws within each block: every rho* near 1 apart from the block mixture's own correlation.
            Assert.That(forward.All(r => r > 0.5 && r < 1.5));
        }

        [Test]
        public void SecondOrderZCurveIsStrandInvariantAndCatchesNeighbourStructure()
        {
            var sequence = Structured(3);
            var forward = GenomeSignatures.SecondOrder(GenomeSignatures.DinucleotideCounts(sequence));
            var reverse = GenomeSignatures.SecondOrder(GenomeSignatures.DinucleotideCounts(ReverseComplement(sequence)));
            Assert.That(forward, Is.EqualTo(reverse).Within(1e-12));

            // Independent bases: no correlation between neighbours, whatever the composition.
            var random = GenomeSignatures.SecondOrder(GenomeSignatures.DinucleotideCounts(Generator.GetRandomDNA(400_000, seed: 5)));
            Assert.That(random.All(r => Math.Abs(r) < 0.01), string.Join(", ", random));

            // ATAT...: x (purine +1 / pyrimidine -1) alternates, so corr(x(n), x(n + 1)) = -1; y alternates too (A amino, T keto).
            var alternating = GenomeSignatures.SecondOrder(GenomeSignatures.DinucleotideCounts(string.Concat(Enumerable.Repeat("AT", 5000))));
            Assert.That(alternating[0], Is.EqualTo(-1).Within(1e-3));
            Assert.That(alternating[4], Is.EqualTo(-1).Within(1e-3));
        }

        [Test]
        public void SlidingDistanceFindsAnEmbeddedSeriesOnEitherStrand()
        {
            var rng = new Random(9);
            var longer = Enumerable.Range(0, 30).Select(_ => Enumerable.Range(0, 9).Select(_ => rng.NextDouble()).ToArray()).ToArray();
            var inside = longer.Skip(7).Take(10).ToArray();
            Assert.That(GenomeSignatures.SlidingDistance(inside, longer), Is.EqualTo(0));
            Assert.That(GenomeSignatures.SlidingDistance(longer, inside.Reverse().ToArray()), Is.EqualTo(0));
            Assert.That(GenomeSignatures.SlidingDistance(longer, longer.Select(w => w.Select(v => v + 0.5).ToArray()).ToArray()),
                Is.EqualTo(Math.Sqrt(9 * 0.25)).Within(1e-12));
            Assert.That(GenomeSignatures.SlidingDistance(Array.Empty<double[]>(), longer), Is.EqualTo(double.PositiveInfinity));
        }

        [Test]
        public void ProfileDerivedSignaturesMatchTheExactOnes()
        {
            var sequence = Structured(2);
            var fine = new ZCurveWindowedProfile(sequence);

            var exactVector = GenomeSignatures.NaturalVector(sequence);
            var derivedVector = GenomeSignatures.NaturalVector(fine, sequence.Length);
            for (int i = 0; i < 12; i++)
                Assert.That(derivedVector[i], Is.EqualTo(exactVector[i]).Within(Math.Abs(exactVector[i]) * 1e-3 + 1e-6), $"Natural Vector [{i}]");

            Assert.That(GenomeSignatures.ZCurveMean(fine), Is.EqualTo(GenomeSignatures.ZCurveMean(sequence)).Within(1e-6));

            var exactCurve = GenomeSignatures.ZhengCurve(sequence);
            var derivedCurve = GenomeSignatures.ZhengCurve(fine);
            for (int i = 0; i < exactCurve.Length; i++)
                Assert.That(derivedCurve[i], Is.EqualTo(exactCurve[i]).Within(0.01), $"Zheng curve [{i}]");
        }
    }
}
