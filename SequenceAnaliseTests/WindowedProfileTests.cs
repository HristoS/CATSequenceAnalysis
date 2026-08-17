using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class WindowedProfileTests
    {
        // Mutate a fraction of bases deterministically, to simulate a diverged homolog.
        private static string Mutate(string s, double rate, int seed)
        {
            var rng = new Random(seed);
            char[] bases = { 'A', 'C', 'G', 'T' };
            var chars = s.ToCharArray();
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < rate) chars[i] = bases[rng.Next(bases.Length)];
            return new string(chars);
        }

        [Test]
        public void WindowCountMatchesSequenceLength()
        {
            var profile = new WindowedProfile(Generator.GetRandomDNA(5000, seed: 1), windowSize: 500, k: 4);
            Assert.That(profile.Windows.Count, Is.EqualTo(10));
        }

        [Test]
        public void ShorterThanOneWindowYieldsNoWindows()
        {
            var profile = new WindowedProfile(Generator.GetRandomDNA(100, seed: 1), windowSize: 500, k: 4);
            Assert.That(profile.Windows.Count, Is.EqualTo(0));
            // Comparing against anything returns 0 rather than throwing.
            var other = new WindowedProfile(Generator.GetRandomDNA(5000, seed: 2), windowSize: 500, k: 4);
            Assert.That(WindowedProfile.Compare(profile, other), Is.EqualTo(0));
        }

        [Test]
        public void IdenticalSequencesCompareAsOne()
        {
            var dna = Generator.GetRandomDNA(5000, seed: 42);
            var a = new WindowedProfile(dna);
            var b = new WindowedProfile(dna);
            Assert.That(WindowedProfile.Compare(a, b), Is.EqualTo(1.0).Within(1e-9));
        }

        [Test]
        public void ExactFragmentLocalizesToCorrectOffsetWithPerfectScore()
        {
            var reference = new WindowedProfile(Generator.GetRandomDNA(5000, seed: 7), windowSize: 500, k: 4);
            // A fragment cut exactly on window boundaries [1500, 3500) == reference windows 3..6.
            var fragment = new WindowedProfile(reference.DnaString.Substring(1500, 2000), windowSize: 500, k: 4);

            var (score, offset) = WindowedProfile.BestAlignment(fragment, reference);
            Assert.That(offset, Is.EqualTo(3), "fragment starting at base 1500 should localize to window 3");
            Assert.That(score, Is.EqualTo(1.0).Within(1e-9), "an exact fragment's windows match perfectly");
        }

        [Test]
        public void TrueFragmentScoresHigherThanUnrelatedSequence()
        {
            var reference = new WindowedProfile(Generator.GetRandomDNA(5000, seed: 7), windowSize: 500, k: 4);
            var trueFragment = new WindowedProfile(reference.DnaString.Substring(1500, 2000));
            var unrelated = new WindowedProfile(Generator.GetRandomDNA(2000, seed: 999));

            double fragScore = WindowedProfile.Compare(reference, trueFragment);
            double unrelatedScore = WindowedProfile.Compare(reference, unrelated);

            Assert.That(fragScore, Is.GreaterThan(unrelatedScore),
                $"true fragment ({fragScore:F4}) should beat unrelated ({unrelatedScore:F4})");
            Assert.That(fragScore, Is.GreaterThan(0.99), "an exact fragment should score near 1");
        }

        [Test]
        public void DivergedFragmentStillLocalizesButScoresLower()
        {
            var reference = new WindowedProfile(Generator.GetRandomDNA(5000, seed: 7), windowSize: 500, k: 4);
            var exact = reference.DnaString.Substring(1500, 2000);
            var diverged = new WindowedProfile(Mutate(exact, rate: 0.15, seed: 3)); // ~85% identity

            var (score, offset) = WindowedProfile.BestAlignment(diverged, reference);
            Assert.That(offset, Is.EqualTo(3), "a diverged fragment should still localize to its true region");
            Assert.That(score, Is.InRange(0.80, 0.98), "a diverged fragment scores below an exact one but well above unrelated");
        }

        [Test]
        public void AmbiguityCharactersDoNotThrow()
        {
            // Interspersed 'N' runs should be skipped, not crash the k-mer roller.
            var dna = Generator.GetRandomDNA(1000, seed: 5) + new string('N', 50) + Generator.GetRandomDNA(1000, seed: 6);
            var profile = new WindowedProfile(dna);
            Assert.That(profile.Windows.Count, Is.GreaterThan(0));
        }
    }
}
