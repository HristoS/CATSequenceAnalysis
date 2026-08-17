using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class ZCurveWindowedProfileTests
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

        private static string BiasedDNA(char[] weightedBases, int length, int seed)
        {
            var rng = new Random(seed);
            return new string(Enumerable.Range(0, length).Select(_ => weightedBases[rng.Next(weightedBases.Length)]).ToArray());
        }

        [Test]
        public void IdenticalSequencesCompareAsOne()
        {
            var dna = Generator.GetRandomDNA(5000, seed: 42);
            var a = new ZCurveWindowedProfile(dna);
            var b = new ZCurveWindowedProfile(dna);
            Assert.That(ZCurveWindowedProfile.Compare(a, b), Is.EqualTo(1.0).Within(1e-9));
        }

        [Test]
        public void ExactFragmentLocalizesToCorrectOffsetWithPerfectScore()
        {
            var reference = new ZCurveWindowedProfile(Generator.GetRandomDNA(5000, seed: 7), windowSize: 500);
            // A fragment cut exactly on window boundaries [1500, 3500) == reference windows 3..6.
            // DnaString is non-null here: this profile came from the string constructor, not FromPrecomputedWindows.
            var fragment = new ZCurveWindowedProfile(reference.DnaString!.Substring(1500, 2000), windowSize: 500);

            var (score, offset) = ZCurveWindowedProfile.BestAlignment(fragment, reference);
            Assert.That(offset, Is.EqualTo(3), "fragment starting at base 1500 should localize to window 3");
            Assert.That(score, Is.EqualTo(1.0).Within(1e-9), "an exact fragment's windows match perfectly");
        }

        [Test]
        public void FixedWindowSizeAvoidsTheLargeNCollapse()
        {
            // The whole point of switching from fixed checkpoint *count* to fixed window *size*:
            // each window's own noise floor stays constant regardless of total sequence length, so
            // unrelated-sequence scores should stay roughly flat as length grows, unlike
            // ZCurveProfile.CompareTrajectory (which drifted from 0.91 to 0.99 over the same
            // length range in ZCurveProfileTests).
            var shortScore = ZCurveWindowedProfile.Compare(
                new ZCurveWindowedProfile(Generator.GetRandomDNA(2000, seed: 1)),
                new ZCurveWindowedProfile(Generator.GetRandomDNA(2000, seed: 2)));
            var longScore = ZCurveWindowedProfile.Compare(
                new ZCurveWindowedProfile(Generator.GetRandomDNA(200_000, seed: 1)),
                new ZCurveWindowedProfile(Generator.GetRandomDNA(200_000, seed: 2)));

            TestContext.Out.WriteLine($"unrelated Compare at length 2000: {shortScore}");
            TestContext.Out.WriteLine($"unrelated Compare at length 200000: {longScore}");
            TestContext.Out.WriteLine($"drift: {longScore - shortScore}");

            Assert.That(Math.Abs(longScore - shortScore), Is.LessThan(0.05),
                "unlike the fixed-checkpoint-count trajectory, a fixed window size should not let unrelated scores drift toward 1 as length grows");
        }

        [Test]
        public void TrajectoryDetectsRegionalStructureThatTheEndpointCannotSee()
        {
            char[] atRich = { 'A', 'A', 'T', 'T', 'A', 'G' }; // ~67% AT
            char[] gcRich = { 'G', 'C', 'G', 'C', 'A', 'T' }; // ~67% GC
            var structured = BiasedDNA(atRich, 2500, seed: 10) + BiasedDNA(gcRich, 2500, seed: 11);
            var shuffled = Shuffle(structured, seed: 20);

            var original = new ZCurveWindowedProfile(structured, windowSize: 500);
            var reordered = new ZCurveWindowedProfile(shuffled, windowSize: 500);

            var score = ZCurveWindowedProfile.Compare(original, reordered);
            TestContext.Out.WriteLine($"structured vs shuffled Compare: {score}");

            Assert.That(score, Is.LessThan(0.99),
                "windowed Z-curve should also notice the regional AT->GC trend was destroyed by shuffling");
        }

        [Test]
        public void NearDuplicateBeatsUnrelatedAtTheSameLength()
        {
            var dna = Generator.GetRandomDNA(5000, seed: 7);
            var nearDuplicate = Mutate(dna, rate: 0.03, seed: 11); // 97% identity
            var unrelated = Generator.GetRandomDNA(5000, seed: 8);

            var x = new ZCurveWindowedProfile(dna);
            var duplicateScore = ZCurveWindowedProfile.Compare(x, new ZCurveWindowedProfile(nearDuplicate));
            var unrelatedScore = ZCurveWindowedProfile.Compare(x, new ZCurveWindowedProfile(unrelated));

            TestContext.Out.WriteLine($"97%-similar Compare: {duplicateScore}");
            TestContext.Out.WriteLine($"unrelated (same length) Compare: {unrelatedScore}");
            TestContext.Out.WriteLine($"margin: {duplicateScore - unrelatedScore}");

            Assert.That(duplicateScore, Is.GreaterThan(unrelatedScore));
        }

        [Test]
        public void RehydratedFromPrecomputedWindowsMatchesFreshlyBuiltProfile()
        {
            // Answers a reviewer question about the article: is profile construction genuinely a
            // separate, cacheable step, or just conceptually separable? This proves it end-to-end -
            // a profile rebuilt only from stored window vectors (no original sequence) scores
            // identically to the profile it was cached from.
            var dna = Generator.GetRandomDNA(5000, seed: 123);
            var fresh = new ZCurveWindowedProfile(dna, windowSize: 500);

            var rehydrated = ZCurveWindowedProfile.FromPrecomputedWindows(fresh.WindowSize, fresh.Windows);

            Assert.That(rehydrated.DnaString, Is.Null, "rehydrated-without-a-string profile should not fabricate one");
            Assert.That(rehydrated.WindowSize, Is.EqualTo(fresh.WindowSize));
            Assert.That(rehydrated.Windows, Is.EqualTo(fresh.Windows), "window vectors must round-trip tuple-for-tuple");

            var other = new ZCurveWindowedProfile(Generator.GetRandomDNA(5000, seed: 456), windowSize: 500);

            // Exact equality (no tolerance): Compare/BestAlignment only ever read WindowSize and
            // Windows, so this is the same arithmetic over the same inputs, not two independently
            // constructed profiles that merely ought to agree approximately.
            Assert.That(ZCurveWindowedProfile.Compare(rehydrated, other),
                Is.EqualTo(ZCurveWindowedProfile.Compare(fresh, other)));

            var (freshScore, freshOffset) = ZCurveWindowedProfile.BestAlignment(fresh, other);
            var (rehydratedScore, rehydratedOffset) = ZCurveWindowedProfile.BestAlignment(rehydrated, other);
            Assert.That(rehydratedScore, Is.EqualTo(freshScore));
            Assert.That(rehydratedOffset, Is.EqualTo(freshOffset));
        }

        [Test]
        public void FromPrecomputedWindowsValidatesLengthAgainstSuppliedDnaString()
        {
            var dna = Generator.GetRandomDNA(5000, seed: 1);
            var fresh = new ZCurveWindowedProfile(dna, windowSize: 500);

            Assert.That(
                () => ZCurveWindowedProfile.FromPrecomputedWindows(fresh.WindowSize, fresh.Windows, dnaString: dna.Substring(0, 4999)),
                Throws.ArgumentException);

            Assert.That(
                () => ZCurveWindowedProfile.FromPrecomputedWindows(fresh.WindowSize, fresh.Windows, dnaString: dna),
                Throws.Nothing);
        }
    }
}
