using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class ZCurveProfileTests
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

        // Samples from a weighted base pool to build a sequence with a skewed composition.
        private static string BiasedDNA(char[] weightedBases, int length, int seed)
        {
            var rng = new Random(seed);
            return new string(Enumerable.Range(0, length).Select(_ => weightedBases[rng.Next(weightedBases.Length)]).ToArray());
        }

        [Test]
        public void IdenticalSequencesCompareAsOne()
        {
            var dna = Generator.GetRandomDNA(5000, seed: 42);
            var a = new ZCurveProfile(dna);
            var b = new ZCurveProfile(dna);
            Assert.That(ZCurveProfile.Compare(a, b), Is.EqualTo(1.0).Within(1e-9));
        }

        [Test]
        public void ShuffledSameCompositionSequenceHitsTheCompositionCeiling()
        {
            // Shuffling preserves the exact base multiset, so the Z-curve endpoint - a linear
            // recombination of base counts - must be identical even though every base moved.
            // This is the direct, exact demonstration of the ceiling documented on the class.
            var dna = Generator.GetRandomDNA(5000, seed: 7);
            var shuffled = Shuffle(dna, seed: 99);

            var original = new ZCurveProfile(dna);
            var reordered = new ZCurveProfile(shuffled);

            TestContext.Out.WriteLine($"original point: {original.Point}");
            TestContext.Out.WriteLine($"shuffled point: {reordered.Point}");
            TestContext.Out.WriteLine($"Compare: {ZCurveProfile.Compare(original, reordered)}");

            Assert.That(ZCurveProfile.Compare(original, reordered), Is.EqualTo(1.0).Within(1e-9),
                "endpoint-only Z-curve cannot see re-ordering that preserves composition");
        }

        [Test]
        public void UnrelatedRandomSequencesOfSameLengthScoreNearOne()
        {
            // Two independently generated random sequences converge to ~25% each base, so their
            // Z-curve endpoints land close together too - the same false-positive shape reported
            // for the original CAT method on real unrelated same-composition sequences (phiX174 vs
            // lambda phage, ~0.99) is expected here for the same underlying reason.
            var x = new ZCurveProfile(Generator.GetRandomDNA(50_000, seed: 1));
            var y = new ZCurveProfile(Generator.GetRandomDNA(50_000, seed: 2));

            var score = ZCurveProfile.Compare(x, y);
            TestContext.Out.WriteLine($"unrelated random (length 50000) Compare: {score}");

            Assert.That(score, Is.GreaterThan(0.98));
        }

        [Test]
        public void NearDuplicateBarelyBeatsUnrelatedAtTheSameLength()
        {
            // 97%-similar should score higher than something truly unrelated. Measured: at length
            // 1000 the near-duplicate scores ~0.988 vs ~0.920 for unrelated-same-length - a
            // comfortable ~0.068 margin here. But per UnrelatedRandomSequencesOfSameLengthScoreNearOne,
            // two *unrelated* random sequences score higher the longer they are (0.977 at length
            // 5000, 0.983 at length 50000) - unrelated pairs drift toward 1 as length grows, since
            // the endpoint is a length-normalized sum concentrating by the law of large numbers.
            // That means this margin should be expected to shrink at longer lengths, the same
            // large-N collapse the original CAT method's single point showed.
            var dna = Generator.GetRandomDNA(1000, seed: 7);
            var nearDuplicate = Mutate(dna, rate: 0.03, seed: 11); // 97% identity
            var unrelated = Generator.GetRandomDNA(1000, seed: 8);

            var x = new ZCurveProfile(dna);
            var duplicateScore = ZCurveProfile.Compare(x, new ZCurveProfile(nearDuplicate));
            var unrelatedScore = ZCurveProfile.Compare(x, new ZCurveProfile(unrelated));

            TestContext.Out.WriteLine($"97%-similar Compare: {duplicateScore}");
            TestContext.Out.WriteLine($"unrelated (same length) Compare: {unrelatedScore}");
            TestContext.Out.WriteLine($"margin: {duplicateScore - unrelatedScore}");

            Assert.That(duplicateScore, Is.GreaterThan(unrelatedScore),
                "a near-duplicate should still edge out something unrelated, even if the margin is thin");
        }

        [Test]
        public void DifferentCompositionSequencesAreMeasurablyFartherApart()
        {
            // Sanity check that the axes carry real signal at all: force a GC-rich sequence
            // (skewed base frequencies) against a uniform-random one of matched length.
            char[] gcRich = { 'G', 'C', 'G', 'C', 'A', 'T' }; // ~67% GC
            var gcRichSeq = BiasedDNA(gcRich, 5000, seed: 3);
            var uniform = Generator.GetRandomDNA(5000, seed: 4);

            var x = new ZCurveProfile(gcRichSeq);
            var y = new ZCurveProfile(uniform);
            var z = new ZCurveProfile(Generator.GetRandomDNA(5000, seed: 1));

            var skewedVsUniform = ZCurveProfile.Compare(x, y);
            var uniformVsUniform = ZCurveProfile.Compare(y, z);

            TestContext.Out.WriteLine($"GC-rich vs uniform Compare: {skewedVsUniform}");
            TestContext.Out.WriteLine($"uniform vs uniform Compare: {uniformVsUniform}");

            Assert.That(skewedVsUniform, Is.LessThan(uniformVsUniform),
                "a real composition difference should separate the two endpoints further than two same-composition sequences do");
        }

        [Test]
        public void TrajectoryDetectsRegionalStructureThatTheEndpointCannotSee()
        {
            // Build a sequence with a real regional trend: an AT-rich first half followed by a
            // GC-rich second half (the kind of isochore/GC-skew structure real chromosomes
            // actually have). Shuffling it destroys that regional trend while leaving the overall
            // composition - and therefore the endpoint - completely unchanged, since shuffling is
            // just a permutation of the same multiset of bases.
            char[] atRich = { 'A', 'A', 'T', 'T', 'A', 'G' }; // ~67% AT
            char[] gcRich = { 'G', 'C', 'G', 'C', 'A', 'T' }; // ~67% GC
            var structured = BiasedDNA(atRich, 2500, seed: 10) + BiasedDNA(gcRich, 2500, seed: 11);
            var shuffled = Shuffle(structured, seed: 20);

            const int checkpointCount = 20;
            var original = new ZCurveProfile(structured, checkpointCount);
            var reordered = new ZCurveProfile(shuffled, checkpointCount);

            var endpointScore = ZCurveProfile.Compare(original, reordered);
            var trajectoryScore = ZCurveProfile.CompareTrajectory(original, reordered);

            TestContext.Out.WriteLine($"endpoint Compare (blind to order): {endpointScore}");
            TestContext.Out.WriteLine($"trajectory CompareTrajectory ({checkpointCount} checkpoints): {trajectoryScore}");

            Assert.That(endpointScore, Is.EqualTo(1.0).Within(1e-9),
                "endpoint still can't see the shuffle - same ceiling as before");
            Assert.That(trajectoryScore, Is.LessThan(endpointScore),
                "the trajectory should notice the regional AT->GC trend was destroyed, unlike the endpoint");
        }

        [Test]
        public void TrajectoryWithFixedCheckpointCountStillCollapsesAtLargeEnoughLength()
        {
            // A fixed *count* of checkpoints (unlike WindowedProfile's fixed window *size*) means
            // each checkpoint's own span - and therefore its own noise floor - still grows with
            // total sequence length. So the same large-N collapse the endpoint has should still
            // show up here too, just at a length pushed out roughly checkpointCount-fold further.
            const int checkpointCount = 20;
            var shortScore = ZCurveProfile.CompareTrajectory(
                new ZCurveProfile(Generator.GetRandomDNA(2000, seed: 1), checkpointCount),
                new ZCurveProfile(Generator.GetRandomDNA(2000, seed: 2), checkpointCount));
            var longScore = ZCurveProfile.CompareTrajectory(
                new ZCurveProfile(Generator.GetRandomDNA(200_000, seed: 1), checkpointCount),
                new ZCurveProfile(Generator.GetRandomDNA(200_000, seed: 2), checkpointCount));

            TestContext.Out.WriteLine($"unrelated CompareTrajectory at length 2000: {shortScore}");
            TestContext.Out.WriteLine($"unrelated CompareTrajectory at length 200000: {longScore}");

            Assert.That(longScore, Is.GreaterThan(shortScore),
                "fixed checkpoint count delays but does not remove the large-N collapse - window size, not window count, would need to stay fixed to avoid it");
        }

        [Test]
        public void TrajectoryMarginForNearDuplicateVsUnrelated()
        {
            // Same experiment as NearDuplicateBarelyBeatsUnrelatedAtTheSameLength, but scored with
            // CompareTrajectory (20 checkpoints) instead of the endpoint-only Compare, to see
            // whether positional resolution helps or hurts the margin on sequences with no real
            // regional structure of their own (plain i.i.d. random DNA).
            const int checkpointCount = 20;
            var dna = Generator.GetRandomDNA(1000, seed: 7);
            var nearDuplicate = Mutate(dna, rate: 0.03, seed: 11); // 97% identity
            var unrelated = Generator.GetRandomDNA(1000, seed: 8);

            var x = new ZCurveProfile(dna, checkpointCount);
            var duplicateScore = ZCurveProfile.CompareTrajectory(x, new ZCurveProfile(nearDuplicate, checkpointCount));
            var unrelatedScore = ZCurveProfile.CompareTrajectory(x, new ZCurveProfile(unrelated, checkpointCount));

            TestContext.Out.WriteLine($"97%-similar CompareTrajectory: {duplicateScore}");
            TestContext.Out.WriteLine($"unrelated (same length) CompareTrajectory: {unrelatedScore}");
            TestContext.Out.WriteLine($"margin: {duplicateScore - unrelatedScore}");

            Assert.That(duplicateScore, Is.GreaterThan(unrelatedScore),
                "a near-duplicate should still edge out something unrelated under the trajectory metric too");
        }
    }
}
