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
            // Answers a likely reviewer question: is profile construction genuinely a
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

        // 10,000 is an exact multiple of 2,000; the others leave a trailing partial coarse window
        // that both paths must drop identically.
        [TestCase(10_000)]
        [TestCase(10_499)]
        [TestCase(11_999)]
        [TestCase(1_999)]
        public void DownsampleMatchesProfileBuiltAtTheCoarserWindowSize(int length)
        {
            var dna = Generator.GetRandomDNA(length, seed: length);
            var fine = new ZCurveWindowedProfile(dna, windowSize: 500);
            var direct = new ZCurveWindowedProfile(dna, windowSize: 2000);

            var downsampled = fine.Downsample(4);

            Assert.That(downsampled.WindowSize, Is.EqualTo(direct.WindowSize));
            Assert.That(downsampled.WindowCount, Is.EqualTo(direct.WindowCount));
            for (int i = 0; i < direct.WindowCount; i++)
            {
                Assert.That(downsampled.Windows[i].X, Is.EqualTo(direct.Windows[i].X).Within(1e-6));
                Assert.That(downsampled.Windows[i].Y, Is.EqualTo(direct.Windows[i].Y).Within(1e-6));
                Assert.That(downsampled.Windows[i].Z, Is.EqualTo(direct.Windows[i].Z).Within(1e-6));
            }
        }

        [Test]
        public void DownsampleWithPhaseMatchesProfileOfTheShiftedSequence()
        {
            var dna = Generator.GetRandomDNA(20_000, seed: 3);
            var fine = new ZCurveWindowedProfile(dna, windowSize: 500);

            for (int phase = 0; phase < 4; phase++)
            {
                var shifted = new ZCurveWindowedProfile(dna.Substring(phase * 500), windowSize: 2000);
                var downsampled = fine.Downsample(4, phase);

                Assert.That(downsampled.WindowCount, Is.EqualTo(shifted.WindowCount), $"phase {phase}");
                for (int i = 0; i < shifted.WindowCount; i++)
                    Assert.That(downsampled.Windows[i].X, Is.EqualTo(shifted.Windows[i].X).Within(1e-6), $"phase {phase}, window {i}");
            }
        }

        private static ZCurveProfileLadder LadderOf(int length, int seed) =>
            new(new ZCurveWindowedProfile(Generator.GetRandomDNA(length, seed), windowSize: 500));

        // Default ladder: 500 bp / 2 kbp / 8 kbp / 32 kbp, minimum 8 windows for the shorter sequence.
        [TestCase(10_000, 12_000, 1)]      // virus-sized: 20 windows at 500 bp, only 5 at 2 kbp
        [TestCase(40_000, 50_000, 4)]      // 20 windows at 2 kbp, 5 at 8 kbp
        [TestCase(600_000, 900_000, 64)]   // bacterium-sized: 18 windows at 32 kbp
        [TestCase(10_000, 600_000, 1)]     // mixed: the shorter (virus-sized) sequence decides
        [TestCase(2_000, 600_000, 1)]      // fewer than 8 windows anywhere: falls back to the fine level
        public void LadderPicksCoarsestLevelKeepingEnoughWindows(int lengthX, int lengthY, int expectedFactor)
        {
            var x = LadderOf(lengthX, seed: 1);
            var y = LadderOf(lengthY, seed: 2);

            Assert.That(ZCurveProfileLadder.SelectFactor(x, y), Is.EqualTo(expectedFactor));
            Assert.That(ZCurveProfileLadder.Compare(x, y),
                Is.EqualTo(ZCurveWindowedProfile.Compare(x.Level(expectedFactor), y.Level(expectedFactor))));
        }

        [Test]
        public void LadderReturnsNullRatherThanZeroWhenNotComparable()
        {
            var tooShort = LadderOf(499, seed: 1);
            var genome = LadderOf(50_000, seed: 2);

            Assert.That(ZCurveWindowedProfile.Compare(tooShort.Fine, genome.Fine), Is.EqualTo(0),
                "the plain Compare reports 0 - indistinguishable from 'unrelated' - which is why the ladder exists");
            Assert.That(ZCurveProfileLadder.SelectFactor(tooShort, genome), Is.Null);
            Assert.That(ZCurveProfileLadder.Compare(tooShort, genome), Is.Null);
        }

        [TestCase(10_000, 1)]       // 5 windows at 2 kbp
        [TestCase(40_000, 4)]       // 20 at 2 kbp, 5 at 8 kbp
        [TestCase(600_000, 64)]     // 18 at 32 kbp
        [TestCase(2_000, 1)]        // 4 windows: falls back to the fine level
        [TestCase(499, null)]       // no fine window
        public void QueryFactorDependsOnTheQueryAlone(int length, int? expectedFactor)
        {
            Assert.That(ZCurveProfileLadder.SelectQueryFactor(LadderOf(length, seed: 1)), Is.EqualTo(expectedFactor));
        }

        [TestCase(600_000, 0.5, 64)]   // 18 windows at 32 kbp: 9 >= 8
        [TestCase(500_000, 0.5, 16)]   // 15 windows at 32 kbp: 7.5 < 8
        [TestCase(500_000, 1.0, 64)]   // ratio 1 is SelectQueryFactor
        [TestCase(40_000, 0.5, 4)]     // 20 windows at 2 kbp: 10 >= 8
        [TestCase(10_000, 0.5, 1)]     // 20 fine windows: 10 >= 8
        [TestCase(6_000, 0.5, 1)]      // 12 fine windows: 6 < 8 anywhere, falls back to the fine level
        [TestCase(499, 0.5, null)]
        public void BandFactorLeavesTheShortestAllowedCandidateEnoughWindows(int length, double ratio, int? expectedFactor)
        {
            var query = LadderOf(length, seed: 1);
            Assert.That(ZCurveProfileLadder.SelectBandFactor(query, ratio), Is.EqualTo(expectedFactor));

            if (expectedFactor is int f && f > 1)
            {
                var shortest = LadderOf((int)Math.Ceiling(ratio * length), seed: 2);
                Assert.That(shortest.Level(f).WindowCount, Is.GreaterThanOrEqualTo(ZCurveProfileLadder.DefaultMinWindows));
            }
        }

        [Test]
        public void SearchScoresEveryRankedCandidateAtTheQueryLevel()
        {
            var query = LadderOf(600_000, seed: 1);
            var longer = LadderOf(900_000, seed: 2);
            var shorterButEnough = LadderOf(300_000, seed: 3);   // 9 windows at 32 kbp

            var a = ZCurveProfileLadder.CompareForQuery(query, longer);
            var b = ZCurveProfileLadder.CompareForQuery(query, shorterButEnough);

            Assert.That((a.List, a.Factor), Is.EqualTo((SearchList.Ranked, 64)));
            Assert.That((b.List, b.Factor), Is.EqualTo((SearchList.Ranked, 64)));
            Assert.That(b.Score, Is.EqualTo(ZCurveWindowedProfile.Compare(query.Level(64), shorterButEnough.Level(64))));
        }

        [Test]
        public void SearchSendsCandidatesTooShortForTheQueryLevelToTheShortList()
        {
            var query = LadderOf(600_000, seed: 1);
            var virusSized = LadderOf(10_000, seed: 2);

            var match = ZCurveProfileLadder.CompareForQuery(query, virusSized);

            Assert.That((match.List, match.Factor), Is.EqualTo((SearchList.ShortList, 1)));
            Assert.That(match.Score, Is.EqualTo(ZCurveProfileLadder.Compare(query, virusSized)),
                "the short list is scored at the per-pair level");
        }

        [Test]
        public void SearchByAShortQueryRanksCandidatesWithAtLeastTheQuerysWindows()
        {
            var query = LadderOf(2_000, seed: 1);                    // 4 fine windows
            var genome = LadderOf(50_000, seed: 2);
            var shorterThanQuery = LadderOf(1_000, seed: 3);         // 2 fine windows

            Assert.That(ZCurveProfileLadder.CompareForQuery(query, genome).List, Is.EqualTo(SearchList.Ranked));
            Assert.That(ZCurveProfileLadder.CompareForQuery(query, shorterThanQuery).List, Is.EqualTo(SearchList.ShortList));
        }

        [Test]
        public void SearchReturnsNotComparableWithoutAFineWindowOnEitherSide()
        {
            var genome = LadderOf(50_000, seed: 1);
            var tooShort = LadderOf(499, seed: 2);

            Assert.That(ZCurveProfileLadder.CompareForQuery(tooShort, genome), Is.EqualTo(new SearchMatch(SearchList.NotComparable, null, null)));
            Assert.That(ZCurveProfileLadder.CompareForQuery(genome, tooShort), Is.EqualTo(new SearchMatch(SearchList.NotComparable, null, null)));
        }

        // Random DNA with a replication-like GC skew: G-rich from 0 to the middle (the leading strand from an origin at 0),
        // C-rich after it, so the origin is at 0 and the terminus in the middle.
        private static string Skewed(int length, int seed)
        {
            var chars = Generator.GetRandomDNA(length, seed).ToCharArray();
            var rng = new Random(seed + 1);
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < 0.1)
                    chars[i] = i < length / 2 ? (chars[i] == 'C' ? 'G' : chars[i]) : (chars[i] == 'G' ? 'C' : chars[i]);
            return new string(chars);
        }

        private static int CircularDistance(int a, int b, int count) => Math.Min(Math.Abs(a - b), count - Math.Abs(a - b));

        [Test]
        public void GcSkewExtremesFindTheOriginAndTerminusWhereverTheRingWasCut()
        {
            var dna = Skewed(200_000, seed: 1);
            var cut = new ZCurveWindowedProfile(dna[30_000..] + dna[..30_000]);   // 400 windows; the origin now sits at window 340
            var reverse = new ZCurveWindowedProfile(ReverseComplement(dna[30_000..] + dna[..30_000]));

            var skew = cut.GcSkewExtremes();
            var reverseSkew = reverse.GcSkewExtremes();

            Assert.That(CircularDistance(skew.OriginWindow, 340, 400), Is.LessThanOrEqualTo(2));
            Assert.That(CircularDistance(skew.TerminusWindow, 140, 400), Is.LessThanOrEqualTo(2));
            Assert.That(skew.AmplitudeBp, Is.GreaterThan(1_000));
            Assert.That(CircularDistance(reverseSkew.OriginWindow, 400 - 340, 400), Is.LessThanOrEqualTo(2), "the other strand finds the same origin");
        }

        [Test]
        public void RotatingToTheOriginLaysOutCopiesCutElsewhereAlike()
        {
            var dna = Skewed(200_000, seed: 2);
            var original = new ZCurveWindowedProfile(dna);
            var copy = new ZCurveWindowedProfile(dna[100_000..] + dna[..100_000]);

            var rotated = ZCurveWindowedProfile.Compare(original.RotateTo(original.GcSkewExtremes().OriginWindow),
                copy.RotateTo(copy.GcSkewExtremes().OriginWindow));

            Assert.That(rotated, Is.GreaterThan(0.999));
            Assert.That(rotated, Is.GreaterThan(ZCurveWindowedProfile.Compare(original, copy)));
        }

        [Test]
        public void RotateToWrapsTheWindowsRound()
        {
            var profile = ZCurveWindowedProfile.FromPrecomputedWindows(500, new[] { (0.1, 0.0, 0.0), (0.2, 0.0, 0.0), (0.3, 0.0, 0.0), (0.4, 0.0, 0.0) });

            Assert.That(profile.RotateTo(1).Windows.Select(w => Math.Round(w.X, 3)), Is.EqualTo(new[] { 0.2, 0.3, 0.4, 0.1 }));
            Assert.That(profile.RotateTo(0).Windows, Is.EqualTo(profile.Windows));
            Assert.That(() => profile.RotateTo(4), Throws.InstanceOf<ArgumentOutOfRangeException>());
        }

        private static string ReverseComplement(string dna) =>
            new(dna.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', var o => o }).ToArray());

        // Localization: the coarsest level leaving the fragment 4 windows, i.e. 500 bp below 8 kbp, 2 kbp below 32 kbp, ...
        [TestCase(7_999, 1)]
        [TestCase(8_000, 4)]
        [TestCase(31_999, 4)]
        [TestCase(32_000, 16)]
        [TestCase(127_999, 16)]
        [TestCase(128_000, 64)]
        public void FragmentFactorKeepsFourCoarseWindows(int length, int expectedFactor)
        {
            var fragment = LadderOf(length, seed: 1);
            Assert.That(ZCurveProfileLadder.SelectQueryFactor(fragment, ZCurveProfileLadder.DefaultMinFragmentWindows), Is.EqualTo(expectedFactor));
        }

        [Test]
        public void RefineAlignmentMatchesTheFullSearchNearItsOffset()
        {
            var dna = Generator.GetRandomDNA(200_000, seed: 1);
            var reference = new ZCurveWindowedProfile(dna);
            var fragment = new ZCurveWindowedProfile(dna.Substring(50_250, 20_000));   // off the window grid
            var (score, offset) = ZCurveWindowedProfile.BestAlignment(fragment, reference);

            Assert.That(ZCurveWindowedProfile.RefineAlignment(fragment, reference, offset + 3, radius: 4), Is.EqualTo((score, offset)));
            Assert.That(ZCurveWindowedProfile.RefineAlignment(fragment, reference, offset - 4, radius: 4), Is.EqualTo((score, offset)));
        }

        [Test]
        public void RefineAlignmentClipsItsNeighbourhoodToTheReference()
        {
            var dna = Generator.GetRandomDNA(100_000, seed: 2);
            var reference = new ZCurveWindowedProfile(dna);
            var start = new ZCurveWindowedProfile(dna.Substring(0, 10_000));
            var end = new ZCurveWindowedProfile(dna.Substring(90_000));

            Assert.That(ZCurveWindowedProfile.RefineAlignment(start, reference, -3, radius: 4), Is.EqualTo((1.0, 0)));
            Assert.That(ZCurveWindowedProfile.RefineAlignment(end, reference, 182, radius: 4), Is.EqualTo((1.0, 180)));
            Assert.That(ZCurveWindowedProfile.RefineAlignment(start, reference, 300, radius: 4), Is.EqualTo((double.NaN, -1)));
            Assert.That(() => ZCurveWindowedProfile.RefineAlignment(start, reference, 0, radius: -1), Throws.InstanceOf<ArgumentOutOfRangeException>());
        }

        [Test]
        public void LadderRejectsMismatchedOrInvalidFactors()
        {
            var fine = new ZCurveWindowedProfile(Generator.GetRandomDNA(10_000, seed: 1), windowSize: 500);

            Assert.That(() => new ZCurveProfileLadder(fine, new[] { 4, 16 }), Throws.ArgumentException);
            Assert.That(() => new ZCurveProfileLadder(fine, new[] { 1, 16, 4 }), Throws.ArgumentException);
            Assert.That(() => ZCurveProfileLadder.SelectFactor(new ZCurveProfileLadder(fine), new ZCurveProfileLadder(fine, new[] { 1, 2 })),
                Throws.ArgumentException);
            Assert.That(() => ZCurveProfileLadder.CompareForQuery(new ZCurveProfileLadder(fine), new ZCurveProfileLadder(fine, new[] { 1, 2 })),
                Throws.ArgumentException);
        }

        [Test]
        public void BinaryFormatRoundTripsExactlyWithoutTheDnaString()
        {
            var fresh = new ZCurveWindowedProfile(Generator.GetRandomDNA(12_345, seed: 5), windowSize: 500);

            var bytes = fresh.ToBinary();
            var restored = ZCurveWindowedProfile.FromBinary(bytes);

            Assert.That(bytes.Length, Is.EqualTo(16 + 12 * fresh.WindowCount));
            Assert.That(restored.DnaString, Is.Null, "the stored format must never carry the sequence");
            Assert.That(restored.WindowSize, Is.EqualTo(fresh.WindowSize));
            Assert.That(restored.Windows, Is.EqualTo(fresh.Windows));
            Assert.That(ZCurveWindowedProfile.FromBinary(fresh.Downsample(4).ToBinary()).Windows,
                Is.EqualTo(fresh.Downsample(4).Windows));
        }

        [Test]
        public void BinaryFormatRejectsCorruptInput()
        {
            var bytes = new ZCurveWindowedProfile(Generator.GetRandomDNA(5_000, seed: 6), windowSize: 500).ToBinary();

            Assert.That(() => ZCurveWindowedProfile.FromBinary(bytes.AsSpan(0, bytes.Length - 1)), Throws.TypeOf<FormatException>());
            var wrongMagic = (byte[])bytes.Clone();
            wrongMagic[0] ^= 0xFF;
            Assert.That(() => ZCurveWindowedProfile.FromBinary(wrongMagic), Throws.TypeOf<FormatException>());
        }

        [Test]
        public void ZCurveDistanceIsTwiceTheBaseCompositionDistance()
        {
            // The Z-curve axes are three orthogonal Hadamard rows of norm 2 applied to (a,c,g,t), and
            // the difference of two composition vectors has no component along the fourth row
            // (1,1,1,1). So the Euclidean distance between two windows is exactly twice the distance
            // between their base frequencies - any rotation of this frame would change nothing.
            var x = Generator.GetRandomDNA(500, seed: 10);
            var y = BiasedDNA(new[] { 'G', 'G', 'C', 'A', 'T' }, 500, seed: 11);
            var zx = new ZCurveWindowedProfile(x, windowSize: 500).Windows[0];
            var zy = new ZCurveWindowedProfile(y, windowSize: 500).Windows[0];

            double[] Composition(string s) => "ACGT".Select(b => s.Count(c => c == b) / (double)s.Length).ToArray();
            var cx = Composition(x);
            var cy = Composition(y);

            var zDistance = Math.Sqrt(Math.Pow(zx.X - zy.X, 2) + Math.Pow(zx.Y - zy.Y, 2) + Math.Pow(zx.Z - zy.Z, 2));
            var compositionDistance = Math.Sqrt(cx.Zip(cy, (a, b) => (a - b) * (a - b)).Sum());

            Assert.That(zDistance, Is.EqualTo(2 * compositionDistance).Within(1e-6));
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
