using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class ScoreNullTests
    {
        // 0.000, 0.001, ..., 0.999: a cell whose p-values can be read off by hand.
        private static ScoreNullCell Uniform(int windowsBucket = 6, int offsetsBucket = 6, int count = 1_000) =>
            ScoreNullCell.FromScores(500, windowsBucket, offsetsBucket, Enumerable.Range(0, count).Select(i => i / (double)count));

        [TestCase(1, 0)]
        [TestCase(2, 2)]
        [TestCase(3, 3)]
        [TestCase(8, 6)]
        [TestCase(11, 6)]
        [TestCase(12, 7)]
        public void BucketsAreHalfPowersOfTwo(int count, int expected)
        {
            Assert.That(ScoreNull.Bucket(count), Is.EqualTo(expected));
        }

        [Test]
        public void ShapeIsTheShorterWindowCountAndTheOffsetsTried()
        {
            Assert.That(ScoreNull.Shape(145, 100), Is.EqualTo((100, 46)));
            Assert.That(ScoreNull.Shape(20, 20), Is.EqualTo((20, 1)));
        }

        [Test]
        public void FromScoresSortsTheSample()
        {
            var cell = ScoreNullCell.FromScores(500, 6, 6, new[] { 0.7, 0.2, 0.5 });

            Assert.That(cell.SortedScores, Is.EqualTo(new[] { 0.2f, 0.5f, 0.7f }));
            Assert.That(() => ScoreNullCell.FromScores(500, 6, 6, Array.Empty<double>()), Throws.ArgumentException);
        }

        [Test]
        public void PValueIsEmpiricalWhileTenScoresReachItAndABoundPastThat()
        {
            var cell = Uniform();

            Assert.That(cell.PValue(0.5), Is.EqualTo(new PValue(Math.Log10(501.0 / 1001), PValueKind.Empirical)));
            Assert.That(cell.PValue(0.990), Is.EqualTo(new PValue(Math.Log10(11.0 / 1001), PValueKind.Empirical)));
            Assert.That(cell.PValue(0.991), Is.EqualTo(new PValue(Math.Log10(10.0 / 1001), PValueKind.Bound)));
            Assert.That(cell.PValue(1.5), Is.EqualTo(new PValue(Math.Log10(1.0 / 1001), PValueKind.Bound)),
                "past the sample, the sample's own bound and no extrapolation");
        }

        [Test]
        public void MissingCellFallsBackToTheNearestOneOfTheSameWindowSize()
        {
            var near = Uniform(windowsBucket: 6, offsetsBucket: 6, count: 200);
            var far = Uniform(windowsBucket: 12, offsetsBucket: 6, count: 300);
            var null_ = new ScoreNull(new[] { near, far });

            Assert.That(null_.PValue(500, 10, 7, 0.85), Is.EqualTo(near.PValue(0.85)), "windows bucket 6, offsets bucket 5");
            Assert.That(null_.PValue(2000, 10, 7, 0.85), Is.Null, "no cell at that window size");
            Assert.That(null_.PValue(500, 0, 7, 0.85), Is.Null);
        }
    }
}
