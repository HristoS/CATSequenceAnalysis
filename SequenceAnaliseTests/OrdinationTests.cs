using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class OrdinationTests
    {
        [Test]
        public void SymmetricEigenRecoversAKnownDecomposition()
        {
            // [[2, 1, 0], [1, 2, 0], [0, 0, 5]]: eigenvalues 1, 3 and 5.
            var (values, vectors) = Ordination.SymmetricEigen(new double[,] { { 2, 1, 0 }, { 1, 2, 0 }, { 0, 0, 5 } });
            Assert.That(values.Order(), Is.EqualTo(new[] { 1.0, 3, 5 }).Within(1e-12));
            int three = Array.FindIndex(values, v => Math.Abs(v - 3) < 1e-9);
            Assert.That(Math.Abs(vectors[0, three]), Is.EqualTo(Math.Sqrt(0.5)).Within(1e-12));
            Assert.That(vectors[0, three] * vectors[1, three], Is.GreaterThan(0));
        }

        [Test]
        public void PcaFindsTheDirectionOfSpreadAndItsShare()
        {
            // Points spread along (1, 1, 0) with a little noise along (1, -1, 0) and none along z.
            var rng = new Random(4);
            var rows = Enumerable.Range(0, 500).Select(_ =>
            {
                double along = rng.NextDouble() * 10 - 5, across = (rng.NextDouble() - 0.5) * 0.2;
                return new[] { along + across + 7, along - across - 3, 2.0 };
            }).ToArray();

            var pca = Ordination.Pca(rows, 2);

            Assert.That(pca.Loadings[0][0], Is.EqualTo(Math.Sqrt(0.5)).Within(1e-3));
            Assert.That(pca.Loadings[0][1], Is.EqualTo(Math.Sqrt(0.5)).Within(1e-3));
            Assert.That(pca.ExplainedVariance[0], Is.GreaterThan(0.999));
            Assert.That(pca.Means, Is.EqualTo(new[] { rows.Average(r => r[0]), rows.Average(r => r[1]), 2.0 }).Within(1e-9));
            // The first score is the position along the line: sqrt(2) times "along", up to the mean.
            Assert.That(pca.Scores.Select(s => s[0]).Average(), Is.EqualTo(0).Within(1e-9));
        }
    }
}
