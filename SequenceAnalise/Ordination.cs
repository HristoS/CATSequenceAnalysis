namespace SequenceAnalyses
{
    /// <summary>Principal components of a set of feature vectors: scores per row, loadings per component, variance explained.</summary>
    /// <param name="Loadings">Loadings[c][f]: feature f's weight in component c, each component a unit vector signed so its largest weight is positive.</param>
    public sealed record PcaResult(double[][] Scores, double[][] Loadings, double[] ExplainedVariance, double[] Means);

    public static class Ordination
    {
        /// <summary>
        /// PCA of the rows as given (centred, not rescaled: use it where the features share one unit, and standardize first where
        /// they don't). Components are ordered by variance; ExplainedVariance is each one's share of the total.
        /// </summary>
        public static PcaResult Pca(double[][] rows, int components)
        {
            int n = rows.Length, width = rows[0].Length;
            var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
            var covariance = new double[width, width];
            foreach (var row in rows)
                for (int i = 0; i < width; i++)
                    for (int j = 0; j < width; j++)
                        covariance[i, j] += (row[i] - means[i]) * (row[j] - means[j]) / n;

            var (values, vectors) = SymmetricEigen(covariance);
            double total = values.Sum();
            var order = Enumerable.Range(0, width).OrderByDescending(c => values[c]).Take(components).ToArray();
            var loadings = order.Select(c =>
            {
                var v = Enumerable.Range(0, width).Select(r => vectors[r, c]).ToArray();
                double sign = Math.Sign(v.MaxBy(Math.Abs));
                return v.Select(e => e * (sign == 0 ? 1 : sign)).ToArray();
            }).ToArray();
            var scores = rows.Select(r => loadings.Select(l => l.Select((w, f) => w * (r[f] - means[f])).Sum()).ToArray()).ToArray();
            return new PcaResult(scores, loadings, order.Select(c => total == 0 ? 0 : values[c] / total).ToArray(), means);
        }

        /// <summary>Cyclic Jacobi rotations for a symmetric matrix: eigenvalues, and eigenvectors as the columns of Vectors.</summary>
        public static (double[] Values, double[,] Vectors) SymmetricEigen(double[,] matrix)
        {
            int size = matrix.GetLength(0);
            var a = (double[,])matrix.Clone();
            var v = new double[size, size];
            for (int i = 0; i < size; i++) v[i, i] = 1;
            for (int sweep = 0; sweep < 100; sweep++)
            {
                double off = 0;
                for (int p = 0; p < size; p++)
                    for (int q = p + 1; q < size; q++)
                        off += Math.Abs(a[p, q]);
                if (off < 1e-30) break;
                for (int p = 0; p < size; p++)
                    for (int q = p + 1; q < size; q++)
                    {
                        if (a[p, q] == 0) continue;
                        double theta = (a[q, q] - a[p, p]) / (2 * a[p, q]);
                        double t = Math.Sign(theta == 0 ? 1 : theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                        double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
                        for (int k = 0; k < size; k++)
                        {
                            double akp = a[k, p], akq = a[k, q];
                            a[k, p] = c * akp - s * akq;
                            a[k, q] = s * akp + c * akq;
                        }
                        for (int k = 0; k < size; k++)
                        {
                            double apk = a[p, k], aqk = a[q, k];
                            a[p, k] = c * apk - s * aqk;
                            a[q, k] = s * apk + c * aqk;
                        }
                        for (int k = 0; k < size; k++)
                        {
                            double vkp = v[k, p], vkq = v[k, q];
                            v[k, p] = c * vkp - s * vkq;
                            v[k, q] = s * vkp + c * vkq;
                        }
                    }
            }
            return (Enumerable.Range(0, size).Select(i => a[i, i]).ToArray(), v);
        }
    }
}
