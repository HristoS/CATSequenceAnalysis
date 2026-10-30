namespace SequenceAnalyses
{
    /// <summary>
    /// Fixed-length whole-genome signatures, the composition baselines a Z-curve tree is judged against (P5, R1's list), plus the
    /// Z-curve's own window summary. Each one that the platform needs is computable from what it stores: the dinucleotide counts
    /// (D8) or the 500 bp profile. A window's base counts follow exactly from its Z-curve averages, because A, C, G and T map to
    /// (1, 1, -1), (-1, 1, 1), (1, -1, 1) and (-1, -1, -1): with W bases, A = W (1 + x + y - z) / 4 and so on. The string versions
    /// are the exact definitions, used to check the profile-derived ones.
    /// </summary>
    public static class GenomeSignatures
    {
        private const string bases = "ACGT";

        /// <summary>
        /// Per axis: |mean| for x and y (they change sign with the strand), the mean for z (GC, the same on both strands), the SD
        /// across windows and the lag-1 autocorrelation across windows. 9 values.
        /// </summary>
        public static double[] WindowSummary(ZCurveWindowedProfile level)
        {
            var result = new List<double>();
            foreach (var (values, signed) in new[] { (level.Windows.Select(w => w.X).ToArray(), false),
                (level.Windows.Select(w => w.Y).ToArray(), false), (level.Windows.Select(w => w.Z).ToArray(), true) })
            {
                double mean = values.Average();
                double variance = values.Sum(v => (v - mean) * (v - mean)) / values.Length;
                double lag1 = 0;
                for (int i = 1; i < values.Length; i++) lag1 += (values[i] - mean) * (values[i - 1] - mean);
                result.Add(signed ? mean : Math.Abs(mean));
                result.Add(Math.Sqrt(variance));
                result.Add(variance == 0 ? 0 : lag1 / ((values.Length - 1) * variance));
            }
            return result.ToArray();
        }

        /// <summary>The whole-genome mean (x, y, z): the one-point "global" Z-curve, a function of base composition only.</summary>
        public static double[] ZCurveMean(ZCurveWindowedProfile fine) =>
            new[] { fine.Windows.Average(w => w.X), fine.Windows.Average(w => w.Y), fine.Windows.Average(w => w.Z) };

        public static double[] ZCurveMean(string sequence)
        {
            double x = 0, y = 0, z = 0;
            long count = 0;
            foreach (char c in sequence)
            {
                var (ax, ay, az) = ZCurveProfile.Axes(c);
                if (ax == 0 && ay == 0 && az == 0) continue;
                (x, y, z, count) = (x + ax, y + ay, z + az, count + 1);
            }
            return new[] { x / count, y / count, z / count };
        }

        /// <summary>
        /// Zheng-style geometric summary of the cumulative Z-curve (after Zheng et al. 2005, whose exact normalization is unverified):
        /// the curve's points divided by the genome length, then their centre (3 values) and the three principal directions of the
        /// point cloud, largest first, each signed so its largest component is positive (9 values). From the profile the points are
        /// the window boundaries, every 500 bp.
        /// </summary>
        public static double[] ZhengCurve(ZCurveWindowedProfile fine)
        {
            var moments = new CurveMoments();
            double x = 0, y = 0, z = 0;
            double length = (double)fine.WindowCount * fine.WindowSize;
            foreach (var w in fine.Windows)
            {
                (x, y, z) = (x + w.X * fine.WindowSize, y + w.Y * fine.WindowSize, z + w.Z * fine.WindowSize);
                moments.Add(x / length, y / length, z / length);
            }
            return moments.Summary();
        }

        public static double[] ZhengCurve(string sequence)
        {
            var moments = new CurveMoments();
            double x = 0, y = 0, z = 0, length = sequence.Length;
            foreach (char c in sequence)
            {
                var (ax, ay, az) = ZCurveProfile.Axes(c);
                (x, y, z) = (x + ax, y + ay, z + az);
                moments.Add(x / length, y / length, z / length);
            }
            return moments.Summary();
        }

        /// <summary>
        /// The Natural Vector (Deng et al. 2011; Yu et al. 2013) in A, C, G, T order: counts n_k, mean positions mu_k, and
        /// D_k = sum (p - mu_k)^2 / (n_k * N). 12 values. From the profile each window's counts are solved from its averages and
        /// placed at the window's centre, plus the within-window spread of a uniform placement, (W^2 - 1) / 12 per base; N is the
        /// genome length, so the dropped partial last window is the only other loss.
        /// </summary>
        public static double[] NaturalVector(ZCurveWindowedProfile fine, long lengthBp)
        {
            int size = fine.WindowSize;
            var windows = fine.Windows;
            var counts = new double[windows.Count, 4];
            for (int w = 0; w < windows.Count; w++)
            {
                var (x, y, z) = windows[w];
                counts[w, 0] = Math.Max(0, size * (1 + x + y - z) / 4);
                counts[w, 1] = Math.Max(0, size * (1 - x + y + z) / 4);
                counts[w, 2] = Math.Max(0, size * (1 + x - y + z) / 4);
                counts[w, 3] = Math.Max(0, size * (1 - x - y - z) / 4);
            }
            var result = new double[12];
            for (int k = 0; k < 4; k++)
            {
                double n = 0, sum = 0;
                for (int w = 0; w < windows.Count; w++)
                {
                    n += counts[w, k];
                    sum += counts[w, k] * Centre(w);
                }
                double mu = n == 0 ? 0 : sum / n, spread = 0;
                for (int w = 0; w < windows.Count; w++)
                    spread += counts[w, k] * ((Centre(w) - mu) * (Centre(w) - mu) + (size * (double)size - 1) / 12);
                (result[k], result[4 + k], result[8 + k]) = (n, mu, n == 0 ? 0 : spread / (n * lengthBp));
            }
            return result;

            double Centre(int w) => (double)w * size + (size + 1) / 2.0;   // 1-based positions, as the definition uses
        }

        public static double[] NaturalVector(string sequence)
        {
            var n = new double[4];
            var sum = new double[4];
            var squares = new double[4];
            for (int i = 0; i < sequence.Length; i++)
            {
                int k = bases.IndexOf(char.ToUpperInvariant(sequence[i]));
                if (k < 0) continue;
                double p = i + 1;
                (n[k], sum[k], squares[k]) = (n[k] + 1, sum[k] + p, squares[k] + p * p);
            }
            var result = new double[12];
            for (int k = 0; k < 4; k++)
            {
                double mu = n[k] == 0 ? 0 : sum[k] / n[k];
                double spread = squares[k] - n[k] * mu * mu;
                (result[k], result[4 + k], result[8 + k]) = (n[k], mu, n[k] == 0 ? 0 : spread / (n[k] * sequence.Length));
            }
            return result;
        }

        /// <summary>Counts of adjacent base pairs, index 4 * first + second in A, C, G, T; an ambiguity code breaks the pair (as SequenceFeatureReader).</summary>
        public static long[] DinucleotideCounts(string sequence)
        {
            var counts = new long[16];
            int previous = -1;
            foreach (char c in sequence)
            {
                int k = bases.IndexOf(char.ToUpperInvariant(c));
                if (k >= 0 && previous >= 0) counts[4 * previous + k]++;
                previous = k;
            }
            return counts;
        }

        /// <summary>
        /// Karlin's symmetrized dinucleotide relative abundance rho*(XY) = f*(XY) / (f*(X) f*(Y)), over the sequence and its reverse
        /// complement (Karlin and Burge 1995), so it's the same for either strand. 16 values, index 4 * first + second.
        /// </summary>
        public static double[] KarlinRho(long[] dinucleotideCounts)
        {
            var both = new double[16];
            for (int first = 0; first < 4; first++)
                for (int second = 0; second < 4; second++)
                    both[4 * first + second] = dinucleotideCounts[4 * first + second] + dinucleotideCounts[4 * (3 - second) + (3 - first)];
            double total = both.Sum();
            var single = new double[4];
            for (int i = 0; i < 16; i++) single[i / 4] += both[i] / total;   // A, C, G, T are 0..3, so 3 - k is the complement
            var rho = new double[16];
            for (int i = 0; i < 16; i++)
                rho[i] = single[i / 4] * single[i % 4] == 0 ? 0 : both[i] / total / (single[i / 4] * single[i % 4]);
            return rho;
        }

        /// <summary>Karlin's delta*: the mean absolute difference of rho* over the 16 dinucleotides.</summary>
        public static double KarlinDelta(double[] rhoA, double[] rhoB) => rhoA.Zip(rhoB, (a, b) => Math.Abs(a - b)).Average();

        private sealed class CurveMoments
        {
            private readonly double[] sum = new double[3];
            private readonly double[,] products = new double[3, 3];
            private long count;

            public void Add(double x, double y, double z)
            {
                var point = new[] { x, y, z };
                for (int i = 0; i < 3; i++)
                {
                    sum[i] += point[i];
                    for (int j = 0; j < 3; j++) products[i, j] += point[i] * point[j];
                }
                count++;
            }

            public double[] Summary()
            {
                var centre = sum.Select(s => s / count).ToArray();
                var covariance = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        covariance[i, j] = products[i, j] / count - centre[i] * centre[j];
                var (values, vectors) = SymmetricEigen(covariance);
                var result = new List<double>(centre);
                foreach (int column in Enumerable.Range(0, 3).OrderByDescending(c => values[c]))
                {
                    var v = Enumerable.Range(0, 3).Select(r => vectors[r, column]).ToArray();
                    double sign = Math.Sign(v.MaxBy(Math.Abs));
                    result.AddRange(v.Select(e => e * (sign == 0 ? 1 : sign)));
                }
                return result.ToArray();
            }
        }

        // Cyclic Jacobi rotations for a symmetric 3x3 matrix: eigenvalues on the diagonal, eigenvectors as the columns of v.
        private static (double[] Values, double[,] Vectors) SymmetricEigen(double[,] matrix)
        {
            var a = (double[,])matrix.Clone();
            var v = new double[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            for (int sweep = 0; sweep < 50; sweep++)
            {
                double off = Math.Abs(a[0, 1]) + Math.Abs(a[0, 2]) + Math.Abs(a[1, 2]);
                if (off < 1e-30) break;
                foreach (var (p, q) in new[] { (0, 1), (0, 2), (1, 2) })
                {
                    if (a[p, q] == 0) continue;
                    double theta = (a[q, q] - a[p, p]) / (2 * a[p, q]);
                    double t = Math.Sign(theta == 0 ? 1 : theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                    double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
                    for (int k = 0; k < 3; k++)
                    {
                        double akp = a[k, p], akq = a[k, q];
                        a[k, p] = c * akp - s * akq;
                        a[k, q] = s * akp + c * akq;
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        double apk = a[p, k], aqk = a[q, k];
                        a[p, k] = c * apk - s * aqk;
                        a[q, k] = s * apk + c * aqk;
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        double vkp = v[k, p], vkq = v[k, q];
                        v[k, p] = c * vkp - s * vkq;
                        v[k, q] = s * vkp + c * vkq;
                    }
                }
            }
            return (new[] { a[0, 0], a[1, 1], a[2, 2] }, v);
        }
    }
}
