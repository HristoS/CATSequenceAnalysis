namespace SequenceAnalyses
{
    /// <summary>How a <see cref="PValue"/> was obtained.</summary>
    public enum PValueKind
    {
        /// <summary>At least 10 null scores reach the score: (k + 1) / (N + 1).</summary>
        Empirical,
        /// <summary>Fewer than 10 do: (k + 1) / (N + 1) is only an upper bound, the true p is likely smaller.</summary>
        Bound,
    }

    public readonly record struct PValue(double Log10, PValueKind Kind);

    /// <summary>
    /// One cell of a <see cref="ScoreNull"/>: the best-alignment scores of unrelated pairs that share a level
    /// (window size), a bucket of the shorter profile's window count and a bucket of the number of offsets
    /// tried, sorted ascending.
    ///
    /// There is deliberately no fitted tail beyond the sample. On the local catalogue's null, held-out halves of the
    /// full cells showed a whole-cell Gumbel and an exponential top-decile tail both far too conservative (0 and 2
    /// exceedances where 600 were expected at p = 1e-3), and a generalized Pareto too liberal (2.1x at 1e-3, 13.7x
    /// at 1e-4): results/S36_TailCalibration.csv. Past the sample, p is reported as the sample's bound.
    /// </summary>
    public sealed record ScoreNullCell(int WindowSize, int WindowsBucket, int OffsetsBucket, float[] SortedScores)
    {
        private const int minEmpiricalExceedances = 10;

        public static ScoreNullCell FromScores(int windowSize, int windowsBucket, int offsetsBucket, IEnumerable<double> scores)
        {
            var sorted = scores.Select(s => (float)s).Order().ToArray();
            if (sorted.Length == 0) throw new ArgumentException("A cell needs at least one score.", nameof(scores));
            return new ScoreNullCell(windowSize, windowsBucket, offsetsBucket, sorted);
        }

        /// <summary>
        /// P(null score >= <paramref name="score"/>) = (k + 1) / (N + 1), with k the null scores reaching it. With fewer
        /// than 10 such scores that is only an upper bound, and the kind says so. A bound depends on N, so it isn't
        /// comparable across cells of different sizes: it's reported, not ranked by (results/S36_CatalogueRanking.csv:
        /// ranking by it put an identical sequence at 416).
        /// </summary>
        public PValue PValue(double score)
        {
            int n = SortedScores.Length;
            int exceedances = n - LowerBound(SortedScores, (float)score);
            double log10 = Math.Log10((exceedances + 1.0) / (n + 1.0));
            return new PValue(log10, exceedances >= minEmpiricalExceedances ? PValueKind.Empirical : PValueKind.Bound);
        }

        // First index whose value is >= target.
        private static int LowerBound(float[] sorted, float target)
        {
            int low = 0, high = sorted.Length;
            while (low < high)
            {
                int middle = (low + high) / 2;
                if (sorted[middle] < target) low = middle + 1;
                else high = middle;
            }
            return low;
        }
    }

    /// <summary>
    /// Empirical null for <see cref="ZCurveWindowedProfile.BestAlignment"/> scores (platform R5, P3 step S3.6).
    /// A best-alignment score is a maximum over offsets of an average over the shorter profile's windows, so
    /// what an unrelated pair scores depends on the level, on that window count (fewer windows, more spread)
    /// and on how many offsets were tried (more tries, higher maximum). A p-value from the cell with the same
    /// shape accounts for all three. Cells are keyed by window size and by two buckets of half-powers of two.
    /// </summary>
    public sealed class ScoreNull
    {
        private const int minNeighbourSamples = 100;

        private readonly Dictionary<(int WindowSize, int WindowsBucket, int OffsetsBucket), ScoreNullCell> cells;

        public IReadOnlyCollection<ScoreNullCell> Cells => cells.Values;

        public ScoreNull(IEnumerable<ScoreNullCell> cells) =>
            this.cells = cells.ToDictionary(c => (c.WindowSize, c.WindowsBucket, c.OffsetsBucket));

        /// <summary>floor(2 log2 count): buckets 1, 1.41, 2, 2.83, 4, ... wide.</summary>
        public static int Bucket(int count) =>
            count >= 1 ? (int)Math.Floor(2 * Math.Log2(count)) : throw new ArgumentOutOfRangeException(nameof(count));

        /// <summary>The shorter profile's window count and the offsets BestAlignment tries for two window counts.</summary>
        public static (int ShorterWindows, int Offsets) Shape(int windowsA, int windowsB) =>
            (Math.Min(windowsA, windowsB), Math.Abs(windowsA - windowsB) + 1);

        /// <summary>
        /// The p-value of <paramref name="score"/> for an alignment of that shape, from its own cell or, when that cell
        /// is missing, the nearest cell of the same window size with at least 100 scores; null when there is none.
        /// </summary>
        public PValue? PValue(int windowSize, int shorterWindows, int offsets, double score)
        {
            if (shorterWindows < 1) return null;
            int windowsBucket = Bucket(shorterWindows), offsetsBucket = Bucket(offsets);
            // A neighbour with the same window count is preferred: the window count sets the spread of the average,
            // the offset count only how far into the tail the maximum reaches.
            if (!cells.TryGetValue((windowSize, windowsBucket, offsetsBucket), out var cell))
                cell = cells.Values
                    .Where(c => c.WindowSize == windowSize && c.SortedScores.Length >= minNeighbourSamples)
                    .MinBy(c => Math.Abs(c.WindowsBucket - windowsBucket) * 4 + Math.Abs(c.OffsetsBucket - offsetsBucket));
            return cell?.PValue(score);
        }
    }
}
