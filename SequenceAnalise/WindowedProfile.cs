using System;
using System.Collections.Generic;

namespace SequenceAnalyses
{
    /// <summary>
    /// Alignment-free positional profile. Splits a DNA sequence into fixed-size windows and
    /// represents each window by its k-mer frequency vector. Two profiles are compared by sliding
    /// the shorter window-vector along the longer one to find the best-matching contiguous region,
    /// which both scores similarity and localizes where a fragment sits within a larger sequence.
    ///
    /// Complements <see cref="CATProfile"/>: CATProfile compresses a whole sequence to a single
    /// point — fast and O(n), but composition-only, so it cannot separate unrelated sequences that
    /// share base composition (empirically phiX174 vs lambda phage score ~0.99). WindowedProfile
    /// keeps per-window k-mer detail, which restores that discrimination: on real sequences it
    /// localizes exact fragments to the correct window offset, detects diverged fragments down to
    /// ~70% identity, and separates unrelated sequences (cosine ~0.5–0.67) from true fragments
    /// (~0.77–1.0). Cost: O(n) build, and memory O((n / windowSize) * distinct-k-mers-per-window),
    /// so it is heavier than a single CATProfile point — suited to reads and small/medium genomes
    /// or read-vs-genome localization rather than whole-large-genome-vs-genome at once.
    /// </summary>
    public sealed class WindowedProfile
    {
        public const int DefaultWindowSize = 500;
        public const int DefaultK = 4;

        /// <summary>
        /// One window's k-mer content, stored sparsely (k-mer id → count) so memory is bounded by
        /// the window size, not by 4^k. The precomputed L2 norm makes cosine O(distinct k-mers).
        /// </summary>
        public sealed class Window
        {
            internal readonly Dictionary<int, int> Counts;
            internal readonly double Norm;

            internal Window(Dictionary<int, int> counts)
            {
                Counts = counts;
                double sumSquares = 0;
                foreach (var c in counts.Values) sumSquares += (double)c * c;
                Norm = Math.Sqrt(sumSquares);
            }

            /// <summary>Number of distinct k-mers observed in this window.</summary>
            public int DistinctKmers => Counts.Count;
        }

        public string DnaString { get; }

        public int WindowSize { get; }

        public int K { get; }

        /// <summary>Ordered per-window k-mer profiles; window i covers bases [i*WindowSize, (i+1)*WindowSize).</summary>
        public IReadOnlyList<Window> Windows { get; }

        public WindowedProfile(string dnaString, int windowSize = DefaultWindowSize, int k = DefaultK)
        {
            if (dnaString is null) throw new ArgumentNullException(nameof(dnaString));
            if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));
            // 2*k must stay within an int rolling hash; k also can't exceed the window it lives in.
            if (k < 1 || k > 15 || k > windowSize) throw new ArgumentOutOfRangeException(nameof(k));

            DnaString = dnaString;
            WindowSize = windowSize;
            K = k;

            var windows = new List<Window>(dnaString.Length / windowSize);
            for (int start = 0; start + windowSize <= dnaString.Length; start += windowSize)
                windows.Add(BuildWindow(dnaString, start, windowSize, k));
            Windows = windows;
        }

        private static Window BuildWindow(string seq, int start, int length, int k)
        {
            var counts = new Dictionary<int, int>();
            int mask = (1 << (2 * k)) - 1;
            int rolling = 0, valid = 0, end = start + length;

            for (int i = start; i < end; i++)
            {
                int baseIndex = BaseIndex(seq[i]);
                if (baseIndex < 0)
                {
                    // Ambiguity code (N etc.) breaks the current k-mer; restart the rolling window.
                    valid = 0;
                    rolling = 0;
                    continue;
                }

                rolling = ((rolling << 2) | baseIndex) & mask;
                if (++valid >= k)
                    counts[rolling] = counts.TryGetValue(rolling, out var c) ? c + 1 : 1;
            }

            return new Window(counts);
        }

        private static int BaseIndex(char c) => c switch
        {
            'A' or 'a' => 0,
            'C' or 'c' => 1,
            'G' or 'g' => 2,
            'T' or 't' => 3,
            _ => -1
        };

        private static double Cosine(Window a, Window b)
        {
            if (a.Norm == 0 || b.Norm == 0) return 0;

            // Walk the smaller map and probe the larger, so the dot product is O(min distinct k-mers).
            var (small, large) = a.Counts.Count <= b.Counts.Count ? (a.Counts, b.Counts) : (b.Counts, a.Counts);
            double dot = 0;
            foreach (var kv in small)
                if (large.TryGetValue(kv.Key, out var otherCount))
                    dot += (double)kv.Value * otherCount;

            return dot / (a.Norm * b.Norm);
        }

        /// <summary>
        /// Slides <paramref name="query"/> (which must have no more windows than
        /// <paramref name="reference"/>) along the reference and returns the best average
        /// per-window cosine similarity together with the window offset where that alignment
        /// starts. <c>Offset * WindowSize</c> is the approximate base position of the match in the
        /// reference. Returns (NaN, -1) when the query is empty or longer than the reference.
        /// </summary>
        public static (double Score, int Offset) BestAlignment(WindowedProfile query, WindowedProfile reference)
        {
            if (query.WindowSize != reference.WindowSize || query.K != reference.K)
                throw new ArgumentException("Profiles must share the same window size and k.");

            var q = query.Windows;
            var r = reference.Windows;
            if (q.Count == 0 || r.Count < q.Count) return (double.NaN, -1);

            double best = double.NegativeInfinity;
            int bestOffset = -1;
            for (int offset = 0; offset + q.Count <= r.Count; offset++)
            {
                double sum = 0;
                for (int j = 0; j < q.Count; j++)
                    sum += Cosine(q[j], r[offset + j]);

                double score = sum / q.Count;
                if (score > best)
                {
                    best = score;
                    bestOffset = offset;
                }
            }

            return (best, bestOffset);
        }

        /// <summary>
        /// Symmetric similarity in [0, 1]: treats the profile with fewer windows as the query and
        /// localizes it within the other, returning the best-alignment score (0 if either is empty).
        /// </summary>
        public static double Compare(WindowedProfile x, WindowedProfile y)
        {
            var (shorter, longer) = x.Windows.Count <= y.Windows.Count ? (x, y) : (y, x);
            var (score, _) = BestAlignment(shorter, longer);
            return double.IsNaN(score) ? 0 : score;
        }
    }
}
