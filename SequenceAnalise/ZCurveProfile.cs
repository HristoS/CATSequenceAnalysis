using System;

namespace SequenceAnalyses
{
    /// <summary>
    /// Z-curve profile: maps each base onto three independent biochemical dichotomies - purine
    /// (A,G) vs pyrimidine (C,T), amino (A,C) vs keto (G,T), and strong/3-H-bond (G,C) vs weak/
    /// 2-H-bond (A,T) - then walks the sequence once, summing each axis. Each axis's running
    /// average (sum / bases-seen-so-far) is sampled at <see cref="CheckpointCount"/> evenly-spaced
    /// positions, giving a fixed-size trajectory - a coarse sketch of the walk's shape - in one
    /// real, shared 3D frame derived directly from the bases' own chemistry, rather than
    /// trilaterated against an arbitrary periodic benchmark string the way
    /// <see cref="CATProfile.Point"/> is.
    ///
    /// With CheckpointCount = 1 (the default), the single checkpoint is the whole-sequence
    /// endpoint - mathematically just a linear recombination of base composition (#A,#C,#G,#T), so
    /// it cannot separate two sequences that share composition but differ in base order (verified:
    /// shuffling a sequence leaves this endpoint exactly unchanged). That's the same ceiling
    /// documented for CATProfile's single-point compression in CLAUDE.md. Raising CheckpointCount
    /// trades that blindness for sensitivity to *regional* composition trends along the sequence
    /// (e.g. an AT-rich half followed by a GC-rich half reads differently from either half alone,
    /// and differently from a global shuffle of the same two halves) - still O(n) to build and
    /// O(CheckpointCount) to compare, without WindowedProfile's O(n/windowSize) k-mer memory. It
    /// remains blind to reordering *within* what a single checkpoint's span covers, and, like any
    /// running average, each checkpoint's own noise still shrinks as that span grows - so more
    /// checkpoints (finer positional resolution) trade off against noisier individual checkpoints.
    /// </summary>
    public sealed class ZCurveProfile
    {
        public const int DefaultCheckpointCount = 1;

        public string DnaString { get; }

        public int CheckpointCount { get; }

        /// <summary>
        /// Running-average (x,y,z) sampled at <see cref="CheckpointCount"/> evenly-spaced
        /// positions: Checkpoints[i] is the length-normalized cumulative sum over the prefix
        /// ending at base index (i+1) * DnaString.Length / CheckpointCount (integer division).
        /// </summary>
        public IReadOnlyList<(double X, double Y, double Z)> Checkpoints { get; }

        /// <summary>The whole-sequence endpoint - equivalent to the last checkpoint.</summary>
        public (double X, double Y, double Z) Point => Checkpoints[CheckpointCount - 1];

        public ZCurveProfile(string dnaString, int checkpointCount = DefaultCheckpointCount)
        {
            if (dnaString is null) throw new ArgumentNullException(nameof(dnaString));
            if (checkpointCount < 1 || checkpointCount > dnaString.Length)
                throw new ArgumentOutOfRangeException(nameof(checkpointCount));

            DnaString = dnaString;
            CheckpointCount = checkpointCount;

            int length = dnaString.Length;
            var checkpoints = new (double X, double Y, double Z)[checkpointCount];
            long sumX = 0, sumY = 0, sumZ = 0;
            int checkpointIndex = 0;
            long nextBoundary = (long)(checkpointIndex + 1) * length / checkpointCount;

            for (int i = 0; i < length; i++)
            {
                var (x, y, z) = Axes(dnaString[i]);
                sumX += x;
                sumY += y;
                sumZ += z;

                if (i + 1 == nextBoundary)
                {
                    double count = i + 1;
                    checkpoints[checkpointIndex] = (sumX / count, sumY / count, sumZ / count);
                    checkpointIndex++;
                    if (checkpointIndex < checkpointCount)
                        nextBoundary = (long)(checkpointIndex + 1) * length / checkpointCount;
                }
            }

            Checkpoints = checkpoints;
        }

        // Purine(+1)/pyrimidine(-1), amino(+1)/keto(-1), strong 3-H-bond(+1)/weak 2-H-bond(-1).
        // Internal so ZCurveWindowedProfile can reuse the same mapping.
        internal static (int X, int Y, int Z) Axes(char @base) => @base switch
        {
            'A' or 'a' => (1, 1, -1),
            'C' or 'c' => (-1, 1, 1),
            'G' or 'g' => (1, -1, 1),
            'T' or 't' => (-1, -1, -1),
            _ => (0, 0, 0), // ambiguity codes (N etc.) contribute nothing to any axis
        };

        /// <summary>
        /// Endpoint-only similarity in the same "1 - Euclidean distance" convention as
        /// <see cref="CATProfile.CompareTrilateration"/>: 1 for identical endpoints, falling off
        /// as the two sequences' composition-derived coordinates diverge. Ignores every checkpoint
        /// but the last; use <see cref="CompareTrajectory"/> to use all of them.
        /// </summary>
        public static double Compare(ZCurveProfile x, ZCurveProfile y)
        {
            var dx = x.Point.X - y.Point.X;
            var dy = x.Point.Y - y.Point.Y;
            var dz = x.Point.Z - y.Point.Z;
            return 1 - Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        /// <summary>
        /// Average per-checkpoint Euclidean distance, expressed as "1 - average distance" to match
        /// <see cref="Compare"/>'s convention. Both profiles must share the same CheckpointCount;
        /// since checkpoint positions are defined by *fraction* of each sequence's own length, two
        /// different-length sequences still compare corresponding relative positions (checkpoint i
        /// of a 1000bp sequence against checkpoint i of a 10000bp sequence lines up i/CheckpointCount
        /// of the way through each), not corresponding absolute base positions.
        /// </summary>
        public static double CompareTrajectory(ZCurveProfile x, ZCurveProfile y)
        {
            if (x.CheckpointCount != y.CheckpointCount)
                throw new ArgumentException("Profiles must share the same checkpoint count.");

            double sumDistance = 0;
            for (int i = 0; i < x.CheckpointCount; i++)
            {
                var dx = x.Checkpoints[i].X - y.Checkpoints[i].X;
                var dy = x.Checkpoints[i].Y - y.Checkpoints[i].Y;
                var dz = x.Checkpoints[i].Z - y.Checkpoints[i].Z;
                sumDistance += Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }

            return 1 - sumDistance / x.CheckpointCount;
        }
    }

    /// <summary>
    /// Windowed Z-curve profile: structured exactly like <see cref="WindowedProfile"/> (fixed
    /// window *size*, sliding-window alignment for comparison), but each window's feature is a
    /// local Z-curve axis average - purine/pyrimidine, amino/keto, strong/weak-H-bond, computed
    /// from only that window's own bases - instead of a k-mer frequency vector.
    ///
    /// Exists because <see cref="ZCurveProfile"/>'s fixed-*count* checkpoints still collapse at
    /// large N: with checkpoint count fixed, each checkpoint's span (and its own noise floor)
    /// grows with total sequence length, so unrelated sequences drift toward looking identical the
    /// longer they get (confirmed empirically in ZCurveProfileTests -
    /// TrajectoryWithFixedCheckpointCountStillCollapsesAtLargeEnoughLength). Fixing window size
    /// instead - the same choice WindowedProfile already makes for k-mers - keeps each window's
    /// noise floor constant regardless of total length, so window count grows with the sequence
    /// instead of window span growing with it. Far cheaper per window than WindowedProfile's k-mer
    /// counts (3 doubles vs a sparse dictionary), at the cost of the coarser signal a 3-axis
    /// composition average carries versus full k-mer content.
    /// </summary>
    public sealed class ZCurveWindowedProfile
    {
        public const int DefaultWindowSize = WindowedProfile.DefaultWindowSize;

        // Nullable: a profile rehydrated via FromPrecomputedWindows may not have the original
        // string available (that's the point - see FromPrecomputedWindows below). Compare and
        // BestAlignment never read this field, only WindowSize and Windows.
        public string? DnaString { get; }

        public int WindowSize { get; }

        /// <summary>Ordered per-window local Z-curve axis averages; window i covers bases [i*WindowSize, (i+1)*WindowSize).</summary>
        public IReadOnlyList<(double X, double Y, double Z)> Windows { get; }

        public ZCurveWindowedProfile(string dnaString, int windowSize = DefaultWindowSize)
        {
            if (dnaString is null) throw new ArgumentNullException(nameof(dnaString));
            if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));

            DnaString = dnaString;
            WindowSize = windowSize;

            var windows = new List<(double X, double Y, double Z)>(dnaString.Length / windowSize);
            for (int start = 0; start + windowSize <= dnaString.Length; start += windowSize)
            {
                long sumX = 0, sumY = 0, sumZ = 0;
                for (int i = start; i < start + windowSize; i++)
                {
                    var (x, y, z) = ZCurveProfile.Axes(dnaString[i]);
                    sumX += x;
                    sumY += y;
                    sumZ += z;
                }
                windows.Add((sumX / (double)windowSize, sumY / (double)windowSize, sumZ / (double)windowSize));
            }
            Windows = windows;
        }

        private ZCurveWindowedProfile(string? dnaString, int windowSize, IReadOnlyList<(double X, double Y, double Z)> windows)
        {
            DnaString = dnaString;
            WindowSize = windowSize;
            Windows = windows;
        }

        /// <summary>
        /// Rehydrates a profile directly from previously-computed window vectors, without
        /// re-deriving them from the raw DNA string - lets a cached profile be reloaded and
        /// compared without ever touching the original sequence again (verified identical to a
        /// freshly-built profile in ZCurveWindowedProfileTests).
        ///
        /// <paramref name="dnaString"/> is optional; when supplied, its length is cross-checked
        /// against <paramref name="windows"/>.Count, otherwise that check is skipped and becomes
        /// the caller's responsibility.
        /// </summary>
        public static ZCurveWindowedProfile FromPrecomputedWindows(
            int windowSize,
            IReadOnlyList<(double X, double Y, double Z)> windows,
            string? dnaString = null)
        {
            if (windows is null) throw new ArgumentNullException(nameof(windows));
            if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));
            if (dnaString is not null)
            {
                int expectedWindows = dnaString.Length / windowSize;
                if (windows.Count != expectedWindows)
                    throw new ArgumentException(
                        $"windows.Count ({windows.Count}) does not match the {expectedWindows} windows " +
                        $"implied by dnaString.Length ({dnaString.Length}) and windowSize ({windowSize}).",
                        nameof(windows));
            }

            // Defensive copy: unlike the string-based constructor (which always builds its own
            // fresh List), this factory would otherwise share storage with a caller-owned list
            // that could be mutated later out from under the profile.
            return new ZCurveWindowedProfile(dnaString, windowSize, windows.ToList());
        }

        /// <summary>
        /// Builds a profile by reading bases from <paramref name="reader"/> in fixed-size chunks,
        /// never holding more than one buffer's worth of raw sequence in memory - unlike the
        /// string-based constructor, which needs the whole sequence materialized as one
        /// <see cref="string"/> first. Peak memory is O(windowCount) plus a small fixed buffer, not
        /// O(sequence length), so this scales past .NET's ~1.07 billion character string ceiling
        /// (see CLAUDE.md); the source never has to exist as one in-memory string at all.
        ///
        /// A trailing partial window (fewer than <paramref name="windowSize"/> bases left) is
        /// dropped, matching the string-based constructor's boundary exactly - confirmed
        /// bit-for-bit identical output for the same sequence
        /// (ZCurveWindowedProfileTests.FromStreamMatchesStringConstructor).
        /// </summary>
        public static ZCurveWindowedProfile FromStream(TextReader reader, int windowSize = DefaultWindowSize)
        {
            if (reader is null) throw new ArgumentNullException(nameof(reader));
            if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));

            var windows = new List<(double X, double Y, double Z)>();
            long sumX = 0, sumY = 0, sumZ = 0;
            int inWindow = 0;

            var buffer = new char[65536];
            int read;
            while ((read = reader.Read(buffer, 0, buffer.Length)) > 0)
            {
                for (int i = 0; i < read; i++)
                {
                    var (x, y, z) = ZCurveProfile.Axes(buffer[i]);
                    sumX += x;
                    sumY += y;
                    sumZ += z;
                    inWindow++;

                    if (inWindow == windowSize)
                    {
                        windows.Add((sumX / (double)windowSize, sumY / (double)windowSize, sumZ / (double)windowSize));
                        sumX = sumY = sumZ = 0;
                        inWindow = 0;
                    }
                }
            }

            return new ZCurveWindowedProfile(null, windowSize, windows);
        }

        private static double Distance((double X, double Y, double Z) a, (double X, double Y, double Z) b)
        {
            var dx = a.X - b.X;
            var dy = a.Y - b.Y;
            var dz = a.Z - b.Z;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        /// <summary>
        /// Slides <paramref name="query"/> (which must have no more windows than
        /// <paramref name="reference"/>) along the reference and returns the best average
        /// per-window score together with the window offset where that alignment starts, mirroring
        /// <see cref="WindowedProfile.BestAlignment"/>. Returns (NaN, -1) when the query is empty
        /// or longer than the reference.
        /// </summary>
        public static (double Score, int Offset) BestAlignment(ZCurveWindowedProfile query, ZCurveWindowedProfile reference)
        {
            if (query.WindowSize != reference.WindowSize)
                throw new ArgumentException("Profiles must share the same window size.");

            var q = query.Windows;
            var r = reference.Windows;
            if (q.Count == 0 || r.Count < q.Count) return (double.NaN, -1);

            double best = double.NegativeInfinity;
            int bestOffset = -1;
            for (int offset = 0; offset + q.Count <= r.Count; offset++)
            {
                double sumDistance = 0;
                for (int j = 0; j < q.Count; j++)
                    sumDistance += Distance(q[j], r[offset + j]);

                // Same "1 - average distance" convention as ZCurveProfile.Compare.
                double score = 1 - sumDistance / q.Count;
                if (score > best)
                {
                    best = score;
                    bestOffset = offset;
                }
            }

            return (best, bestOffset);
        }

        /// <summary>
        /// Symmetric similarity: treats the profile with fewer windows as the query and localizes
        /// it within the other, returning the best-alignment score (0 if either is empty).
        /// </summary>
        public static double Compare(ZCurveWindowedProfile x, ZCurveWindowedProfile y)
        {
            var (shorter, longer) = x.Windows.Count <= y.Windows.Count ? (x, y) : (y, x);
            var (score, _) = BestAlignment(shorter, longer);
            return double.IsNaN(score) ? 0 : score;
        }
    }
}
