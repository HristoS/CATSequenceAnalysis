using System;
using System.Buffers.Binary;
using System.Linq;
using System.Runtime.Intrinsics;

namespace SequenceAnalyses
{
    /// <summary>
    /// Z-curve profile: maps each base onto three independent biochemical dichotomies - purine
    /// (A,G) vs pyrimidine (C,T), amino (A,C) vs keto (G,T), and strong/3-H-bond (G,C) vs weak/
    /// 2-H-bond (A,T) - then walks the sequence once, summing each axis. Each axis's running
    /// average (sum / bases-seen-so-far) is sampled at <see cref="CheckpointCount"/> evenly-spaced
    /// positions, giving a fixed-size trajectory - a coarse sketch of the walk's shape - in one
    /// real, shared 3D frame derived directly from the bases' own chemistry, rather than
    /// trilaterated against an arbitrary periodic benchmark string the way the original CAT
    /// method was.
    ///
    /// With CheckpointCount = 1 (the default), the single checkpoint is the whole-sequence
    /// endpoint - mathematically just a linear recombination of base composition (#A,#C,#G,#T), so
    /// it cannot separate two sequences that share composition but differ in base order (verified:
    /// shuffling a sequence leaves this endpoint exactly unchanged). That's the same ceiling
    /// the original CAT method's single-point compression hit. Raising CheckpointCount
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
        /// Endpoint-only similarity as "1 - Euclidean distance": 1 for identical endpoints, falling off
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

        // Structure-of-arrays, and float rather than double. Both choices exist for BestAlignment,
        // which indexes these millions of times per call when comparing whole genomes of very
        // different lengths:
        //  - Separate per-axis arrays mean 8 consecutive windows' X values are contiguous in memory,
        //    so a vector register can be filled by one direct load. An array of (X,Y,Z) tuples
        //    interleaves the axes, forcing a gather - assembling a register from scattered scalars,
        //    which measured *slower* than plain scalar code when tried that way.
        //  - float halves the bytes per window (12 vs 24), so twice as many windows fit in cache on
        //    the streaming reads this loop does, and one Vector256 holds 8 lanes instead of 4.
        //    Precision is ample: every stored value is a mean of windowSize terms each in {-1,0,+1},
        //    so it lies in [-1,1] with only 2*windowSize+1 distinct possible values (1001 at the
        //    default 500) - about 10 bits, against float's 24-bit mantissa.
        private readonly float[] xs;
        private readonly float[] ys;
        private readonly float[] zs;

        // Materialized on first access and cached: the public shape stays
        // IReadOnlyList<(double,double,double)> so every existing caller of Windows.Count/Windows[i]
        // is unaffected by the layout change above. BestAlignment/Compare never touch this - they
        // read the three arrays directly - so the hot path never pays for it.
        private IReadOnlyList<(double X, double Y, double Z)>? windowsView;

        /// <summary>
        /// Number of windows, without materializing <see cref="Windows"/> - prefer this over
        /// <c>Windows.Count</c>, which builds the whole widened view just to read a length (about
        /// 149 MB for a genome-scale profile).
        /// </summary>
        public int WindowCount => xs.Length;

        /// <summary>Ordered per-window local Z-curve axis averages; window i covers bases [i*WindowSize, (i+1)*WindowSize).</summary>
        public IReadOnlyList<(double X, double Y, double Z)> Windows
        {
            get
            {
                if (windowsView is null)
                {
                    var view = new (double X, double Y, double Z)[xs.Length];
                    for (int i = 0; i < view.Length; i++)
                        view[i] = (xs[i], ys[i], zs[i]);
                    windowsView = view;
                }
                return windowsView;
            }
        }

        public ZCurveWindowedProfile(string dnaString, int windowSize = DefaultWindowSize)
        {
            if (dnaString is null) throw new ArgumentNullException(nameof(dnaString));
            if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize));

            DnaString = dnaString;
            WindowSize = windowSize;

            int count = dnaString.Length / windowSize;
            var xs = new float[count];
            var ys = new float[count];
            var zs = new float[count];
            for (int w = 0, start = 0; w < count; w++, start += windowSize)
            {
                long sumX = 0, sumY = 0, sumZ = 0;
                for (int i = start; i < start + windowSize; i++)
                {
                    var (x, y, z) = ZCurveProfile.Axes(dnaString[i]);
                    sumX += x;
                    sumY += y;
                    sumZ += z;
                }
                // Divide in double, store in float: the division itself stays exact-as-before, only
                // the stored result is narrowed.
                xs[w] = (float)(sumX / (double)windowSize);
                ys[w] = (float)(sumY / (double)windowSize);
                zs[w] = (float)(sumZ / (double)windowSize);
            }
            this.xs = xs;
            this.ys = ys;
            this.zs = zs;
        }

        private ZCurveWindowedProfile(string? dnaString, int windowSize, float[] xs, float[] ys, float[] zs)
        {
            DnaString = dnaString;
            WindowSize = windowSize;
            this.xs = xs;
            this.ys = ys;
            this.zs = zs;
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

            // Always copies into the profile's own storage, never aliases the caller's collection -
            // which could otherwise be mutated later out from under the profile.
            int count = windows.Count;
            var xs = new float[count];
            var ys = new float[count];
            var zs = new float[count];
            for (int i = 0; i < count; i++)
            {
                var (x, y, z) = windows[i];
                xs[i] = (float)x;
                ys[i] = (float)y;
                zs[i] = (float)z;
            }
            return new ZCurveWindowedProfile(dnaString, windowSize, xs, ys, zs);
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

            var xs = new List<float>();
            var ys = new List<float>();
            var zs = new List<float>();
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
                        // Same divide-in-double-store-in-float as the string constructor, so the two
                        // paths stay bit-for-bit identical (FromStreamMatchesStringConstructor).
                        xs.Add((float)(sumX / (double)windowSize));
                        ys.Add((float)(sumY / (double)windowSize));
                        zs.Add((float)(sumZ / (double)windowSize));
                        sumX = sumY = sumZ = 0;
                        inWindow = 0;
                    }
                }
            }

            return new ZCurveWindowedProfile(null, windowSize, xs.ToArray(), ys.ToArray(), zs.ToArray());
        }

        /// <summary>
        /// Coarser profile derived from this one without the original sequence: each output window
        /// is the mean of <paramref name="factor"/> consecutive windows, so the result has
        /// WindowSize * factor. Every window holds exactly WindowSize bases (the trailing partial
        /// window is always dropped), so a mean of equal-sized means is the mean over all their
        /// bases - this reproduces <c>new ZCurveWindowedProfile(dna, WindowSize * factor)</c> up to
        /// float rounding, including the window count (a trailing partial group is dropped too).
        ///
        /// <paramref name="phase"/> skips that many leading windows first. Coarse windows can only
        /// start at multiples of <paramref name="factor"/> fine windows, so a fragment that begins
        /// mid-group lines up with none of them; building each phase 0..factor-1 of a short query
        /// covers every alignment for coarse-to-fine search.
        /// </summary>
        private const uint binaryMagic = 0x5057435A; // "ZCWP" read as a little-endian uint32
        private const ushort binaryFormatVersion = 1;
        private const int binaryHeaderSize = 16;

        /// <summary>
        /// Serializes the windows (never the DNA string) in the storage format of
        /// docs/database-decision.md §4.3: "ZCWP" | uint16 version | uint16 reserved | int32
        /// windowSize | int32 windowCount | float32 xs[n] | ys[n] | zs[n], all little-endian. The
        /// per-axis blocks match the in-memory layout, so a reader can fetch just windows [a, b) of
        /// a stored object with three byte-range reads.
        /// </summary>
        public byte[] ToBinary()
        {
            int n = xs.Length;
            var bytes = new byte[binaryHeaderSize + 12 * n];
            var span = bytes.AsSpan();
            BinaryPrimitives.WriteUInt32LittleEndian(span, binaryMagic);
            BinaryPrimitives.WriteUInt16LittleEndian(span[4..], binaryFormatVersion);
            BinaryPrimitives.WriteInt32LittleEndian(span[8..], WindowSize);
            BinaryPrimitives.WriteInt32LittleEndian(span[12..], n);
            WriteAxis(span.Slice(binaryHeaderSize, 4 * n), xs);
            WriteAxis(span.Slice(binaryHeaderSize + 4 * n, 4 * n), ys);
            WriteAxis(span.Slice(binaryHeaderSize + 8 * n, 4 * n), zs);
            return bytes;
        }

        /// <summary>Inverse of <see cref="ToBinary"/>; the result carries no DNA string.</summary>
        public static ZCurveWindowedProfile FromBinary(ReadOnlySpan<byte> bytes)
        {
            if (bytes.Length < binaryHeaderSize || BinaryPrimitives.ReadUInt32LittleEndian(bytes) != binaryMagic)
                throw new FormatException("Not a ZCWP profile.");
            var version = BinaryPrimitives.ReadUInt16LittleEndian(bytes[4..]);
            if (version != binaryFormatVersion)
                throw new FormatException($"Unsupported ZCWP format version {version}.");

            int windowSize = BinaryPrimitives.ReadInt32LittleEndian(bytes[8..]);
            int n = BinaryPrimitives.ReadInt32LittleEndian(bytes[12..]);
            if (windowSize < 1 || n < 0 || bytes.Length != binaryHeaderSize + 12L * n)
                throw new FormatException("ZCWP header does not match the payload length.");

            return new ZCurveWindowedProfile(null, windowSize,
                ReadAxis(bytes.Slice(binaryHeaderSize, 4 * n)),
                ReadAxis(bytes.Slice(binaryHeaderSize + 4 * n, 4 * n)),
                ReadAxis(bytes.Slice(binaryHeaderSize + 8 * n, 4 * n)));
        }

        private static void WriteAxis(Span<byte> destination, float[] values)
        {
            for (int i = 0; i < values.Length; i++)
                BinaryPrimitives.WriteSingleLittleEndian(destination[(4 * i)..], values[i]);
        }

        private static float[] ReadAxis(ReadOnlySpan<byte> source)
        {
            var values = new float[source.Length / 4];
            for (int i = 0; i < values.Length; i++)
                values[i] = BinaryPrimitives.ReadSingleLittleEndian(source[(4 * i)..]);
            return values;
        }

        public ZCurveWindowedProfile Downsample(int factor, int phase = 0)
        {
            if (factor < 1) throw new ArgumentOutOfRangeException(nameof(factor));
            if (phase < 0 || phase >= factor) throw new ArgumentOutOfRangeException(nameof(phase));

            int count = Math.Max(0, (xs.Length - phase) / factor);
            var coarseXs = new float[count];
            var coarseYs = new float[count];
            var coarseZs = new float[count];
            for (int w = 0, start = phase; w < count; w++, start += factor)
            {
                double sumX = 0, sumY = 0, sumZ = 0;
                for (int i = start; i < start + factor; i++)
                {
                    sumX += xs[i];
                    sumY += ys[i];
                    sumZ += zs[i];
                }
                coarseXs[w] = (float)(sumX / factor);
                coarseYs[w] = (float)(sumY / factor);
                coarseZs[w] = (float)(sumZ / factor);
            }

            // A phase-shifted profile no longer starts at base 0 of DnaString, so it doesn't carry it.
            return new ZCurveWindowedProfile(phase == 0 ? DnaString : null, WindowSize * factor, coarseXs, coarseYs, coarseZs);
        }

        // Euclidean distance, sqrt included. Four cheaper alternatives (squared Euclidean, Chebyshev,
        // Manhattan, and one sqrt of the summed distance) were each measured against the real-genome
        // corpus and all reverted: they either compress the near-duplicate-vs-unrelated spacing the
        // method depends on, or - in Chebyshev's case - run slower despite fewer operations. All of
        // them correlate highly with this one, which is precisely why correlation alone was not a
        // sufficient test. See DistanceMetricComparisonBenchmark for the harness that measured them.
        private static float Distance(float ax, float ay, float az, float bx, float by, float bz)
        {
            var dx = ax - bx;
            var dy = ay - by;
            var dz = az - bz;
            return MathF.Sqrt(dx * dx + dy * dy + dz * dz);
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

            // The measured hot path - up to O(windowCount^2) per call for two whole genomes of very
            // different lengths. Reads the three per-axis arrays directly rather than the public
            // Windows property, and takes one of two vectorized shapes depending on how many windows
            // the query has; both keep every load contiguous. Roughly 3-4x faster than the scalar,
            // interleaved-layout implementation it replaced, on both whole-genome comparison and
            // short-fragment localization.
            //
            // Getting here took one failed attempt worth recording: keeping the interleaved
            // (X,Y,Z)[] layout and vectorizing only the sqrt - building each vector with
            // Vector256.Create from four separately-computed scalars - measured materially *slower*
            // than scalar code, because that gather plus a per-batch horizontal reduce cost more than
            // the batched sqrt saved. Splitting the axes into separate arrays is what removes the
            // gather; without that layout change, vectorizing here is a net loss.
            //
            // A second one: replacing Vector256.Sqrt with an rsqrt estimate + one Newton step (+ a
            // zero mask, since rsqrt(0) is inf) measured 0.72-0.96x - slower on every shape, 0.725x
            // over all 435 real-genome pairs (i7-10510U, Release) - despite negligible score drift.
            // The sqrt runs on its own divider unit alongside the sub/mul/add work; the replacement
            // puts ~7 extra ops on the FP ports that are already the bottleneck.
            var qx = query.xs;
            var qy = query.ys;
            var qz = query.zs;
            var rx = reference.xs;
            var ry = reference.ys;
            var rz = reference.zs;
            int m = qx.Length;
            int n = rx.Length;
            if (m == 0 || n < m) return (double.NaN, -1);

            int width = Vector256<float>.Count;
            double best = double.NegativeInfinity;
            int bestOffset = -1;

            if (Vector256.IsHardwareAccelerated && m >= width)
            {
                // Wide query: vectorize across j. Each lane handles a different window of the same
                // offset, so the per-offset total needs one horizontal reduce at the end.
                int vectorizedLength = m - (m % width);
                for (int offset = 0; offset + m <= n; offset++)
                {
                    var acc = Vector256<float>.Zero;
                    int j = 0;
                    for (; j < vectorizedLength; j += width)
                    {
                        var dx = Vector256.LoadUnsafe(ref qx[j]) - Vector256.LoadUnsafe(ref rx[offset + j]);
                        var dy = Vector256.LoadUnsafe(ref qy[j]) - Vector256.LoadUnsafe(ref ry[offset + j]);
                        var dz = Vector256.LoadUnsafe(ref qz[j]) - Vector256.LoadUnsafe(ref rz[offset + j]);
                        acc += Vector256.Sqrt(dx * dx + dy * dy + dz * dz);
                    }
                    double sumDistance = Vector256.Sum(acc);
                    for (; j < m; j++)
                        sumDistance += Distance(qx[j], qy[j], qz[j], rx[offset + j], ry[offset + j], rz[offset + j]);

                    // Same "1 - average distance" convention as ZCurveProfile.Compare.
                    double score = 1 - sumDistance / m;
                    if (score > best)
                    {
                        best = score;
                        bestOffset = offset;
                    }
                }
            }
            else
            {
                // Narrow query (fewer windows than one vector holds - e.g. localizing a 2,000 bp
                // fragment, 4 windows at the default size): vectorizing across j would leave most
                // lanes idle, so vectorize across *offsets* instead. Each lane accumulates a
                // different offset's running total, while the query window is broadcast to all lanes.
                // The loads stay contiguous for the same reason as the wide path - offsets o..o+7 read
                // rx[o+j..o+7+j], which is adjacent memory - so this avoids the gather that made an
                // earlier attempt slower, and avoids the cache penalty of walking three separate
                // arrays scalar-wise (three cache lines per window instead of one, which measured
                // materially slower here than the old interleaved layout - this path is not a
                // micro-optimization, without it the narrow-query case regresses).
                int lastOffset = n - m;
                int offsetCount = lastOffset + 1;
                int vectorizedOffsets = Vector256.IsHardwareAccelerated ? offsetCount - (offsetCount % width) : 0;

                int offset = 0;
                for (; offset < vectorizedOffsets; offset += width)
                {
                    var acc = Vector256<float>.Zero;
                    for (int j = 0; j < m; j++)
                    {
                        var dx = Vector256.Create(qx[j]) - Vector256.LoadUnsafe(ref rx[offset + j]);
                        var dy = Vector256.Create(qy[j]) - Vector256.LoadUnsafe(ref ry[offset + j]);
                        var dz = Vector256.Create(qz[j]) - Vector256.LoadUnsafe(ref rz[offset + j]);
                        acc += Vector256.Sqrt(dx * dx + dy * dy + dz * dz);
                    }

                    // Lane k holds offset+k's total. Scanned in increasing k so ties resolve to the
                    // lowest offset, matching the scalar loop's "first best wins" behavior. Each lane
                    // also accumulated j in order, so a lane's total is bit-identical to what the
                    // scalar path computes for that offset.
                    for (int k = 0; k < width; k++)
                    {
                        double score = 1 - (double)acc.GetElement(k) / m;
                        if (score > best)
                        {
                            best = score;
                            bestOffset = offset + k;
                        }
                    }
                }

                for (; offset <= lastOffset; offset++)
                {
                    double sumDistance = 0;
                    for (int j = 0; j < m; j++)
                        sumDistance += Distance(qx[j], qy[j], qz[j], rx[offset + j], ry[offset + j], rz[offset + j]);

                    double score = 1 - sumDistance / m;
                    if (score > best)
                    {
                        best = score;
                        bestOffset = offset;
                    }
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
            var (shorter, longer) = x.xs.Length <= y.xs.Length ? (x, y) : (y, x);
            var (score, _) = BestAlignment(shorter, longer);
            return double.IsNaN(score) ? 0 : score;
        }
    }

    /// <summary>
    /// A fine <see cref="ZCurveWindowedProfile"/> plus coarser levels derived from it with
    /// <see cref="ZCurveWindowedProfile.Downsample"/> (default 500 bp, 2 kbp, 8 kbp, 32 kbp), so one
    /// stored profile can be compared cheaply against both a 5 kbp virus and a 5 Mbp bacterium.
    ///
    /// Ranking quality depends on how many windows the *shorter* sequence keeps, not on window size
    /// itself (MultiResolutionBenchmark.WindowCountSweep): bacterial retrieval AUC showed no loss
    /// down to about 4 windows, viral same-genus AUC lost about 0.005 at ~9 windows and 0.03 at ~2,
    /// and zero windows is a hard failure - <see cref="ZCurveWindowedProfile.Compare"/> then returns
    /// 0, indistinguishable from "unrelated". <see cref="Compare"/> therefore picks the coarsest level
    /// that still leaves the shorter sequence <see cref="DefaultMinWindows"/> windows, and returns null
    /// rather than 0 when the pair isn't comparable at any level.
    /// </summary>
    public sealed class ZCurveProfileLadder
    {
        public const int DefaultMinWindows = 8;

        public static readonly IReadOnlyList<int> DefaultFactors = new[] { 1, 4, 16, 64 };

        private readonly int[] factors;
        private readonly ZCurveWindowedProfile[] levels;

        public IReadOnlyList<int> Factors => factors;

        public ZCurveWindowedProfile Fine => levels[0];

        /// <param name="factors">Ascending downsample factors; the first must be 1 (the fine profile itself).</param>
        public ZCurveProfileLadder(ZCurveWindowedProfile fine, IReadOnlyList<int>? factors = null)
        {
            if (fine is null) throw new ArgumentNullException(nameof(fine));
            var chosen = (factors ?? DefaultFactors).ToArray();
            if (chosen.Length == 0 || chosen[0] != 1)
                throw new ArgumentException("The first factor must be 1.", nameof(factors));
            for (int i = 1; i < chosen.Length; i++)
                if (chosen[i] <= chosen[i - 1])
                    throw new ArgumentException("Factors must be strictly ascending.", nameof(factors));

            this.factors = chosen;
            levels = chosen.Select(f => f == 1 ? fine : fine.Downsample(f)).ToArray();
        }

        public ZCurveWindowedProfile Level(int factor)
        {
            int index = Array.IndexOf(factors, factor);
            if (index < 0) throw new ArgumentOutOfRangeException(nameof(factor), $"This ladder has no level with factor {factor}.");
            return levels[index];
        }

        /// <summary>
        /// The coarsest factor at which the shorter of the two sequences still has at least
        /// <paramref name="minWindows"/> windows; the finest factor if no level reaches that but the
        /// fine level has at least one window; null if even the fine level leaves it with none.
        /// </summary>
        public static int? SelectFactor(ZCurveProfileLadder x, ZCurveProfileLadder y, int minWindows = DefaultMinWindows)
        {
            if (minWindows < 1) throw new ArgumentOutOfRangeException(nameof(minWindows));
            RequireCompatible(x, y);

            for (int i = x.factors.Length - 1; i >= 0; i--)
                if (Math.Min(x.levels[i].WindowCount, y.levels[i].WindowCount) >= minWindows)
                    return x.factors[i];

            return Math.Min(x.Fine.WindowCount, y.Fine.WindowCount) >= 1 ? 1 : null;
        }

        private static void RequireCompatible(ZCurveProfileLadder x, ZCurveProfileLadder y)
        {
            if (x.Fine.WindowSize != y.Fine.WindowSize || !x.factors.AsSpan().SequenceEqual(y.factors))
                throw new ArgumentException("Ladders must share the same fine window size and factors.");
        }

        /// <summary>
        /// <see cref="ZCurveWindowedProfile.Compare"/> at the level <see cref="SelectFactor"/> picks,
        /// or null when the pair isn't comparable at all (a sequence shorter than one fine window).
        /// </summary>
        public static double? Compare(ZCurveProfileLadder x, ZCurveProfileLadder y, int minWindows = DefaultMinWindows)
        {
            var factor = SelectFactor(x, y, minWindows);
            return factor is int f ? ZCurveWindowedProfile.Compare(x.Level(f), y.Level(f)) : null;
        }

        /// <summary>
        /// The one level a search scores every candidate at: the coarsest factor at which the query
        /// alone keeps at least <paramref name="minWindows"/> windows; the finest factor if no level
        /// reaches that but the fine level has a window; null if the query has no fine window.
        /// </summary>
        public static int? SelectQueryFactor(ZCurveProfileLadder query, int minWindows = DefaultMinWindows)
        {
            if (minWindows < 1) throw new ArgumentOutOfRangeException(nameof(minWindows));

            for (int i = query.factors.Length - 1; i >= 0; i--)
                if (query.levels[i].WindowCount >= minWindows)
                    return query.factors[i];

            return query.Fine.WindowCount >= 1 ? 1 : null;
        }

        /// <summary>
        /// The level for a search that only compares the query with candidates at least
        /// <paramref name="lengthRatio"/> times its length: the coarsest factor at which
        /// lengthRatio × the query's window count still reaches <paramref name="minWindows"/>, so every such
        /// candidate keeps that many windows there too (it has at least floor(lengthRatio × query windows)).
        /// The finest factor if no level reaches that but the query has a fine window; null otherwise.
        /// Comparing whole sequences of very different lengths measures containment, not similarity: the
        /// shorter one slides along the longer and keeps its best stretch (results/S34_LengthGate.csv).
        /// </summary>
        public static int? SelectBandFactor(ZCurveProfileLadder query, double lengthRatio, int minWindows = DefaultMinWindows)
        {
            if (lengthRatio <= 0 || lengthRatio > 1) throw new ArgumentOutOfRangeException(nameof(lengthRatio));
            if (minWindows < 1) throw new ArgumentOutOfRangeException(nameof(minWindows));

            for (int i = query.factors.Length - 1; i >= 0; i--)
                if (lengthRatio * query.levels[i].WindowCount >= minWindows)
                    return query.factors[i];

            return query.Fine.WindowCount >= 1 ? 1 : null;
        }

        /// <summary>
        /// Scores <paramref name="candidate"/> for a search by <paramref name="query"/>. Candidates
        /// that keep at least <paramref name="minWindows"/> windows at <see cref="SelectQueryFactor"/>'s
        /// level (or, when the query itself falls back below that, at least as many as the query) are
        /// <see cref="SearchList.Ranked"/> at that level, so one query's ranking never mixes levels.
        /// Shorter candidates go to <see cref="SearchList.ShortList"/>, scored at the per-pair level
        /// <see cref="SelectFactor"/> picks; those scores can come from different levels, so the short
        /// list is reported separately and isn't a ranking on one scale.
        /// </summary>
        public static SearchMatch CompareForQuery(ZCurveProfileLadder query, ZCurveProfileLadder candidate, int minWindows = DefaultMinWindows)
        {
            RequireCompatible(query, candidate);
            if (SelectQueryFactor(query, minWindows) is not int queryFactor)
                return new SearchMatch(SearchList.NotComparable, null, null);

            var queryLevel = query.Level(queryFactor);
            var candidateLevel = candidate.Level(queryFactor);
            if (candidateLevel.WindowCount >= Math.Min(minWindows, queryLevel.WindowCount))
                return new SearchMatch(SearchList.Ranked, queryFactor, ZCurveWindowedProfile.Compare(queryLevel, candidateLevel));

            if (SelectFactor(query, candidate, minWindows) is not int pairFactor)
                return new SearchMatch(SearchList.NotComparable, null, null);

            return new SearchMatch(SearchList.ShortList, pairFactor,
                ZCurveWindowedProfile.Compare(query.Level(pairFactor), candidate.Level(pairFactor)));
        }
    }

    public enum SearchList { Ranked, ShortList, NotComparable }

    /// <summary>Factor and Score are null exactly when List is <see cref="SearchList.NotComparable"/>.</summary>
    public readonly record struct SearchMatch(SearchList List, int? Factor, double? Score);
}
