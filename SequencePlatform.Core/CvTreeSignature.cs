using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// A genome's amino-acid CVTree vector (Hao's composition vector of the 4-mers of its annotated proteins, with the 3- and 2-mers
    /// derived from the stored 4-mer counts, <see cref="ProteinKmers.Shorter"/>), feature-hashed to <see cref="Dimensions"/> and stored
    /// as unit length in 16-bit fixed point: the relatives signal for genomes with a gene layer (results/README.md §4: same-genus share
    /// of the 10 nearest +13 points over log rho* for bacteria and archaea, +16 for viruses). Hashing to 8,192 keeps the full vector's
    /// accuracy (results/CvTree_Compression.csv at commit f5cb655). Fixed point rather than half floats because .NET has no vectorized
    /// half conversion, and a unit vector's components are within [-1, 1], so by Cauchy-Schwarz no int32 partial sum of a dot product
    /// can pass <see cref="scale"/>^2 and overflow.
    /// </summary>
    public static class CvTreeSignature
    {
        public const int Dimensions = 8192;
        public const string Description = "amino-acid CVTree (4-mers of the annotated proteins, hashed to 8,192 dimensions), (1 - cosine) / 2";
        public const int K = 4;
        private const int batch = 1000;
        private const double scale = short.MaxValue;

        /// <summary>The vector of a record's stored 4-mer counts; null when they hold no word at all.</summary>
        public static short[]? FromFourMers(int[] counts)
        {
            var hashed = Hash(Composition(counts));
            double norm = Math.Sqrt(hashed.Sum(x => x * x));
            if (norm == 0) return null;
            var vector = new short[Dimensions];
            for (int i = 0; i < Dimensions; i++) vector[i] = (short)Math.Round(hashed[i] / norm * scale);
            return vector;
        }

        /// <summary>(1 - cosine) / 2: 0 for the same direction, 1 for opposite ones.</summary>
        public static double Distance(short[] a, short[] b) => Math.Clamp((1 - Dot(a, b) / (scale * scale)) / 2, 0, 1);

        /// <summary>
        /// The vectors of these records, from their stored 4-mer counts, built on <paramref name="parallelism"/> threads in batches so the
        /// payloads never sit in memory at once; records without counts are left out.
        /// </summary>
        public static async Task<Dictionary<long, short[]>> ReadAsync(NpgsqlDataSource dataSource, IReadOnlyCollection<long> ids,
            int parallelism, CancellationToken cancellationToken = default)
        {
            var vectors = new Dictionary<long, short[]>();
            foreach (var chunk in ids.Chunk(batch))
            {
                var payloads = new List<(long Id, byte[] Payload)>();
                await using (var command = dataSource.CreateCommand(
                    "select sequence_id, payload from protein_kmer where k = $1 and status = $2 and sequence_id = any($3)"))
                {
                    command.Parameters.AddWithValue((short)K);
                    command.Parameters.AddWithValue(ProteinKmers.Done);
                    command.Parameters.AddWithValue(chunk);
                    await using var reader = await command.ExecuteReaderAsync(cancellationToken);
                    while (await reader.ReadAsync(cancellationToken)) payloads.Add((reader.GetInt64(0), reader.GetFieldValue<byte[]>(1)));
                }
                foreach (var (id, vector) in payloads.AsParallel().WithDegreeOfParallelism(Math.Max(1, parallelism)).WithCancellation(cancellationToken)
                    .Select(p => (p.Id, FromFourMers(ProteinKmers.Decode(p.Payload, K)))).ToList())
                    if (vector is not null) vectors[id] = vector;
            }
            return vectors;
        }

        // AVX2's multiply-add of adjacent pairs (vpmaddwd: 16 products per instruction, summed in pairs into int32) is what keeps a
        // relatives query over 28,000 viruses under 100 ms; the portable path widens first and is the fallback.
        private static long Dot(short[] a, short[] b)
        {
            if (Avx2.IsSupported && a.Length % (2 * Vector256<short>.Count) == 0)
            {
                ref short left = ref MemoryMarshal.GetArrayDataReference(a);
                ref short right = ref MemoryMarshal.GetArrayDataReference(b);
                var first = Vector256<int>.Zero;
                var second = Vector256<int>.Zero;
                for (int j = 0; j < a.Length; j += 2 * Vector256<short>.Count)
                {
                    first = Avx2.Add(first, Avx2.MultiplyAddAdjacent(Vector256.LoadUnsafe(ref left, (nuint)j), Vector256.LoadUnsafe(ref right, (nuint)j)));
                    int next = j + Vector256<short>.Count;
                    second = Avx2.Add(second, Avx2.MultiplyAddAdjacent(Vector256.LoadUnsafe(ref left, (nuint)next), Vector256.LoadUnsafe(ref right, (nuint)next)));
                }
                return Vector256.Sum(Avx2.Add(first, second));
            }

            var sum = Vector<int>.Zero;
            int i = 0;
            for (; i <= a.Length - Vector<short>.Count; i += Vector<short>.Count)
            {
                Vector.Widen(new Vector<short>(a, i), out var aLow, out var aHigh);
                Vector.Widen(new Vector<short>(b, i), out var bLow, out var bHigh);
                sum += aLow * bLow + aHigh * bHigh;
            }
            long total = Vector.Sum(sum);
            for (; i < a.Length; i++) total += a[i] * b[i];
            return total;
        }

        // a(w) = (f(w) - f0(w)) / f0(w), f0(w1..w4) = f(w1..w3) f(w2..w4) / f(w2w3), with the shorter frequencies from the 4-mer counts.
        private static double[] Composition(int[] words)
        {
            var shorterCounts = ProteinKmers.Shorter(words);
            var shortestCounts = ProteinKmers.Shorter(shorterCounts);
            double total = words.Sum(w => (double)w), shorterTotal = shorterCounts.Sum(w => (double)w), shortestTotal = shortestCounts.Sum(w => (double)w);
            var vector = new double[words.Length];
            if (total == 0) return vector;
            int letters = AminoAcidKmers.Alphabet.Length, suffixSpace = shorterCounts.Length, middleSpace = shortestCounts.Length;
            for (int w = 0; w < words.Length; w++)
            {
                int prefix = w / letters, suffix = w % suffixSpace, middle = prefix % middleSpace;
                double shortest = shortestCounts[middle] / shortestTotal;
                double expected = shortest > 0 ? (shorterCounts[prefix] / shorterTotal) * (shorterCounts[suffix] / shorterTotal) / shortest : 0;
                vector[w] = expected > 0 ? (words[w] / total - expected) / expected : 0;
            }
            return vector;
        }

        // Each word to one bucket with a fixed pseudo-random sign (splitmix64 finalizer), as cvtree_compression_check.cs at commit f5cb655.
        private static double[] Hash(double[] vector)
        {
            var hashed = new double[Dimensions];
            for (int w = 0; w < vector.Length; w++)
            {
                ulong h = (ulong)w * 0x9E3779B97F4A7C15UL;
                h ^= h >> 29;
                h *= 0xBF58476D1CE4E5B9UL;
                h ^= h >> 32;
                hashed[(int)(h % Dimensions)] += (h & (1UL << 63)) == 0 ? vector[w] : -vector[w];
            }
            return hashed;
        }
    }
}
