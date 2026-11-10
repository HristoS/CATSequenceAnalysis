using System.Security.Cryptography;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// A record's gene profile (S7.S, migration 0017): <see cref="ZCurveGeneProfile.ToBinary"/> rows in gene_profile, one per length bucket
    /// (<see cref="BucketOf"/>), so a protein search reads only the buckets it compares with (<see cref="BucketsFor"/>). Before 0017 the
    /// whole profile was one profile_level row of method 'zcurve-cds'; 0018 removed the last of them.
    /// </summary>
    public static class GeneProfiles
    {
        public const int LongBucket = 0;
        private const string over = " over codon positions 1 and 2 in a window";

        /// <summary>The 9 window features in stored order (3 * order + axis), with what each measures.</summary>
        public static readonly IReadOnlyList<(string Name, string Meaning)> Features =
        [
            ("x mean", "purine (A, G) minus pyrimidine (C, T)" + over),
            ("y mean", "amino (A, C) minus keto (G, T)" + over),
            ("z mean", "strong (C, G) minus weak (A, T)" + over),
            ("x trend", "how x changes from the window's 5' end to its 3' end (Legendre P1, scaled)"),
            ("y trend", "how y changes from the window's 5' end to its 3' end (Legendre P1, scaled)"),
            ("z trend", "how z changes from the window's 5' end to its 3' end (Legendre P1, scaled)"),
            ("x curvature", "whether x bulges or dips in the window's middle (Legendre P2, scaled)"),
            ("y curvature", "whether y bulges or dips in the window's middle (Legendre P2, scaled)"),
            ("z curvature", "whether z bulges or dips in the window's middle (Legendre P2, scaled)"),
        ];

        public const int MaxBucket = 25;
        // Gene rows have no object-storage path, unlike genome levels (ProfileLevels.MaxInlineBytes, 1 MB). The largest vertebrate
        // proteomes put 1-1.5 MB of long genes in one bucket (2026-11-09, S7.T), so gene rows get their own limit (your choice).
        public const int MaxRowBytes = 4 << 20;
        private const double bucketRatio = 1.1;

        /// <summary>The bucket of a gene of this many residues: <see cref="LongBucket"/> from 500 aa, else 1 + floor(log1.1(L / 50)).</summary>
        public static int BucketOf(int residues) => residues >= ZCurveGeneProfile.LongResidues ? LongBucket
            : 1 + (int)Math.Floor(Math.Log(residues / (double)ZCurveGeneProfile.MinResidues, bucketRatio));

        /// <summary>The buckets holding every gene within <paramref name="band"/> of a protein of this length, as a search compares it.</summary>
        public static int[] BucketsFor(int residues, double band)
        {
            if (residues >= ZCurveGeneProfile.LongResidues) return new[] { LongBucket };
            int low = BucketOf(Math.Max(ZCurveGeneProfile.MinResidues, (int)Math.Floor(residues * (1 - band))));
            int high = BucketOf(Math.Min(ZCurveGeneProfile.LongResidues - 1, (int)Math.Ceiling(residues * (1 + band))));
            return Enumerable.Range(low, high - low + 1).ToArray();
        }

        /// <summary>
        /// Replaces the record's gene profile inside the caller's transaction: a row per bucket that has genes, and always a
        /// <see cref="LongBucket"/> row (it marks the record done).
        /// </summary>
        /// <returns>false, with nothing written, when a bucket's payload is over <see cref="MaxRowBytes"/>.</returns>
        public static async Task<bool> WriteAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId, ZCurveGeneProfile profile)
        {
            var buckets = Enumerable.Range(0, profile.Genes.Count).GroupBy(g => BucketOf(profile.Genes[g].Residues))
                .ToDictionary(b => b.Key, b => profile.Select(b));
            if (!buckets.ContainsKey(LongBucket)) buckets[LongBucket] = profile.Select(Array.Empty<int>());
            var payloads = buckets.ToDictionary(b => b.Key, b => b.Value.ToBinary());
            if (payloads.Values.Any(p => p.Length > MaxRowBytes)) return false;

            await using (var delete = new NpgsqlCommand("delete from gene_profile where sequence_id = $1", connection, transaction))
            {
                delete.Parameters.AddWithValue(sequenceId);
                await delete.ExecuteNonQueryAsync();
            }
            foreach (var (bucket, part) in buckets.OrderBy(b => b.Key))
            {
                await using var insert = new NpgsqlCommand(
                    "insert into gene_profile (sequence_id, bucket, format_version, gene_count, window_count, payload, byte_size, sha256) " +
                    "values ($1, $2, $3, $4, $5, $6, $7, $8)", connection, transaction);
                var p = insert.Parameters;
                p.AddWithValue(sequenceId);
                p.AddWithValue((short)bucket);
                p.AddWithValue((short)part.BinaryFormatVersion);
                p.AddWithValue(part.Genes.Count);
                p.AddWithValue(part.WindowCount);
                p.AddWithValue(payloads[bucket]);
                p.AddWithValue(payloads[bucket].Length);
                p.AddWithValue(SHA256.HashData(payloads[bucket]));
                await insert.ExecuteNonQueryAsync();
            }
            return true;
        }

        /// <summary>
        /// Records a whole-assembly record's gene layer in its annotations, as ncbi.gene_layer: where the genes came from and, under
        /// "sequences", the names of the sequences its genes lie on, in the order <see cref="ZCurveGeneProfile.Gene.Chromosome"/> numbers them.
        /// </summary>
        public static async Task AnnotateLayerAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId, string layerJson)
        {
            await using var update = new NpgsqlCommand(
                "update sequence_record set annotations = annotations || " +
                "jsonb_build_object('ncbi', coalesce(annotations -> 'ncbi', '{}'::jsonb) || jsonb_build_object('gene_layer', $2::jsonb)) where id = $1",
                connection, transaction);
            update.Parameters.AddWithValue(sequenceId);
            update.Parameters.AddWithValue(layerJson);
            if (await update.ExecuteNonQueryAsync() != 1) throw new InvalidOperationException($"No record {sequenceId}.");
        }

        /// <summary>For each of these records that has them, the names of the sequences its genes lie on (see <see cref="AnnotateLayerAsync"/>).</summary>
        public static async Task<Dictionary<long, string[]>> SequenceNamesAsync(NpgsqlDataSource dataSource, IReadOnlyCollection<long> ids,
            CancellationToken cancellationToken = default)
        {
            var names = new Dictionary<long, string[]>();
            if (ids.Count == 0) return names;
            await using var command = dataSource.CreateCommand(
                "select id, annotations #>> '{ncbi,gene_layer,sequences}' from sequence_record " +
                "where id = any($1) and annotations #> '{ncbi,gene_layer,sequences}' is not null");
            command.Parameters.AddWithValue(ids.ToArray());
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                names[reader.GetInt64(0)] = JsonSerializer.Deserialize<string[]>(reader.GetString(1)) ?? [];
            return names;
        }

        /// <summary>
        /// Every gene row checked against its stored SHA-256, the check <see cref="ProteinSearch"/> no longer makes on each search
        /// (it cost as much as decoding, results/SearchBreakdown_Before.csv). <see cref="ReadAsync"/> still checks the rows it reads.
        /// </summary>
        /// <returns>Rows checked, and the (record, bucket) of every row that doesn't match.</returns>
        public static async Task<(long Rows, List<(long SequenceId, int Bucket)> Mismatches)> VerifyAsync(NpgsqlDataSource dataSource,
            Action<long>? progress = null, CancellationToken cancellationToken = default)
        {
            var mismatches = new List<(long, int)>();
            long rows = 0;
            await using var command = dataSource.CreateCommand("select sequence_id, bucket, sha256, payload from gene_profile");
            await using var reader = await command.ExecuteReaderAsync(System.Data.CommandBehavior.SequentialAccess, cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
            {
                long sequenceId = reader.GetInt64(0);
                int bucket = reader.GetInt16(1);
                var sha256 = reader.GetFieldValue<byte[]>(2);
                if (!SHA256.HashData(reader.GetFieldValue<byte[]>(3)).AsSpan().SequenceEqual(sha256)) mismatches.Add((sequenceId, bucket));
                if (++rows % 100_000 == 0) progress?.Invoke(rows);
            }
            return (rows, mismatches);
        }

        /// <summary>The record's whole gene profile, its buckets combined; null when it has none.</summary>
        public static async Task<ZCurveGeneProfile?> ReadAsync(NpgsqlDataSource dataSource, long sequenceId)
        {
            await using var command = dataSource.CreateCommand("select payload, sha256 from gene_profile where sequence_id = $1 order by bucket");
            command.Parameters.AddWithValue(sequenceId);
            var parts = new List<ZCurveGeneProfile>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync())
            {
                var payload = reader.GetFieldValue<byte[]>(0);
                if (!SHA256.HashData(payload).AsSpan().SequenceEqual(reader.GetFieldValue<byte[]>(1)))
                    throw new InvalidDataException($"Gene profile of record {sequenceId} doesn't match its SHA-256.");
                parts.Add(ZCurveGeneProfile.FromBinary(payload));
            }
            return parts.Count == 0 ? null : ZCurveGeneProfile.Combine(parts);
        }
    }
}
