using System.Security.Cryptography;
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
        public const int MaxBucket = 25;
        private const double bucketRatio = 1.1;
        private const short formatVersion = 2; // the ZCGP header version ToBinary writes

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
        /// <returns>false, with nothing written, when a bucket's payload is over <see cref="ProfileLevels.MaxInlineBytes"/>.</returns>
        public static async Task<bool> WriteAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId, ZCurveGeneProfile profile)
        {
            var buckets = Enumerable.Range(0, profile.Genes.Count).GroupBy(g => BucketOf(profile.Genes[g].Residues))
                .ToDictionary(b => b.Key, b => profile.Select(b));
            if (!buckets.ContainsKey(LongBucket)) buckets[LongBucket] = profile.Select(Array.Empty<int>());
            var payloads = buckets.ToDictionary(b => b.Key, b => b.Value.ToBinary());
            if (payloads.Values.Any(p => p.Length > ProfileLevels.MaxInlineBytes)) return false;

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
                p.AddWithValue(formatVersion);
                p.AddWithValue(part.Genes.Count);
                p.AddWithValue(part.WindowCount);
                p.AddWithValue(payloads[bucket]);
                p.AddWithValue(payloads[bucket].Length);
                p.AddWithValue(SHA256.HashData(payloads[bucket]));
                await insert.ExecuteNonQueryAsync();
            }
            return true;
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
