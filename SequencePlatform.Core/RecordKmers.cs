using Npgsql;
using NpgsqlTypes;

namespace SequencePlatform.Core
{
    /// <summary>
    /// A record's 6-mer counts (migration 0013, D8 widened 2026-11-04), as <see cref="SequenceFeatureReader"/> counts them. New records
    /// get them at ingest with their profile; older ones by re-downloading (Ingest backfill-kmers), which also records why a record
    /// couldn't be counted, so a rerun skips it.
    /// </summary>
    public static class RecordKmers
    {
        public const int K = SequenceFeatureReader.HexamerLength;
        public const string Done = "done";
        public const string Mismatch = "mismatch";
        public const string Missing = "missing";

        /// <param name="counts">Required for <see cref="Done"/>, null otherwise.</param>
        public static async Task WriteAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId, string status, int[]? counts)
        {
            if ((status == Done) != (counts is not null))
                throw new ArgumentException($"Counts go with status '{Done}' and only with it.", nameof(counts));
            if (counts is not null && counts.Length != SequenceFeatureReader.HexamerCount)
                throw new ArgumentException($"Expected {SequenceFeatureReader.HexamerCount} counts, got {counts.Length}.", nameof(counts));

            await using var command = new NpgsqlCommand(
                "insert into record_kmer (sequence_id, k, counts, status) values ($1, $2, $3, $4) " +
                "on conflict (sequence_id, k) do update set counts = excluded.counts, status = excluded.status, checked_at = now()",
                connection, transaction);
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue((short)K);
            command.Parameters.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Array | NpgsqlDbType.Integer, Value = (object?)counts ?? DBNull.Value });
            command.Parameters.AddWithValue(status);
            await command.ExecuteNonQueryAsync();
        }

        /// <summary>The record's counts, or null when it has none stored.</summary>
        public static async Task<int[]?> ReadAsync(NpgsqlDataSource dataSource, long sequenceId)
        {
            await using var command = dataSource.CreateCommand("select counts from record_kmer where sequence_id = $1 and k = $2 and status = $3");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue((short)K);
            command.Parameters.AddWithValue(Done);
            return await command.ExecuteScalarAsync() is int[] counts ? counts : null;
        }
    }
}
