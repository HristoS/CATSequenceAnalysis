using System.Data;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// profile_level I/O (docs/database-decision.md rev 2, section 5.1): one row per ladder level, the
    /// ZCWP bytes of <see cref="ZCurveWindowedProfile.ToBinary"/> inline up to 1 MB. Larger levels belong
    /// in object storage, which doesn't exist yet, so they're refused on write and on read.
    /// </summary>
    public static class ProfileLevels
    {
        public const int MaxInlineBytes = 1 << 20;
        private const short formatVersion = 1; // the ZCWP header version ToBinary writes

        /// <summary>
        /// Inserts one level inside the caller's transaction, so a genome's records and levels commit together.
        /// </summary>
        /// <returns>false, with nothing written, when the level is over <see cref="MaxInlineBytes"/>.</returns>
        public static async Task<bool> InsertAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId,
            string method, string methodVersion, int factor, ZCurveWindowedProfile level)
        {
            if (factor < 1 || level.WindowSize % factor != 0)
                throw new ArgumentException($"Window size {level.WindowSize} isn't a multiple of factor {factor}.", nameof(factor));

            var payload = level.ToBinary();
            if (payload.Length > MaxInlineBytes) return false;

            await using var insert = new NpgsqlCommand(
                "insert into profile_level (sequence_id, method, method_version, fine_window_size, factor, window_count, " +
                "format_version, payload, byte_size, sha256) values ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10)",
                connection, transaction);
            var p = insert.Parameters;
            p.AddWithValue(sequenceId);
            p.AddWithValue(method);
            p.AddWithValue(methodVersion);
            p.AddWithValue(level.WindowSize / factor);
            p.AddWithValue(factor);
            p.AddWithValue(level.WindowCount);
            p.AddWithValue(formatVersion);
            p.AddWithValue(payload);
            p.AddWithValue((long)payload.Length);
            p.AddWithValue(SHA256.HashData(payload));
            await insert.ExecuteNonQueryAsync();
            return true;
        }

        /// <returns>null when the sequence has no level of that window size.</returns>
        public static async Task<ZCurveWindowedProfile?> ReadAsync(NpgsqlDataSource dataSource, long sequenceId,
            string method, string methodVersion, int windowSize)
        {
            await using var command = dataSource.CreateCommand(
                "select sequence_id, payload, object_key, sha256 from profile_level " +
                "where sequence_id = $1 and method = $2 and method_version = $3 and window_size = $4");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(methodVersion);
            command.Parameters.AddWithValue(windowSize);

            await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess);
            return await reader.ReadAsync() ? (await ReadLevelAsync(reader)).Level : null;
        }

        /// <summary>
        /// Every level of one method, version and window size, in sequence id order: the worker's cache
        /// load (profile_level_load_idx). Rows are read one at a time, never buffered as a whole.
        /// </summary>
        public static async IAsyncEnumerable<(long SequenceId, ZCurveWindowedProfile Level)> StreamAsync(
            NpgsqlDataSource dataSource, string method, string methodVersion, int windowSize,
            [EnumeratorCancellation] CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select sequence_id, payload, object_key, sha256 from profile_level " +
                "where method = $1 and method_version = $2 and window_size = $3 order by sequence_id");
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(methodVersion);
            command.Parameters.AddWithValue(windowSize);

            await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess, cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                yield return await ReadLevelAsync(reader, cancellationToken);
        }

        /// <summary>
        /// The levels of one method, version and window size whose profile_level.id is above
        /// <paramref name="afterLevelId"/>, in id order: what <see cref="ProfileCache"/> adds on a refresh.
        /// </summary>
        public static async IAsyncEnumerable<(long LevelId, long SequenceId, ZCurveWindowedProfile Level)> StreamAfterAsync(
            NpgsqlDataSource dataSource, string method, string methodVersion, int windowSize, long afterLevelId,
            [EnumeratorCancellation] CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select sequence_id, payload, object_key, sha256, id from profile_level " +
                "where method = $1 and method_version = $2 and window_size = $3 and id > $4 order by id");
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(methodVersion);
            command.Parameters.AddWithValue(windowSize);
            command.Parameters.AddWithValue(afterLevelId);

            await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess, cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
            {
                var (sequenceId, level) = await ReadLevelAsync(reader, cancellationToken);
                yield return (reader.GetInt64(4), sequenceId, level);
            }
        }

        /// <summary>How many levels of one method, version and window size have an id up to <paramref name="maxLevelId"/>.</summary>
        public static async Task<long> CountUpToAsync(NpgsqlDataSource dataSource, string method, string methodVersion, int windowSize,
            long maxLevelId, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select count(*) from profile_level where method = $1 and method_version = $2 and window_size = $3 and id <= $4");
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(methodVersion);
            command.Parameters.AddWithValue(windowSize);
            command.Parameters.AddWithValue(maxLevelId);
            return (long)(await command.ExecuteScalarAsync(cancellationToken))!;
        }

        // Columns in select order (SequentialAccess): sequence_id, payload, object_key, sha256.
        private static async Task<(long SequenceId, ZCurveWindowedProfile Level)> ReadLevelAsync(
            NpgsqlDataReader reader, CancellationToken cancellationToken = default)
        {
            var sequenceId = reader.GetInt64(0);
            var payload = await reader.IsDBNullAsync(1, cancellationToken) ? null : await reader.GetFieldValueAsync<byte[]>(1, cancellationToken);
            var objectKey = await reader.IsDBNullAsync(2, cancellationToken) ? null : reader.GetString(2);
            var sha256 = await reader.GetFieldValueAsync<byte[]>(3, cancellationToken);

            if (payload is null)
                throw new NotSupportedException($"Level of sequence {sequenceId} is in object storage ('{objectKey}'), which isn't implemented yet.");
            if (!SHA256.HashData(payload).AsSpan().SequenceEqual(sha256))
                throw new InvalidDataException($"Level of sequence {sequenceId} doesn't match its stored SHA-256.");
            return (sequenceId, ZCurveWindowedProfile.FromBinary(payload));
        }
    }
}
