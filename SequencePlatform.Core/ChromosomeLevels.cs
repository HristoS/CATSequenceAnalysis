using System.Data;
using System.Security.Cryptography;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// chromosome_level I/O (migration 0022): a eukaryote's windowed Z-curve per nuclear sequence of at least <see cref="MinLengthBp"/>,
    /// at 500 bp and 32 kbp, as ZCWP version 2 (float16). For drawing a chromosome's path; nothing searches these rows. The 32 kbp level is
    /// the 500 bp one downsampled, so the two agree up to float16 rounding.
    /// </summary>
    public static class ChromosomeLevels
    {
        /// <summary>Shorter sequences (unplaced scaffolds, mostly) aren't kept; how many and how much of the genome is recorded instead.</summary>
        public const long MinLengthBp = 100_000;
        public const int CoarseFactor = 64;
        public static readonly int[] WindowSizes = [RecordStore.FineWindowSize, RecordStore.FineWindowSize * CoarseFactor];
        private const short formatVersion = 2;

        /// <summary>One kept sequence: its FASTA id, length and 500 bp level.</summary>
        public sealed record Chromosome(string Name, long LengthBp, ZCurveWindowedProfile Fine);

        /// <summary>One stored level of one sequence, as listed.</summary>
        public sealed record LevelInfo(int Index, string Name, long LengthBp, int WindowSize, int WindowCount, long ByteSize);

        /// <summary>What a genome's sequences came to, as recorded under annotations' ncbi.chromosome_levels.</summary>
        public sealed record Summary(int Kept, long KeptBp, int LeftOut, long LeftOutBp, string Status = Stored)
        {
            public string ToJson() => System.Text.Json.JsonSerializer.Serialize(new Dictionary<string, object>
            {
                ["status"] = Status, ["sequences"] = Kept, ["sequences_bp"] = KeptBp, ["left_out"] = LeftOut, ["left_out_bp"] = LeftOutBp,
                ["min_length_bp"] = MinLengthBp, ["window_sizes"] = WindowSizes,
            });
        }

        public const string Stored = "stored";
        public const string Mismatch = "mismatch";

        /// <summary>Both levels of each sequence, numbered in list order, inside the caller's transaction.</summary>
        public static async Task InsertAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId,
            IReadOnlyList<Chromosome> chromosomes)
        {
            for (int index = 0; index < chromosomes.Count; index++)
            {
                var c = chromosomes[index];
                foreach (var level in new[] { c.Fine, c.Fine.Downsample(CoarseFactor) })
                {
                    if (level.WindowCount == 0) continue;
                    var payload = level.ToBinary(half: true);
                    await using var insert = new NpgsqlCommand(
                        "insert into chromosome_level (sequence_id, sequence_index, sequence_name, length_bp, method, method_version, window_size, " +
                        "window_count, format_version, payload, byte_size, sha256) values ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12)",
                        connection, transaction);
                    var p = insert.Parameters;
                    p.AddWithValue(sequenceId);
                    p.AddWithValue(index);
                    p.AddWithValue(c.Name);
                    p.AddWithValue(c.LengthBp);
                    p.AddWithValue(RecordStore.Method);
                    p.AddWithValue(RecordStore.MethodVersion);
                    p.AddWithValue(level.WindowSize);
                    p.AddWithValue(level.WindowCount);
                    p.AddWithValue(formatVersion);
                    p.AddWithValue(payload);
                    p.AddWithValue((long)payload.Length);
                    p.AddWithValue(SHA256.HashData(payload));
                    await insert.ExecuteNonQueryAsync();
                }
            }
        }

        /// <summary>Records the genome's <see cref="Summary"/> as annotations' ncbi.chromosome_levels.</summary>
        public static async Task AnnotateAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId, Summary summary)
        {
            await using var update = new NpgsqlCommand(
                "update sequence_record set annotations = annotations || " +
                "jsonb_build_object('ncbi', coalesce(annotations -> 'ncbi', '{}'::jsonb) || jsonb_build_object('chromosome_levels', $2::jsonb)) where id = $1",
                connection, transaction);
            update.Parameters.AddWithValue(sequenceId);
            update.Parameters.AddWithValue(summary.ToJson());
            if (await update.ExecuteNonQueryAsync() != 1) throw new InvalidOperationException($"No record {sequenceId}.");
        }

        /// <summary>Every stored level of one record, by sequence then window size; empty for a record without any.</summary>
        public static async Task<List<LevelInfo>> ListAsync(NpgsqlDataSource dataSource, long sequenceId, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select sequence_index, sequence_name, length_bp, window_size, window_count, byte_size from chromosome_level " +
                "where sequence_id = $1 and method = $2 and method_version = $3 order by sequence_index, window_size");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue(RecordStore.Method);
            command.Parameters.AddWithValue(RecordStore.MethodVersion);
            var levels = new List<LevelInfo>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                levels.Add(new LevelInfo(reader.GetInt32(0), reader.GetString(1), reader.GetInt64(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetInt64(5)));
            return levels;
        }

        /// <returns>null when the record has no such sequence or no level of that window size.</returns>
        public static async Task<ZCurveWindowedProfile?> ReadAsync(NpgsqlDataSource dataSource, long sequenceId, int index, int windowSize,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select payload, sha256 from chromosome_level " +
                "where sequence_id = $1 and sequence_index = $2 and method = $3 and method_version = $4 and window_size = $5");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue(index);
            command.Parameters.AddWithValue(RecordStore.Method);
            command.Parameters.AddWithValue(RecordStore.MethodVersion);
            command.Parameters.AddWithValue(windowSize);
            await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess, cancellationToken);
            if (!await reader.ReadAsync(cancellationToken)) return null;
            var payload = await reader.GetFieldValueAsync<byte[]>(0, cancellationToken);
            var sha256 = await reader.GetFieldValueAsync<byte[]>(1, cancellationToken);
            if (!SHA256.HashData(payload).AsSpan().SequenceEqual(sha256))
                throw new InvalidDataException($"Chromosome level {index} of sequence {sequenceId} doesn't match its stored SHA-256.");
            return ZCurveWindowedProfile.FromBinary(payload);
        }
    }
}
