using System.Globalization;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>score_null I/O (migration 0007): one row per <see cref="ScoreNullCell"/>, and a CSV export of the same.</summary>
    public static class ScoreNullStore
    {
        public const string Genome = "genome";
        public const string Fragment = "fragment";
        public const string Permuted = "permuted";

        /// <summary>Replaces every cell of this method, version and source in one transaction.</summary>
        public static async Task SaveAsync(NpgsqlDataSource dataSource, string method, string methodVersion, string source,
            IReadOnlyCollection<ScoreNullCell> cells)
        {
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var delete = new NpgsqlCommand(
                "delete from score_null where method = $1 and method_version = $2 and source = $3", connection, transaction))
            {
                delete.Parameters.AddWithValue(method);
                delete.Parameters.AddWithValue(methodVersion);
                delete.Parameters.AddWithValue(source);
                await delete.ExecuteNonQueryAsync();
            }

            foreach (var cell in cells)
            {
                await using var insert = new NpgsqlCommand(
                    "insert into score_null (method, method_version, source, window_size, windows_bucket, offsets_bucket, scores) " +
                    "values ($1, $2, $3, $4, $5, $6, $7)", connection, transaction);
                var p = insert.Parameters;
                p.AddWithValue(method);
                p.AddWithValue(methodVersion);
                p.AddWithValue(source);
                p.AddWithValue(cell.WindowSize);
                p.AddWithValue((short)cell.WindowsBucket);
                p.AddWithValue((short)cell.OffsetsBucket);
                p.AddWithValue(cell.SortedScores);
                await insert.ExecuteNonQueryAsync();
            }
            await transaction.CommitAsync();
        }

        /// <returns>null when no cell of this method, version and source has been built.</returns>
        public static async Task<ScoreNull?> LoadAsync(NpgsqlDataSource dataSource, string method, string methodVersion, string source,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select window_size, windows_bucket, offsets_bucket, scores from score_null " +
                "where method = $1 and method_version = $2 and source = $3");
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(methodVersion);
            command.Parameters.AddWithValue(source);

            var cells = new List<ScoreNullCell>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                cells.Add(new ScoreNullCell(reader.GetInt32(0), reader.GetInt16(1), reader.GetInt16(2), reader.GetFieldValue<float[]>(3)));
            return cells.Count == 0 ? null : new ScoreNull(cells);
        }

        /// <summary>
        /// Source, WindowSize, WindowsBucket, OffsetsBucket, SampleCount, then the sorted scores space-separated: what the
        /// analysis tests read without a database (results/S36_ScoreNull.csv).
        /// </summary>
        public static void WriteCsv(TextWriter writer, IEnumerable<(string Source, ScoreNullCell Cell)> cells)
        {
            writer.WriteLine("Source, WindowSize, WindowsBucket, OffsetsBucket, SampleCount, Scores");
            foreach (var (source, c) in cells)
                writer.WriteLine(string.Join(", ", source, c.WindowSize, c.WindowsBucket, c.OffsetsBucket, c.SortedScores.Length,
                    string.Join(' ', c.SortedScores.Select(s => s.ToString("R", CultureInfo.InvariantCulture)))));
        }
    }
}
