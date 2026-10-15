using Npgsql;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// crawl_item reads and writes (migration 0003). Final states are done, duplicate and too_large; a run
    /// never re-requests those. failed keeps its attempt count and last error and is retried only with
    /// --retry-failed, so one broken genome can't make every batch slower.
    /// </summary>
    internal static class CrawlState
    {
        public const string Done = "done";
        public const string Failed = "failed";
        public const string Duplicate = "duplicate";
        public const string TooLarge = "too_large";
        private const int maxErrorLength = 2000;

        public static bool IsFinal(string status) => status is Done or Duplicate or TooLarge;

        // Registers newly listed items as pending and returns every listed item's current status.
        public static async Task<Dictionary<string, string>> RegisterAsync(NpgsqlDataSource dataSource, string provider, string? group,
            IReadOnlyCollection<string> externalIds)
        {
            await using (var insert = dataSource.CreateCommand(
                "insert into crawl_item (provider, external_id, item_group) select $1, unnest($2::text[]), $3 " +
                "on conflict (provider, external_id) do nothing"))
            {
                insert.Parameters.AddWithValue(provider);
                insert.Parameters.AddWithValue(externalIds.ToArray());
                insert.Parameters.AddWithValue((object?)group ?? DBNull.Value);
                await insert.ExecuteNonQueryAsync();
            }

            var statuses = new Dictionary<string, string>(StringComparer.Ordinal);
            await using var select = dataSource.CreateCommand(
                "select external_id, status from crawl_item where provider = $1 and external_id = any($2)");
            select.Parameters.AddWithValue(provider);
            select.Parameters.AddWithValue(externalIds.ToArray());
            await using var reader = await select.ExecuteReaderAsync();
            while (await reader.ReadAsync()) statuses[reader.GetString(0)] = reader.GetString(1);
            return statuses;
        }

        /// <param name="group">null retries failed items of every group.</param>
        public static async Task<List<string>> FailedAsync(NpgsqlDataSource dataSource, string provider, string? group, int limit)
        {
            await using var command = dataSource.CreateCommand(
                "select external_id from crawl_item where provider = $1 and status = 'failed' and ($2::text is null or item_group = $2) " +
                "order by external_id limit $3");
            command.Parameters.AddWithValue(provider);
            command.Parameters.AddWithValue((object?)group ?? DBNull.Value);
            command.Parameters.AddWithValue(limit);
            var ids = new List<string>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) ids.Add(reader.GetString(0));
            return ids;
        }

        // For outcomes that store nothing (duplicate, failed). Each counts as an attempt.
        public static async Task MarkAsync(NpgsqlDataSource dataSource, string provider, string externalId, string status, string? error = null)
        {
            await using var command = dataSource.CreateCommand(
                "update crawl_item set status = $3, attempts = attempts + 1, last_error = $4, updated_at = now() " +
                "where provider = $1 and external_id = $2");
            command.Parameters.AddWithValue(provider);
            command.Parameters.AddWithValue(externalId);
            command.Parameters.AddWithValue(status);
            command.Parameters.AddWithValue(error is null ? DBNull.Value : error.Length <= maxErrorLength ? error : error[..maxErrorLength]);
            await command.ExecuteNonQueryAsync();
        }

        // Inside the genome's store transaction, so the catalogue rows and the done mark commit together.
        public static async Task MarkStoredAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, string provider,
            string externalId, long sourceId, bool tooLarge)
        {
            await using var command = new NpgsqlCommand(
                "update crawl_item set status = $3, attempts = attempts + 1, last_error = null, source_id = $4, updated_at = now() " +
                "where provider = $1 and external_id = $2", connection, transaction);
            command.Parameters.AddWithValue(provider);
            command.Parameters.AddWithValue(externalId);
            command.Parameters.AddWithValue(tooLarge ? TooLarge : Done);
            command.Parameters.AddWithValue(sourceId);
            if (await command.ExecuteNonQueryAsync() != 1)
                throw new InvalidOperationException($"No crawl_item for {provider} {externalId}; register it before ingesting.");
        }
    }
}
