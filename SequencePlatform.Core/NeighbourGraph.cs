using Npgsql;
using NpgsqlTypes;

namespace SequencePlatform.Core
{
    /// <summary>
    /// neighbour_list (migration 0008, P3 step S3.4b): every record's own find-similar top K at one length ratio,
    /// excluding its own assembly, found by the same <see cref="SimilaritySearch.SearchRecordAsync"/> a user would run.
    /// <see cref="SimilaritySearch"/> uses it to re-rank a query's top hits by shared neighbours.
    /// </summary>
    public static class NeighbourGraph
    {
        /// <summary>Searches every record with a stored fine level, then replaces this ratio's rows in one transaction.</summary>
        /// <returns>The number of records listed.</returns>
        public static async Task<int> BuildAsync(NpgsqlDataSource dataSource, double lengthRatio, int neighbours, Action<string> log,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion, CancellationToken cancellationToken = default)
        {
            var ids = new List<long>();
            await using (var command = dataSource.CreateCommand(
                // Public records only (D3): a private upload neither gets a list nor, through the search's own filter, joins one.
                "select r.id from sequence_record r join sequence_source s on s.id = r.source_id " +
                "where s.owner_id is null and exists (select 1 from profile_level p where p.sequence_id = r.id " +
                "and p.method = $1 and p.method_version = $2 and p.window_size = $3) order by r.id"))
            {
                command.Parameters.AddWithValue(method);
                command.Parameters.AddWithValue(methodVersion);
                command.Parameters.AddWithValue(RecordStore.FineWindowSize);
                await using var reader = await command.ExecuteReaderAsync(cancellationToken);
                while (await reader.ReadAsync(cancellationToken)) ids.Add(reader.GetInt64(0));
            }

            var search = new SimilaritySearch(dataSource, method: method, methodVersion: methodVersion);
            await search.LoadAllAsync(cancellationToken);   // concurrent searches below then only read the caches

            var lists = new System.Collections.Concurrent.ConcurrentDictionary<long, long[]>();
            int done = 0;
            await Parallel.ForEachAsync(ids, new ParallelOptions { MaxDegreeOfParallelism = 4, CancellationToken = cancellationToken },
                async (id, token) =>
                {
                    var result = await search.SearchRecordAsync(id, lengthRatio: lengthRatio, top: neighbours, useNeighbours: false,
                        cancellationToken: token);
                    lists[id] = result.Ranked.Select(h => h.SequenceId).ToArray();
                    if (Interlocked.Increment(ref done) % 2_000 == 0) log($"{done:N0} of {ids.Count:N0} records searched");
                });

            await SaveAsync(dataSource, method, methodVersion, lengthRatio, lists, cancellationToken);
            return lists.Count;
        }

        private static async Task SaveAsync(NpgsqlDataSource dataSource, string method, string methodVersion, double lengthRatio,
            IReadOnlyDictionary<long, long[]> lists, CancellationToken cancellationToken)
        {
            await using var connection = await dataSource.OpenConnectionAsync(cancellationToken);
            await using var transaction = await connection.BeginTransactionAsync(cancellationToken);
            await using (var delete = new NpgsqlCommand(
                "delete from neighbour_list where method = $1 and method_version = $2 and length_ratio = $3", connection, transaction))
            {
                delete.Parameters.AddWithValue(method);
                delete.Parameters.AddWithValue(methodVersion);
                delete.Parameters.AddWithValue((float)lengthRatio);
                await delete.ExecuteNonQueryAsync(cancellationToken);
            }

            await using (var import = await connection.BeginBinaryImportAsync(
                "copy neighbour_list (method, method_version, length_ratio, sequence_id, neighbour_ids) from stdin (format binary)", cancellationToken))
            {
                foreach (var (sequenceId, neighbourIds) in lists)
                {
                    await import.StartRowAsync(cancellationToken);
                    await import.WriteAsync(method, NpgsqlDbType.Text, cancellationToken);
                    await import.WriteAsync(methodVersion, NpgsqlDbType.Text, cancellationToken);
                    await import.WriteAsync((float)lengthRatio, NpgsqlDbType.Real, cancellationToken);
                    await import.WriteAsync(sequenceId, NpgsqlDbType.Bigint, cancellationToken);
                    await import.WriteAsync(neighbourIds, NpgsqlDbType.Array | NpgsqlDbType.Bigint, cancellationToken);
                }
                await import.CompleteAsync(cancellationToken);
            }
            await transaction.CommitAsync(cancellationToken);
        }

        /// <returns>Each listed record's neighbour ids, best first; empty when this ratio's graph hasn't been built.</returns>
        public static async Task<IReadOnlyDictionary<long, long[]>> LoadAsync(NpgsqlDataSource dataSource, double lengthRatio,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select sequence_id, neighbour_ids from neighbour_list where method = $1 and method_version = $2 and length_ratio = $3");
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(methodVersion);
            command.Parameters.AddWithValue((float)lengthRatio);

            var lists = new Dictionary<long, long[]>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                lists[reader.GetInt64(0)] = reader.GetFieldValue<long[]>(1);
            return lists;
        }
    }
}
