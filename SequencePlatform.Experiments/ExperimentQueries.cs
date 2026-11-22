using System.Text.Json;
using Npgsql;
using NpgsqlTypes;

namespace SequencePlatform.Experiments
{
    /// <summary>One stored query of an experiment: its record, the fields the experiment needs, and its profile. Never a sequence.</summary>
    public sealed record ExperimentQuery(long SequenceId, JsonElement Fields, double[] Profile);

    /// <summary>
    /// experiment_query I/O (migration 0024): the query set of an experiment whose queries came from outside the database, stored once
    /// (by the Ingest's seed-experiment-queries, as the owner) so that it reruns from the database alone.
    /// </summary>
    public static class ExperimentQueries
    {
        /// <summary>Replaces the experiment's whole query set, in order, in one transaction.</summary>
        public static async Task ReplaceAsync(NpgsqlDataSource dataSource, string experiment, IReadOnlyList<ExperimentQuery> queries)
        {
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var delete = new NpgsqlCommand("delete from experiment_query where experiment = $1", connection, transaction))
            {
                delete.Parameters.AddWithValue(experiment);
                await delete.ExecuteNonQueryAsync();
            }
            for (int i = 0; i < queries.Count; i++)
            {
                await using var insert = new NpgsqlCommand(
                    "insert into experiment_query (experiment, position, sequence_id, fields, profile) values ($1, $2, $3, $4, $5)", connection, transaction);
                insert.Parameters.AddWithValue(experiment);
                insert.Parameters.AddWithValue(i);
                insert.Parameters.AddWithValue(queries[i].SequenceId);
                insert.Parameters.Add(new NpgsqlParameter { Value = queries[i].Fields.GetRawText(), NpgsqlDbType = NpgsqlDbType.Jsonb });
                insert.Parameters.AddWithValue(queries[i].Profile);
                await insert.ExecuteNonQueryAsync();
            }
            await transaction.CommitAsync();
        }

        /// <returns>The experiment's queries in their stored order; empty when none are stored.</returns>
        public static async Task<IReadOnlyList<ExperimentQuery>> ReadAsync(NpgsqlDataSource dataSource, string experiment,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select sequence_id, fields::text, profile from experiment_query where experiment = $1 order by position");
            command.Parameters.AddWithValue(experiment);
            var queries = new List<ExperimentQuery>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
            {
                using var fields = JsonDocument.Parse(reader.GetString(1));
                queries.Add(new ExperimentQuery(reader.GetInt64(0), fields.RootElement.Clone(), reader.GetFieldValue<double[]>(2)));
            }
            return queries;
        }
    }
}
