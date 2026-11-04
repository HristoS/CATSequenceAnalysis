using Npgsql;

namespace SequencePlatform.Core
{
    /// <summary>A P9 experiment as published: what it asked, how it came out, and its summary numbers as JSON.</summary>
    /// <param name="Verdict">pass, fail, or context (an experiment with no pass rule).</param>
    /// <param name="Snapshot">The snapshot database it was computed on.</param>
    public sealed record ExperimentRun(long Id, string Name, string Title, string Question, string Verdict, string Snapshot, string Script,
        string SummaryJson, DateTime PublishedAt, int Findings);

    /// <summary>What an experiment found about one record.</summary>
    /// <param name="Kind">What sort of finding (an experiment may have several), such as "cpg" or "outlier".</param>
    public sealed record ExperimentFinding(long SequenceId, string Kind, double? Score, string PayloadJson);

    /// <summary>A finding together with the experiment it belongs to, for a record's own page.</summary>
    public sealed record RecordFinding(string Experiment, string Title, string Verdict, string Kind, double? Score, string PayloadJson);

    /// <summary>
    /// experiment_run and experiment_finding I/O (migration 0012). Experiments are computed on a snapshot of the catalogue; record ids
    /// are the live database's, because the snapshot is a dump of it. A finding about a record that has since gone, or that isn't public,
    /// is skipped when published.
    /// </summary>
    public static class ExperimentStore
    {
        public const string Pass = "pass";
        public const string Fail = "fail";
        public const string Context = "context";

        private const string runSelect =
            "select r.id, r.name, r.title, r.question, r.verdict, r.snapshot, r.script, r.summary::text, r.published_at, " +
            "(select count(*) from experiment_finding f where f.run_id = r.id)::int from experiment_run r";

        /// <summary>Replaces the run of this name, findings included, in one transaction.</summary>
        /// <returns>The findings stored, and those skipped because their record isn't a public record of this database.</returns>
        public static async Task<(int Stored, int Skipped)> PublishAsync(NpgsqlDataSource dataSource, string name, string title, string question,
            string verdict, string snapshot, string script, string summaryJson, IReadOnlyCollection<ExperimentFinding> findings)
        {
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var delete = new NpgsqlCommand("delete from experiment_run where name = $1", connection, transaction))
            {
                delete.Parameters.AddWithValue(name);
                await delete.ExecuteNonQueryAsync();
            }

            long runId;
            await using (var insert = new NpgsqlCommand(
                "insert into experiment_run (name, title, question, verdict, snapshot, script, summary) values ($1, $2, $3, $4, $5, $6, $7::jsonb) returning id",
                connection, transaction))
            {
                foreach (var value in new[] { name, title, question, verdict, snapshot, script, summaryJson }) insert.Parameters.AddWithValue(value);
                runId = (long)(await insert.ExecuteScalarAsync())!;
            }

            int stored;
            await using (var insert = new NpgsqlCommand(
                "insert into experiment_finding (run_id, sequence_id, kind, score, payload) " +
                "select $1, f.sequence_id, f.kind, f.score, f.payload::jsonb from unnest($2::bigint[], $3::text[], $4::float8[], $5::text[]) " +
                "  as f (sequence_id, kind, score, payload) " +
                "where exists (select 1 from sequence_record r join sequence_source s on s.id = r.source_id where r.id = f.sequence_id and s.owner_id is null)",
                connection, transaction))
            {
                var p = insert.Parameters;
                p.AddWithValue(runId);
                p.AddWithValue(findings.Select(f => f.SequenceId).ToArray());
                p.AddWithValue(findings.Select(f => f.Kind).ToArray());
                p.AddWithValue(findings.Select(f => f.Score).ToArray());
                p.AddWithValue(findings.Select(f => f.PayloadJson).ToArray());
                stored = await insert.ExecuteNonQueryAsync();
            }
            await transaction.CommitAsync();
            return (stored, findings.Count - stored);
        }

        public static async Task<IReadOnlyList<ExperimentRun>> ListAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(runSelect + " order by r.published_at desc");
            return await ReadRunsAsync(command, cancellationToken);
        }

        /// <returns>null when no run of this name has been published.</returns>
        public static async Task<ExperimentRun?> GetAsync(NpgsqlDataSource dataSource, string name, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(runSelect + " where r.name = $1");
            command.Parameters.AddWithValue(name);
            return (await ReadRunsAsync(command, cancellationToken)).SingleOrDefault();
        }

        /// <summary>A page of a run's findings, optionally of one kind, by score (descending when asked) then record id.</summary>
        /// <returns>The page and how many findings match in all.</returns>
        public static async Task<(IReadOnlyList<ExperimentFinding> Items, long Total)> FindingsAsync(NpgsqlDataSource dataSource, long runId,
            string? kind, bool descending, int offset, int limit, CancellationToken cancellationToken = default)
        {
            const string where = "where run_id = $1 and ($2::text is null or kind = $2)";
            var order = descending ? "score desc nulls last" : "score asc nulls last";
            NpgsqlCommand Command(string sql)
            {
                var command = dataSource.CreateCommand(sql);
                command.Parameters.AddWithValue(runId);
                command.Parameters.Add(new NpgsqlParameter { Value = (object?)kind ?? DBNull.Value, NpgsqlDbType = NpgsqlTypes.NpgsqlDbType.Text });
                return command;
            }

            long total;
            await using (var count = Command($"select count(*) from experiment_finding {where}"))
                total = (long)(await count.ExecuteScalarAsync(cancellationToken))!;

            var items = new List<ExperimentFinding>();
            await using var page = Command($"select sequence_id, kind, score, payload::text from experiment_finding {where} order by {order}, sequence_id offset $3 limit $4");
            page.Parameters.AddWithValue(offset);
            page.Parameters.AddWithValue(limit);
            await using var reader = await page.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                items.Add(new ExperimentFinding(reader.GetInt64(0), reader.GetString(1), reader.IsDBNull(2) ? null : reader.GetDouble(2), reader.GetString(3)));
            return (items, total);
        }

        /// <summary>Every finding about one record, newest experiment first. Visibility is the caller's to check.</summary>
        public static async Task<IReadOnlyList<RecordFinding>> ForRecordAsync(NpgsqlDataSource dataSource, long sequenceId,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select r.name, r.title, r.verdict, f.kind, f.score, f.payload::text from experiment_finding f " +
                "join experiment_run r on r.id = f.run_id where f.sequence_id = $1 order by r.published_at desc, f.kind");
            command.Parameters.AddWithValue(sequenceId);
            var findings = new List<RecordFinding>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                findings.Add(new RecordFinding(reader.GetString(0), reader.GetString(1), reader.GetString(2), reader.GetString(3),
                    reader.IsDBNull(4) ? null : reader.GetDouble(4), reader.GetString(5)));
            return findings;
        }

        private static async Task<IReadOnlyList<ExperimentRun>> ReadRunsAsync(NpgsqlCommand command, CancellationToken cancellationToken)
        {
            var runs = new List<ExperimentRun>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                runs.Add(new ExperimentRun(reader.GetInt64(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetString(4),
                    reader.GetString(5), reader.GetString(6), reader.GetString(7), reader.GetDateTime(8), reader.GetInt32(9)));
            return runs;
        }
    }
}
