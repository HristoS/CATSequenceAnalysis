using System.IO.Compression;
using Npgsql;

namespace SequencePlatform.Core
{
    /// <summary>An experiment as published: what it asked, how it came out, and its summary numbers as JSON.</summary>
    /// <param name="Verdict">pass, fail, or context (an experiment with no pass rule).</param>
    /// <param name="Snapshot">The database it was computed on (a frozen snapshot, the live catalogue, or the bundled genomes).</param>
    /// <param name="Details">Which article it backs, where it was computed and how to rerun it (migration 0023).</param>
    /// <param name="Files">Its result files by name, with their uncompressed sizes.</param>
    public sealed record ExperimentRun(long Id, string Name, string Title, string Question, string Verdict, string Snapshot, string Script,
        string SummaryJson, DateTime PublishedAt, int Findings, ExperimentDetails Details, IReadOnlyList<ExperimentFileInfo> Files);

    /// <param name="Article">The article section it backs, such as "A-genome", "A-gene", "B" or "P9"; null when none.</param>
    /// <param name="Origin"><see cref="ExperimentStore.ArticleOrigin"/> for the run an article cites, <see cref="ExperimentStore.ServerOrigin"/>
    /// for a rerun on this database.</param>
    /// <param name="Commit">The commit of the code that computed it.</param>
    /// <param name="Machine">Where it ran, for timings.</param>
    /// <param name="Command">The command that reproduces it from the repository root.</param>
    /// <param name="Runnable">Whether it can be rerun from this database alone; experiments that need DNA never can.</param>
    public sealed record ExperimentDetails(string? Article = null, string Origin = ExperimentStore.ArticleOrigin, string? Commit = null,
        string? Machine = null, string? Command = null, bool Runnable = false);

    /// <summary>A result file to publish: a table of numbers, never sequences.</summary>
    public sealed record ExperimentFile(string Name, byte[] Content);

    public sealed record ExperimentFileInfo(string Name, long ByteSize);

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
        public const string ArticleOrigin = "article";
        public const string ServerOrigin = "server";

        private const string runSelect =
            "select r.id, r.name, r.title, r.question, r.verdict, r.snapshot, r.script, r.summary::text, r.published_at, " +
            "(select count(*) from experiment_finding f where f.run_id = r.id)::int, r.article, r.origin, r.source_commit, r.machine, r.command, " +
            "r.runnable, array(select e.name from experiment_file e where e.run_id = r.id order by e.name), " +
            "array(select e.byte_size from experiment_file e where e.run_id = r.id order by e.name) from experiment_run r";

        /// <summary>Replaces the run of this name and origin, findings and files included, in one transaction.</summary>
        /// <param name="details">Defaults to an article run with nothing else known.</param>
        /// <param name="files">Stored gzip-compressed; names are letters, digits, '_', '.' and '-'.</param>
        /// <returns>The findings stored, and those skipped because their record isn't a public record of this database.</returns>
        public static async Task<(int Stored, int Skipped)> PublishAsync(NpgsqlDataSource dataSource, string name, string title, string question,
            string verdict, string snapshot, string script, string summaryJson, IReadOnlyCollection<ExperimentFinding> findings,
            ExperimentDetails? details = null, IReadOnlyCollection<ExperimentFile>? files = null)
        {
            details ??= new ExperimentDetails();
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var delete = new NpgsqlCommand("delete from experiment_run where name = $1 and origin = $2", connection, transaction))
            {
                delete.Parameters.AddWithValue(name);
                delete.Parameters.AddWithValue(details.Origin);
                await delete.ExecuteNonQueryAsync();
            }

            long runId;
            await using (var insert = new NpgsqlCommand(
                "insert into experiment_run (name, title, question, verdict, snapshot, script, summary, article, origin, source_commit, machine, command, runnable) " +
                "values ($1, $2, $3, $4, $5, $6, $7::jsonb, $8, $9, $10, $11, $12, $13) returning id",
                connection, transaction))
            {
                foreach (var value in new[] { name, title, question, verdict, snapshot, script, summaryJson }) insert.Parameters.AddWithValue(value);
                foreach (var value in new[] { details.Article, details.Origin, details.Commit, details.Machine, details.Command })
                    insert.Parameters.Add(new NpgsqlParameter { Value = (object?)value ?? DBNull.Value, NpgsqlDbType = NpgsqlTypes.NpgsqlDbType.Text });
                insert.Parameters.AddWithValue(details.Runnable);
                runId = (long)(await insert.ExecuteScalarAsync())!;
            }

            foreach (var file in files ?? [])
            {
                await using var insert = new NpgsqlCommand("insert into experiment_file (run_id, name, byte_size, content) values ($1, $2, $3, $4)",
                    connection, transaction);
                insert.Parameters.AddWithValue(runId);
                insert.Parameters.AddWithValue(file.Name);
                insert.Parameters.AddWithValue((long)file.Content.Length);
                insert.Parameters.AddWithValue(Compress(file.Content));
                await insert.ExecuteNonQueryAsync();
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

        /// <summary>Every run of one origin, newest first.</summary>
        public static async Task<IReadOnlyList<ExperimentRun>> ListAsync(NpgsqlDataSource dataSource, string origin = ArticleOrigin,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(runSelect + " where r.origin = $1 order by r.published_at desc");
            command.Parameters.AddWithValue(origin);
            return await ReadRunsAsync(command, cancellationToken);
        }

        /// <returns>null when no run of this name and origin has been published.</returns>
        public static async Task<ExperimentRun?> GetAsync(NpgsqlDataSource dataSource, string name, string origin = ArticleOrigin,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(runSelect + " where r.name = $1 and r.origin = $2");
            command.Parameters.AddWithValue(name);
            command.Parameters.AddWithValue(origin);
            return (await ReadRunsAsync(command, cancellationToken)).SingleOrDefault();
        }

        /// <summary>One result file of a run, uncompressed.</summary>
        /// <returns>null when the run has no file of this name.</returns>
        public static async Task<byte[]?> FileAsync(NpgsqlDataSource dataSource, long runId, string fileName, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand("select content from experiment_file where run_id = $1 and name = $2");
            command.Parameters.AddWithValue(runId);
            command.Parameters.AddWithValue(fileName);
            return await command.ExecuteScalarAsync(cancellationToken) is byte[] content ? Decompress(content) : null;
        }

        private static byte[] Compress(byte[] content)
        {
            using var output = new MemoryStream();
            using (var gzip = new GZipStream(output, CompressionLevel.Optimal, leaveOpen: true)) gzip.Write(content);
            return output.ToArray();
        }

        private static byte[] Decompress(byte[] content)
        {
            using var gzip = new GZipStream(new MemoryStream(content), CompressionMode.Decompress);
            using var output = new MemoryStream();
            gzip.CopyTo(output);
            return output.ToArray();
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
                "join experiment_run r on r.id = f.run_id where f.sequence_id = $1 and r.origin = 'article' order by r.published_at desc, f.kind");
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
            {
                string? Text(int i) => reader.IsDBNull(i) ? null : reader.GetString(i);
                var names = reader.GetFieldValue<string[]>(16);
                var sizes = reader.GetFieldValue<long[]>(17);
                runs.Add(new ExperimentRun(reader.GetInt64(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetString(4),
                    reader.GetString(5), reader.GetString(6), reader.GetString(7), reader.GetDateTime(8), reader.GetInt32(9),
                    new ExperimentDetails(Text(10), reader.GetString(11), Text(12), Text(13), Text(14), reader.GetBoolean(15)),
                    names.Zip(sizes, (n, s) => new ExperimentFileInfo(n, s)).ToList()));
            }
            return runs;
        }
    }
}
