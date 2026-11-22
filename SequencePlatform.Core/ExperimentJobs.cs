using Npgsql;

namespace SequencePlatform.Core
{
    public enum ExperimentJobStatus { Queued, Running, Done, Failed }

    /// <summary>A rerun of one experiment on this database, as queued by an admin.</summary>
    /// <param name="RunId">The server run it published; null until done, or once a later rerun replaced that run.</param>
    public sealed record ExperimentJob(Guid Id, string Experiment, string RequestedBy, ExperimentJobStatus Status, string? Error, int Attempts,
        DateTime CreatedAt, DateTime? StartedAt, DateTime? CompletedAt, long? RunId);

    /// <summary>
    /// experiment_job I/O (migration 0024). PostgreSQL is the queue, as for searches (<see cref="SearchLog"/>): <see cref="EnqueueAsync"/>
    /// notifies <see cref="Channel"/>, <see cref="ClaimAsync"/> takes the oldest queued job with <c>for update skip locked</c>, and
    /// <see cref="RequeueStuckAsync"/> returns jobs whose worker stopped mid-run. One queued or running job per experiment at a time.
    /// </summary>
    public static class ExperimentJobs
    {
        public const string Channel = "experiment_queued";
        // An experiment runs for minutes and uses every core: one that brought its worker down twice isn't tried a third time.
        public const int MaxAttempts = 2;

        private const string uniqueViolation = "23505";
        private const string jobSelect =
            "select j.public_id, j.experiment, u.name, j.status, j.error, j.attempts, j.created_at, j.started_at, j.completed_at, j.run_id " +
            "from experiment_job j join app_user u on u.id = j.requested_by";

        /// <returns>The new job's id; null when this experiment already has a queued or running job.</returns>
        public static async Task<Guid?> EnqueueAsync(NpgsqlDataSource dataSource, string experiment, long requestedBy,
            CancellationToken cancellationToken = default)
        {
            await using var connection = await dataSource.OpenConnectionAsync(cancellationToken);
            Guid id;
            try
            {
                await using var insert = new NpgsqlCommand("insert into experiment_job (experiment, requested_by) values ($1, $2) returning public_id", connection);
                insert.Parameters.AddWithValue(experiment);
                insert.Parameters.AddWithValue(requestedBy);
                id = (Guid)(await insert.ExecuteScalarAsync(cancellationToken))!;
            }
            catch (PostgresException e) when (e.SqlState == uniqueViolation && e.ConstraintName == "experiment_job_active")
            {
                return null;
            }
            await using var notify = new NpgsqlCommand("select pg_notify($1, $2)", connection);
            notify.Parameters.AddWithValue(Channel);
            notify.Parameters.AddWithValue(id.ToString());
            await notify.ExecuteNonQueryAsync(cancellationToken);
            return id;
        }

        /// <summary>Marks the oldest queued job running and counts the attempt.</summary>
        /// <returns>Its id and experiment; null when nothing is queued.</returns>
        public static async Task<(Guid Id, string Experiment)?> ClaimAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "update experiment_job set status = 'running', started_at = now(), attempts = attempts + 1 " +
                "where id = (select id from experiment_job where status = 'queued' order by created_at, id limit 1 for update skip locked) " +
                "returning public_id, experiment");
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            return await reader.ReadAsync(cancellationToken) ? (reader.GetGuid(0), reader.GetString(1)) : null;
        }

        /// <summary>
        /// Jobs running since before <paramref name="startedBefore"/>, whose worker stopped: back to queued, or failed once they have
        /// been claimed <see cref="MaxAttempts"/> times.
        /// </summary>
        public static async Task<(int Requeued, int Failed)> RequeueStuckAsync(NpgsqlDataSource dataSource, DateTimeOffset startedBefore,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "with stuck as (select id, attempts from experiment_job where status = 'running' and started_at < $1 for update skip locked), " +
                "failed as (update experiment_job j set status = 'failed', completed_at = now(), " +
                "error = 'The experiment stopped its worker ' || s.attempts || ' times; it was not retried.' " +
                "from stuck s where j.id = s.id and s.attempts >= $2 returning j.id), " +
                "requeued as (update experiment_job j set status = 'queued', started_at = null from stuck s where j.id = s.id and s.attempts < $2 returning j.id) " +
                "select (select count(*)::integer from requeued), (select count(*)::integer from failed)");
            command.Parameters.AddWithValue(startedBefore);
            command.Parameters.AddWithValue((short)MaxAttempts);
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            await reader.ReadAsync(cancellationToken);
            return (reader.GetInt32(0), reader.GetInt32(1));
        }

        /// <param name="runId">The server run the job published.</param>
        public static async Task<bool> FinishAsync(NpgsqlDataSource dataSource, Guid id, long runId, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "update experiment_job set status = 'done', completed_at = now(), run_id = $2 where public_id = $1 and status = 'running'");
            command.Parameters.AddWithValue(id);
            command.Parameters.AddWithValue(runId);
            return await command.ExecuteNonQueryAsync(cancellationToken) == 1;
        }

        public static async Task<bool> FailAsync(NpgsqlDataSource dataSource, Guid id, string error, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "update experiment_job set status = 'failed', error = $2, completed_at = now() where public_id = $1 and status in ('queued', 'running')");
            command.Parameters.AddWithValue(id);
            command.Parameters.AddWithValue(error);
            return await command.ExecuteNonQueryAsync(cancellationToken) == 1;
        }

        /// <returns>null when no job has this id.</returns>
        public static async Task<ExperimentJob?> GetAsync(NpgsqlDataSource dataSource, Guid id, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(jobSelect + " where j.public_id = $1");
            command.Parameters.AddWithValue(id);
            return (await ReadJobsAsync(command, cancellationToken)).SingleOrDefault();
        }

        /// <summary>An experiment's latest jobs, newest first.</summary>
        public static async Task<IReadOnlyList<ExperimentJob>> ListAsync(NpgsqlDataSource dataSource, string experiment, int limit,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(jobSelect + " where j.experiment = $1 order by j.created_at desc, j.id desc limit $2");
            command.Parameters.AddWithValue(experiment);
            command.Parameters.AddWithValue(limit);
            return await ReadJobsAsync(command, cancellationToken);
        }

        private static async Task<IReadOnlyList<ExperimentJob>> ReadJobsAsync(NpgsqlCommand command, CancellationToken cancellationToken)
        {
            var jobs = new List<ExperimentJob>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                jobs.Add(new ExperimentJob(reader.GetGuid(0), reader.GetString(1), reader.GetString(2),
                    Enum.Parse<ExperimentJobStatus>(reader.GetString(3), ignoreCase: true), reader.IsDBNull(4) ? null : reader.GetString(4),
                    reader.GetInt16(5), reader.GetDateTime(6), reader.IsDBNull(7) ? null : reader.GetDateTime(7),
                    reader.IsDBNull(8) ? null : reader.GetDateTime(8), reader.IsDBNull(9) ? null : reader.GetInt64(9)));
            return jobs;
        }

        /// <summary>Whether a search is queued or running: the experiment worker waits for those, so a rerun never delays a user's search.</summary>
        public static async Task<bool> SearchesPendingAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand("select exists (select 1 from search_run where status in ('queued', 'running'))");
            return (bool)(await command.ExecuteScalarAsync(cancellationToken))!;
        }
    }
}
