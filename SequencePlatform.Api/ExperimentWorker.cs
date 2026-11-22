using System.Diagnostics;
using Microsoft.Extensions.Options;
using Npgsql;
using SequencePlatform.Core;
using SequencePlatform.Experiments;

namespace SequencePlatform.Api
{
    /// <summary>The experiment worker's settings, from the "ExperimentWorker" configuration section.</summary>
    public sealed class ExperimentWorkerOptions
    {
        /// <summary>Off by default: exactly one process runs it (the AppHost's "experiments", or the API's local http profile).</summary>
        public bool Enabled { get; set; }

        /// <summary>The cores an experiment may use; 0 means all.</summary>
        public int Parallelism { get; set; }

        /// <summary>Start a rerun only while no search is queued or running; tests turn it off, as other fixtures leave searches queued.</summary>
        public bool WaitForSearches { get; set; } = true;

        /// <summary>How often to look for queued jobs when no notification arrives, and to look again while searches are pending.</summary>
        public TimeSpan PollInterval { get; set; } = TimeSpan.FromSeconds(30);

        /// <summary>A job running longer than this lost its worker. protein-homolog scores for 8 min on 8 cores (2026-11-21), so hours.</summary>
        public TimeSpan StuckAfter { get; set; } = TimeSpan.FromHours(3);
    }

    /// <summary>
    /// Runs queued experiment reruns (migration 0024), one at a time, as a third process beside the API and the search worker. A rerun uses
    /// every core for minutes, so it starts only while no search is queued or running; a search queued during a rerun waits for cores but
    /// not for the rerun. It wakes on <see cref="ExperimentJobs.Channel"/> and also polls. It assumes it is the only experiment worker: a job
    /// still running when it starts lost its worker and is requeued, and after <see cref="ExperimentJobs.MaxAttempts"/> claims it is failed.
    /// </summary>
    internal sealed class ExperimentWorker(NpgsqlDataSource dataSource, IEnumerable<IExperiment> experiments, IOptions<ExperimentWorkerOptions> options,
        ILogger<ExperimentWorker> logger) : BackgroundService
    {
        private readonly SemaphoreSlim wake = new(0, 1);

        protected override async Task ExecuteAsync(CancellationToken stoppingToken)
        {
            var settings = options.Value;
            if (!settings.Enabled) return;

            var listener = ListenAsync(settings.PollInterval, stoppingToken);
            bool recovered = false;
            var nextRecovery = DateTimeOffset.MinValue;
            while (!stoppingToken.IsCancellationRequested)
            {
                try
                {
                    if (DateTimeOffset.UtcNow >= nextRecovery)
                    {
                        var (requeued, failed) = await ExperimentJobs.RequeueStuckAsync(dataSource,
                            recovered ? DateTimeOffset.UtcNow - settings.StuckAfter : DateTimeOffset.UtcNow, stoppingToken);
                        if (requeued + failed > 0) logger.LogWarning("Found {Requeued} interrupted experiment(s) requeued and {Failed} failed", requeued, failed);
                        recovered = true;
                        nextRecovery = DateTimeOffset.UtcNow + settings.PollInterval;
                    }
                    while (!stoppingToken.IsCancellationRequested && !(settings.WaitForSearches && await ExperimentJobs.SearchesPendingAsync(dataSource, stoppingToken))
                        && await ExperimentJobs.ClaimAsync(dataSource, stoppingToken) is { } job)
                        await RunAsync(job.Id, job.Experiment, settings, stoppingToken);
                }
                catch (OperationCanceledException) when (stoppingToken.IsCancellationRequested)
                {
                    break;
                }
                catch (Exception e)
                {
                    logger.LogError(e, "The experiment worker's loop failed; retrying in {Interval}", settings.PollInterval);
                }

                try { await wake.WaitAsync(settings.PollInterval, stoppingToken); }
                catch (OperationCanceledException) { break; }
            }
            await listener;
        }

        private async Task RunAsync(Guid id, string name, ExperimentWorkerOptions settings, CancellationToken cancellationToken)
        {
            var watch = Stopwatch.StartNew();
            try
            {
                if (ExperimentRunner.Find(experiments, name) is not { } experiment)
                    throw new ExperimentNotRunnableException($"This build has no database experiment named {name}.");
                var context = new ExperimentContext(settings.Parallelism > 0 ? settings.Parallelism : Environment.ProcessorCount,
                    message => logger.LogInformation("Experiment {Experiment}: {Progress}", name, message));
                var (runId, result) = await ExperimentRunner.RunAndPublishAsync(dataSource, experiment, context, cancellationToken);
                await ExperimentJobs.FinishAsync(dataSource, id, runId, CancellationToken.None);
                logger.LogInformation("Experiment {Experiment} ({JobId}): {Verdict}, server run {RunId}, in {Seconds:F0} s", name, id, result.Verdict, runId,
                    watch.Elapsed.TotalSeconds);
            }
            catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
            {
                throw;
            }
            catch (Exception e)
            {
                logger.LogError(e, "Experiment {Experiment} ({JobId}) failed after {Seconds:F0} s", name, id, watch.Elapsed.TotalSeconds);
                // The runner's own checks (no stored query set, no runnable article run) are written for the admin; anything else stays in the log.
                var error = e is ExperimentNotRunnableException ? e.Message : "The experiment failed; the server log has the details.";
                await ExperimentJobs.FailAsync(dataSource, id, error, CancellationToken.None);
            }
        }

        // A dedicated connection waiting for notifications; reconnects after a failure. Polling covers the gaps.
        private async Task ListenAsync(TimeSpan retryAfter, CancellationToken cancellationToken)
        {
            while (!cancellationToken.IsCancellationRequested)
            {
                try
                {
                    await using var connection = await dataSource.OpenConnectionAsync(cancellationToken);
                    connection.Notification += (_, _) => Wake();
                    await using (var listen = new NpgsqlCommand($"listen {ExperimentJobs.Channel}", connection))
                        await listen.ExecuteNonQueryAsync(cancellationToken);
                    Wake();
                    while (true) await connection.WaitAsync(cancellationToken);
                }
                catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
                {
                    return;
                }
                catch (Exception e)
                {
                    logger.LogWarning(e, "Lost the {Channel} listener; polling until it reconnects", ExperimentJobs.Channel);
                    try { await Task.Delay(retryAfter, cancellationToken); }
                    catch (OperationCanceledException) { return; }
                }
            }
        }

        private void Wake()
        {
            try { if (wake.CurrentCount == 0) wake.Release(); }
            catch (SemaphoreFullException) { }
        }
    }
}
