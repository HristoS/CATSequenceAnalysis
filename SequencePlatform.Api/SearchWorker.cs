using System.Diagnostics;
using Microsoft.Extensions.Options;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>The worker's settings, from the "SearchWorker" configuration section.</summary>
    public sealed class SearchWorkerOptions
    {
        public bool Enabled { get; set; } = true;

        /// <summary>The stored profiles to search; tests point it at their own.</summary>
        public string MethodVersion { get; set; } = RecordStore.MethodVersion;

        /// <summary>
        /// false (the default, and the server's): localization streams its level from the database, 0.85-9.6 s per fragment from 128
        /// down to 4 kbp. true (local runs: the API's http launch profile and the AppHost in run mode): every level in memory, 0.52-1.4 s,
        /// for about 0.8 GB and a load on first use (results/SearchBreakdown_Before.csv).
        /// </summary>
        public bool CacheLocalization { get; set; }

        /// <summary>How often to look for queued jobs when no notification arrives.</summary>
        public TimeSpan PollInterval { get; set; } = TimeSpan.FromSeconds(30);

        /// <summary>A job running longer than this lost its worker (another process, before this one started).</summary>
        public TimeSpan StuckAfter { get; set; } = TimeSpan.FromMinutes(10);

        /// <summary>How often to apply D7's expiry and reload the score null and neighbour graphs (rebuilt by the Ingest).</summary>
        public TimeSpan MaintenanceInterval { get; set; } = TimeSpan.FromMinutes(10);
    }

    /// <summary>
    /// Runs queued searches (P4 S4.3), one at a time: each search already uses every core. It wakes on the <see cref="SearchLog.Channel"/>
    /// notification a new job sends and also polls, so a lost listener only delays jobs. It assumes it is the only worker (the API is
    /// one process until P6): a job still running when it starts lost its worker and is requeued, and after
    /// <see cref="SearchLog.MaxAttempts"/> claims it is failed instead. A job interrupted by shutdown stays running and is requeued
    /// on the next start.
    ///
    /// Find-similar keeps its band levels in memory (24 MB at ratio 0.5) and refreshes them incrementally before every job, so a
    /// record stored a moment ago is a candidate with its profile. A protein search streams the gene rows from the database every
    /// time (about 20 s on the local catalogue).
    /// </summary>
    internal sealed class SearchWorker(NpgsqlDataSource dataSource, IOptions<SearchWorkerOptions> options, ILogger<SearchWorker> logger)
        : BackgroundService
    {
        private readonly SemaphoreSlim wake = new(0, 1);

        protected override async Task ExecuteAsync(CancellationToken stoppingToken)
        {
            var settings = options.Value;
            if (!settings.Enabled) return;

            var similar = new SimilaritySearch(dataSource, useCache: true, method: RecordStore.Method, methodVersion: settings.MethodVersion);
            var localization = new Localization(dataSource, settings.CacheLocalization, RecordStore.Method, settings.MethodVersion);
            var listener = ListenAsync(settings.PollInterval, stoppingToken);
            bool recovered = false;
            var nextMaintenance = DateTimeOffset.MinValue;

            while (!stoppingToken.IsCancellationRequested)
            {
                try
                {
                    if (!recovered)
                    {
                        var (requeued, failed) = await SearchLog.RequeueStuckAsync(dataSource, DateTimeOffset.UtcNow, stoppingToken);
                        if (requeued + failed > 0) logger.LogWarning("Found {Requeued} interrupted search(es) requeued and {Failed} failed", requeued, failed);
                        recovered = true;
                    }
                    if (DateTimeOffset.UtcNow >= nextMaintenance)
                    {
                        await MaintainAsync(settings, similar, localization, stoppingToken);
                        nextMaintenance = DateTimeOffset.UtcNow + settings.MaintenanceInterval;
                    }
                    while (!stoppingToken.IsCancellationRequested && await SearchLog.ClaimAsync(dataSource, stoppingToken) is Guid id)
                        await RunAsync(id, similar, localization, stoppingToken);
                }
                catch (OperationCanceledException) when (stoppingToken.IsCancellationRequested)
                {
                    break;
                }
                catch (Exception e)
                {
                    // The database is unreachable or similar: try again at the next wake-up rather than stopping the API.
                    logger.LogError(e, "The search worker's loop failed; retrying in {Interval}", settings.PollInterval);
                }

                try { await wake.WaitAsync(settings.PollInterval, stoppingToken); }
                catch (OperationCanceledException) { break; }
            }
            await listener;
        }

        private async Task MaintainAsync(SearchWorkerOptions settings, SimilaritySearch similar, Localization localization, CancellationToken cancellationToken)
        {
            var expired = await SearchLog.ExpireDueAsync(dataSource, cancellationToken);
            var (requeued, failed) = await SearchLog.RequeueStuckAsync(dataSource, DateTimeOffset.UtcNow - settings.StuckAfter, cancellationToken);
            await similar.RefreshAsync(cancellationToken);
            await localization.RefreshAsync(cancellationToken);
            if (expired + requeued + failed > 0)
                logger.LogInformation("Maintenance: {Expired} search(es) expired (D7), {Requeued} stuck requeued, {Failed} failed", expired, requeued, failed);
        }

        private async Task RunAsync(Guid id, SimilaritySearch similar, Localization localization, CancellationToken cancellationToken)
        {
            var watch = Stopwatch.StartNew();
            try
            {
                // Null when the search was deleted after it was claimed.
                if (await SearchLog.GetAsync(dataSource, id, cancellationToken) is not { } run) return;

                var request = run.Request;
                bool saved;
                if (request.Kind == SearchKind.Protein)
                {
                    if (await SearchLog.ProteinQueryAsync(dataSource, id, cancellationToken) is not { } protein) return;
                    var result = await ProteinSearch.SearchAsync(dataSource, protein, request.Filter, request.Top, cancellationToken);
                    saved = await SearchLog.SaveAsync(dataSource, id, result, cancellationToken);
                }
                else if (await SearchLog.QueryAsync(dataSource, id, cancellationToken) is not { } query)
                {
                    return;
                }
                else if (request.Kind == SearchKind.Similar)
                {
                    await similar.RefreshProfilesAsync(cancellationToken);
                    var result = await similar.SearchAsync(query, request.Filter, request.LengthRatio ?? SimilaritySearch.DefaultLengthRatio,
                        request.Top, cancellationToken);
                    saved = await SearchLog.SaveAsync(dataSource, id, result, cancellationToken);
                }
                else
                {
                    await localization.RefreshProfilesAsync(cancellationToken);
                    var result = request.WithinSequenceId is long record
                        ? await localization.LocalizeInRecordAsync(query, record, cancellationToken)
                        : await localization.LocalizeAsync(query, request.Filter, request.Top, cancellationToken);
                    saved = await SearchLog.SaveAsync(dataSource, id, result, cancellationToken);
                }
                logger.LogInformation("Search {SearchId} ({Kind}, {Length}): {Outcome} in {Milliseconds} ms", id, request.Kind,
                    run.QueryLengthAa is int aa ? $"{aa} aa" : $"{run.QueryLengthBp} bp",
                    saved ? "done" : "deleted while it ran, result dropped", watch.ElapsedMilliseconds);
            }
            catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
            {
                throw;
            }
            catch (Exception e)
            {
                logger.LogError(e, "Search {SearchId} failed after {Milliseconds} ms", id, watch.ElapsedMilliseconds);
                // Core's argument checks are written for the caller (a record without a fine level, say); anything else stays in the log.
                var error = e is ArgumentException a
                    ? (a.ParamName is null ? a.Message : a.Message.Replace($" (Parameter '{a.ParamName}')", ""))
                    : "The search failed; the server log has the details.";
                await SearchLog.FailAsync(dataSource, id, error, CancellationToken.None);
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
                    await using (var listen = new NpgsqlCommand($"listen {SearchLog.Channel}", connection))
                        await listen.ExecuteNonQueryAsync(cancellationToken);
                    Wake(); // whatever was queued while nobody listened
                    while (true) await connection.WaitAsync(cancellationToken);
                }
                catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested)
                {
                    return;
                }
                catch (Exception e)
                {
                    logger.LogWarning(e, "Lost the {Channel} listener; polling until it reconnects", SearchLog.Channel);
                    try { await Task.Delay(retryAfter, cancellationToken); }
                    catch (OperationCanceledException) { return; }
                }
            }
        }

        private void Wake()
        {
            try { if (wake.CurrentCount == 0) wake.Release(); }
            catch (SemaphoreFullException) { } // a wake-up is already pending
        }
    }
}
