using System.Net.Http.Headers;
using Microsoft.AspNetCore.Hosting;
using Microsoft.AspNetCore.Mvc.Testing;
using Microsoft.AspNetCore.TestHost;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.DependencyInjection.Extensions;
using Microsoft.Extensions.Logging;
using Npgsql;
using SequencePlatform.Api;
using SequencePlatform.Core;
using SequencePlatform.Experiments;

namespace SequencePlatform.Tests
{
    /// <summary>The API in process over the throwaway database, connected as <see cref="Database.AppRole"/>, and clients holding a new user's key.</summary>
    internal static class ApiHost
    {
        /// <param name="workerMethodVersion">Runs the search worker over the profiles stored under this method version; null: no worker,
        /// so a test's queued jobs stay queued until it runs them itself.</param>
        /// <param name="logs">Also receives every log entry, at every level.</param>
        /// <param name="experiments">Runs the experiment worker, over these experiments in place of the real ones; null: no experiment worker.</param>
        // Any public type of the API assembly names its entry point; Program would be ambiguous with the Ingest's.
        public static WebApplicationFactory<SequencePageResponse> Create(string? workerMethodVersion = null, LogCapture? logs = null,
            IReadOnlyList<IExperiment>? experiments = null) =>
            new WebApplicationFactory<SequencePageResponse>().WithWebHostBuilder(host =>
            {
                if (logs is not null)
                    host.ConfigureLogging(logging => logging.AddProvider(logs).AddFilter<LogCapture>(null, LogLevel.Trace));
                host.ConfigureTestServices(services =>
                {
                    services.RemoveAll<NpgsqlDataSource>();
                    // As the hosted stack connects: the least-privilege role, not the owner the tests seed with.
                    services.AddSingleton(ThrowawayDatabase.AppDataSource);
                    services.Configure<SearchWorkerOptions>(worker =>
                    {
                        worker.Enabled = workerMethodVersion is not null;
                        worker.MethodVersion = workerMethodVersion ?? worker.MethodVersion;
                        worker.PollInterval = TimeSpan.FromSeconds(1);
                    });
                    services.Configure<ExperimentWorkerOptions>(worker =>
                    {
                        worker.Enabled = experiments is not null;
                        worker.WaitForSearches = false;
                        worker.PollInterval = TimeSpan.FromSeconds(1);
                    });
                    if (experiments is not null)
                    {
                        services.RemoveAll<IExperiment>();
                        foreach (var experiment in experiments) services.AddSingleton(experiment);
                    }
                });
            });

        public static async Task<(AppUser User, HttpClient Client)> NewClientAsync(this WebApplicationFactory<SequencePageResponse> factory, string name)
        {
            var (user, apiKey) = await Users.CreateAsync(ThrowawayDatabase.DataSource, $"{name}-{Guid.NewGuid():N}");
            var client = factory.CreateClient();
            client.DefaultRequestHeaders.Authorization = new AuthenticationHeaderValue("Bearer", apiKey);
            return (user, client);
        }
    }
}
