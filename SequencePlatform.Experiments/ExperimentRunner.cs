using System.Reflection;
using System.Runtime.InteropServices;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// Reruns an experiment and publishes the result as its <see cref="ExperimentStore.ServerOrigin"/> run, replacing the previous one, beside
    /// the article run it repeats. Used by the experiment worker and by the Ingest's run-experiment.
    /// </summary>
    public static class ExperimentRunner
    {
        /// <summary>Every experiment that reruns from the database alone, by name.</summary>
        public static IReadOnlyList<IExperiment> All { get; } = [new GcStratum(), new ProteinHomolog()];

        public static IExperiment? Find(IEnumerable<IExperiment> experiments, string name) => experiments.FirstOrDefault(e => e.Name == name);

        /// <returns>The server run's id, and what it published.</returns>
        /// <exception cref="ExperimentNotRunnableException">No article run of that name is published, or it isn't marked runnable.</exception>
        public static async Task<(long RunId, ExperimentResult Result)> RunAndPublishAsync(NpgsqlDataSource dataSource, IExperiment experiment, ExperimentContext context,
            CancellationToken cancellationToken)
        {
            var article = await ExperimentStore.GetAsync(dataSource, experiment.Name, cancellationToken: cancellationToken);
            if (article is null || !article.Details.Runnable)
                throw new ExperimentNotRunnableException($"{experiment.Name} has no runnable article run to repeat: publish it with publish-experiments first.");

            var result = await experiment.RunAsync(dataSource, context, cancellationToken);
            var details = new ExperimentDetails(article.Details.Article, ExperimentStore.ServerOrigin, SourceCommit(),
                $"{Environment.ProcessorCount} cores, {RuntimeInformation.OSDescription}",
                $"dotnet run --project SequencePlatform.Ingest -- run-experiment {experiment.Name} --publish", Runnable: true);
            await ExperimentStore.PublishAsync(dataSource, experiment.Name, article.Title, article.Question, result.Verdict, result.Snapshot, result.Script,
                result.SummaryJson, result.Findings, details, result.Files);
            return ((await ExperimentStore.GetAsync(dataSource, experiment.Name, ExperimentStore.ServerOrigin, cancellationToken))!.Id, result);
        }

        // The SDK stamps the commit after '+' when it builds from a git checkout; a build without one says nothing.
        private static string? SourceCommit()
        {
            var version = typeof(ExperimentRunner).Assembly.GetCustomAttribute<AssemblyInformationalVersionAttribute>()?.InformationalVersion;
            int plus = version?.IndexOf('+') ?? -1;
            return plus < 0 ? null : version![(plus + 1)..][..Math.Min(7, version.Length - plus - 1)];
        }
    }
}
