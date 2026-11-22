using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// An article experiment that reruns from the database alone (Phase 1 of the experiments page): no DNA and no download. Its title,
    /// question and article section are the article run's; a rerun is published beside it as the <see cref="ExperimentStore.ServerOrigin"/> run.
    /// </summary>
    public interface IExperiment
    {
        /// <summary>The name its article run is published under (results/experiments.json).</summary>
        string Name { get; }

        Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken);
    }

    /// <param name="Parallelism">The cores it may use.</param>
    /// <param name="Progress">Where it reports what it is doing, for the log.</param>
    public sealed record ExperimentContext(int Parallelism, Action<string> Progress)
    {
        public static ExperimentContext Default { get; } = new(Environment.ProcessorCount, _ => { });
    }

    /// <summary>What a run found, ready for <see cref="ExperimentStore.PublishAsync"/>.</summary>
    /// <param name="Verdict">pass, fail or context, by the rule the article run fixed.</param>
    /// <param name="Snapshot">The database it ran on, as of when.</param>
    /// <param name="Script">The code that computed it.</param>
    /// <param name="SummaryJson">At least keyNumbers, with the article run's labels where it reports the same number, and the rule's outcome.</param>
    public sealed record ExperimentResult(string Verdict, string Snapshot, string Script, string SummaryJson, IReadOnlyCollection<ExperimentFile> Files,
        IReadOnlyCollection<ExperimentFinding> Findings);

    /// <summary>Why an experiment can't run here, written for the admin who queued it.</summary>
    public sealed class ExperimentNotRunnableException(string message) : InvalidOperationException(message);

    internal static class Snapshots
    {
        /// <summary>The live database as of now: its largest record id, which says which records a rerun could see.</summary>
        public static async Task<string> LiveAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken)
        {
            await using var command = dataSource.CreateCommand("select coalesce(max(id), 0) from sequence_record");
            var maxId = (long)(await command.ExecuteScalarAsync(cancellationToken))!;
            return $"live database, {DateTime.UtcNow:yyyy-MM-dd}, records ≤ {maxId}";
        }
    }
}
