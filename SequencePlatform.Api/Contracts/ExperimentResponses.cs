using System.Text.Json;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>A published experiment.</summary>
    /// <param name="Verdict">pass, fail, or context (no pass rule).</param>
    /// <param name="Snapshot">What it was computed on: a frozen copy of the catalogue, the live catalogue, or the bundled genomes.</param>
    /// <param name="Summary">The experiment's own numbers: its rule, key figures against their baselines, and any tables its page draws.</param>
    /// <param name="Findings">How many records it says something about.</param>
    /// <param name="Article">The article section it backs (A-genome, A-gene, B, P9), or null.</param>
    /// <param name="Origin">article (the run an article cites, computed offline) or server (a rerun on this database).</param>
    /// <param name="Command">How to reproduce it from the repository root.</param>
    /// <param name="Runnable">Whether it can be rerun from the database alone; experiments that need DNA can't.</param>
    /// <param name="Files">Its result tables, for download.</param>
    public sealed record ExperimentResponse(string Name, string Title, string Question, string Verdict, string Snapshot, string Script,
        DateTime PublishedAt, int Findings, JsonElement Summary, string Href, string? Article, string Origin, string? Commit, string? Machine,
        string? Command, bool Runnable, IReadOnlyList<ExperimentFileResponse> Files)
    {
        internal static ExperimentResponse From(ExperimentRun run)
        {
            var href = $"/api/experiments/{Uri.EscapeDataString(run.Name)}";
            var query = run.Details.Origin == ExperimentStore.ArticleOrigin ? "" : $"?origin={run.Details.Origin}";
            return new(run.Name, run.Title, run.Question, run.Verdict, run.Snapshot, run.Script, run.PublishedAt, run.Findings, Json(run.SummaryJson),
                href + query, run.Details.Article, run.Details.Origin, run.Details.Commit, run.Details.Machine, run.Details.Command, run.Details.Runnable,
                run.Files.Select(f => new ExperimentFileResponse(f.Name, f.ByteSize, $"{href}/files/{Uri.EscapeDataString(f.Name)}{query}")).ToList());
        }

        internal static JsonElement Json(string text)
        {
            using var document = JsonDocument.Parse(text);
            return document.RootElement.Clone();
        }
    }

    /// <summary>One result table of an experiment run.</summary>
    /// <param name="ByteSize">Uncompressed.</param>
    public sealed record ExperimentFileResponse(string Name, long ByteSize, string Href);

    /// <summary>Every published experiment, newest first.</summary>
    public sealed record ExperimentListResponse(IReadOnlyList<ExperimentResponse> Items);

    /// <summary>One finding: the record it's about, and what the experiment found.</summary>
    public sealed record ExperimentFindingResponse(SequenceSummaryResponse Sequence, string Kind, double? Score, JsonElement Payload, string Href);

    /// <summary>An experiment and a page of its findings.</summary>
    public sealed record ExperimentPageResponse(ExperimentResponse Experiment, long Total, int Offset, IReadOnlyList<ExperimentFindingResponse> Items);

    /// <summary>One experiment's finding about a record, for the record's page.</summary>
    public sealed record RecordFindingResponse(string Experiment, string Title, string Verdict, string Kind, double? Score, JsonElement Payload, string Href)
    {
        internal static RecordFindingResponse From(RecordFinding f) =>
            new(f.Experiment, f.Title, f.Verdict, f.Kind, f.Score, ExperimentResponse.Json(f.PayloadJson), $"/api/experiments/{f.Experiment}");
    }

    /// <summary>What the published experiments say about one record; empty when nothing.</summary>
    public sealed record RecordFindingsResponse(long Id, IReadOnlyList<RecordFindingResponse> Items);

    /// <summary>A rerun of an experiment on this database (migration 0024).</summary>
    /// <param name="RequestedBy">The admin who queued it.</param>
    /// <param name="RunHref">The server run it published, while that run is still the latest; null otherwise.</param>
    public sealed record ExperimentJobResponse(Guid Id, string Experiment, string RequestedBy, ExperimentJobStatus Status, string? Error, int Attempts,
        DateTime CreatedAt, DateTime? StartedAt, DateTime? CompletedAt, string? RunHref, string Href)
    {
        internal static ExperimentJobResponse From(ExperimentJob job) => new(job.Id, job.Experiment, job.RequestedBy, job.Status, job.Error, job.Attempts,
            job.CreatedAt, job.StartedAt, job.CompletedAt,
            job.RunId is null ? null : $"/api/experiments/{Uri.EscapeDataString(job.Experiment)}?origin={ExperimentStore.ServerOrigin}",
            $"/api/experiments/jobs/{job.Id}");
    }

    /// <summary>Whether an experiment reruns here, whether the caller may queue it, and its latest reruns.</summary>
    /// <param name="Rerunnable">Its article run is marked runnable and this build has its code.</param>
    /// <param name="CanRun">Rerunnable, and the caller is an admin.</param>
    public sealed record ExperimentRerunsResponse(string Experiment, bool Rerunnable, bool CanRun, IReadOnlyList<ExperimentJobResponse> Jobs);
}
