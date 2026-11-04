using System.Text.Json;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>A published P9 experiment.</summary>
    /// <param name="Verdict">pass, fail, or context (no pass rule).</param>
    /// <param name="Snapshot">The frozen copy of the catalogue it was computed on.</param>
    /// <param name="Summary">The experiment's own numbers: its rule, key figures against their baselines, and any tables its page draws.</param>
    /// <param name="Findings">How many records it says something about.</param>
    public sealed record ExperimentResponse(string Name, string Title, string Question, string Verdict, string Snapshot, string Script,
        DateTime PublishedAt, int Findings, JsonElement Summary, string Href)
    {
        internal static ExperimentResponse From(ExperimentRun run) => new(run.Name, run.Title, run.Question, run.Verdict, run.Snapshot, run.Script,
            run.PublishedAt, run.Findings, Json(run.SummaryJson), $"/api/experiments/{run.Name}");

        internal static JsonElement Json(string text)
        {
            using var document = JsonDocument.Parse(text);
            return document.RootElement.Clone();
        }
    }

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
}
