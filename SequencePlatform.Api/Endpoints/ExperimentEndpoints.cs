using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// The P9 experiments as published by their scripts (<see cref="ExperimentStore"/>): computed offline on a frozen snapshot of the
    /// catalogue, read-only here. Findings are about public records only, so every caller sees the same results.
    /// </summary>
    internal static class ExperimentEndpoints
    {
        public const int MaxPage = 1000;

        public static void MapExperiments(this WebApplication app)
        {
            var experiments = app.MapGroup("/api/experiments").WithTags("Experiments");

            experiments.MapGet("/", ListAsync)
                .WithName("ListExperiments")
                .WithSummary("Every published experiment, newest first, with its verdict and summary numbers");

            experiments.MapGet("/{name}", GetAsync)
                .WithName("GetExperiment")
                .WithSummary("One experiment and a page of its findings")
                .WithDescription("kind keeps one sort of finding; order is score (ascending, default) or -score; limit is 1-1000, default 100.")
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status404NotFound);

            app.MapGroup("/api/sequences").WithTags("Experiments").MapGet("/{id:long}/findings", RecordAsync)
                .WithName("GetRecordFindings")
                .WithSummary("What the published experiments say about one record")
                .ProducesProblem(StatusCodes.Status404NotFound);
        }

        private static async Task<Ok<ExperimentListResponse>> ListAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken) =>
            TypedResults.Ok(new ExperimentListResponse((await ExperimentStore.ListAsync(dataSource, cancellationToken)).Select(ExperimentResponse.From).ToList()));

        private static async Task<Results<Ok<ExperimentPageResponse>, ProblemHttpResult>> GetAsync(NpgsqlDataSource dataSource, string name,
            string? kind, string order = "score", int offset = 0, int limit = 100, CancellationToken cancellationToken = default)
        {
            if (order is not ("score" or "-score")) throw new ArgumentException("order is score or -score.");
            if (offset < 0 || limit < 1 || limit > MaxPage) throw new ArgumentException($"offset must be 0 or more and limit 1 to {MaxPage}.");
            var run = await ExperimentStore.GetAsync(dataSource, name, cancellationToken);
            if (run is null) return TypedResults.Problem($"No experiment '{name}'.", statusCode: StatusCodes.Status404NotFound);

            var (findings, total) = await ExperimentStore.FindingsAsync(dataSource, run.Id, kind, order == "-score", offset, limit, cancellationToken);
            var records = await Catalogue.ByIdsAsync(dataSource, findings.Select(f => f.SequenceId).Distinct().ToList(), cancellationToken);
            return TypedResults.Ok(new ExperimentPageResponse(ExperimentResponse.From(run), total, offset, findings
                .Where(f => records.ContainsKey(f.SequenceId))
                .Select(f => new ExperimentFindingResponse(SequenceSummaryResponse.From(records[f.SequenceId]), f.Kind, f.Score,
                    ExperimentResponse.Json(f.PayloadJson), $"/api/sequences/{f.SequenceId}"))
                .ToList()));
        }

        private static async Task<Results<Ok<RecordFindingsResponse>, ProblemHttpResult>> RecordAsync(NpgsqlDataSource dataSource,
            ClaimsPrincipal user, long id, CancellationToken cancellationToken)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var findings = await ExperimentStore.ForRecordAsync(dataSource, id, cancellationToken);
            return TypedResults.Ok(new RecordFindingsResponse(id, findings.Select(RecordFindingResponse.From).ToList()));
        }
    }
}
