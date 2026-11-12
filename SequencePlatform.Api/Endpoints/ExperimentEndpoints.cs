using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// The experiments the articles state, as published by their scripts and the Ingest's publish-experiments (<see cref="ExperimentStore"/>):
    /// computed offline (on a frozen snapshot of the catalogue, or on the bundled genomes for the library ones), read-only here, with their
    /// result tables for download. Findings are about public records only, so every caller sees the same results.
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

            experiments.MapGet("/{name}/files/{file}", FileAsync)
                .WithName("GetExperimentFile")
                .WithSummary("Download one result table of an experiment")
                .WithDescription("origin is article (default: the run an article cites) or server (the latest rerun on this database).")
                .Produces(StatusCodes.Status200OK, contentType: "text/csv")
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status404NotFound);

            app.MapGroup("/api/sequences").WithTags("Experiments").MapGet("/{id:long}/findings", RecordAsync)
                .WithName("GetRecordFindings")
                .WithSummary("What the published experiments say about one record")
                .ProducesProblem(StatusCodes.Status404NotFound);
        }

        private static async Task<Ok<ExperimentListResponse>> ListAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken) =>
            TypedResults.Ok(new ExperimentListResponse((await ExperimentStore.ListAsync(dataSource, cancellationToken: cancellationToken)).Select(ExperimentResponse.From).ToList()));

        private static string CheckOrigin(string origin) => origin is ExperimentStore.ArticleOrigin or ExperimentStore.ServerOrigin
            ? origin
            : throw new ArgumentException($"origin is {ExperimentStore.ArticleOrigin} or {ExperimentStore.ServerOrigin}.");

        private static async Task<Results<FileContentHttpResult, ProblemHttpResult>> FileAsync(NpgsqlDataSource dataSource, string name, string file,
            string origin = ExperimentStore.ArticleOrigin, CancellationToken cancellationToken = default)
        {
            var run = await ExperimentStore.GetAsync(dataSource, name, CheckOrigin(origin), cancellationToken);
            var content = run is null ? null : await ExperimentStore.FileAsync(dataSource, run.Id, file, cancellationToken);
            if (content is null) return TypedResults.Problem($"No file '{file}' in experiment '{name}' ({origin}).", statusCode: StatusCodes.Status404NotFound);
            var contentType = file.EndsWith(".csv", StringComparison.OrdinalIgnoreCase) ? "text/csv; charset=utf-8" : "text/plain; charset=utf-8";
            return TypedResults.File(content, contentType, fileDownloadName: file);
        }

        private static async Task<Results<Ok<ExperimentPageResponse>, ProblemHttpResult>> GetAsync(NpgsqlDataSource dataSource, string name,
            string? kind, string order = "score", int offset = 0, int limit = 100, string origin = ExperimentStore.ArticleOrigin,
            CancellationToken cancellationToken = default)
        {
            if (order is not ("score" or "-score")) throw new ArgumentException("order is score or -score.");
            if (offset < 0 || limit < 1 || limit > MaxPage) throw new ArgumentException($"offset must be 0 or more and limit 1 to {MaxPage}.");
            var run = await ExperimentStore.GetAsync(dataSource, name, CheckOrigin(origin), cancellationToken);
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
