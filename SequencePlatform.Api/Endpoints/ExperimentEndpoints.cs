using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequencePlatform.Core;
using SequencePlatform.Experiments;

namespace SequencePlatform.Api
{
    /// <summary>
    /// The experiments the articles state, as published by their scripts and the Ingest's publish-experiments (<see cref="ExperimentStore"/>):
    /// computed offline (on a frozen snapshot of the catalogue, or on the bundled genomes for the library ones), read-only here, with their
    /// result tables for download. Findings are about public records only, so every caller sees the same results. An admin may rerun an
    /// experiment that needs only the database (migration 0024): the experiment worker publishes it as the server run beside the article's.
    /// </summary>
    internal static class ExperimentEndpoints
    {
        public const int MaxPage = 1000;
        public const string AdminPolicy = "admin";

        private const int jobsShown = 10;

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

            experiments.MapGet("/{name}/runs", RerunsAsync)
                .WithName("ListExperimentReruns")
                .WithSummary("Whether an experiment reruns on this database, whether the caller may queue it, and its latest reruns");

            experiments.MapPost("/{name}/runs", RerunAsync)
                .RequireAuthorization(AdminPolicy)
                .WithName("RerunExperiment")
                .WithSummary("Queue a rerun of a database experiment (admins only)")
                .WithDescription("Returns 202 and the job's URL. The experiment worker publishes the result as the experiment's server run " +
                    "(GET /api/experiments/{name}?origin=server), replacing the previous one. 409 while a rerun of it is queued or running.")
                .ProducesProblem(StatusCodes.Status403Forbidden)
                .ProducesProblem(StatusCodes.Status404NotFound)
                .ProducesProblem(StatusCodes.Status409Conflict);

            experiments.MapGet("/jobs/{id:guid}", JobAsync)
                .WithName("GetExperimentRerun")
                .WithSummary("One experiment rerun's status")
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

        // Rerunnable: the article run says it needs only the database, and this build has the code to run it.
        private static async Task<bool> RerunnableAsync(NpgsqlDataSource dataSource, IEnumerable<IExperiment> experiments, string name,
            CancellationToken cancellationToken) =>
            ExperimentRunner.Find(experiments, name) is not null
            && await ExperimentStore.GetAsync(dataSource, name, cancellationToken: cancellationToken) is { Details.Runnable: true };

        private static async Task<Ok<ExperimentRerunsResponse>> RerunsAsync(NpgsqlDataSource dataSource, IEnumerable<IExperiment> experiments,
            ClaimsPrincipal user, string name, CancellationToken cancellationToken)
        {
            bool rerunnable = await RerunnableAsync(dataSource, experiments, name, cancellationToken);
            var jobs = await ExperimentJobs.ListAsync(dataSource, name, jobsShown, cancellationToken);
            return TypedResults.Ok(new ExperimentRerunsResponse(name, rerunnable, rerunnable && user.IsInRole(ApiKeyAuthentication.AdminRole),
                jobs.Select(ExperimentJobResponse.From).ToList()));
        }

        private static async Task<Results<Accepted<ExperimentJobResponse>, ProblemHttpResult>> RerunAsync(NpgsqlDataSource dataSource,
            IEnumerable<IExperiment> experiments, ClaimsPrincipal user, string name, CancellationToken cancellationToken)
        {
            if (!await RerunnableAsync(dataSource, experiments, name, cancellationToken))
                return TypedResults.Problem($"'{name}' doesn't rerun from the database.", statusCode: StatusCodes.Status404NotFound);
            if (await ExperimentJobs.EnqueueAsync(dataSource, name, user.UserId(), cancellationToken) is not Guid id)
                return TypedResults.Problem($"A rerun of '{name}' is already queued or running.", statusCode: StatusCodes.Status409Conflict);
            var job = ExperimentJobResponse.From((await ExperimentJobs.GetAsync(dataSource, id, cancellationToken))!);
            return TypedResults.Accepted(job.Href, job);
        }

        private static async Task<Results<Ok<ExperimentJobResponse>, ProblemHttpResult>> JobAsync(NpgsqlDataSource dataSource, Guid id,
            CancellationToken cancellationToken) =>
            await ExperimentJobs.GetAsync(dataSource, id, cancellationToken) is { } job
                ? TypedResults.Ok(ExperimentJobResponse.From(job))
                : TypedResults.Problem($"No experiment rerun {id}.", statusCode: StatusCodes.Status404NotFound);

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
