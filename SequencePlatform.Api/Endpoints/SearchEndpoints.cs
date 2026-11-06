using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// Find-similar, localize and protein search (P7) as asynchronous jobs (P4 S4.2c): POST profiles the query, checks every parameter, queues a
    /// <see cref="SearchLog"/> job and answers 202 with its address; GET reports its status and, once done, pages through its hit
    /// lists; DELETE drops its results early (D7). The query's DNA lives only for the POST: the job keeps both strands' profiles and
    /// hashes. A search belongs to whoever ran it; anyone else gets a 404, and its candidates are the public catalogue plus the
    /// owner's uploads (<see cref="CatalogueFilter.ViewerId"/>). The worker that runs queued jobs is S4.3.
    /// </summary>
    internal static class SearchEndpoints
    {
        public const int MaxTop = 1000;
        public const int MaxPageSize = 1000;

        private const int defaultSimilarTop = 50;
        private const int defaultLocalizeTop = 20;

        public static void MapSearches(this WebApplication app)
        {
            var searches = app.MapGroup("/api/searches").WithTags("Searches");

            searches.MapPost("/similar", SimilarAsync)
                .WithName("FindSimilar")
                .WithSummary("Queue a find-similar search")
                .WithDescription("The body is the query: one FASTA record (text/x-fasta), or its 500 bp ZCWP level (application/x-zcwp), " +
                    "optionally followed directly by the reverse strand's. Candidates are the records of length [ratio × query, query / ratio] " +
                    "(ratio 0.25-1, default 0.5) matching the catalogue filters. 202 with the search's address.")
                .Accepts<Stream>(RequestBody.FastaContentType, RequestBody.ProfileContentType)
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status415UnsupportedMediaType);

            searches.MapPost("/localize", LocalizeAsync)
                .WithName("Localize")
                .WithSummary("Queue a localization: where in the catalogue, or in one record, a fragment comes from")
                .WithDescription("The body is the fragment, as for find-similar. Over the catalogue it must be at least 4 kbp, and the " +
                    "candidates are records at least as long; with record, it is placed inside that record only, and the filters don't apply. " +
                    "top is 1-100. 202 with the search's address.")
                .Accepts<Stream>(RequestBody.FastaContentType, RequestBody.ProfileContentType)
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status404NotFound)
                .ProducesProblem(StatusCodes.Status415UnsupportedMediaType);

            searches.MapPost("/protein", ProteinAsync)
                .WithName("FindGenes")
                .WithSummary("Queue a protein search: which stored genes encode this protein")
                .WithDescription($"The body is one protein: a FASTA record (text/x-fasta) or the bare letters (text/plain), at least " +
                    $"{ProteinSearch.MinLengthAa} aa. It is read with genetic code 11 (default) or 4, and scored from the protein alone against " +
                    "the genes BV-BRC annotates on the records matching the catalogue filters: from 500 aa the genes of 500 aa and up, below that " +
                    "the genes within 10% of its length. " +
                    "top is 1-1000, default 50. 202 with the search's address.")
                .Accepts<Stream>(RequestBody.FastaContentType, "text/plain")
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status415UnsupportedMediaType);

            searches.MapGet("/{id:guid}", GetAsync)
                .WithName("GetSearch")
                .WithSummary("A search's status and, once done, one page of one hit list")
                .WithDescription("list is Ranked (default), Identical or NotComparable; page with after (the previous page's nextAfter) and limit.");

            searches.MapGet("/{id:guid}/query", QueryWindowsAsync)
                .WithName("GetProteinQueryWindows")
                .WithSummary("A protein search's query as windows, to draw over a hit gene")
                .WithDescription("The expected windows the protein was compared with (mean, trend and curvature of x, y and z over codon positions 1 " +
                    "and 2), in the gene windows' format; never the protein. 404 for another user's search, a DNA search, or once the search expired.");

            searches.MapDelete("/{id:guid}", DeleteAsync)
                .WithName("DeleteSearch")
                .WithSummary("Drop a search's query profile and results now rather than after 30 days")
                .WithDescription("The search stays listed as Expired with its parameters, the query's length and hash, and its summary.");
        }

        /// <param name="ratio">Candidates are records of length [ratio × query length, query length / ratio]; 0.25 to 1, default 0.5.</param>
        /// <param name="top">Ranked hits kept, 1 to 1000; default 50.</param>
        /// <param name="lengthBp">ZCWP only: the query's exact length; by default windows × 500.</param>
        private static async Task<Results<Accepted<SearchAcceptedResponse>, ProblemHttpResult>> SimilarAsync(HttpRequest request,
            NpgsqlDataSource dataSource, ClaimsPrincipal user, ILoggerFactory loggers, string? kind, string? rank, string? taxon, long? minLength,
            long? maxLength, string? provider, string? name, double ratio = SimilaritySearch.DefaultLengthRatio, int top = defaultSimilarTop,
            long? lengthBp = null, CancellationToken cancellationToken = default)
        {
            if (ratio < SimilaritySearch.SmallestLengthRatio || ratio > 1)
                throw new ArgumentException($"ratio must be {SimilaritySearch.SmallestLengthRatio} to 1.");
            var filter = Filter(user, kind, rank, taxon, minLength, maxLength, provider, name);
            return await EnqueueAsync(request, dataSource, user, loggers, new SearchRequest(SearchKind.Similar, filter, CheckTop(top, MaxTop), ratio),
                lengthBp, null, cancellationToken);
        }

        /// <param name="record">Place the fragment inside this record only (yours or public); the filters then don't apply.</param>
        /// <param name="top">Placements kept, 1 to 100; default 20.</param>
        /// <param name="lengthBp">ZCWP only: the fragment's exact length; by default windows × 500.</param>
        private static async Task<Results<Accepted<SearchAcceptedResponse>, ProblemHttpResult>> LocalizeAsync(HttpRequest request,
            NpgsqlDataSource dataSource, ClaimsPrincipal user, ILoggerFactory loggers, string? kind, string? rank, string? taxon, long? minLength,
            long? maxLength, string? provider, string? name, long? record = null, int top = defaultLocalizeTop, long? lengthBp = null,
            CancellationToken cancellationToken = default)
        {
            var filter = Filter(user, kind, rank, taxon, minLength, maxLength, provider, name);
            if (record is long id && !await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            return await EnqueueAsync(request, dataSource, user, loggers,
                new SearchRequest(SearchKind.Localize, filter, CheckTop(top, Localization.RefineDepth), WithinSequenceId: record),
                lengthBp, query =>
                {
                    if (record is null && query.LengthBp < Localization.MinCatalogueLengthBp)
                        throw new ArgumentException($"A fragment under {Localization.MinCatalogueLengthBp:N0} bp ({query.LengthBp:N0} here) can't be " +
                            "localized over the catalogue: name a record to place it in.");
                }, cancellationToken);
        }

        /// <param name="code">The genetic code the protein is read with: 11 (default) or 4 (Mycoplasma and relatives, TGA = Trp).</param>
        /// <param name="top">Genes kept, 1 to 1000; default 50.</param>
        private static async Task<Results<Accepted<SearchAcceptedResponse>, ProblemHttpResult>> ProteinAsync(HttpRequest request,
            NpgsqlDataSource dataSource, ClaimsPrincipal user, ILoggerFactory loggers, string? kind, string? rank, string? taxon, long? minLength,
            long? maxLength, string? provider, string? name, int code = ProteinQuery.StandardCode, int top = defaultSimilarTop,
            CancellationToken cancellationToken = default)
        {
            var filter = Filter(user, kind, rank, taxon, minLength, maxLength, provider, name);
            var search = new SearchRequest(SearchKind.Protein, filter, CheckTop(top, MaxTop), GeneticCode: code);
            if (RequestBody.Format(request) != BodyFormat.Fasta) return RequestBody.Unsupported();
            var query = ProteinQuery.FromProtein(await RequestBody.ReadProteinAsync(request, cancellationToken), code);

            var id = await SearchLog.EnqueueAsync(dataSource, search, query, cancellationToken);
            // Hash and length only: the protein never reaches a log.
            loggers.CreateLogger(typeof(SearchEndpoints).FullName!).LogInformation(
                "User {UserId} queued {Kind} search {SearchId}: {Residues} aa, code {Code}, sha256 {Sha256}",
                user.UserId(), search.Kind, id, query.LengthAa, code, Convert.ToHexStringLower(query.Sha256!));

            var href = $"/api/searches/{id}";
            return TypedResults.Accepted(href, new SearchAcceptedResponse(id, SearchStatus.Queued, href));
        }

        private static CatalogueFilter Filter(ClaimsPrincipal user, string? kind, string? rank, string? taxon, long? minLength, long? maxLength,
            string? provider, string? name)
        {
            var filter = new CatalogueFilter(kind, rank, taxon, minLength, maxLength, provider, name, ViewerId: user.UserId());
            Catalogue.Validate(filter);
            return filter;
        }

        private static int CheckTop(int top, int max) =>
            top >= 1 && top <= max ? top : throw new ArgumentException($"top must be 1 to {max}.");

        private static async Task<Results<Accepted<SearchAcceptedResponse>, ProblemHttpResult>> EnqueueAsync(HttpRequest request,
            NpgsqlDataSource dataSource, ClaimsPrincipal user, ILoggerFactory loggers, SearchRequest search, long? lengthBp,
            Action<QueryProfile>? checkQuery, CancellationToken cancellationToken)
        {
            var format = RequestBody.Format(request);
            if (format == BodyFormat.Unsupported) return RequestBody.Unsupported();
            var query = await RequestBody.ReadQueryAsync(request, format, lengthBp, cancellationToken);
            checkQuery?.Invoke(query);

            var id = await SearchLog.EnqueueAsync(dataSource, search, query, cancellationToken);
            // Hash and length only: the query never reaches a log (S4.4).
            loggers.CreateLogger(typeof(SearchEndpoints).FullName!).LogInformation(
                "User {UserId} queued {Kind} search {SearchId}: {Bases} bp, {Strands} strand(s), sha256 {Sha256}",
                user.UserId(), search.Kind, id, query.LengthBp, query.Reverse is null ? 1 : 2,
                query.Sha256 is null ? "none (profile)" : Convert.ToHexStringLower(query.Sha256));

            var href = $"/api/searches/{id}";
            return TypedResults.Accepted(href, new SearchAcceptedResponse(id, SearchStatus.Queued, href));
        }

        /// <param name="list">Ranked (default), Identical or NotComparable.</param>
        /// <param name="after">The previous page's nextAfter; 0 for the first page.</param>
        /// <param name="limit">Hits per page, 1 to 1000; default 50.</param>
        private static async Task<Results<Ok<SearchRunResponse>, NotFound>> GetAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user, Guid id,
            HitList list = HitList.Ranked, int after = 0, int limit = 50, CancellationToken cancellationToken = default)
        {
            if (limit < 1 || limit > MaxPageSize) throw new ArgumentException($"limit must be 1 to {MaxPageSize}.");
            if (await OwnRunAsync(dataSource, user, id, cancellationToken) is not { } run) return TypedResults.NotFound();
            if (run.Status != SearchStatus.Done) return TypedResults.Ok(SearchRunResponse.From(run, null));

            // One row past the page tells whether another page exists.
            var hits = await SearchLog.HitsAsync(dataSource, id, list, after, limit + 1, cancellationToken);
            var page = hits.Take(limit).ToList();
            var records = await Catalogue.ByIdsAsync(dataSource, page.Select(h => h.SequenceId).ToHashSet(), cancellationToken);
            var products = await ProteinFamilies.ProductsAsync(dataSource, page.Where(h => h.FamilyId is not null).Select(h => h.FamilyId!.Value).ToHashSet(),
                cancellationToken);
            var items = page.Select(h => new SearchHitResponse(h.Position, h.SequenceId,
                records.TryGetValue(h.SequenceId, out var r) ? SequenceSummaryResponse.From(r) : null, h.Strand, h.Score, h.CoarseScore, h.OffsetBp,
                h.P?.Log10, h.P?.Kind, h.Log10E, h.SharedNeighbours, h.Rotated, h.GeneStart, h.GeneEnd,
                h.FamilyId is int family ? ProteinFamilies.Format(family) : null,
                h.FamilyId is int known && products.TryGetValue(known, out var product) ? product : null)).ToList();
            return TypedResults.Ok(SearchRunResponse.From(run,
                new SearchHitPageResponse(list, items, hits.Count > limit ? page[^1].Position : null)));
        }

        private static async Task<Results<Ok<ProteinQueryWindowsResponse>, NotFound>> QueryWindowsAsync(NpgsqlDataSource dataSource,
            ClaimsPrincipal user, Guid id, CancellationToken cancellationToken)
        {
            if (await OwnRunAsync(dataSource, user, id, cancellationToken) is null) return TypedResults.NotFound();
            if (await SearchLog.ProteinQueryAsync(dataSource, id, cancellationToken) is not { } query) return TypedResults.NotFound();
            const int features = ZCurveGeneProfile.FeatureCount;
            return TypedResults.Ok(new ProteinQueryWindowsResponse(id, query.LengthAa, ZCurveGeneProfile.WindowOf(query.LengthAa),
                Enumerable.Range(0, query.Axes.Length / features).Select(w => (IReadOnlyList<double>)query.Axes[(w * features)..((w + 1) * features)]).ToList()));
        }

        private static async Task<Results<NoContent, NotFound>> DeleteAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user, Guid id,
            CancellationToken cancellationToken)
        {
            if (await OwnRunAsync(dataSource, user, id, cancellationToken) is null) return TypedResults.NotFound();
            await SearchLog.ExpireAsync(dataSource, id, cancellationToken);
            return TypedResults.NoContent();
        }

        // Someone else's search, or one run from the Ingest (no owner), is a 404 like one that doesn't exist.
        private static async Task<SearchRun?> OwnRunAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user, Guid id, CancellationToken cancellationToken) =>
            await SearchLog.GetAsync(dataSource, id, cancellationToken) is { } run && run.Request.Filter.ViewerId == user.UserId() ? run : null;
    }
}
