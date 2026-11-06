using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// The catalogue (P4 S4.2a). Every query passes the caller as the catalogue's viewer, so a list, a record or a profile is
    /// either public or the caller's own upload; someone else's upload is a 404, the same as a record that doesn't exist.
    /// </summary>
    internal static class SequenceEndpoints
    {
        public static void MapSequences(this WebApplication app)
        {
            var sequences = app.MapGroup("/api/sequences").WithTags("Sequences");

            sequences.MapGet("/", ListAsync)
                .WithName("ListSequences")
                .WithSummary("List the catalogue")
                .WithDescription("Keyset-paged in id order: pass the previous page's nextAfter as after. Filters combine with AND; " +
                    "rank and taxon go together (rank: superkingdom, realm, kingdom, phylum, class, order, family, genus, species). " +
                    "kind is a stored kind or a group from /api/sequences/facets (prokaryote: bacterium and archaeon).")
                .ProducesValidationProblem();

            sequences.MapGet("/facets", FacetsAsync)
                .WithName("GetSequenceFacets")
                .WithSummary("Filter choices: each kind, kind group and provider with its record count")
                .WithDescription("Counts cover what the caller can see: the public catalogue and the caller's own uploads.");

            sequences.MapGet("/{id:long}", GetAsync)
                .WithName("GetSequence")
                .WithSummary("One record: metadata, lineages, annotations, stored levels and provenance");

            sequences.MapGet("/{id:long}/profile", GetProfileAsync)
                .WithName("GetSequenceProfile")
                .WithSummary("One stored profile level as ZCWP bytes")
                .WithDescription("The binary format of docs/database-decision.md section 4.3: 'ZCWP', version, window size, window count, " +
                    "then the x, y and z window means as little-endian float32 blocks. windowSize is 500, 2000, 8000 or 32000 (default 500).")
                .Produces(StatusCodes.Status200OK, contentType: RequestBody.ProfileContentType);

            sequences.MapPost("/", SequenceUpload.UploadAsync)
                .WithName("UploadSequences")
                .WithSummary("Upload FASTA or a ZCWP profile as your own private records")
                .WithDescription("The body is the file: text/x-fasta (one or more records) or application/x-zcwp (the 500 bp level, with " +
                    "lengthBp). Metadata goes in the query string. Each record is profiled as it streams in and only its profile, hashes and " +
                    "metadata are stored, never the DNA. 201 with the new records; Location is the first one.")
                .Accepts<Stream>(RequestBody.FastaContentType, RequestBody.ProfileContentType)
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status413PayloadTooLarge)
                .ProducesProblem(StatusCodes.Status415UnsupportedMediaType);
        }

        private static async Task<Ok<SequencePageResponse>> ListAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            string? kind, string? rank, string? taxon, long? minLength, long? maxLength, string? provider, string? name,
            long after = 0, int limit = 50, CancellationToken cancellationToken = default)
        {
            var filter = new CatalogueFilter(kind, rank, taxon, minLength, maxLength, provider, name, ViewerId: user.UserId());
            var page = await Catalogue.ListAsync(dataSource, filter, after, limit, cancellationToken);
            return TypedResults.Ok(new SequencePageResponse(page.Records.Select(SequenceSummaryResponse.From).ToList(), page.NextAfter));
        }

        private static async Task<Ok<SequenceFacetsResponse>> FacetsAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            CancellationToken cancellationToken) =>
            TypedResults.Ok(SequenceFacetsResponse.From(await Catalogue.FacetsAsync(dataSource, user.UserId(), cancellationToken)));

        private static async Task<Results<Ok<SequenceDetailResponse>, NotFound>> GetAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            long id, CancellationToken cancellationToken)
        {
            var detail = await Catalogue.DetailAsync(dataSource, id, user.UserId(), cancellationToken);
            return detail is null ? TypedResults.NotFound() : TypedResults.Ok(SequenceDetailResponse.From(detail));
        }

        private static async Task<Results<FileContentHttpResult, ProblemHttpResult>> GetProfileAsync(NpgsqlDataSource dataSource,
            ClaimsPrincipal user, long id, int windowSize = RecordStore.FineWindowSize, CancellationToken cancellationToken = default)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var level = await ProfileLevels.ReadAsync(dataSource, id, RecordStore.Method, RecordStore.MethodVersion, windowSize);
            return level is null
                ? TypedResults.Problem($"Sequence {id} has no stored level of {windowSize} bp windows.", statusCode: StatusCodes.Status404NotFound)
                : TypedResults.File(level.ToBinary(), RequestBody.ProfileContentType, $"{id}-{windowSize}.zcwp");
        }
    }
}
