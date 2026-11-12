using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// A eukaryote's chromosome levels (migration 0022): each nuclear sequence of at least 100 kbp at 500 bp and 32 kbp windows, for drawing
    /// one chromosome's Z-curve path. Visibility is the catalogue's.
    /// </summary>
    internal static class ChromosomeEndpoints
    {
        public static void MapChromosomes(this WebApplication app)
        {
            var sequences = app.MapGroup("/api/sequences").WithTags("Chromosomes");

            sequences.MapGet("/{id:long}/chromosomes", ListAsync)
                .WithName("GetChromosomes")
                .WithSummary("A eukaryote's stored chromosome levels: each nuclear sequence of at least 100 kbp, its length and window sizes")
                .WithDescription("In the genome file's order. Records without chromosome levels (prokaryotes, viruses, uploads, eukaryotes not yet " +
                    "backfilled) answer with none; shorter sequences aren't kept, and annotations' ncbi.chromosome_levels says how many and how " +
                    "much of the genome they were.")
                .ProducesProblem(StatusCodes.Status404NotFound);

            sequences.MapGet("/{id:long}/chromosomes/{index:int}/profile", ProfileAsync)
                .WithName("GetChromosomeProfile")
                .WithSummary("One chromosome level as ZCWP bytes")
                .WithDescription("The same format as /api/sequences/{id}/profile (float32 window means). windowSize is 500 or 32000 (default 500).")
                .Produces(StatusCodes.Status200OK, contentType: RequestBody.ProfileContentType)
                .ProducesProblem(StatusCodes.Status404NotFound);
        }

        private static async Task<Results<Ok<ChromosomeListResponse>, ProblemHttpResult>> ListAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            long id, CancellationToken cancellationToken)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var levels = await ChromosomeLevels.ListAsync(dataSource, id, cancellationToken);
            return TypedResults.Ok(new ChromosomeListResponse(id, levels.GroupBy(l => l.Index).Select(g => new ChromosomeResponse(g.Key, g.First().Name,
                g.First().LengthBp, g.Select(l => new ChromosomeLevelResponse(l.WindowSize, l.WindowCount, l.ByteSize,
                    $"/api/sequences/{id}/chromosomes/{g.Key}/profile?windowSize={l.WindowSize}")).ToList())).ToList()));
        }

        private static async Task<Results<FileContentHttpResult, ProblemHttpResult>> ProfileAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            long id, int index, int windowSize = RecordStore.FineWindowSize, CancellationToken cancellationToken = default)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var level = await ChromosomeLevels.ReadAsync(dataSource, id, index, windowSize, cancellationToken);
            return level is null
                ? TypedResults.Problem($"Sequence {id} has no stored chromosome {index} at {windowSize} bp windows.", statusCode: StatusCodes.Status404NotFound)
                : TypedResults.File(level.ToBinary(), RequestBody.ProfileContentType, $"{id}-chromosome{index}-{windowSize}.zcwp");
        }
    }
}
