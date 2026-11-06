using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// A record's gene layer (P7): its genes with their BV-BRC families, its mean gene window, and one gene's windows for drawing its
    /// Z-curve. Read from that record's gene_profile rows (about 100 KB for a bacterium); visibility is the catalogue's.
    /// </summary>
    internal static class GeneEndpoints
    {
        public static void MapGenes(this WebApplication app)
        {
            var sequences = app.MapGroup("/api/sequences").WithTags("Genes");

            sequences.MapGet("/{id:long}/genes", ListAsync)
                .WithName("GetGenes")
                .WithSummary("A record's stored genes and its mean gene window")
                .WithDescription("Genes of at least 50 aa that BV-BRC annotates on the record, in order of start, with family and product. Records " +
                    "without a gene layer (uploads, NCBI eukaryotes) answer with no genes.")
                .ProducesProblem(StatusCodes.Status404NotFound);

            sequences.MapGet("/{id:long}/genes/{index:int}", WindowsAsync)
                .WithName("GetGeneWindows")
                .WithSummary("One gene's stored windows: mean, trend and curvature of x, y and z over codon positions 1 and 2")
                .ProducesProblem(StatusCodes.Status404NotFound);
        }

        private static async Task<Results<Ok<GeneListResponse>, ProblemHttpResult>> ListAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            long id, CancellationToken cancellationToken)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var profile = await GeneProfiles.ReadAsync(dataSource, id);
            if (profile is null) return TypedResults.Ok(new GeneListResponse(id, 0, 0, [], []));

            var products = await ProteinFamilies.ProductsAsync(dataSource,
                profile.Genes.Where(g => g.Family != ZCurveGeneProfile.NoFamily).Select(g => g.Family).Distinct().ToList(), cancellationToken);
            var mean = profile.MeanWindow();
            return TypedResults.Ok(new GeneListResponse(id, profile.Genes.Count, profile.WindowCount,
                GeneProfiles.Features.Select((f, i) => new SignatureValueResponse(f.Name, mean[i], f.Meaning, null)).ToList(),
                profile.Genes.Select((g, i) => Gene(g, i, products)).ToList()));
        }

        private static async Task<Results<Ok<GeneWindowsResponse>, ProblemHttpResult>> WindowsAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            long id, int index, CancellationToken cancellationToken)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var profile = await GeneProfiles.ReadAsync(dataSource, id);
            if (profile is null || index < 0 || index >= profile.Genes.Count)
                return TypedResults.Problem($"Sequence {id} has no gene {index}.", statusCode: StatusCodes.Status404NotFound);

            var gene = profile.Genes[index];
            var products = gene.Family == ZCurveGeneProfile.NoFamily ? [] : await ProteinFamilies.ProductsAsync(dataSource, [gene.Family], cancellationToken);
            var windows = profile.WindowsOf(index);
            return TypedResults.Ok(new GeneWindowsResponse(id, Gene(gene, index, products), ZCurveGeneProfile.WindowOf(gene.Residues),
                Enumerable.Range(0, gene.WindowCount)
                    .Select(w => (IReadOnlyList<double>)windows.AsSpan(w * ZCurveGeneProfile.FeatureCount, ZCurveGeneProfile.FeatureCount).ToArray()).ToList()));
        }

        private static GeneResponse Gene(ZCurveGeneProfile.Gene g, int index, Dictionary<int, string> products) =>
            new(index, g.Start, g.End, g.Reverse ? Strand.Reverse : Strand.Forward, g.Residues, g.WindowCount,
                g.Family == ZCurveGeneProfile.NoFamily ? null : ProteinFamilies.Format(g.Family), products.GetValueOrDefault(g.Family));
    }
}
