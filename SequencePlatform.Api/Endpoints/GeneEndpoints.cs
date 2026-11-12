using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Microsoft.AspNetCore.Mvc;
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
                .WithDescription("Genes of at least 50 aa that BV-BRC annotates on the record, with family and product, or that NCBI annotates on " +
                    "a eukaryote's nuclear sequences (no family), in order of sequence and start. Records without a gene layer (uploads, " +
                    "unannotated eukaryotes) answer with no genes.")
                .ProducesProblem(StatusCodes.Status404NotFound);

            sequences.MapGet("/{id:long}/genes/path", PathAsync)
                .WithName("GetGenePath")
                .WithSummary("The record's Z-curve path built from its genes alone, in chromosome order")
                .WithDescription("origin=true starts at the replication origin when the record has one (as the 1-origin levels do), else at 0. " +
                    "sequence=N draws only the genes on the record's sequence N (an index into /genes' sequences; a eukaryote's chromosomes " +
                    "are drawn one at a time, never joined); without it the sequences follow one another. " +
                    "Each gene window is one step: window positions times its means of x, y, z over codon positions 1 and 2; reverse-strand " +
                    "genes are flipped onto the forward strand, and intergenic DNA is left out. Records without genes answer with no points.")
                .ProducesProblem(StatusCodes.Status404NotFound);

            sequences.MapGet("/{id:long}/genes/{index:int}", WindowsAsync)
                .WithName("GetGeneWindows")
                .WithSummary("One gene's stored windows: mean, trend and curvature of x, y and z over codon positions 1 and 2")
                .ProducesProblem(StatusCodes.Status404NotFound);

            sequences.MapGet("/{id:long}/genes/{index:int}/matches", MatchesAsync)
                .WithName("GetGeneMatches")
                .WithSummary("A gene's best matches among another record's genes, by the Z-curve gene score, flagged when reciprocal")
                .WithDescription("in=ID is the record searched; top (1–20, default 5) how many matches. Genes from 500 aa are compared with that " +
                    "record's long genes, shorter ones with its genes within ±10% of their length. A reciprocal best hit (the match's own best " +
                    "match back is this gene) is how orthologs are usually paired; between same-genus eukaryotes 93–99.8% of them carry the same " +
                    "protein name (results/GeneMatch_Check.csv). For records without families, such as NCBI eukaryotes. With in=ID equal to " +
                    "the gene's own record, the gene itself is left out, so the matches are its closest paralogs.")
                .ProducesValidationProblem()
                .ProducesProblem(StatusCodes.Status404NotFound);
        }

        private static async Task<Results<Ok<GeneListResponse>, ProblemHttpResult>> ListAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            long id, CancellationToken cancellationToken)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken))
                return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var profile = await GeneProfiles.ReadAsync(dataSource, id);
            if (profile is null) return TypedResults.Ok(new GeneListResponse(id, 0, 0, [], [], []));

            var products = await ProteinFamilies.ProductsAsync(dataSource,
                profile.Genes.Where(g => g.Family != ZCurveGeneProfile.NoFamily).Select(g => g.Family).Distinct().ToList(), cancellationToken);
            var names = await GeneProfiles.SequenceNamesAsync(dataSource, [id], cancellationToken);
            var mean = profile.MeanWindow();
            return TypedResults.Ok(new GeneListResponse(id, profile.Genes.Count, profile.WindowCount,
                GeneProfiles.Features.Select((f, i) => new SignatureValueResponse(f.Name, mean[i], f.Meaning, null)).ToList(),
                profile.Genes.Select((g, i) => Gene(g, i, products)).ToList(), names.GetValueOrDefault(id) ?? []));
        }

        private static async Task<Results<Ok<GenePathResponse>, ProblemHttpResult>> PathAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user,
            long id, bool origin = false, int? sequence = null, CancellationToken cancellationToken = default)
        {
            var detail = await Catalogue.DetailAsync(dataSource, id, user.UserId(), cancellationToken);
            if (detail is null) return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var profile = await GeneProfiles.ReadAsync(dataSource, id);
            if (profile is null) return TypedResults.Ok(new GenePathResponse(id, 0, 0, [], [], [], sequence));
            if (sequence is { } s && !profile.Genes.Any(g => g.Chromosome == s))
                return TypedResults.Problem($"Sequence {id} has no genes on its sequence {s}.", statusCode: StatusCodes.Status404NotFound);

            // A replication origin is a single circular chromosome's; one sequence of a eukaryote is drawn from its start.
            long lengthBp = detail.Record.LengthBp;
            long startBp = origin && sequence is null ? await OriginRotation.OriginBpAsync(dataSource, id, lengthBp) ?? 0 : 0;
            var path = GenePath.Build(profile, lengthBp, startBp, sequence);
            var coordinates = new double[3 * path.Points.Count];
            for (int i = 0; i < path.Points.Count; i++)
                (coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]) =
                    (Math.Round(path.Points[i].X, 1), Math.Round(path.Points[i].Y, 1), Math.Round(path.Points[i].Z, 1));
            return TypedResults.Ok(new GenePathResponse(id, path.StartBp, path.Genes.Count, coordinates,
                path.Genes.Select(g => g.Gene).ToList(), path.Genes.Select(g => g.FirstPoint).ToList(), sequence));
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

        private static async Task<Results<Ok<GeneMatchListResponse>, ValidationProblem, ProblemHttpResult>> MatchesAsync(NpgsqlDataSource dataSource,
            ClaimsPrincipal user, long id, int index, [FromQuery(Name = "in")] long? other, int top = 5, CancellationToken cancellationToken = default)
        {
            if (other is not { } target)
                return TypedResults.ValidationProblem(new Dictionary<string, string[]> { ["in"] = ["Give the record to search, as in=ID."] });
            if (top is < 1 or > 20) return TypedResults.ValidationProblem(new Dictionary<string, string[]> { ["top"] = ["top is 1 to 20."] });
            foreach (var record in new[] { id, target })
                if (!await Catalogue.IsVisibleAsync(dataSource, record, user.UserId(), cancellationToken))
                    return TypedResults.Problem($"No sequence {record}.", statusCode: StatusCodes.Status404NotFound);
            var profile = await GeneProfiles.ReadAsync(dataSource, id);
            if (profile is null || index < 0 || index >= profile.Genes.Count)
                return TypedResults.Problem($"Sequence {id} has no gene {index}.", statusCode: StatusCodes.Status404NotFound);
            var searched = target == id ? profile : await GeneProfiles.ReadAsync(dataSource, target);
            if (searched is null) return TypedResults.Problem($"Sequence {target} has no stored genes.", statusCode: StatusCodes.Status404NotFound);

            var (matches, compared) = GeneMatching.Find(profile, index, searched, top, sameRecord: target == id);
            var families = matches.Select(m => searched.Genes[m.Gene].Family).Append(profile.Genes[index].Family)
                .Where(f => f != ZCurveGeneProfile.NoFamily).Distinct().ToList();
            var products = families.Count == 0 ? new Dictionary<int, string>() : await ProteinFamilies.ProductsAsync(dataSource, families, cancellationToken);
            return TypedResults.Ok(new GeneMatchListResponse(id, Gene(profile.Genes[index], index, products), target, compared,
                matches.Select(m => new GeneMatchResponse(Gene(searched.Genes[m.Gene], m.Gene, products), Math.Round(m.Score, 6), m.Reciprocal)).ToList()));
        }

        private static GeneResponse Gene(ZCurveGeneProfile.Gene g, int index, Dictionary<int, string> products) =>
            new(index, g.Start, g.End, g.Reverse ? Strand.Reverse : Strand.Forward, g.Residues, g.WindowCount,
                g.Family == ZCurveGeneProfile.NoFamily ? null : ProteinFamilies.Format(g.Family), products.GetValueOrDefault(g.Family), g.Chromosome);
    }
}
