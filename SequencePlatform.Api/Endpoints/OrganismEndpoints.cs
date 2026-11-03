using System.Globalization;
using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// Organizing organisms (P5 S5.5): nearest relatives, the signature behind them in either view, and trees, all from the
    /// dinucleotide signature (<see cref="SignatureIndex"/>). They answer at once: a query reads the in-memory index, a tree of up to
    /// <see cref="SignatureIndex.MaxTreeLeaves"/> leaves takes well under a second. Visibility is the catalogue's: public records
    /// plus the caller's own uploads; anything else is a 404, as if it didn't exist.
    /// </summary>
    internal static class OrganismEndpoints
    {
        public const int MaxRelatives = 100;

        public static void MapOrganisms(this WebApplication app)
        {
            var sequences = app.MapGroup("/api/sequences").WithTags("Organisms");

            sequences.MapGet("/{id:long}/relatives", RelativesAsync)
                .WithName("GetRelatives")
                .WithSummary("The genomes nearest to a record by dinucleotide signature")
                .WithDescription("Relatives come from the record's group: bacteria and archaea together, or viruses. One record stands for each " +
                    "catalogue genome (its largest chromosome); plasmids and profile-only uploads have no relatives. k is 1-100, default 10.")
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status404NotFound);

            sequences.MapGet("/{id:long}/signature", SignatureAsync)
                .WithName("GetSignature")
                .WithSummary("A record's dinucleotide signature, as Z-curve coefficients or as Karlin's per-dinucleotide values")
                .WithDescription("view is zcurve (default) or karlin: the same 16 numbers in two coordinate systems.")
                .ProducesProblem(StatusCodes.Status400BadRequest)
                .ProducesProblem(StatusCodes.Status404NotFound);

            app.MapGroup("/api/trees").WithTags("Organisms").MapGet("/", TreeAsync)
                .WithName("GetTree")
                .WithSummary("A tree of genomes from their dinucleotide signatures")
                .WithDescription("Leaves are either ids (comma-separated record ids) or every genome of one taxon (rank and taxon, as the catalogue " +
                    "filter). 3 to 500 leaves, all bacteria and archaea or all viruses. method is upgma (default) or nj. Newick labels are record ids.")
                .ProducesProblem(StatusCodes.Status400BadRequest);

            app.MapGroup("/api/landscape").WithTags("Organisms").MapGet("/", LandscapeAsync)
                .WithName("GetLandscape")
                .WithSummary("The public catalogue on two axes, for a map")
                .WithDescription("PCA of the dinucleotide signature. group is cellular (default), viral or eukaryotic, one genome per species, or all: every " +
                    "public genome on one set of axes. Two axes keep phyla apart but not genera, so this is a picture; nearest relatives come from the " +
                    "full signature. Private uploads are never on it.")
                .ProducesProblem(StatusCodes.Status400BadRequest);
        }

        private static async Task<Ok<LandscapeResponse>> LandscapeAsync(Landscape landscape, string group = SignatureIndex.Cellular,
            CancellationToken cancellationToken = default) =>
            TypedResults.Ok(LandscapeResponse.From(await landscape.GetAsync(group, cancellationToken)));

        private static async Task<Results<Ok<RelativesResponse>, ProblemHttpResult>> RelativesAsync(NpgsqlDataSource dataSource,
            SignatureIndex index, ClaimsPrincipal user, long id, int k = 10, CancellationToken cancellationToken = default)
        {
            if (k < 1 || k > MaxRelatives) throw new ArgumentException($"k must be 1 to {MaxRelatives}.");
            var query = await QueryAsync(dataSource, user, id, cancellationToken);
            if (query is null) return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var (group, counts) = query.Value;

            await index.RefreshAsync(cancellationToken);
            var relatives = index.Nearest(SignatureIndex.Signature(counts), group, user.UserId(), k, exclude: id);
            var records = await Catalogue.ByIdsAsync(dataSource, relatives.Select(r => r.Id).ToList(), cancellationToken);
            return TypedResults.Ok(new RelativesResponse(id, group, SignatureIndex.Description, relatives
                .Where(r => records.ContainsKey(r.Id))
                .Select(r => new RelativeResponse(r.Distance, SequenceSummaryResponse.From(records[r.Id]), $"/api/sequences/{r.Id}"))
                .ToList()));
        }

        private static async Task<Results<Ok<SignatureResponse>, ProblemHttpResult>> SignatureAsync(NpgsqlDataSource dataSource,
            ClaimsPrincipal user, long id, string view = "zcurve", CancellationToken cancellationToken = default)
        {
            var query = await QueryAsync(dataSource, user, id, cancellationToken);
            if (query is null) return TypedResults.Problem($"No sequence {id}.", statusCode: StatusCodes.Status404NotFound);
            var values = SignatureIndex.Views(query.Value.Counts, view);
            return TypedResults.Ok(new SignatureResponse(id, view.ToLowerInvariant(), SignatureIndex.Description,
                values.Select(SignatureValueResponse.From).ToList()));
        }

        // The record's group and counts; null when the caller may not see it. A plasmid or a profile-only upload is a 400.
        private static async Task<(string Group, long[] Counts)?> QueryAsync(NpgsqlDataSource dataSource, ClaimsPrincipal user, long id,
            CancellationToken cancellationToken)
        {
            if (!await Catalogue.IsVisibleAsync(dataSource, id, user.UserId(), cancellationToken)) return null;
            var (kind, counts) = (await SignatureIndex.ReadCountsAsync(dataSource, [id], cancellationToken))[id];
            return (Checked(id, kind, counts), counts!);
        }

        private static string Checked(long id, string kind, long[]? counts) =>
            SignatureIndex.GroupOf(kind) is not string group
                ? throw new ArgumentException($"Sequence {id} is a {kind}: signatures are compared per genome, and plasmids are left out.")
                : counts is null
                    ? throw new ArgumentException($"Sequence {id} was uploaded as a profile, so it has no dinucleotide counts to compare.")
                    : group;

        private static async Task<Ok<TreeResponse>> TreeAsync(NpgsqlDataSource dataSource, SignatureIndex index, ClaimsPrincipal user,
            string? ids, string? rank, string? taxon, string method = "upgma", CancellationToken cancellationToken = default)
        {
            List<long> leaves;
            if (ids is not null)
            {
                if (rank is not null || taxon is not null) throw new ArgumentException("Give ids or rank and taxon, not both.");
                leaves = ids.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                    .Select(text => long.TryParse(text, NumberStyles.None, CultureInfo.InvariantCulture, out long id)
                        ? id : throw new ArgumentException($"'{text}' is not a record id."))
                    .Distinct().ToList();
                CheckSize(leaves.Count);
                var visible = await Catalogue.VisibleAsync(dataSource, leaves, user.UserId(), cancellationToken);
                var missing = leaves.Where(id => !visible.Contains(id)).ToList();
                if (missing.Count > 0) throw new ArgumentException($"No sequence {string.Join(", ", missing)}.");
            }
            else
            {
                if (rank is null || taxon is null) throw new ArgumentException("Give ids, or rank and taxon.");
                var filter = new CatalogueFilter(Rank: rank, Taxon: taxon, ViewerId: user.UserId());
                Catalogue.Validate(filter);
                await index.RefreshAsync(cancellationToken);
                var matching = await Catalogue.IdsAsync(dataSource, filter, cancellationToken);
                leaves = matching.Where(index.Entries.ContainsKey).Order().ToList();
                CheckSize(leaves.Count);
            }

            var rows = await SignatureIndex.ReadCountsAsync(dataSource, leaves, cancellationToken);
            var groups = leaves.Select(id => Checked(id, rows[id].Kind, rows[id].Counts)).Distinct().ToList();
            if (groups.Count > 1) throw new ArgumentException("A tree holds bacteria and archaea, or viruses, not both.");
            var tree = SignatureIndex.Tree(
                leaves.Select(id => (id.ToString(CultureInfo.InvariantCulture), SignatureIndex.Signature(rows[id].Counts!))).ToList(), method);
            var records = await Catalogue.ByIdsAsync(dataSource, leaves, cancellationToken);
            return TypedResults.Ok(new TreeResponse(method.ToLowerInvariant(), SignatureIndex.Description, tree.ToNewick(),
                leaves.Select(id => SequenceSummaryResponse.From(records[id])).ToList()));
        }

        private static void CheckSize(int leaves)
        {
            if (leaves < 3 || leaves > SignatureIndex.MaxTreeLeaves)
                throw new ArgumentException($"A tree needs 3 to {SignatureIndex.MaxTreeLeaves} genomes; this one has {leaves}.");
        }
    }
}
