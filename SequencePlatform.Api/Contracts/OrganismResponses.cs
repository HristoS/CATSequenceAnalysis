using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>One relative: a genome and how far its dinucleotide signature is from the query's.</summary>
    public sealed record RelativeResponse(double Distance, SequenceSummaryResponse Sequence, string Href);

    /// <summary>A record's nearest genomes by dinucleotide signature, nearest first.</summary>
    /// <param name="Group">cellular (bacteria and archaea) or viral: relatives come from the query's own group.</param>
    /// <param name="Signature">What was compared.</param>
    public sealed record RelativesResponse(long Id, string Group, string Signature, IReadOnlyList<RelativeResponse> Items);

    /// <summary>One value of a signature view.</summary>
    /// <param name="Ratio">Karlin view only: rho* itself, observed / expected (below 1: the pair is avoided).</param>
    public sealed record SignatureValueResponse(string Name, double Value, string Meaning, double? Ratio)
    {
        internal static SignatureValueResponse From(SignatureValue v) => new(v.Name, v.Value, v.Meaning, v.Ratio);
    }

    /// <summary>
    /// A record's dinucleotide signature. The two views are the same 16 numbers in two coordinate systems: zcurve gives the
    /// coefficients on the Z-curve axes of neighbouring bases (the complete second-order Z-curve), karlin the per-dinucleotide log rho*.
    /// </summary>
    public sealed record SignatureResponse(long Id, string View, string Signature, IReadOnlyList<SignatureValueResponse> Values);

    /// <summary>A tree of genomes built from their dinucleotide signatures.</summary>
    /// <param name="Newick">Leaf labels are record ids.</param>
    /// <param name="Leaves">The records at the leaves, in id order.</param>
    public sealed record TreeResponse(string Method, string Signature, string Newick, IReadOnlyList<SequenceSummaryResponse> Leaves);

    /// <summary>One dinucleotide's weight on a map axis.</summary>
    public sealed record LandscapeLoadingResponse(string Dinucleotide, double Weight);

    /// <summary>One axis of the map.</summary>
    /// <param name="ExplainedVariance">The axis's share of the signatures' total variance.</param>
    /// <param name="TopLoadings">The dinucleotides that weigh most on it, heaviest first.</param>
    public sealed record LandscapeAxisResponse(double ExplainedVariance, IReadOnlyList<LandscapeLoadingResponse> TopLoadings);

    /// <summary>
    /// The public catalogue on two axes, one genome per species. Parallel arrays with one entry per genome, because a map holds thousands;
    /// a taxon is an index into its list, -1 when the genome has none.
    /// </summary>
    /// <param name="Group">all, cellular (bacteria and archaea) or viral.</param>
    /// <param name="Signature">What was mapped.</param>
    public sealed record LandscapeResponse(string Group, string Signature, long[] Ids, string[] Names, double[] X, double[] Y, double[] Gc, long[] LengthBp,
        string[] Domains, int[] Domain, string[] Phyla, int[] Phylum, string[] Families, int[] Family, string[] Genera, int[] Genus,
        IReadOnlyList<LandscapeAxisResponse> Axes)
    {
        private const int coordinateDecimals = 5;
        private const int gcDecimals = 4;

        internal static LandscapeResponse From(LandscapeMap map)
        {
            var points = map.Points;
            var (domains, domain) = Indexed(points.Select(p => p.Domain));
            var (phyla, phylum) = Indexed(points.Select(p => p.Phylum));
            var (families, family) = Indexed(points.Select(p => p.Family));
            var (genera, genus) = Indexed(points.Select(p => p.Genus));
            return new LandscapeResponse(map.Group, SignatureIndex.Description,
                points.Select(p => p.Id).ToArray(), points.Select(p => p.Name).ToArray(),
                points.Select(p => Math.Round(p.X, coordinateDecimals)).ToArray(), points.Select(p => Math.Round(p.Y, coordinateDecimals)).ToArray(),
                points.Select(p => Math.Round(p.Gc, gcDecimals)).ToArray(), points.Select(p => p.LengthBp).ToArray(),
                domains, domain, phyla, phylum, families, family, genera, genus,
                map.Axes.Select(a => new LandscapeAxisResponse(a.ExplainedVariance,
                    a.TopLoadings.Select(l => new LandscapeLoadingResponse(l.Dinucleotide, l.Weight)).ToList())).ToList());
        }

        private static (string[] Names, int[] Indexes) Indexed(IEnumerable<string?> taxa)
        {
            var values = taxa.ToList();
            var names = values.OfType<string>().Distinct().Order(StringComparer.Ordinal).ToArray();
            var positions = names.Select((name, i) => (name, i)).ToDictionary(x => x.name, x => x.i);
            return (names, values.Select(t => t is null ? -1 : positions[t]).ToArray());
        }
    }
}
