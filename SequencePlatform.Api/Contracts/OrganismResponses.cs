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
}
