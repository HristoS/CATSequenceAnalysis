using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>One stored gene of a record (the gene layer, P7): where it is, its protein's length and its BV-BRC family. No bases.</summary>
    /// <param name="Index">Its position in the record's gene list, in order of start; the key of <c>/genes/{index}</c>.</param>
    /// <param name="StartBp">0-based on the record's stored strand.</param>
    /// <param name="EndBp">Exclusive; the stop codon is included when the annotation has one.</param>
    /// <param name="Residues">The protein's length in amino acids, without the stop.</param>
    /// <param name="Windows">How many windows it has: 5 for 50–499 aa, one per 333 codon positions 1 and 2 from 500 aa.</param>
    public sealed record GeneResponse(int Index, int StartBp, int EndBp, Strand Strand, int Residues, int Windows, string? FamilyId, string? Product);

    /// <summary>A record's genes and its mean gene window: the 9 window features averaged over every window of every gene.</summary>
    /// <param name="MeanWindow">Mean, trend and curvature of x, y and z over codon positions 1 and 2; in the dinucleotide signature's format.</param>
    public sealed record GeneListResponse(long Id, int GeneCount, int WindowCount, IReadOnlyList<SignatureValueResponse> MeanWindow,
        IReadOnlyList<GeneResponse> Genes);

    /// <summary>One gene's stored windows, from which its Z-curve over codon positions 1 and 2 can be redrawn without any bases.</summary>
    /// <param name="WindowPositions">Codon positions 1 and 2 per window.</param>
    /// <param name="Windows">Per window, 9 values in <see cref="GeneProfiles.Features"/> order: x, y, z mean, then trend, then curvature,
    /// each a Legendre moment scaled by sqrt(2k + 1). Inside a window of W positions, the running sum of an axis at t in [-1, 1] is
    /// W / 2 * (m0 (t + 1) + sqrt(3) m1 (t^2 - 1) / 2 + sqrt(5) m2 (t^3 - t) / 2).</param>
    public sealed record GeneWindowsResponse(long Id, GeneResponse Gene, int WindowPositions, IReadOnlyList<IReadOnlyList<double>> Windows);

    /// <summary>A protein search's query as the windows it was compared with: the protein's expected Z-curve, never its letters.</summary>
    /// <param name="Windows">Per window, the same 9 values as <see cref="GeneWindowsResponse.Windows"/>.</param>
    public sealed record ProteinQueryWindowsResponse(Guid SearchId, int LengthAa, int WindowPositions, IReadOnlyList<IReadOnlyList<double>> Windows);
}
