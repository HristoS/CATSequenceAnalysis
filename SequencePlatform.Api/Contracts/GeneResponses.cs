using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>One stored gene of a record (the gene layer, P7): where it is, its protein's length and its BV-BRC family. No bases.</summary>
    /// <param name="Index">Its position in the record's gene list, in order of start; the key of <c>/genes/{index}</c>.</param>
    /// <param name="StartBp">0-based on the record's stored strand.</param>
    /// <param name="EndBp">Exclusive; the stop codon is included when the annotation has one.</param>
    /// <param name="Residues">The protein's length in amino acids, without the stop.</param>
    /// <param name="Windows">How many windows it has: 5 for 50–499 aa, one per 333 codon positions 1 and 2 from 500 aa.</param>
    /// <param name="Sequence">The record's sequence the gene lies on, an index into <see cref="GeneListResponse.Sequences"/>; always 0 for a
    /// single-sequence record. StartBp and EndBp are on that sequence.</param>
    public sealed record GeneResponse(int Index, int StartBp, int EndBp, Strand Strand, int Residues, int Windows, string? FamilyId, string? Product,
        int Sequence = 0);

    /// <summary>A record's genes and its mean gene window: the 9 window features averaged over every window of every gene.</summary>
    /// <param name="MeanWindow">Mean, trend and curvature of x, y and z over codon positions 1 and 2; in the dinucleotide signature's format.</param>
    /// <param name="Sequences">For an NCBI eukaryote, the accessions of the sequences its genes lie on, in assembly order; empty for a
    /// single-sequence record.</param>
    public sealed record GeneListResponse(long Id, int GeneCount, int WindowCount, IReadOnlyList<SignatureValueResponse> MeanWindow,
        IReadOnlyList<GeneResponse> Genes, IReadOnlyList<string> Sequences);

    /// <summary>One gene's stored windows, from which its Z-curve over codon positions 1 and 2 can be redrawn without any bases.</summary>
    /// <param name="WindowPositions">Codon positions 1 and 2 per window.</param>
    /// <param name="Windows">Per window, 9 values in <see cref="GeneProfiles.Features"/> order: x, y, z mean, then trend, then curvature,
    /// each a Legendre moment scaled by sqrt(2k + 1). Inside a window of W positions, the running sum of an axis at t in [-1, 1] is
    /// W / 2 * (m0 (t + 1) + sqrt(3) m1 (t^2 - 1) / 2 + sqrt(5) m2 (t^3 - t) / 2).</param>
    public sealed record GeneWindowsResponse(long Id, GeneResponse Gene, int WindowPositions, IReadOnlyList<IReadOnlyList<double>> Windows);

    /// <summary>
    /// A record's Z-curve path built from its gene layer alone (<see cref="GenePath"/>): genes in chromosome order from <c>StartBp</c>,
    /// each window one step of window positions times its means over codon positions 1 and 2, reverse-strand genes flipped onto the
    /// forward strand; intergenic DNA left out.
    /// </summary>
    /// <param name="Coordinates">The points flattened, x, y, z per point, rounded to 0.1 base; the first point is the axes' origin.</param>
    /// <param name="GeneIndex">The genes in path order, by their index in <c>/genes</c>.</param>
    /// <param name="FirstPoint">Where each of those genes' stretch starts; it ends where the next one starts.</param>
    /// <param name="Sequence">The one sequence drawn, when asked for; absent when the path runs over all of them.</param>
    /// <param name="GeneCount">Genes on the path: the record's, or the asked sequence's.</param>
    public sealed record GenePathResponse(long Id, long StartBp, int GeneCount, IReadOnlyList<double> Coordinates, IReadOnlyList<int> GeneIndex,
        IReadOnlyList<int> FirstPoint, int? Sequence = null);

    /// <summary>One of a gene's best matches among another record's genes (<see cref="GeneMatching"/>).</summary>
    /// <param name="Score">1 - the mean per-window distance between the two genes' stored windows; 1 is identical.</param>
    /// <param name="Reciprocal">A reciprocal best hit: the first match, and its own best match back in the first record is the queried gene.</param>
    public sealed record GeneMatchResponse(GeneResponse Gene, double Score, bool Reciprocal);

    /// <summary>A gene's best matches in record <c>In</c>, best first.</summary>
    /// <param name="Compared">How many of <c>In</c>'s genes were comparable: from 500 aa its long genes, below that those within ±10% of the length.</param>
    public sealed record GeneMatchListResponse(long Id, GeneResponse Gene, long In, int Compared, IReadOnlyList<GeneMatchResponse> Matches);

    /// <summary>A protein search's query as the windows it was compared with: the protein's expected Z-curve, never its letters.</summary>
    /// <param name="Windows">Per window, the same 9 values as <see cref="GeneWindowsResponse.Windows"/>.</param>
    public sealed record ProteinQueryWindowsResponse(Guid SearchId, int LengthAa, int WindowPositions, IReadOnlyList<IReadOnlyList<double>> Windows);
}
