namespace SequencePlatform.Api
{
    /// <summary>One stored level of a chromosome.</summary>
    /// <param name="Href">Its ZCWP bytes.</param>
    public sealed record ChromosomeLevelResponse(int WindowSize, int WindowCount, long ByteSize, string Href);

    /// <summary>One kept nuclear sequence of a eukaryote: its index (the key of <c>/chromosomes/{index}</c>), FASTA id and length.</summary>
    /// <param name="Name">The sequence's accession.version, as the genome file and the gene layer's sequences name it.</param>
    public sealed record ChromosomeResponse(int Index, string Name, long LengthBp, IReadOnlyList<ChromosomeLevelResponse> Levels);

    /// <summary>A eukaryote's stored chromosome levels, in the genome file's order; empty for any other record.</summary>
    public sealed record ChromosomeListResponse(long Id, IReadOnlyList<ChromosomeResponse> Chromosomes);
}
