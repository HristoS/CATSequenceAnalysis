using System.Text.Json;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>A taxonomic lineage, domain down to species; ranks the source doesn't give are absent.</summary>
    public sealed record LineageResponse(string? Superkingdom, string? Realm, string? Kingdom, string? Phylum, string? Class,
        string? Order, string? Family, string? Genus, string? Species)
    {
        internal static LineageResponse From(Lineage l) => new(l.Superkingdom, l.Realm, l.Kingdom, l.Phylum, l.Class, l.Order, l.Family, l.Genus, l.Species);
    }

    /// <summary>One catalogue record, as listed. A record holds no DNA: only its Z-curve profile levels are stored.</summary>
    /// <param name="Accession">The source's versioned accession; null for an upload.</param>
    /// <param name="RecordId">The FASTA record id within the source.</param>
    /// <param name="ProviderLineage">The lineage as the provider gave it; the detail adds NCBI's and GTDB's current ones.</param>
    public sealed record SequenceSummaryResponse(long Id, string Provider, string? Accession, string RecordId, string? Name, string Kind,
        long LengthBp, float? GcFraction, int? NcbiTaxid, LineageResponse ProviderLineage)
    {
        internal static SequenceSummaryResponse From(CatalogueRecord r) =>
            new(r.Id, r.Provider, r.SourceAccession, r.RecordId, r.Name, r.Kind, r.LengthBp, r.GcFraction, r.NcbiTaxid, LineageResponse.From(r.Lineage));
    }

    /// <summary>A page of the catalogue, in id order.</summary>
    /// <param name="NextAfter">Pass as <c>after</c> for the next page; absent on the last page.</param>
    public sealed record SequencePageResponse(IReadOnlyList<SequenceSummaryResponse> Items, long? NextAfter);

    /// <summary>Where a record's profile came from. The downloaded file was deleted after profiling; its hash and size remain.</summary>
    /// <param name="Url">Absent for an upload.</param>
    /// <param name="RawSha256">SHA-256 of the downloaded or uploaded file, hex.</param>
    public sealed record ProvenanceResponse(string? Url, DateTimeOffset RetrievedAt, string RawSha256, long RawByteSize, string? License);

    /// <summary>NCBI's current name and lineage of the record's taxid.</summary>
    public sealed record NcbiTaxonResponse(string? Name, LineageResponse Lineage);

    /// <summary>GTDB's lineage of the record's assembly; GTDB has no realm or kingdom, and its domain is given as superkingdom.</summary>
    public sealed record GtdbTaxonResponse(string Release, string Accession, LineageResponse Lineage);

    /// <summary>One stored profile level.</summary>
    /// <param name="MethodVersion">"1" as profiled; "1-origin" for a circular chromosome laid out from its replication origin.</param>
    /// <param name="Href">Where to download it as ZCWP bytes; absent for rotated levels, which only searches use.</param>
    public sealed record ProfileLevelResponse(string MethodVersion, int WindowSize, int WindowCount, long ByteSize, string? Href);

    /// <summary>Everything the platform knows about one record. It never stores the DNA.</summary>
    /// <param name="SequenceSha256">SHA-256 of the normalized sequence (uppercase letters only), hex: equal hashes mean identical sequences.</param>
    /// <param name="Owned">The record is your own upload; uploads are visible only to their owner.</param>
    /// <param name="Annotations">Free-form provider metadata (host, isolation source, assembly accession, ...).</param>
    public sealed record SequenceDetailResponse(SequenceSummaryResponse Sequence, string DataClass, string? SequenceSha256, bool Owned,
        ProvenanceResponse Provenance, NcbiTaxonResponse? Ncbi, GtdbTaxonResponse? Gtdb, IReadOnlyList<ProfileLevelResponse> Levels,
        JsonElement Annotations)
    {
        internal static SequenceDetailResponse From(CatalogueDetail d)
        {
            using var annotations = JsonDocument.Parse(d.AnnotationsJson);
            return new SequenceDetailResponse(
                SequenceSummaryResponse.From(d.Record), d.DataClass, d.SequenceSha256 is null ? null : Convert.ToHexStringLower(d.SequenceSha256), d.Owned,
                new ProvenanceResponse(d.SourceUrl, d.RetrievedAt, Convert.ToHexStringLower(d.RawSha256), d.RawByteSize, d.License),
                d.NcbiLineage is null ? null : new NcbiTaxonResponse(d.NcbiName, LineageResponse.From(d.NcbiLineage)),
                d.Gtdb is null ? null : new GtdbTaxonResponse(d.Gtdb.Release, d.Gtdb.Accession, LineageResponse.From(d.Gtdb.Lineage)),
                d.Levels.Select(l => new ProfileLevelResponse(l.MethodVersion, l.WindowSize, l.WindowCount, l.ByteSize,
                    l.MethodVersion == RecordStore.MethodVersion ? $"/api/sequences/{d.Record.Id}/profile?windowSize={l.WindowSize}" : null)).ToList(),
                annotations.RootElement.Clone());
        }
    }

    /// <summary>One record of an upload.</summary>
    /// <param name="RecordId">Its FASTA header id, or the recordId given with a profile.</param>
    /// <param name="WindowCount">Its 500 bp windows; 0 for a record shorter than one window, which is stored but can't be searched.</param>
    /// <param name="Href">The record in the catalogue, visible to you only.</param>
    public sealed record UploadedRecordResponse(long Id, string RecordId, long LengthBp, int WindowCount, string Href);

    /// <summary>The records an upload created, in file order.</summary>
    public sealed record UploadResponse(IReadOnlyList<UploadedRecordResponse> Records);
}
