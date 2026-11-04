using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequencePlatform.Core;

// The genomes every P9 experiment reads (S9.0), from the frozen snapshot so the catalogue growing can't move a result. The selection
// is SignatureIndex's, over public records only: the largest chromosome of each bacterium or archaeon, every virus and every eukaryote;
// no plasmids, no uploads. A record described as partial is marked Excluded and left out of every experiment. Most viral records
// under 1 kbp are complete segments or viroids, and "gene, complete cds" usually covers a whole small genome, so neither is excluded.
internal sealed record CorpusGenome(long Id, string Group, string Kind, string Name, long LengthBp, double Gc, long[] Counts,
    double[] Signature, Lineage Provider, Lineage? Ncbi, Lineage? Gtdb, IReadOnlyDictionary<string, string> Fields, bool Excluded)
{
    /// <summary>A BV-BRC annotation (host_name, isolation_source, checkm_contamination, ...), or the contig's "description"; null if empty.</summary>
    public string? Field(string name) => Fields.TryGetValue(name, out var value) && value != "" ? value : null;
}

internal static class ExperimentCorpus
{
    public const string SnapshotDatabase = "sequence_platform_p9";
    public const string SnapshotFile = "sequence_platform_2026-11-03_P9.dump";

    private static readonly Regex partial = new(@"\bpartial\b", RegexOptions.IgnoreCase | RegexOptions.Compiled);

    /// <summary>The usual connection settings pointed at the snapshot database; it must have been restored from <see cref="SnapshotFile"/>.</summary>
    public static NpgsqlDataSource OpenSnapshot()
    {
        var builder = new NpgsqlConnectionStringBuilder(ConnectionString.Resolve()) { Database = SnapshotDatabase };
        return NpgsqlDataSource.Create(builder.ConnectionString);
    }

    public static async Task<List<CorpusGenome>> LoadAsync(NpgsqlDataSource dataSource)
    {
        var genomes = new List<CorpusGenome>();
        await using var command = dataSource.CreateCommand(
            "select r.id, r.kind, coalesce(r.name, r.record_id), r.length_bp, r.dinucleotide_counts, r.annotations::text, " +
            "r.tax_superkingdom, r.tax_realm, r.tax_kingdom, r.tax_phylum, r.tax_class, r.tax_order, r.tax_family, r.tax_genus, r.tax_species, " +
            "t.taxid is not null, t.tax_superkingdom, t.tax_realm, t.tax_kingdom, t.tax_phylum, t.tax_class, t.tax_order, t.tax_family, t.tax_genus, t.tax_species, " +
            "g.assembly_key is not null, g.tax_domain, g.tax_phylum, g.tax_class, g.tax_order, g.tax_family, g.tax_genus, g.tax_species " +
            "from (select r.*, row_number() over (partition by r.source_id, r.kind in ('bacterium', 'archaeon') order by r.length_bp desc, r.id) as n " +
            "      from sequence_record r join sequence_source s on s.id = r.source_id" +
            "      where s.owner_id is null and r.dinucleotide_counts is not null and r.kind in ('bacterium', 'archaeon', 'virus', 'eukaryote')) r " +
            "left join taxon t on t.taxid = r.ncbi_taxid " +
            "left join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
            "  and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
            "where r.kind = 'virus' or r.n = 1 order by r.id");
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync())
        {
            string? Text(int i) => reader.IsDBNull(i) || reader.GetString(i) == "" ? null : reader.GetString(i);
            Lineage Full(int i) => new(Text(i), Text(i + 1), Text(i + 2), Text(i + 3), Text(i + 4), Text(i + 5), Text(i + 6), Text(i + 7), Text(i + 8));
            var kind = reader.GetString(1);
            var counts = reader.GetFieldValue<long[]>(4);
            var fields = Flatten(Text(5));
            var gtdb = reader.GetBoolean(25)
                ? new Lineage(Text(26), null, null, Text(27), Text(28), Text(29), Text(30), Text(31), Text(32)) : null;
            // GC from the counts (first base C or G: indexes 4..11), not gc_fraction, which some kinds don't store; one definition for all.
            double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();
            genomes.Add(new CorpusGenome(reader.GetInt64(0), SignatureIndex.GroupOf(kind)!, kind, reader.GetString(2), reader.GetInt64(3),
                gc, counts, SignatureIndex.Signature(counts), Full(6), reader.GetBoolean(15) ? Full(16) : null, gtdb, fields,
                fields.TryGetValue("description", out var description) && partial.IsMatch(description)));
        }
        return genomes;
    }

    // BV-BRC's fields as strings, plus the contig description; the experiments only ever read these as text.
    private static Dictionary<string, string> Flatten(string? annotations)
    {
        var fields = new Dictionary<string, string>(StringComparer.Ordinal);
        if (annotations is null) return fields;
        using var document = JsonDocument.Parse(annotations);
        if (document.RootElement.TryGetProperty("bvbrc", out var bvbrc))
            foreach (var property in bvbrc.EnumerateObject())
                fields[property.Name] = property.Value.ValueKind == JsonValueKind.String ? property.Value.GetString()! : property.Value.GetRawText();
        if (document.RootElement.TryGetProperty("contig", out var contig) && contig.TryGetProperty("description", out var description)
            && description.ValueKind == JsonValueKind.String)
            fields["description"] = description.GetString()!;
        return fields;
    }
}
