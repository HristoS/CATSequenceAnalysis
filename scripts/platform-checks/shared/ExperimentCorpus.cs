using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequencePlatform.Core;

// The genomes every P9 experiment reads. Round 1 (S9.0) reads the frozen snapshot database; round 2 (S9.9, D-P9d) reads the live one
// cut off at RoundTwoMaxId and refuses to run unless the fingerprint of everything round 2 reads still matches. The selection
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
    public const long RoundTwoMaxId = 61_841;
    public const string RoundTwoArchiveFile = "sequence_platform_2026-11-12_P9b.dump";

    // Per-table hashes taken on 2026-11-12 (D-P9d). A mismatch means something wrote to data round 2 reads: stop, don't update these.
    private static readonly (string Table, string Hash)[] roundTwoExpected =
    [
        ("gene_profile", "1eb6e9538996100b8639bb1ea578e86e"),
        ("gtdb_taxon", "98f21921effdcb4289d36ec411e20f6f"),
        ("profile_level", "0b01cead1a7fb7016e22ca3d8011f62c"),
        ("protein_family", "7037cb9f7574bc6400b11ba1c097bf0a"),
        ("protein_kmer", "022b7e620f27714664a80ea082e8904c"),
        ("record_kmer", "cfbf87d97682469d3ef96e0be152177c"),
        ("sequence_record", "a4f60ef121ff51013df843db7ec21b7f"),
        ("taxon", "2468728558554be561b61a310b50b884"),
    ];

    // Every row round 2 can read, cut off by record id where the table is per record. Arrays and payloads are hashed in binary form
    // (6-mer counts are about 1 GB; text rendering doubled the time). Rows with nothing stored still count, through coalesce.
    private const string fingerprintSql =
        "select 'sequence_record', md5(string_agg(md5(row(r.id, r.source_id, r.record_id, r.name, r.kind, r.data_class, r.length_bp, " +
        "  r.gc_fraction, r.ncbi_taxid, r.tax_superkingdom, r.tax_realm, r.tax_kingdom, r.tax_phylum, r.tax_class, r.tax_order, r.tax_family, " +
        "  r.tax_genus, r.tax_species, r.annotations, r.sequence_sha256, r.dinucleotide_counts, s.owner_id)::text), '' order by r.id)) " +
        "  from sequence_record r join sequence_source s on s.id = r.source_id where r.id <= $1 " +
        "union all select 'record_kmer', md5(string_agg(md5(row(sequence_id, k, status)::text || coalesce(md5(array_send(counts)), '-')), '' " +
        "  order by sequence_id, k)) from record_kmer where sequence_id <= $1 " +
        "union all select 'protein_kmer', md5(string_agg(md5(row(sequence_id, k, format_version, gene_count, status, md5(payload))::text), '' " +
        "  order by sequence_id, k)) from protein_kmer where sequence_id <= $1 " +
        "union all select 'gene_profile', md5(string_agg(md5(row(sequence_id, bucket, format_version, gene_count, window_count, sha256)::text), '' " +
        "  order by sequence_id, bucket)) from gene_profile where sequence_id <= $1 " +
        "union all select 'profile_level', md5(string_agg(md5(row(sequence_id, method, method_version, fine_window_size, factor, window_count, " +
        "  format_version, sha256)::text), '' order by sequence_id, method, method_version, fine_window_size, factor)) " +
        "  from profile_level where sequence_id <= $1 " +
        "union all select 'protein_family', md5(string_agg(md5(row(id, product)::text), '' order by id)) from protein_family " +
        "union all select 'taxon', md5(string_agg(md5(t::text), '' order by taxid)) from taxon t " +
        "union all select 'gtdb_taxon', md5(string_agg(md5(g::text), '' order by assembly_key)) from gtdb_taxon g";

    private static readonly Regex partial = new(@"\bpartial\b", RegexOptions.IgnoreCase | RegexOptions.Compiled);

    /// <summary>The fingerprint every round-2 CSV carries: one hash over the expected table hashes.</summary>
    public static string RoundTwoFingerprint => Combine(roundTwoExpected);

    /// <summary>The live database (connection settings unchanged), for round 2.</summary>
    public static NpgsqlDataSource OpenLive() => NpgsqlDataSource.Create(ConnectionString.Resolve());

    /// <summary>Throws unless every table round 2 reads still hashes as on 2026-11-12. Takes about 40 s; start it beside <see cref="LoadAsync"/>.</summary>
    public static async Task<string> VerifyRoundTwoAsync(NpgsqlDataSource dataSource)
    {
        var actual = await FingerprintAsync(dataSource, RoundTwoMaxId);
        var moved = actual.Where(a => !roundTwoExpected.Contains(a)).Select(a => $"{a.Table} {a.Hash}").ToList();
        if (moved.Count > 0 || actual.Count != roundTwoExpected.Length)
            throw new InvalidOperationException($"Round-2 corpus moved (D-P9d). Tables now: {string.Join(", ", moved)}. " +
                $"Restore {RoundTwoArchiveFile} and run there, or start a new round.");
        return RoundTwoFingerprint;
    }

    public static async Task<List<(string Table, string Hash)>> FingerprintAsync(NpgsqlDataSource dataSource, long maxId)
    {
        var tables = new List<(string Table, string Hash)>();
        await using var command = dataSource.CreateCommand(fingerprintSql);
        command.Parameters.AddWithValue(maxId);
        command.CommandTimeout = 600;
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync()) tables.Add((reader.GetString(0), reader.IsDBNull(1) ? "empty" : reader.GetString(1)));
        return tables.OrderBy(t => t.Table, StringComparer.Ordinal).ToList();
    }

    public static string Combine(IEnumerable<(string Table, string Hash)> tables) =>
        Convert.ToHexStringLower(System.Security.Cryptography.MD5.HashData(System.Text.Encoding.UTF8.GetBytes(
            string.Join(",", tables.OrderBy(t => t.Table, StringComparer.Ordinal).Select(t => t.Table + t.Hash)))));

    /// <summary>The usual connection settings pointed at the snapshot database; it must have been restored from <see cref="SnapshotFile"/>.</summary>
    public static NpgsqlDataSource OpenSnapshot()
    {
        var builder = new NpgsqlConnectionStringBuilder(ConnectionString.Resolve()) { Database = SnapshotDatabase };
        return NpgsqlDataSource.Create(builder.ConnectionString);
    }

    /// <summary>Round 1 reads the snapshot with no cut-off; round 2 passes <see cref="RoundTwoMaxId"/>.</summary>
    public static async Task<List<CorpusGenome>> LoadAsync(NpgsqlDataSource dataSource, long maxId = long.MaxValue)
    {
        var genomes = new List<CorpusGenome>();
        await using var command = dataSource.CreateCommand(
            "select r.id, r.kind, coalesce(r.name, r.record_id), r.length_bp, r.dinucleotide_counts, r.annotations::text, " +
            "r.tax_superkingdom, r.tax_realm, r.tax_kingdom, r.tax_phylum, r.tax_class, r.tax_order, r.tax_family, r.tax_genus, r.tax_species, " +
            "t.taxid is not null, t.tax_superkingdom, t.tax_realm, t.tax_kingdom, t.tax_phylum, t.tax_class, t.tax_order, t.tax_family, t.tax_genus, t.tax_species, " +
            "g.assembly_key is not null, g.tax_domain, g.tax_phylum, g.tax_class, g.tax_order, g.tax_family, g.tax_genus, g.tax_species " +
            "from (select r.*, row_number() over (partition by r.source_id, r.kind in ('bacterium', 'archaeon') order by r.length_bp desc, r.id) as n " +
            "      from sequence_record r join sequence_source s on s.id = r.source_id" +
            "      where r.id <= $1 and s.owner_id is null and r.dinucleotide_counts is not null and r.kind in ('bacterium', 'archaeon', 'virus', 'eukaryote')) r " +
            "left join taxon t on t.taxid = r.ncbi_taxid " +
            "left join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
            "  and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
            "where r.kind = 'virus' or r.n = 1 order by r.id");
        command.Parameters.AddWithValue(maxId);
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
