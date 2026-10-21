using System.IO.Compression;
using System.Text.Json;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <param name="CurrentTaxid">NCBI's id today; differs from <paramref name="QueriedTaxid"/> when it was merged, null when unknown.</param>
    internal sealed record NcbiTaxon(int QueriedTaxid, string Status, int? CurrentTaxid, string? Name, string? Rank, Lineage Lineage);

    internal sealed record GtdbTaxon(string AssemblyKey, string GtdbAccession, string? Domain, string? Phylum, string? Class,
        string? Order, string? Family, string? Genus, string? Species);

    /// <summary>
    /// S2.4 (migration 0006): fills <c>taxon</c> with NCBI's current lineage for every taxid the catalogue references,
    /// and <c>gtdb_taxon</c> with GTDB's lineage for every referenced bacterial or archaeal assembly. Only referenced
    /// taxa are stored: about 16k rows instead of the 2.6M of the full NCBI dump, at the price of fetching the new
    /// taxids after each load (a plain re-run does that). Nothing in sequence_record is changed.
    /// </summary>
    internal sealed class TaxonomyLoader : IDisposable
    {
        public const string Current = "current";
        public const string Merged = "merged";
        public const string Unknown = "unknown";
        private const string datasetsBase = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/";
        private const string gtdbBase = "https://data.gtdb.ecogenomic.org/releases/latest/";
        private const int taxaPerRequest = 200;
        private static readonly string[] gtdbFiles = { "bac120_taxonomy.tsv.gz", "ar53_taxonomy.tsv.gz" };
        private static readonly string[] ncbiRanks = { "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species" };

        private readonly NpgsqlDataSource dataSource;
        // NCBI allows 5 requests a second without an API key; GTDB is three small static files.
        private readonly PoliteHttpClient http;

        public TaxonomyLoader(NpgsqlDataSource dataSource, string? contactEmail)
        {
            this.dataSource = dataSource;
            http = new PoliteHttpClient("NCBI/GTDB", contactEmail, TimeSpan.FromSeconds(1));
        }

        /// <param name="refresh">Re-fetch every referenced taxid, not only those missing from taxon.</param>
        /// <returns>Taxids asked for, and how many came back current, merged and unknown; the rest got no answer.</returns>
        public async Task<(int Asked, int Current, int Merged, int Unknown)> LoadNcbiAsync(bool refresh)
        {
            var taxids = new List<int>();
            await using (var command = dataSource.CreateCommand(
                "select distinct ncbi_taxid from sequence_record where ncbi_taxid is not null " +
                "and ($1 or ncbi_taxid not in (select taxid from taxon)) order by 1"))
            {
                command.Parameters.AddWithValue(refresh);
                await using var reader = await command.ExecuteReaderAsync();
                while (await reader.ReadAsync()) taxids.Add(reader.GetInt32(0));
            }

            int current = 0, merged = 0, unknown = 0;
            for (int i = 0; i < taxids.Count; i += taxaPerRequest)
            {
                var chunk = taxids.Skip(i).Take(taxaPerRequest).ToList();
                using var response = await http.SendAsync(
                    $"{datasetsBase}{string.Join(",", chunk)}/dataset_report?page_size={taxaPerRequest * 2}", "application/json");
                var taxa = ParseDatasetReport(await response.Content.ReadAsStringAsync()).Where(t => chunk.Contains(t.QueriedTaxid)).ToList();
                await UpsertTaxaAsync(dataSource, taxa);
                current += taxa.Count(t => t.Status == Current);
                merged += taxa.Count(t => t.Status == Merged);
                unknown += taxa.Count(t => t.Status == Unknown);
                Console.WriteLine($"  NCBI: {Math.Min(i + taxaPerRequest, taxids.Count),6:N0} of {taxids.Count:N0} taxids");
            }
            return (taxids.Count, current, merged, unknown);
        }

        // The Datasets taxonomy report answers each queried id once: with the taxon it resolves to today (a merged id
        // resolves to the taxon it was merged into, listing the old id in "query"), or with an error for an unknown id.
        internal static IEnumerable<NcbiTaxon> ParseDatasetReport(string json)
        {
            using var document = JsonDocument.Parse(json);
            if (!document.RootElement.TryGetProperty("reports", out var reports)) yield break;
            foreach (var report in reports.EnumerateArray())
            {
                if (!report.TryGetProperty("query", out var query) || query.GetArrayLength() == 0
                    || !int.TryParse(query[0].GetString(), out var queried))
                    continue;
                if (!report.TryGetProperty("taxonomy", out var taxonomy))
                {
                    yield return new NcbiTaxon(queried, Unknown, null, null, null, new Lineage(null, null, null, null, null, null, null, null, null));
                    continue;
                }

                var classification = taxonomy.TryGetProperty("classification", out var c) ? c : default;
                string? Name(string rank) =>
                    classification.ValueKind == JsonValueKind.Object && classification.TryGetProperty(rank, out var r)
                    && r.TryGetProperty("name", out var n) ? n.GetString() : null;
                var ranks = ncbiRanks.Select(Name).ToArray();
                var taxid = taxonomy.GetProperty("tax_id").GetInt32();
                var name = taxonomy.TryGetProperty("current_scientific_name", out var scientific) && scientific.TryGetProperty("name", out var sn)
                    ? sn.GetString() : null;
                var rank = taxonomy.TryGetProperty("rank", out var rk) ? rk.GetString()?.ToLowerInvariant().Replace('_', ' ') : null;
                yield return new NcbiTaxon(queried, taxid == queried ? Current : Merged, taxid, name, rank,
                    new Lineage(Name("domain") ?? Name("acellular_root"), ranks[0], ranks[1], ranks[2], ranks[3], ranks[4], ranks[5], ranks[6], ranks[7]));
            }
        }

        internal static async Task UpsertTaxaAsync(NpgsqlDataSource dataSource, IReadOnlyList<NcbiTaxon> taxa)
        {
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var t in taxa)
            {
                await using var command = new NpgsqlCommand(
                    "insert into taxon (taxid, status, current_taxid, name, rank, tax_superkingdom, tax_realm, tax_kingdom, tax_phylum, " +
                    "tax_class, tax_order, tax_family, tax_genus, tax_species) values ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14) " +
                    "on conflict (taxid) do update set status = excluded.status, current_taxid = excluded.current_taxid, name = excluded.name, " +
                    "rank = excluded.rank, tax_superkingdom = excluded.tax_superkingdom, tax_realm = excluded.tax_realm, " +
                    "tax_kingdom = excluded.tax_kingdom, tax_phylum = excluded.tax_phylum, tax_class = excluded.tax_class, " +
                    "tax_order = excluded.tax_order, tax_family = excluded.tax_family, tax_genus = excluded.tax_genus, " +
                    "tax_species = excluded.tax_species, fetched_at = now()", connection, transaction);
                var l = t.Lineage;
                foreach (var value in new object?[] { t.QueriedTaxid, t.Status, t.CurrentTaxid, t.Name, t.Rank, l.Superkingdom, l.Realm,
                             l.Kingdom, l.Phylum, l.Class, l.Order, l.Family, l.Genus, l.Species })
                    command.Parameters.AddWithValue(value ?? DBNull.Value);
                await command.ExecuteNonQueryAsync();
            }
            await transaction.CommitAsync();
        }

        /// <returns>The GTDB release, the referenced assemblies, and how many of them GTDB lists.</returns>
        public async Task<(string Release, int Referenced, int Matched)> LoadGtdbAsync()
        {
            var referenced = new HashSet<string>(StringComparer.Ordinal);
            await using (var command = dataSource.CreateCommand(
                "select distinct annotations->'bvbrc'->>'assembly_accession' from sequence_record " +
                "where tax_superkingdom in ('Bacteria', 'Archaea') and annotations->'bvbrc'->>'assembly_accession' is not null"))
            await using (var reader = await command.ExecuteReaderAsync())
                while (await reader.ReadAsync())
                    if (AssemblyKey(reader.GetString(0)) is { } key) referenced.Add(key);

            string release;
            using (var version = await http.SendAsync($"{gtdbBase}VERSION.txt", "text/plain"))
                release = (await version.Content.ReadAsStringAsync()).Split('\n')[0].Trim();

            var matched = new List<GtdbTaxon>();
            foreach (var file in gtdbFiles)
            {
                using var response = await http.SendAsync($"{gtdbBase}{file}", "application/gzip", HttpCompletionOption.ResponseHeadersRead);
                await using var gzip = new GZipStream(await response.Content.ReadAsStreamAsync(), CompressionMode.Decompress);
                using var reader = new StreamReader(gzip);
                string? line;
                while ((line = await reader.ReadLineAsync()) != null)
                    if (ParseGtdbLine(line) is { } taxon && referenced.Contains(taxon.AssemblyKey))
                        matched.Add(taxon);
            }

            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var t in matched.DistinctBy(t => t.AssemblyKey))
            {
                await using var command = new NpgsqlCommand(
                    "insert into gtdb_taxon (assembly_key, gtdb_accession, release, tax_domain, tax_phylum, tax_class, tax_order, tax_family, " +
                    "tax_genus, tax_species) values ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10) on conflict (assembly_key) do update set " +
                    "gtdb_accession = excluded.gtdb_accession, release = excluded.release, tax_domain = excluded.tax_domain, " +
                    "tax_phylum = excluded.tax_phylum, tax_class = excluded.tax_class, tax_order = excluded.tax_order, " +
                    "tax_family = excluded.tax_family, tax_genus = excluded.tax_genus, tax_species = excluded.tax_species, fetched_at = now()",
                    connection, transaction);
                foreach (var value in new object?[] { t.AssemblyKey, t.GtdbAccession, release, t.Domain, t.Phylum, t.Class, t.Order, t.Family, t.Genus, t.Species })
                    command.Parameters.AddWithValue(value ?? DBNull.Value);
                await command.ExecuteNonQueryAsync();
            }
            await transaction.CommitAsync();
            return (release, referenced.Count, matched.Select(t => t.AssemblyKey).Distinct().Count());
        }

        // GCA_000005845.2 and GCF_000005845.2 -> 000005845.2. Null for anything that isn't an assembly accession.
        internal static string? AssemblyKey(string accession) =>
            accession.Length > 4 && accession.StartsWith("GC", StringComparison.Ordinal) && accession[3] == '_' ? accession[4..] : null;

        // "RS_GCF_000005845.2<TAB>d__Bacteria;p__Pseudomonadota;...;s__Escherichia coli"; an empty rank ("g__") is null.
        internal static GtdbTaxon? ParseGtdbLine(string line)
        {
            var tab = line.IndexOf('\t');
            if (tab < 4) return null;
            var accession = line[..tab];
            if (AssemblyKey(accession[3..]) is not { } key) return null;

            var ranks = new Dictionary<char, string?>();
            foreach (var part in line[(tab + 1)..].Split(';'))
                if (part.Length >= 3 && part[1] == '_' && part[2] == '_')
                    ranks[part[0]] = part.Length > 3 ? part[3..] : null;
            string? Rank(char prefix) => ranks.GetValueOrDefault(prefix);
            return new GtdbTaxon(key, accession, Rank('d'), Rank('p'), Rank('c'), Rank('o'), Rank('f'), Rank('g'), Rank('s'));
        }

        public void Dispose() => http.Dispose();
    }
}
