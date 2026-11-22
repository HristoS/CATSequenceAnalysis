using System.Globalization;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Experiments;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// seed-experiment-queries (migration 0024): stores the query set of an experiment whose queries come from outside the database, so the
    /// experiment worker can rerun it from the database alone. For protein-homolog: the proteins scripts/platform-checks/protein_homolog_check.cs
    /// draws (the protein_search_check_v2.cs rule: whole BV-BRC proteins, up to 5 per length band from each of 6 bacteria, 2 archaea,
    /// 2 plasmids and 3 viruses, seed 20261105), downloaded, profiled with QueryFromProtein and dropped: only the record, the gene's place,
    /// the protein's length and its expected axes are stored.
    /// </summary>
    internal static class ExperimentQuerySeeder
    {
        private const int seed = 20261105, perBand = 5;
        private static readonly Dictionary<string, int> genomesPerKind = new() { ["bacterium"] = 6, ["archaeon"] = 2, ["plasmid"] = 2, ["virus"] = 3 };
        private static readonly string[] mycoplasmaOrders = ["Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales"];

        public static async Task<int> SeedAsync(NpgsqlDataSource dataSource, string experiment, string? contactEmail)
        {
            if (experiment != ProteinHomolog.ExperimentName)
            {
                Console.WriteLine($"Only {ProteinHomolog.ExperimentName} has a query set to seed.");
                return 1;
            }
            using var http = new PoliteHttpClient("BV-BRC", contactEmail, TimeSpan.FromSeconds(1));
            async Task<JsonElement> GetAsync(string url)
            {
                using var response = await http.SendAsync(url, "application/json");
                using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
                return document.RootElement.Clone();
            }

            var queries = new List<ExperimentQuery>();
            foreach (var (kind, count) in genomesPerKind)
            {
                var samples = new List<(string Genome, long Id, string Record, string? Order)>();
                await using (var command = dataSource.CreateCommand(
                    "select s.accession, r.id, r.record_id, r.tax_order from sequence_record r join sequence_source s on s.id = r.source_id " +
                    "where r.kind = $1 and s.provider = 'BV-BRC' and s.owner_id is null and s.accession in (" +
                    "  select s2.accession from sequence_source s2 join sequence_record r2 on r2.source_id = s2.id " +
                    "  join gene_profile g on g.sequence_id = r2.id and g.bucket = 0 " +
                    "  where r2.kind = $1 group by s2.accession order by md5(s2.accession || $2) limit $3) and exists (" +
                    "  select 1 from gene_profile g where g.sequence_id = r.id)"))
                {
                    command.Parameters.AddWithValue(kind);
                    command.Parameters.AddWithValue(seed.ToString(CultureInfo.InvariantCulture));
                    command.Parameters.AddWithValue(count);
                    await using var reader = await command.ExecuteReaderAsync();
                    while (await reader.ReadAsync())
                        samples.Add((reader.GetString(0), reader.GetInt64(1), reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetString(3)));
                }

                foreach (var genome in samples.GroupBy(s => s.Genome))
                {
                    var records = genome.ToDictionary(s => s.Record);
                    var features = await GetAsync($"https://www.bv-brc.org/api/genome_feature/?and(eq(genome_id,{genome.Key}),eq(feature_type,CDS)," +
                        $"eq(annotation,PATRIC),ge(aa_length,{ZCurveGeneProfile.MinResidues}))&select(accession,start,end,strand,aa_length,aa_sequence_md5)&sort(+feature_id)&limit(25000)");
                    var picker = new Random(seed);
                    var picked = features.EnumerateArray()
                        .Where(f => f.TryGetProperty("aa_sequence_md5", out _) && records.ContainsKey(f.GetProperty("accession").GetString()!))
                        .GroupBy(f => ProteinHomolog.BandOf(f.GetProperty("aa_length").GetInt32()))
                        .SelectMany(band => band.OrderBy(_ => picker.Next()).Take(2 * perBand)).ToList();
                    if (picked.Count == 0) continue;
                    var proteins = new Dictionary<string, string>();
                    foreach (var chunk in picked.Select(f => f.GetProperty("aa_sequence_md5").GetString()!).Distinct().Chunk(50))
                        foreach (var p in (await GetAsync($"https://www.bv-brc.org/api/feature_sequence/?in(md5,({string.Join(",", chunk)}))&select(md5,sequence)&limit(1000)")).EnumerateArray())
                            proteins[p.GetProperty("md5").GetString()!] = p.GetProperty("sequence").GetString()!;
                    var taken = new int[5];
                    var profiles = new Dictionary<long, ZCurveGeneProfile>();
                    foreach (var feature in picked)
                    {
                        var record = records[feature.GetProperty("accession").GetString()!];
                        if (!proteins.TryGetValue(feature.GetProperty("aa_sequence_md5").GetString()!, out var protein) || protein.TrimEnd('*').Contains('*')) continue;
                        protein = protein.TrimEnd('*');
                        if (protein.Length < ZCurveGeneProfile.MinResidues) continue;
                        int band = ProteinHomolog.BandOf(protein.Length);
                        if (taken[band] == perBand) continue;
                        if (!profiles.TryGetValue(record.Id, out var profile)) profiles[record.Id] = profile = (await GeneProfiles.ReadAsync(dataSource, record.Id))!;
                        int start = feature.GetProperty("start").GetInt32() - 1, end = feature.GetProperty("end").GetInt32();
                        bool reverse = feature.GetProperty("strand").GetString() == "-";
                        if (!profile.Genes.Any(g => g.Start == start && g.End == end && g.Reverse == reverse)) continue;
                        var code = record.Order is not null && mycoplasmaOrders.Contains(record.Order) ? ZCurveGeneProfile.mycoplasmaCode : ZCurveGeneProfile.standardCode;
                        queries.Add(ProteinHomolog.ToQuery(record.Id, new ProteinHomolog.QueryFields(kind, protein.Length, start, end, reverse),
                            ZCurveGeneProfile.QueryFromProtein(protein, null, code)));
                        taken[band]++;
                    }
                    Console.WriteLine($"{kind} {genome.Key}: {string.Join(" / ", taken)} proteins by band");
                }
            }
            await ExperimentQueries.ReplaceAsync(dataSource, experiment, queries);
            Console.WriteLine($"Stored {queries.Count} queries for {experiment} (profiles only; the proteins were not kept).");
            return 0;
        }
    }
}
