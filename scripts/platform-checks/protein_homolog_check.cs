#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Protein homolog check (2026-11-12), item 1 of the experiment plan: does the right protein family come first without the query's own
// gene, genome, species or genus? protein_search_check_v2.cs counted a family hit whenever the query's own gene ranked first, so its
// "family top-1 100%" may rest on the self-hit. Queries: the same rule and seed as v2 (whole BV-BRC proteins, up to 5 per length band
// from each of 6 bacteria, 2 archaea, 2 plasmids, 3 viruses, seed 20261105) on today's catalogue, compared as ProteinSearch compares
// them. Per query and exclusion tier, over the genes left: the best hit's PGfam, and whether the query's PGfam is in the top 10. Tiers:
// gene (only the exact gene removed), genome (every record of the query's source), species, genus (records with the same NCBI
// tax_species / tax_genus; a query without that name is left out of the tier). A query counts in a tier only if some gene of its
// family is left there ("reachable"). Viruses have no PGfam and drop out.
// Rule fixed before the run: if family top-1 with the own species removed is under 50% of reachable queries in any band, the protein
// claim becomes "fast lookup of near-identical genes from a protein" and G1 compares against exact/MinHash-style lookup, not BLASTP
// sensitivity; otherwise the identity curve (the full item 1 run) goes ahead.
// Args: output path (default ProteinHomolog_Check.csv).
const int seed = 20261105, perBand = 5, topK = 10;
string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
string[] tierNames = ["gene", "genome", "species", "genus"];
static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 6, ["archaeon"] = 2, ["plasmid"] = 2, ["virus"] = 3 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
var output = args.Length > 0 ? args[0] : "ProteinHomolog_Check.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
using var http = new HttpClient { Timeout = TimeSpan.FromMinutes(5) };
http.DefaultRequestHeaders.UserAgent.ParseAdd("CATSequenceAnalysis-platform-check/1.0");
http.DefaultRequestHeaders.Accept.ParseAdd("application/json");

async Task<JsonElement> GetAsync(string url)
{
    await Task.Delay(1000);
    using var response = await http.GetAsync(url);
    response.EnsureSuccessStatusCode();
    using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
    return document.RootElement.Clone();
}

// Every record's source and names, for the exclusion tiers.
var taxa = new Dictionary<long, Taxon>();
await using (var command = dataSource.CreateCommand("select id, source_id, tax_species, tax_genus from sequence_record"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
        taxa[reader.GetInt64(0)] = new Taxon(reader.GetInt64(1), reader.IsDBNull(2) ? null : reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetString(3));

// Query selection: protein_search_check_v2.cs unchanged.
var queries = new List<Query>();
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
            .GroupBy(f => BandOf(f.GetProperty("aa_length").GetInt32()))
            .SelectMany(band => band.OrderBy(_ => picker.Next()).Take(2 * perBand)).ToList();
        if (picked.Count == 0) continue;
        var proteins = new Dictionary<string, string>();
        foreach (var chunk in picked.Select(f => f.GetProperty("aa_sequence_md5").GetString()!).Distinct().Chunk(50))
            foreach (var p in (await GetAsync($"https://www.bv-brc.org/api/feature_sequence/?in(md5,({string.Join(",", chunk)}))&select(md5,sequence)&limit(1000)")).EnumerateArray())
                proteins[p.GetProperty("md5").GetString()!] = p.GetProperty("sequence").GetString()!;
        var taken = new int[bandNames.Length];
        var profiles = new Dictionary<long, ZCurveGeneProfile>();
        foreach (var feature in picked)
        {
            var record = records[feature.GetProperty("accession").GetString()!];
            if (!proteins.TryGetValue(feature.GetProperty("aa_sequence_md5").GetString()!, out var protein) || protein.TrimEnd('*').Contains('*')) continue;
            protein = protein.TrimEnd('*');
            if (protein.Length < ZCurveGeneProfile.MinResidues) continue;
            int band = BandOf(protein.Length);
            if (taken[band] == perBand) continue;
            if (!profiles.TryGetValue(record.Id, out var profile)) profiles[record.Id] = profile = (await GeneProfiles.ReadAsync(dataSource, record.Id))!;
            int start = feature.GetProperty("start").GetInt32() - 1, end = feature.GetProperty("end").GetInt32();
            bool reverse = feature.GetProperty("strand").GetString() == "-";
            int gene = Enumerable.Range(0, profile.Genes.Count).FirstOrDefault(g => profile.Genes[g].Start == start && profile.Genes[g].End == end
                && profile.Genes[g].Reverse == reverse, -1);
            if (gene < 0) continue;
            var code = record.Order is not null && mycoplasmaOrders.Contains(record.Order) ? ZCurveGeneProfile.mycoplasmaCode : ZCurveGeneProfile.standardCode;
            var axes = ZCurveGeneProfile.QueryFromProtein(protein, null, code);
            queries.Add(new Query(kind, band, protein, record.Id, taxa[record.Id], profile.Genes[gene], axes,
                GeneProfiles.BucketsFor(protein.Length, ProteinSearch.LengthBand).ToHashSet()));
            taken[band]++;
        }
        Console.WriteLine($"{kind} {genome.Key}: {string.Join(" / ", taken)} proteins by band");
    }
}
Console.WriteLine($"{queries.Count} proteins, {queries.Count(q => q.Truth.Family != ZCurveGeneProfile.NoFamily)} with a PGfam.");

// Which tiers a candidate record stays in for a query: tier t keeps it if it isn't excluded at t.
bool Kept(Query query, Taxon candidate, int tier) => tier switch
{
    0 => true,
    1 => candidate.Source != query.Taxon.Source,
    2 => candidate.Source != query.Taxon.Source && (candidate.Species is null || candidate.Species != query.Taxon.Species),
    _ => candidate.Source != query.Taxon.Source && (candidate.Genus is null || candidate.Genus != query.Taxon.Genus),
};
bool TierApplies(Query query, int tier) => tier switch { 2 => query.Taxon.Species is not null, 3 => query.Taxon.Genus is not null, _ => true };
bool Comparable(Query query, ZCurveGeneProfile.Gene gene) => query.Protein.Length >= ZCurveGeneProfile.LongResidues
    ? gene.Residues >= ZCurveGeneProfile.LongResidues
    : gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - query.Protein.Length) <= ProteinSearch.LengthBand * query.Protein.Length;

var scoredQueries = queries.Where(q => q.Truth.Family != ZCurveGeneProfile.NoFamily).ToList();
var needed = scoredQueries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();
var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
var watch = Stopwatch.StartNew();
var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
{
    var tallies = scoredQueries.Select(_ => tierNames.Select(_ => new Tally()).ToArray()).ToArray();
    foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable())
    {
        var candidate = taxa[id];
        var profile = ZCurveGeneProfile.FromBinary(payload);
        for (int q = 0; q < scoredQueries.Count; q++)
        {
            var query = scoredQueries[q];
            if (!query.Buckets.Contains(bucket)) continue;
            for (int g = 0; g < profile.Genes.Count; g++)
            {
                var gene = profile.Genes[g];
                if (!Comparable(query, gene) || (id == query.SequenceId && gene.Start == query.Truth.Start && gene.End == query.Truth.End && gene.Reverse == query.Truth.Reverse)) continue;
                double score = profile.Score(query.Axes, g);
                for (int t = 0; t < tierNames.Length; t++)
                {
                    if (!Kept(query, candidate, t)) continue;
                    var tally = tallies[q][t];
                    if (gene.Family == query.Truth.Family) tally.SameFamily++;
                    if (tally.Top.Count < topK) tally.Top.Enqueue(gene.Family, score);
                    else if (tally.Top.TryPeek(out _, out double lowest) && score > lowest) tally.Top.EnqueueDequeue(gene.Family, score);
                }
            }
        }
    }
    return tallies;
})).ToArray();
await using (var command = dataSource.CreateCommand("select sequence_id, bucket, payload from gene_profile where bucket = any($1) and window_count > 0"))
{
    command.Parameters.AddWithValue(needed);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) rows.Add((reader.GetInt64(0), reader.GetInt16(1), reader.GetFieldValue<byte[]>(2)));
}
rows.CompleteAdding();
var perWorker = await Task.WhenAll(workers);
Console.WriteLine($"Scoring pass: {watch.Elapsed.TotalSeconds:F0} s.");

// Per query and tier: merged top 10 over the workers, and whether any same-family gene was left.
var outcomes = new List<(Query Query, int Tier, bool Reachable, bool Top1, bool Top10)>();
for (int q = 0; q < scoredQueries.Count; q++)
    for (int t = 0; t < tierNames.Length; t++)
    {
        if (!TierApplies(scoredQueries[q], t)) continue;
        var merged = perWorker.SelectMany(w => w[q][t].Top.UnorderedItems).OrderByDescending(e => e.Priority).Take(topK).ToList();
        int family = scoredQueries[q].Truth.Family;
        outcomes.Add((scoredQueries[q], t, perWorker.Sum(w => w[q][t].SameFamily) > 0,
            merged.Count > 0 && merged[0].Element == family, merged.Any(e => e.Element == family)));
    }

string Percent(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
var lines = new List<string> { "Band,Tier,Proteins,Reachable,Family top-1,Family top-10" };
var speciesShares = new Dictionary<int, double>();
foreach (int band in Enumerable.Range(0, bandNames.Length).Append(-1))
    for (int t = 0; t < tierNames.Length; t++)
    {
        var cell = outcomes.Where(o => o.Tier == t && (band < 0 || o.Query.Band == band)).ToList();
        if (cell.Count == 0) continue;
        var reachable = cell.Where(o => o.Reachable).ToList();
        if (t == 2 && band >= 0 && reachable.Count > 0) speciesShares[band] = reachable.Count(o => o.Top1) / (double)reachable.Count;
        lines.Add(string.Join(",", band < 0 ? "all" : bandNames[band], tierNames[t], cell.Count, reachable.Count,
            Percent(reachable.Count(o => o.Top1), reachable.Count), Percent(reachable.Count(o => o.Top10), reachable.Count)));
    }
lines.Add("");
var failing = speciesShares.Where(s => s.Value < 0.5).Select(s => bandNames[s.Key]).ToList();
lines.Add("Rule: species-removed family top-1 under 50% in any band reframes the claim," +
    (failing.Count == 0 ? "no band under 50%: identity curve goes ahead" : "under 50% in " + string.Join(" ", failing) + ": reframe"));
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));

record Taxon(long Source, string? Species, string? Genus);
record Query(string Kind, int Band, string Protein, long SequenceId, Taxon Taxon, ZCurveGeneProfile.Gene Truth, double[] Axes, HashSet<int> Buckets);

sealed class Tally
{
    public int SameFamily;
    public readonly PriorityQueue<int, double> Top = new();
}
