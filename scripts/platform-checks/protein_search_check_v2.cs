#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Protein-search check v2 (2026-11-05), after the S7.S pass: protein_search_check.cs measured the S7.R layer (one profile_level row per
// record, fixed windows); this measures the gene_profile layer at catalogue scale. Queries: whole proteins as BV-BRC serves them
// (genome_feature + feature_sequence, never DNA), up to 5 per length band from each of the same seeded genomes (6 bacteria, 2 archaea,
// 2 plasmids, 3 viruses; seed 20261105), no internal stop; truth is the stored gene with the same record, coordinates and strand. Each
// query is compared as ProteinSearch compares it: from 500 aa with the long genes (bucket 0), below that with the genes within 10% of
// its length (the buckets that band touches). Measured: exact gene top-1 / top-10 and top-1 in the right PGfam, per band. Timing: 3
// queries per band through ProteinSearch.SearchAsync itself, the API's path.
// Rule fixed before the run: if the right-family top-1 at 167-332 aa is under 70%, the API keeps refusing proteins under 333 aa; a band
// whose right-family top-1 is under 90% gets weakSignal.
// Args: output path (default ProteinSearch_CatalogueV2.csv).
const int seed = 20261105, perBand = 5;
string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 6, ["archaeon"] = 2, ["plasmid"] = 2, ["virus"] = 3 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
var output = args.Length > 0 ? args[0] : "ProteinSearch_CatalogueV2.csv";
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
            queries.Add(new Query(kind, band, protein, code == ZCurveGeneProfile.mycoplasmaCode ? 4 : 11, record.Id, profile.Genes[gene], axes,
                profile.Score(axes, gene), GeneProfiles.BucketsFor(protein.Length, ProteinSearch.LengthBand).ToHashSet()));
            taken[band]++;
        }
        Console.WriteLine($"{kind} {genome.Key}: {string.Join(" / ", taken)} proteins by band");
    }
}
Console.WriteLine($"{queries.Count} proteins.");

// One sequential read of every bucket row a query needs, decoding and scoring spread over the cores, as ProteinSearch does per query.
bool Comparable(Query query, ZCurveGeneProfile.Gene gene) => query.Protein.Length >= ZCurveGeneProfile.LongResidues
    ? gene.Residues >= ZCurveGeneProfile.LongResidues
    : gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - query.Protein.Length) <= ProteinSearch.LengthBand * query.Protein.Length;
var needed = queries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();
var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
var watch = Stopwatch.StartNew();
var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
{
    var tallies = queries.Select(_ => new Tally()).ToArray();
    foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable())
    {
        var profile = ZCurveGeneProfile.FromBinary(payload);
        for (int q = 0; q < queries.Count; q++)
        {
            var query = queries[q];
            if (!query.Buckets.Contains(bucket)) continue;
            var tally = tallies[q];
            for (int g = 0; g < profile.Genes.Count; g++)
            {
                var gene = profile.Genes[g];
                if (!Comparable(query, gene) || (id == query.SequenceId && gene.Start == query.Truth.Start && gene.End == query.Truth.End && gene.Reverse == query.Truth.Reverse)) continue;
                double score = profile.Score(query.Axes, g);
                if (score > query.Own) tally.Beaten++;
                if (score > tally.Best) (tally.Best, tally.BestFamily) = (score, gene.Family);
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

var results = queries.Select((query, q) =>
{
    int rank = 1 + perWorker.Sum(w => w[q].Beaten);
    var best = perWorker.Select(w => w[q]).MaxBy(t => t.Best)!;
    bool familyTop = query.Truth.Family != ZCurveGeneProfile.NoFamily && (rank == 1 && query.Own >= best.Best || best.BestFamily == query.Truth.Family);
    return (query, rank, familyTop);
}).ToList();

// The API's path: ProteinSearch.SearchAsync over the whole catalogue, 3 queries per band.
var timings = new Dictionary<int, List<double>>();
foreach (var group in queries.GroupBy(q => q.Band))
    foreach (var query in group.Take(3))
    {
        var timer = Stopwatch.StartNew();
        await ProteinSearch.SearchAsync(dataSource, ProteinQuery.FromProtein(query.Protein, query.GeneticCode), top: 50);
        (timings.TryGetValue(group.Key, out var list) ? list : timings[group.Key] = new List<double>()).Add(timer.Elapsed.TotalSeconds);
    }

string Percent(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
var lines = new List<string> { "Band,Kind,Proteins,Exact top-1,Exact top-10,With PGfam,Family top-1,Median rank,Search seconds (3 queries)" };
var familyShares = new Dictionary<int, double>();
for (int band = 0; band < bandNames.Length; band++)
    foreach (var kind in genomesPerKind.Keys.Append("all"))
    {
        var cell = results.Where(r => r.query.Band == band && (kind == "all" || r.query.Kind == kind)).ToList();
        if (cell.Count == 0) continue;
        var withFamily = cell.Where(r => r.query.Truth.Family != ZCurveGeneProfile.NoFamily).ToList();
        var ranks = cell.Select(r => r.rank).Order().ToList();
        if (kind == "all" && withFamily.Count > 0) familyShares[band] = withFamily.Count(r => r.familyTop) / (double)withFamily.Count;
        lines.Add(string.Join(",", bandNames[band], kind, cell.Count, Percent(cell.Count(r => r.rank == 1), cell.Count), Percent(cell.Count(r => r.rank <= 10), cell.Count),
            withFamily.Count, Percent(withFamily.Count(r => r.familyTop), withFamily.Count), ranks[ranks.Count / 2],
            kind == "all" && timings.TryGetValue(band, out var t) ? string.Join(" / ", t.Select(x => x.ToString("F1", CultureInfo.InvariantCulture))) : ""));
    }
lines.Add("");
double shortFamily = familyShares.GetValueOrDefault(2);
lines.Add(string.Create(CultureInfo.InvariantCulture, $"Rule: 167-332 aa right-family top-1 {100 * shortFamily:F1}% (>= 70%: accept from 50 aa),") +
    (shortFamily >= 0.7 ? "accept from 50 aa" : "keep refusing under 333 aa"));
var weakBands = string.Join(" ", familyShares.Where(f => f.Value < 0.9).Select(f => bandNames[f.Key]));
lines.Add("Rule: bands under 90% right-family top-1 get weakSignal," + (weakBands.Length > 0 ? weakBands : "none"));
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));

record Query(string Kind, int Band, string Protein, int GeneticCode, long SequenceId, ZCurveGeneProfile.Gene Truth, double[] Axes, double Own, HashSet<int> Buckets);

sealed class Tally
{
    public int Beaten;
    public double Best = double.NegativeInfinity;
    public int BestFamily = ZCurveGeneProfile.NoFamily;
}
