#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Protein E-value check (2026-11-06), before building p and E for protein hits. Queries: whole BV-BRC proteins with a PGfam (never
// stored), up to 8 per length band from each of the seeded genomes of protein_search_check_v2.cs. Each is compared as ProteinSearch
// compares it (from 500 aa with the long genes, below with genes within 10% of its length), its own gene left out. Per query: its best
// score against a gene of its own family, and how many genes of other families reach each score (genes with no family are left out of
// both). The proposed E(s) for a band is the mean of those counts over the band's queries.
// Per band: each query's E at its best same-family score from the curve of the other queries (leave one out), and the spread of the
// counts between queries at the score where the band's mean count is 10.
// Rule fixed before the run: no E if, at 167-332 aa, the median query's E at its best same-family hit is above 0.01, or the 90th
// percentile count at that fixed score is more than 100x the 10th (floored at 1).
// Args: output path (default ProteinNull_Check.csv).
const int seed = 20261105, perBand = 8, bins = 2000;
string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 6, ["archaeon"] = 2, ["plasmid"] = 2, ["virus"] = 3 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
var output = args.Length > 0 ? args[0] : "ProteinNull_Check.csv";
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
            $"eq(annotation,PATRIC),ge(aa_length,{ZCurveGeneProfile.MinResidues}))&select(accession,start,end,strand,aa_length,aa_sequence_md5,pgfam_id)&sort(+feature_id)&limit(25000)");
        var picker = new Random(seed);
        var picked = features.EnumerateArray()
            .Where(f => f.TryGetProperty("aa_sequence_md5", out _) && f.TryGetProperty("pgfam_id", out _) && records.ContainsKey(f.GetProperty("accession").GetString()!))
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
            int band = BandOf(protein.Length), family = ProteinFamilies.Parse(feature.GetProperty("pgfam_id").GetString());
            if (taken[band] == perBand || family == ZCurveGeneProfile.NoFamily) continue;
            if (!profiles.TryGetValue(record.Id, out var profile)) profiles[record.Id] = profile = (await GeneProfiles.ReadAsync(dataSource, record.Id))!;
            int start = feature.GetProperty("start").GetInt32() - 1, end = feature.GetProperty("end").GetInt32();
            bool reverse = feature.GetProperty("strand").GetString() == "-";
            int gene = Enumerable.Range(0, profile.Genes.Count).FirstOrDefault(g => profile.Genes[g].Start == start && profile.Genes[g].End == end
                && profile.Genes[g].Reverse == reverse, -1);
            if (gene < 0) continue;
            var code = record.Order is not null && mycoplasmaOrders.Contains(record.Order) ? ZCurveGeneProfile.mycoplasmaCode : ZCurveGeneProfile.standardCode;
            queries.Add(new Query(band, protein.Length, family, record.Id, profile.Genes[gene], ZCurveGeneProfile.QueryFromProtein(protein, null, code),
                GeneProfiles.BucketsFor(protein.Length, ProteinSearch.LengthBand).ToHashSet()));
            taken[band]++;
        }
        Console.WriteLine($"{kind} {genome.Key}: {string.Join(" / ", taken)} proteins by band");
    }
}
Console.WriteLine($"{queries.Count} proteins with a family.");

bool Comparable(Query q, ZCurveGeneProfile.Gene g) => q.LengthAa >= ZCurveGeneProfile.LongResidues
    ? g.Residues >= ZCurveGeneProfile.LongResidues
    : g.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(g.Residues - q.LengthAa) <= ProteinSearch.LengthBand * q.LengthAa;
static int Bin(double score) => Math.Clamp((int)Math.Floor(score * bins), 0, bins - 1);   // scores of 0 to 1; below 0 counts in bin 0

var needed = queries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();
var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
var watch = Stopwatch.StartNew();
var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
{
    var other = queries.Select(_ => new long[bins]).ToArray();
    var bestSame = queries.Select(_ => double.NegativeInfinity).ToArray();
    foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable())
    {
        var profile = ZCurveGeneProfile.FromBinary(payload);
        for (int q = 0; q < queries.Count; q++)
        {
            var query = queries[q];
            if (!query.Buckets.Contains(bucket)) continue;
            for (int g = 0; g < profile.Genes.Count; g++)
            {
                var gene = profile.Genes[g];
                if (gene.Family == ZCurveGeneProfile.NoFamily || !Comparable(query, gene)) continue;
                if (id == query.SequenceId && gene.Start == query.Own.Start && gene.End == query.Own.End && gene.Reverse == query.Own.Reverse) continue;
                double score = profile.Score(query.Axes, g);
                if (gene.Family == query.Family) bestSame[q] = Math.Max(bestSame[q], score);
                else other[q][Bin(score)]++;
            }
        }
    }
    return (other, bestSame);
})).ToArray();
await using (var command = dataSource.CreateCommand("select sequence_id, bucket, payload from gene_profile where bucket = any($1) and window_count > 0"))
{
    command.Parameters.AddWithValue(needed);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) rows.Add((reader.GetInt64(0), reader.GetInt16(1), reader.GetFieldValue<byte[]>(2)));
}
rows.CompleteAdding();
var results = await Task.WhenAll(workers);
Console.WriteLine($"Scoring pass: {watch.Elapsed.TotalSeconds:F0} s.");

// Per query: counts of other-family genes at or above each bin (survival), and its best same-family score.
var survival = queries.Select((_, q) =>
{
    var s = new double[bins + 1];
    for (int b = bins - 1; b >= 0; b--) s[b] = s[b + 1] + results.Sum(r => r.other[q][b]);
    return s;
}).ToArray();
var same = queries.Select((_, q) => results.Max(r => r.bestSame[q])).ToArray();

string F(double v) => v.ToString("G3", CultureInfo.InvariantCulture);
double Quantile(IEnumerable<double> values, double p) { var s = values.Order().ToList(); return s.Count == 0 ? double.NaN : s[(int)Math.Round(p * (s.Count - 1))]; }
var lines = new List<string> { "Band,Proteins,With a same-family gene,Other-family genes per search (mean),Median E at best same-family hit,90% E,Share E < 0.01,Share E < 1," +
    "Score where mean count = 10,Count p10 there,Count p90 there,Spread p90/max(p10 1),Same-family best score median,Other-family best score median" };
var verdict = new Dictionary<int, (double MedianE, double Spread)>();
for (int band = 0; band < bandNames.Length; band++)
{
    var members = Enumerable.Range(0, queries.Count).Where(q => queries[q].Band == band).ToList();
    if (members.Count < 3) continue;
    var withSame = members.Where(q => !double.IsNegativeInfinity(same[q])).ToList();
    double Mean(IEnumerable<int> qs, int b) => qs.Average(q => survival[q][b]);
    var e = withSame.Select(q => Mean(members.Where(o => o != q), Bin(same[q]))).ToList();
    int fixedBin = Enumerable.Range(0, bins).Reverse().FirstOrDefault(b => Mean(members, b) >= 10, 0);
    var counts = members.Select(q => survival[q][fixedBin]).ToList();
    double p10 = Quantile(counts, 0.1), p90 = Quantile(counts, 0.9), spread = p90 / Math.Max(1, p10);
    double bestOther(int q) => Enumerable.Range(0, bins).Reverse().FirstOrDefault(b => survival[q][b] >= 1, 0) / (double)bins;
    verdict[band] = (Quantile(e, 0.5), spread);
    lines.Add(string.Join(",", bandNames[band], members.Count, withSame.Count, F(survival[members[0]][0] > 0 ? members.Average(q => survival[q][0]) : 0),
        F(Quantile(e, 0.5)), F(Quantile(e, 0.9)), (e.Count(x => x < 0.01) / (double)e.Count).ToString("P1", CultureInfo.InvariantCulture),
        (e.Count(x => x < 1) / (double)e.Count).ToString("P1", CultureInfo.InvariantCulture), F(fixedBin / (double)bins), F(p10), F(p90), F(spread),
        F(Quantile(withSame.Select(q => same[q]), 0.5)), F(Quantile(members.Select(bestOther), 0.5))));
}
lines.Add("");
var (medianE, spreadAt) = verdict.GetValueOrDefault(2, (double.NaN, double.NaN));
bool pass = medianE <= 0.01 && spreadAt <= 100;
lines.Add(string.Create(CultureInfo.InvariantCulture, $"Rule (167-332 aa: median E at the best same-family hit <= 0.01 and spread <= 100x),{F(medianE)} and {F(spreadAt)}x: ") + (pass ? "build" : "stop"));
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));

record Query(int Band, int LengthAa, int Family, long SequenceId, ZCurveGeneProfile.Gene Own, double[] Axes, HashSet<int> Buckets);
