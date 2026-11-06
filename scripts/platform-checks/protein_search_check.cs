#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Protein-search check (2026-11-05), before the protein-search API card: with every gene row stored (S7.R), how well does a protein alone
// find its gene among the whole catalogue, and how long does one query take when it streams every gene row from the database (your
// choice)? The 75.7% catalogue figure (results/CdsFollowUp_Search.csv) is extrapolated from 30 genomes and flagged as probably optimistic.
// Queries: proteins as BV-BRC serves them (genome_feature + feature_sequence, never DNA) from a seeded sample of stored genomes per kind,
// at least 500 aa, no internal stop; truth is the stored gene with the same record, coordinates and strand. Query axes from the protein
// alone (uniform codon weights, the record's genetic code). Two arms: the first 500 aa (3 windows, as the CDS sweep) and the whole protein.
// Candidates: every stored gene with at least 3 windows, scored by ZCurveGeneProfile.Score over the common windows.
// Measures: the truth's rank (1 + candidates scoring strictly higher), exact gene top-1 / top-10, and top-1 in the same PGfam (the truth's
// own gene counts; queries whose gene has no PGfam are left out of that measure). Identical or near-identical proteins in sister genomes
// are right answers for a user, so the family measure is the one the rule reads.
// Rule fixed before the run: the API step is worth building if the 500 aa arm's family top-1 is >= 70% over all queries and one streamed
// query takes <= 60 s on this machine.
// Short proteins (second run, 2026-11-05): with a window count w of 1 or 2, the queries are whole proteins of 167-332 aa (one window) or
// 333-499 aa (two), and the candidates every gene with at least w windows. Rule fixed before that run: a window count whose family top-1
// is under 50% puts the API's minimum protein length above it.
// Args: output path (default ProteinSearch_Catalogue.csv), then optionally w.
const int seed = 20261105;
const int proteinsPerGenome = 20;
int shortWindows = args.Length > 1 ? int.Parse(args[1], CultureInfo.InvariantCulture) : 0;
int minCandidateWindows = shortWindows > 0 ? shortWindows : 3;
// A protein of L aa has 2L positions-1-and-2 bases, so floor(2L / 333) windows.
int minProtein = shortWindows > 0 ? (shortWindows * ZCurveGeneProfile.WindowBases + 1) / 2 : 500;
int maxProtein = shortWindows > 0 ? ((shortWindows + 1) * ZCurveGeneProfile.WindowBases + 1) / 2 - 1 : int.MaxValue;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 6, ["archaeon"] = 2, ["plasmid"] = 2, ["virus"] = 3 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
var output = args.Length > 0 ? args[0] : "ProteinSearch_Catalogue.csv";
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
var profiles = new Dictionary<long, ZCurveGeneProfile>();
foreach (var (kind, count) in genomesPerKind)
{
    var samples = new List<(string Genome, long Id, string Record, string? Order)>();
    await using (var command = dataSource.CreateCommand(
        "select s.accession, r.id, r.record_id, r.tax_order from sequence_record r join sequence_source s on s.id = r.source_id " +
        "where r.kind = $1 and s.provider = 'BV-BRC' and s.owner_id is null and s.accession in (" +
        "  select s2.accession from sequence_source s2 join sequence_record r2 on r2.source_id = s2.id " +
        "  join profile_level p on p.sequence_id = r2.id and p.method = $2 and p.window_count >= $5 " +
        "  where r2.kind = $1 group by s2.accession order by md5(s2.accession || $3) limit $4) and exists (" +
        "  select 1 from profile_level p where p.sequence_id = r.id and p.method = $2)"))
    {
        command.Parameters.AddWithValue(kind);
        command.Parameters.AddWithValue(GeneProfiles.Method);
        command.Parameters.AddWithValue(seed.ToString());
        command.Parameters.AddWithValue(count);
        command.Parameters.AddWithValue(minCandidateWindows);
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync())
            samples.Add((reader.GetString(0), reader.GetInt64(1), reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetString(3)));
    }

    foreach (var genome in samples.GroupBy(s => s.Genome))
    {
        var records = genome.ToDictionary(s => s.Record);
        var features = await GetAsync($"https://www.bv-brc.org/api/genome_feature/?and(eq(genome_id,{genome.Key}),eq(feature_type,CDS)," +
            $"eq(annotation,PATRIC),ge(aa_length,{minProtein}){(maxProtein < int.MaxValue ? $",le(aa_length,{maxProtein})" : "")})&select(accession,start,end,strand,aa_sequence_md5)&sort(+feature_id)&limit(25000)");
        var picker = new Random(seed);
        var picked = features.EnumerateArray()
            .Where(f => f.TryGetProperty("aa_sequence_md5", out _) && records.ContainsKey(f.GetProperty("accession").GetString()!))
            .OrderBy(_ => picker.Next()).Take(3 * proteinsPerGenome).ToList();
        if (picked.Count == 0) continue;
        var md5s = string.Join(",", picked.Select(f => f.GetProperty("aa_sequence_md5").GetString()).Distinct());
        var proteins = (await GetAsync($"https://www.bv-brc.org/api/feature_sequence/?in(md5,({md5s}))&select(md5,sequence)&limit(25000)"))
            .EnumerateArray().ToDictionary(p => p.GetProperty("md5").GetString()!, p => p.GetProperty("sequence").GetString()!);
        int taken = 0;
        foreach (var feature in picked)
        {
            if (taken == proteinsPerGenome) break;
            var record = records[feature.GetProperty("accession").GetString()!];
            if (!proteins.TryGetValue(feature.GetProperty("aa_sequence_md5").GetString()!, out var protein) || protein.Contains('*')) continue;
            if (!profiles.TryGetValue(record.Id, out var profile))
                profiles[record.Id] = profile = (await GeneProfiles.ReadAsync(dataSource, record.Id))!;
            int start = feature.GetProperty("start").GetInt32() - 1, end = feature.GetProperty("end").GetInt32();
            bool reverse = feature.GetProperty("strand").GetString() == "-";
            int gene = -1;
            for (int g = 0; g < profile.Genes.Count && gene < 0; g++)
                if (profile.Genes[g].Start == start && profile.Genes[g].End == end && profile.Genes[g].Reverse == reverse) gene = g;
            if (gene < 0 || profile.Genes[gene].WindowCount < minCandidateWindows) continue;
            var code = record.Order is not null && mycoplasmaOrders.Contains(record.Order) ? ZCurveGeneProfile.mycoplasmaCode : ZCurveGeneProfile.standardCode;
            int family = profile.Genes[gene].Family;
            var arms = shortWindows > 0 ? new[] { ($"{shortWindows} window", protein) } : new[] { ("500 aa", protein[..minProtein]), ("whole", protein) };
            foreach (var (arm, text) in arms)
            {
                var axes = ZCurveGeneProfile.QueryFromProtein(text, null, code);
                queries.Add(new Query(kind, arm, record.Id, gene, family, axes, profile.Score(axes, gene)));
            }
            taken++;
        }
        Console.WriteLine($"{kind} {genome.Key}: {taken} proteins");
    }
}
Console.WriteLine($"{queries.Select(q => (q.SequenceId, q.Gene)).Distinct().Count()} proteins, {queries.Count} queries.");

// Streaming as the API would: one sequential reader, decoding and scoring spread over the cores. Pass 1 only decodes (the I/O and
// decode floor of any streamed query); pass 2 scores every query.
async Task<(double Seconds, long Genes, Tally[][] Tallies)> StreamAsync(IReadOnlyList<Query> scored)
{
    var watch = Stopwatch.StartNew();
    var rows = new BlockingCollection<(long Id, byte[] Payload)>(256);
    long genes = 0;
    var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
    {
        var tallies = scored.Select(_ => new Tally()).ToArray();
        long seen = 0;
        foreach (var (id, payload) in rows.GetConsumingEnumerable())
        {
            var profile = ZCurveGeneProfile.FromBinary(payload);
            for (int g = 0; g < profile.Genes.Count; g++)
            {
                if (profile.Genes[g].WindowCount < minCandidateWindows) continue;
                seen++;
                for (int q = 0; q < scored.Count; q++)
                {
                    var query = scored[q];
                    if (id == query.SequenceId && g == query.Gene) continue;
                    double score = profile.Score(query.Axes, g);
                    var tally = tallies[q];
                    if (score > query.Own) tally.Beaten++;
                    if (score > tally.Best) (tally.Best, tally.BestFamily) = (score, profile.Genes[g].Family);
                }
            }
        }
        Interlocked.Add(ref genes, seen);
        return tallies;
    })).ToArray();
    await using (var command = dataSource.CreateCommand("select sequence_id, payload from profile_level where method = $1"))
    {
        command.Parameters.AddWithValue(GeneProfiles.Method);
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync()) rows.Add((reader.GetInt64(0), reader.GetFieldValue<byte[]>(1)));
    }
    rows.CompleteAdding();
    var all = await Task.WhenAll(workers);
    return (watch.Elapsed.TotalSeconds, genes, all);
}

var (decodeSeconds, candidates, _) = await StreamAsync(Array.Empty<Query>());
Console.WriteLine($"Decode-only pass: {decodeSeconds:F1} s, {candidates:N0} candidate genes.");
var (scoreSeconds, _, perWorker) = await StreamAsync(queries);
Console.WriteLine($"Scoring pass: {scoreSeconds:F1} s for {queries.Count} queries.");
double perQuery = decodeSeconds + (scoreSeconds - decodeSeconds) / queries.Count;

var results = queries.Select((query, q) =>
{
    int rank = 1 + perWorker.Sum(w => w[q].Beaten);
    var best = perWorker.Select(w => w[q]).MaxBy(t => t.Best)!;
    bool ownWins = rank == 1 && query.Own >= best.Best;
    bool familyTop = query.Family != ZCurveGeneProfile.NoFamily && (ownWins || best.BestFamily == query.Family);
    return (query, rank, familyTop);
}).ToList();

var lines = new List<string> { "Arm,Kind,Proteins,Exact top-1,Exact top-10,With PGfam,Family top-1,Median rank" };
foreach (var arm in queries.Select(q => q.Arm).Distinct().ToList())
    foreach (var kind in genomesPerKind.Keys.Append("all"))
    {
        var cell = results.Where(r => r.query.Arm == arm && (kind == "all" || r.query.Kind == kind)).ToList();
        if (cell.Count == 0) continue;
        var withFamily = cell.Where(r => r.query.Family != ZCurveGeneProfile.NoFamily).ToList();
        var ranks = cell.Select(r => r.rank).Order().ToList();
        string Percent(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
        lines.Add(string.Join(",", arm, kind, cell.Count, Percent(cell.Count(r => r.rank == 1), cell.Count), Percent(cell.Count(r => r.rank <= 10), cell.Count),
            withFamily.Count, Percent(withFamily.Count(r => r.familyTop), withFamily.Count), ranks[ranks.Count / 2]));
    }
var ruleArm = queries[0].Arm;
var familyCell = results.Where(r => r.query.Arm == ruleArm && r.query.Family != ZCurveGeneProfile.NoFamily).ToList();
double familyShare = 100.0 * familyCell.Count(r => r.familyTop) / familyCell.Count;
lines.Add("");
lines.Add("Measure,Value");
lines.Add($"Candidate genes (>= {minCandidateWindows} windows),{candidates}");
lines.Add(string.Create(CultureInfo.InvariantCulture, $"Decode-only streamed pass (s),{decodeSeconds:F1}"));
lines.Add(string.Create(CultureInfo.InvariantCulture, $"Scoring pass for {queries.Count} queries (s),{scoreSeconds:F1}"));
lines.Add(string.Create(CultureInfo.InvariantCulture, $"One streamed query (s; decode pass + scoring share),{perQuery:F1}"));
lines.Add(string.Create(CultureInfo.InvariantCulture, $"Cores,{Environment.ProcessorCount}"));
string verdict = shortWindows > 0
    ? (familyShare >= 50 ? $"the API accepts proteins of {minProtein} aa" : $"the API's minimum is above {maxProtein} aa")
    : (familyShare >= 70 && perQuery <= 60 ? "pass: build the API step" : "stop");
lines.Add(shortWindows > 0
    ? string.Create(CultureInfo.InvariantCulture, $"Rule ({ruleArm} family top-1 >= 50%),{familyShare:F1}%: {verdict}")
    : string.Create(CultureInfo.InvariantCulture, $"Rule (500 aa family top-1 >= 70% and one query <= 60 s),{familyShare:F1}% and {perQuery:F1} s: {verdict}"));
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));

record Query(string Kind, string Arm, long SequenceId, int Gene, int Family, double[] Axes, double Own);

sealed class Tally
{
    public int Beaten;
    public double Best = double.NegativeInfinity;
    public int BestFamily = ZCurveGeneProfile.NoFamily;
}
