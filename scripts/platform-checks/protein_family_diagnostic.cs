#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Protein family diagnostic (2026-11-12), after protein_homolog_check.cs: with the query's own genus removed, the right PGfam is first for
// 65.9% of proteins and in the top 10 for 79.3% (results/ProteinHomolog_Check.csv). Why do the rest miss, and do two cheap fixes help?
// (1) Weighting the 9 features (P0-P2 x x, y, z) of ZCurveGeneProfile.Score: unweighted (today), P0 only, P0+P1, 1/global variance,
// 1/within-family variance, between/within ratio; the last three learned from a 2% row sample on the tuning half only. (2) Choosing the
// family from the top hits: best hit (today), best hit with a family, majority of the first 3 / 5 / 10 hits with a family.
// Queries: protein_search_check_v2.cs's rule and seed widened to 48 bacteria, 8 archaea, 8 plasmids (viruses have no PGfam). Genomes are
// split into a tuning and a test half by the MD5 of the accession. Measure: family top-1 among genes outside the query's source and NCBI
// genus, over queries with a same-family gene left there. Exact top-1 (own gene first among all genes) is kept per weighting as a guard.
// Diagnostic for today's score: rank of the first same-family gene, score gap to the top hit, and that gene's length against the query.
// Rule fixed before the run: the combination best on the tuning half goes ahead if, on the test half, it gains >= 5 points of family
// top-1 over today's (unweighted, best hit) with a 95% genome-bootstrap CI excluding 0, and loses <= 1 point of exact top-1. Otherwise
// tuning stops and the paper reports today's number.
// --confirm (2026-11-20, card G7): confirmation of the 1/variance weights G4 printed, on new genomes (seed 20261113, G4's 64 skipped),
// with the weights fixed and nothing tuned: two arms (unweighted, 1/variance), every query is a test query. Rule fixed before the run: go
// if exact-gene top-1 gains >= 2 points with the 95% genome-bootstrap CI above 0, and family top-1 (own genus removed, best hit) loses
// <= 1 point.
// --moments (2026-11-22): each moment compared on its own, G4's tuning half only (the test half and the confirmation set stay clean
// for a combined score). Per pair (query against every comparable gene outside its genus), from the 9 stored numbers: per axis and
// pooled over axes, the P0 offset (|mean difference| over the common windows), P0 raw (RMS difference), P0 pattern (RMS difference
// with the offset removed), and for P1 and P2 the RMS difference, the share of windows with the same sign, and Pearson r across
// windows; plus today's weighted score. Per sub-score alone: family top-1 (best hit; genes tied at the best share the credit; genes
// without a family count as in G4) and the AUC of same-family against other-family pairs. Rule fixed before the run: a combined per-moment
// score is checked only if, pooled over axes, P0 pattern beats P0 raw, or P1 / P2 sign or corr beats its diff, by >= 3 points of
// family top-1; otherwise per-moment comparison stops.
// Args: summary path (default ProteinFamily_Diagnostic.csv), per-query path (default ProteinFamily_Diagnostic_Queries.csv),
// [--confirm | --moments].
const int seed = 20261105, confirmSeed = 20261113, perBand = 5, topK = 20, diagnosticK = 200, bootstrapRounds = 2000;
bool confirm = args.Contains("--confirm");
bool perMoment = args.Contains("--moments") && !confirm;
args = args.Where(a => !a.StartsWith("--", StringComparison.Ordinal)).ToArray();
int selectionSeed = confirm ? confirmSeed : seed;
string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 48, ["archaeon"] = 8, ["plasmid"] = 8 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
const int F = ZCurveGeneProfile.FeatureCount;
var summaryPath = args.Length > 0 ? args[0] : "ProteinFamily_Diagnostic.csv";
var queriesPath = args.Length > 1 ? args[1] : "ProteinFamily_Diagnostic_Queries.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
using var http = new HttpClient { Timeout = TimeSpan.FromMinutes(5) };
http.DefaultRequestHeaders.UserAgent.ParseAdd("CATSequenceAnalysis-platform-check/1.0");
http.DefaultRequestHeaders.Accept.ParseAdd("application/json");
var total = Stopwatch.StartNew();

async Task<JsonElement> GetAsync(string url)
{
    await Task.Delay(1000);
    using var response = await http.GetAsync(url);
    response.EnsureSuccessStatusCode();
    using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
    return document.RootElement.Clone();
}

var taxa = new Dictionary<long, Taxon>();
await using (var command = dataSource.CreateCommand("select id, source_id, tax_genus from sequence_record"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
        taxa[reader.GetInt64(0)] = new Taxon(reader.GetInt64(1), reader.IsDBNull(2) ? null : reader.GetString(2));

// Query selection: protein_search_check_v2.cs's rule, more genomes, only proteins with a PGfam and an NCBI genus.
var queries = new List<Query>();
// The genomes of a seed, in md5 order: the same accessions G4's single query picked.
async Task<List<string>> AccessionsAsync(string kind, int forSeed, int limit)
{
    var list = new List<string>();
    await using var command = dataSource.CreateCommand(
        "select s2.accession from sequence_source s2 join sequence_record r2 on r2.source_id = s2.id " +
        "join gene_profile g on g.sequence_id = r2.id and g.bucket = 0 " +
        "where r2.kind = $1 group by s2.accession order by md5(s2.accession || $2) limit $3");
    command.Parameters.AddWithValue(kind);
    command.Parameters.AddWithValue(forSeed.ToString(CultureInfo.InvariantCulture));
    command.Parameters.AddWithValue(limit);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) list.Add(reader.GetString(0));
    return list;
}
foreach (var (kind, count) in genomesPerKind)
{
    var excluded = confirm ? (await AccessionsAsync(kind, seed, count)).ToHashSet() : [];
    var accessions = (await AccessionsAsync(kind, selectionSeed, count + excluded.Count)).Where(a => !excluded.Contains(a)).Take(count).ToArray();
    var samples = new List<(string Genome, long Id, string Record, string? Order)>();
    await using (var command = dataSource.CreateCommand(
        "select s.accession, r.id, r.record_id, r.tax_order from sequence_record r join sequence_source s on s.id = r.source_id " +
        "where r.kind = $1 and s.provider = 'BV-BRC' and s.owner_id is null and s.accession = any($2) and exists (" +
        "  select 1 from gene_profile g where g.sequence_id = r.id)"))
    {
        command.Parameters.AddWithValue(kind);
        command.Parameters.AddWithValue(accessions);
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync())
            samples.Add((reader.GetString(0), reader.GetInt64(1), reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetString(3)));
    }

    foreach (var genome in samples.GroupBy(s => s.Genome).OrderBy(g => Array.IndexOf(accessions, g.Key)))
    {
        bool tuning = !confirm && (MD5.HashData(Encoding.ASCII.GetBytes(genome.Key))[0] & 1) == 0;
        if (perMoment && !tuning) continue;
        var records = genome.ToDictionary(s => s.Record);
        JsonElement features;
        try
        {
            features = await GetAsync($"https://www.bv-brc.org/api/genome_feature/?and(eq(genome_id,{genome.Key}),eq(feature_type,CDS)," +
                $"eq(annotation,PATRIC),ge(aa_length,{ZCurveGeneProfile.MinResidues}))&select(accession,start,end,strand,aa_length,aa_sequence_md5)&sort(+feature_id)&limit(25000)");
        }
        catch (HttpRequestException e) { Console.WriteLine($"{kind} {genome.Key}: skipped ({e.Message})"); continue; }
        var picker = new Random(selectionSeed);
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
            taken[band]++;
            if (profile.Genes[gene].Family == ZCurveGeneProfile.NoFamily || taxa[record.Id].Genus is null) continue;
            var code = record.Order is not null && mycoplasmaOrders.Contains(record.Order) ? ZCurveGeneProfile.mycoplasmaCode : ZCurveGeneProfile.standardCode;
            queries.Add(new Query(genome.Key, tuning, band, protein.Length, record.Id, taxa[record.Id], profile.Genes[gene],
                ZCurveGeneProfile.QueryFromProtein(protein, null, code), profile.WindowsOf(gene),
                GeneProfiles.BucketsFor(protein.Length, ProteinSearch.LengthBand).ToHashSet()));
        }
        Console.WriteLine($"{kind} {genome.Key} ({(tuning ? "tuning" : "test")}): {string.Join(" / ", taken)} proteins by band");
    }
}
Console.WriteLine($"{queries.Count} proteins with a PGfam and a genus: {queries.Count(q => q.Tuning)} tuning, {queries.Count(q => !q.Tuning)} test. " +
    $"{total.Elapsed.TotalMinutes:F1} min.");

bool Comparable(Query query, ZCurveGeneProfile.Gene gene) => query.Residues >= ZCurveGeneProfile.LongResidues
    ? gene.Residues >= ZCurveGeneProfile.LongResidues
    : gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - query.Residues) <= ProteinSearch.LengthBand * query.Residues;
bool OutsideGenus(Query query, Taxon candidate) => candidate.Source != query.Taxon.Source && (candidate.Genus is null || candidate.Genus != query.Taxon.Genus);
bool IsSelf(Query query, long id, ZCurveGeneProfile.Gene gene) =>
    id == query.SequenceId && gene.Start == query.Truth.Start && gene.End == query.Truth.End && gene.Reverse == query.Truth.Reverse;
var byBucket = Enumerable.Range(0, 64).Select(b => queries.Where(q => q.Buckets.Contains(b)).ToArray()).ToArray();
var needed = queries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();

// Per-window squared feature differences over the common windows; returns the window count.
static int Differences(double[] query, double[] gene, double[] squared)
{
    int count = Math.Min(query.Length, gene.Length) / F;
    for (int i = 0; i < count * F; i++) { double d = query[i] - gene[i]; squared[i] = d * d; }
    return count;
}
static double Score(double[] squared, int count, double[] weights)
{
    double sum = 0;
    for (int w = 0; w < count; w++)
    {
        double s = 0;
        for (int f = 0; f < F; f++) s += weights[f] * squared[w * F + f];
        sum += Math.Sqrt(s);
    }
    return 1 - sum / count;
}

async Task<T[]> PassAsync<T>(string sql, Func<T> start, Action<T, long, int, ZCurveGeneProfile> handle)
{
    var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
    var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
    {
        var state = start();
        foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable()) handle(state, id, bucket, ZCurveGeneProfile.FromBinary(payload));
        return state;
    })).ToArray();
    long read = 0;
    await using (var command = dataSource.CreateCommand(sql))
    {
        command.Parameters.AddWithValue(needed);
        command.CommandTimeout = 0;
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync())
        {
            rows.Add((reader.GetInt64(0), reader.GetInt16(1), reader.GetFieldValue<byte[]>(2)));
            if (++read % 20_000 == 0) Console.WriteLine($"  {read:N0} rows, {total.Elapsed.TotalMinutes:F1} min");
        }
    }
    rows.CompleteAdding();
    return await Task.WhenAll(workers);
}

string Vector(double[] v) => string.Join(" ", v.Select(x => x.ToString("G3", CultureInfo.InvariantCulture)));
static double[] Normalized(double[] w) { double s = w.Sum(); return w.Select(x => F * x / s).ToArray(); }

if (perMoment)
{
    const int bins = 4096;
    string[] kindNames = ["offset", "P0 raw", "P0 pattern", "P1 diff", "P1 sign", "P1 corr", "P2 diff", "P2 sign", "P2 corr"];
    string[] slotNames = ["x", "y", "z", "pooled"];
    int S = 4 * kindNames.Length + 1, today = S - 1;
    string Name(int s) => s == today ? "today (weighted 9)" : $"{kindNames[s / 4]} {slotNames[s % 4]}";
    var featureWeights = ZCurveGeneProfile.FeatureWeights.ToArray();

    // Every sub-score is a similarity (distances negated), so one bin order fits all: higher bin = more alike. Distances are binned
    // on sqrt(d / 4), which keeps resolution near 0 where the pairs are.
    int BinOf(int s, double v)
    {
        int kind = s == today ? -1 : s / 4;
        double u = kind is 4 or 7 ? v : kind is 5 or 8 ? (v + 1) / 2 : 1 - Math.Sqrt(Math.Clamp((s == today ? 1 - v : -v) / 4, 0, 1));
        return (int)(Math.Clamp(u, 0, 1) * (bins - 1));
    }

    // Identity check (query against itself) and today's score against the G4 scorer on the own gene.
    var self = new double[S];
    foreach (var query in queries)
    {
        int count = query.Axes.Length / F;
        MomentScores(query.Axes, query.Axes, count, featureWeights, self);
        for (int s = 0; s < S; s++)
        {
            int kind = s == today ? -1 : s / 4;
            double expected = s == today || kind is 4 or 5 or 7 or 8 ? 1 : 0;
            if (kind is 5 or 8 && double.IsNaN(self[s])) continue;
            if (Math.Abs(self[s] - expected) > 1e-9) throw new InvalidOperationException($"Identity check: {Name(s)} = {self[s]}, expected {expected}.");
        }
        var squared = new double[query.TruthWindows.Length];
        int own = Differences(query.Axes, query.TruthWindows, squared);
        MomentScores(query.Axes, query.TruthWindows, own, featureWeights, self);
        if (own > 0 && Math.Abs(self[today] - Score(squared, own, featureWeights)) > 1e-12) throw new InvalidOperationException("Today's score differs from the G4 scorer.");
    }
    Console.WriteLine($"Identity check passed on {queries.Count} queries.");

    var parts = await PassAsync("select sequence_id, bucket, payload from gene_profile where bucket = any($1) and window_count > 0",
        () => new MomentTallies(queries.Count, S, bandNames.Length, bins), (t, id, bucket, profile) =>
        {
            var scores = new double[S];
            var candidate = taxa[id];
            var relevant = byBucket[bucket];
            for (int g = 0; g < profile.Genes.Count; g++)
            {
                var gene = profile.Genes[g];
                double[]? windows = null;
                foreach (var query in relevant)
                {
                    if (!Comparable(query, gene) || IsSelf(query, id, gene) || !OutsideGenus(query, candidate)) continue;
                    windows ??= profile.WindowsOf(g);
                    int count = Math.Min(query.Axes.Length, windows.Length) / F;
                    if (count == 0) continue;
                    MomentScores(query.Axes, windows, count, featureWeights, scores);
                    bool same = gene.Family == query.Truth.Family;
                    if (same) t.SameFamily[query.Index]++;
                    int pairClass = gene.Family == ZCurveGeneProfile.NoFamily ? -1 : same ? 0 : 1;
                    for (int s = 0; s < S; s++)
                    {
                        double v = scores[s];
                        if (double.IsNaN(v)) { t.Undefined[s]++; continue; }
                        int slot = query.Index * S + s;
                        if (v > t.Best[slot]) (t.Best[slot], t.Tied[slot], t.TiedRight[slot]) = (v, 0, 0);
                        if (v == t.Best[slot]) { t.Tied[slot]++; if (same) t.TiedRight[slot]++; }
                        if (pairClass >= 0) t.Histogram[((query.Band * S + s) * 2 + pairClass) * bins + BinOf(s, v)]++;
                    }
                }
            }
        });
    Console.WriteLine($"Moments pass done, {total.Elapsed.TotalMinutes:F1} min.");

    var right = new double[queries.Count, S];
    var ties = new long[queries.Count, S];
    var reach = new bool[queries.Count];
    foreach (var query in queries)
    {
        int q = query.Index;
        reach[q] = parts.Sum(t => t.SameFamily[q]) > 0;
        for (int s = 0; s < S; s++)
        {
            // Genes tied at the best value share the credit: discrete sub-scores (sign shares) tie often, and a tie-break by today's
            // score made the first run's sign rows mostly today's score (2026-11-22).
            double top = parts.Max(t => t.Best[q * S + s]);
            long tied = parts.Where(t => t.Best[q * S + s] == top).Sum(t => t.Tied[q * S + s]);
            long tiedRight = parts.Where(t => t.Best[q * S + s] == top).Sum(t => t.TiedRight[q * S + s]);
            right[q, s] = tied == 0 ? 0 : tiedRight / (double)tied;
            ties[q, s] = tied;
        }
    }
    var histogram = new long[bandNames.Length * S * 2 * bins];
    foreach (var t in parts) for (int i = 0; i < histogram.Length; i++) histogram[i] += t.Histogram[i];
    var undefined = Enumerable.Range(0, S).Select(s => parts.Sum(t => t.Undefined[s])).ToArray();

    (long Positives, long Negatives, double Auc) AucOf(int s, int band)
    {
        var positive = new double[bins];
        var negative = new double[bins];
        foreach (int b in band < 0 ? Enumerable.Range(0, bandNames.Length) : Enumerable.Repeat(band, 1))
            for (int i = 0; i < bins; i++)
            {
                positive[i] += histogram[((b * S + s) * 2) * bins + i];
                negative[i] += histogram[((b * S + s) * 2 + 1) * bins + i];
            }
        double p = positive.Sum(), n = negative.Sum(), below = 0, wins = 0;
        for (int i = 0; i < bins; i++) { wins += positive[i] * (below + 0.5 * negative[i]); below += negative[i]; }
        return ((long)p, (long)n, p == 0 || n == 0 ? double.NaN : wins / (p * n));
    }

    var genomes = queries.GroupBy(q => q.Genome).Select(g => g.ToList()).ToList();
    (double Low, double High) Interval(Func<Query, double> value, int roundsSeed)
    {
        var random = new Random(roundsSeed);
        var draws = new List<double>();
        for (int r = 0; r < bootstrapRounds; r++)
        {
            var resampled = Enumerable.Range(0, genomes.Count).SelectMany(_ => genomes[random.Next(genomes.Count)]).Where(q => reach[q.Index]).ToList();
            if (resampled.Count > 0) draws.Add(resampled.Average(value));
        }
        draws.Sort();
        return (draws[(int)(0.025 * draws.Count)], draws[(int)(0.975 * draws.Count) - 1]);
    }
    string Fixed(double x, int digits) => double.IsNaN(x) ? "" : x.ToString("F" + digits, CultureInfo.InvariantCulture);

    var summary = new List<string> { "Band,Sub-score,Reachable,Family top-1,95% CI low,95% CI high,Same-family pairs,Other-family pairs,AUC,Mean genes tied at best" };
    foreach (int band in Enumerable.Range(0, bandNames.Length).Prepend(-1))
    {
        var cell = queries.Where(q => reach[q.Index] && (band < 0 || q.Band == band)).ToList();
        for (int s = 0; s < S; s++)
        {
            double share = cell.Count == 0 ? double.NaN : cell.Average(q => right[q.Index, s]);
            double meanTied = cell.Count == 0 ? double.NaN : cell.Average(q => (double)ties[q.Index, s]);
            var (lower, upper) = band < 0 ? Interval(q => right[q.Index, s], seed + s) : (double.NaN, double.NaN);
            var (positives, negatives, auc) = AucOf(s, band);
            summary.Add(string.Join(",", band < 0 ? "all" : bandNames[band], Name(s), cell.Count, Fixed(100 * share, 1), Fixed(100 * lower, 1),
                Fixed(100 * upper, 1), positives, negatives, Fixed(auc, 4), Fixed(meanTied, 1)));
        }
    }

    // Rule: pooled over axes, the alternative against its difference-based counterpart, paired genome bootstrap.
    var allReach = queries.Where(q => reach[q.Index]).ToList();
    int Pooled(int kind) => 4 * kind + 3;
    var comparisons = new (int After, int Before)[] { (2, 1), (4, 3), (5, 3), (7, 6), (8, 6) };
    bool go = false;
    summary.Add("");
    foreach (var (after, before) in comparisons)
    {
        double delta = allReach.Average(q => right[q.Index, Pooled(after)] - right[q.Index, Pooled(before)]);
        var (lower, upper) = Interval(q => right[q.Index, Pooled(after)] - right[q.Index, Pooled(before)], seed + 100 + after);
        go |= delta >= 0.03;
        summary.Add($"{Name(Pooled(after))} vs {Name(Pooled(before))},{Points(delta)} points (95% CI {Points(lower)} to {Points(upper)}; {allReach.Count} reachable)");
    }
    summary.Add($"Undefined pairs (corr with < 3 windows or no variance),{string.Join(" ", Enumerable.Range(0, S).Where(s => undefined[s] > 0).Select(s => $"{Name(s)}: {undefined[s]:N0}"))}");
    summary.Add($"Reference,today's weighted score family top-1 {Fixed(100 * allReach.Average(q => right[q.Index, today]), 1)}% (G4 tuning half, 1/variance, best hit: 62.3%)");
    summary.Add("Rule (pooled: pattern vs raw, or sign / corr vs diff, >= +3 points)," + (go ? "go: check a combined per-moment score" : "stop: per-moment comparison adds nothing"));
    summary.Add($"Total,{total.Elapsed.TotalMinutes:F1} min");
    File.WriteAllLines(summaryPath, summary);
    Console.WriteLine(string.Join(Environment.NewLine, summary.Where(l => !bandNames.Any(b => l.StartsWith(b + ",")))));

    var momentQueries = new List<string> { "Genome,Band,Residues,Family,Reachable," + string.Join(",", Enumerable.Range(0, S).Select(Name)) };
    foreach (var q in queries)
        momentQueries.Add(string.Join(",", q.Genome, bandNames[q.Band], q.Residues, q.Truth.Family, reach[q.Index] ? 1 : 0,
            string.Join(",", Enumerable.Range(0, S).Select(s => Fixed(right[q.Index, s], 3)))));
    File.WriteAllLines(queriesPath, momentQueries);
    return;
}

// Sub-scores of one pair over its common windows, all as similarities (distances negated; NaN where undefined), in the order
// 4 * kind + slot (slot x, y, z, pooled over axes), then today's weighted score. Pooled distances weight each axis's squared
// distance by FeatureWeights; pooled sign and corr take the 3 axes' values together.
static void MomentScores(double[] query, double[] gene, int count, double[] weights, double[] scores)
{
    double pooledOffset = 0, pooledRaw = 0, pooledPattern = 0;
    Span<double> pooledDiff = stackalloc double[2];
    Span<double> pooledAgree = stackalloc double[2];
    Span<double> pooledSums = stackalloc double[10];
    for (int axis = 0; axis < 3; axis++)
    {
        double sum = 0, sumSquares = 0;
        for (int w = 0; w < count; w++) { double d = query[w * F + axis] - gene[w * F + axis]; sum += d; sumSquares += d * d; }
        double mean = sum / count, raw = sumSquares / count, pattern = Math.Max(0, raw - mean * mean);
        scores[axis] = -Math.Abs(mean);
        scores[4 + axis] = -Math.Sqrt(raw);
        scores[8 + axis] = -Math.Sqrt(pattern);
        pooledOffset += weights[axis] * mean * mean;
        pooledRaw += weights[axis] * raw;
        pooledPattern += weights[axis] * pattern;
        for (int order = 1; order <= 2; order++)
        {
            int f = 3 * order + axis, kind = order == 1 ? 3 : 6;
            double diff = 0, agree = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
            for (int w = 0; w < count; w++)
            {
                double x = query[w * F + f], y = gene[w * F + f], d = x - y;
                diff += d * d;
                if (Math.Sign(x) == Math.Sign(y)) agree++;
                sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
            }
            scores[4 * kind + axis] = -Math.Sqrt(diff / count);
            scores[4 * (kind + 1) + axis] = agree / count;
            scores[4 * (kind + 2) + axis] = Pearson(count, sx, sy, sxx, syy, sxy);
            pooledDiff[order - 1] += weights[f] * diff / count;
            pooledAgree[order - 1] += agree;
            var sums = pooledSums.Slice(5 * (order - 1), 5);
            sums[0] += sx; sums[1] += sy; sums[2] += sxx; sums[3] += syy; sums[4] += sxy;
        }
    }
    scores[3] = -Math.Sqrt(pooledOffset);
    scores[7] = -Math.Sqrt(pooledRaw);
    scores[11] = -Math.Sqrt(pooledPattern);
    for (int order = 1; order <= 2; order++)
    {
        int kind = order == 1 ? 3 : 6;
        var sums = pooledSums.Slice(5 * (order - 1), 5);
        scores[4 * kind + 3] = -Math.Sqrt(pooledDiff[order - 1]);
        scores[4 * (kind + 1) + 3] = pooledAgree[order - 1] / (3 * count);
        scores[4 * (kind + 2) + 3] = Pearson(3 * count, sums[0], sums[1], sums[2], sums[3], sums[4]);
    }
    double distance = 0;
    for (int w = 0; w < count; w++)
    {
        double s = 0;
        for (int f = 0; f < F; f++) { double d = query[w * F + f] - gene[w * F + f]; s += weights[f] * (d * d); }
        distance += Math.Sqrt(s);
    }
    scores[^1] = 1 - distance / count;
}

static double Pearson(int n, double sx, double sy, double sxx, double syy, double sxy)
{
    if (n < 3) return double.NaN;
    double vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
    return vx <= 1e-15 || vy <= 1e-15 ? double.NaN : (sxy - sx * sy / n) / Math.Sqrt(vx * vy);
}

(string Name, double[] Weights)[] weightings;
if (confirm)
    weightings = [("unweighted", Enumerable.Repeat(1.0, F).ToArray()), ("1/variance (G4 fixed)", [0.741, 0.838, 0.333, 1.17, 1.16, 1.15, 1.19, 1.21, 1.2])];
else
{
// Pre-pass on a 2% row sample: global feature variance over all windows, and within / between-family mean squared differences for
// tuning queries against genes outside their genus.
var sample = await PassAsync(
    "select sequence_id, bucket, payload from gene_profile tablesample system (2) repeatable (20261112) where bucket = any($1) and window_count > 0",
    () => new Moments(), (m, id, bucket, profile) =>
    {
        var squared = new double[64 * F];
        var candidate = taxa[id];
        for (int g = 0; g < profile.Genes.Count; g++)
        {
            var gene = profile.Genes[g];
            var windows = profile.WindowsOf(g);
            for (int i = 0; i < windows.Length; i++) { m.Sum[i % F] += windows[i]; m.SumSquares[i % F] += windows[i] * windows[i]; }
            m.Windows += windows.Length / F;
            foreach (var query in byBucket[bucket])
            {
                if (!query.Tuning || !Comparable(query, gene) || !OutsideGenus(query, candidate)) continue;
                if (squared.Length < windows.Length) squared = new double[windows.Length];
                int count = Differences(query.Axes, windows, squared);
                var (sums, n) = gene.Family == query.Truth.Family ? (m.Within, 0) : (m.Between, 1);
                for (int i = 0; i < count * F; i++) sums[i % F] += squared[i];
                m.Counts[n] += count;
            }
        }
    });
var pooled = new Moments();
foreach (var m in sample) pooled.Add(m);
var variance = Enumerable.Range(0, F).Select(f => pooled.SumSquares[f] / pooled.Windows - Math.Pow(pooled.Sum[f] / pooled.Windows, 2)).ToArray();
var within = Enumerable.Range(0, F).Select(f => pooled.Within[f] / pooled.Counts[0]).ToArray();
var between = Enumerable.Range(0, F).Select(f => pooled.Between[f] / pooled.Counts[1]).ToArray();
Console.WriteLine($"Pre-pass: {pooled.Windows:N0} windows, {pooled.Counts[0]:N0} within / {pooled.Counts[1]:N0} between window pairs, " +
    $"{total.Elapsed.TotalMinutes:F1} min.");
Console.WriteLine($"  variance {Vector(variance)}\n  within {Vector(within)}\n  between {Vector(between)}");

weightings =
[
    ("unweighted", Enumerable.Repeat(1.0, F).ToArray()),
    ("P0", [1, 1, 1, 0, 0, 0, 0, 0, 0]),
    ("P0+P1", [1, 1, 1, 1, 1, 1, 0, 0, 0]),
    ("1/variance", Normalized(variance.Select(v => v > 0 ? 1 / v : 0).ToArray())),
    ("1/within", Normalized(within.Select(v => v > 0 ? 1 / v : 0).ToArray())),
    ("between/within", Normalized(within.Select((v, f) => v > 0 ? between[f] / v : 0).ToArray())),
];
}
int W = weightings.Length;
foreach (var query in queries)
{
    var squared = new double[query.TruthWindows.Length];
    int count = Differences(query.Axes, query.TruthWindows, squared);
    query.Own = weightings.Select(w => count == 0 ? double.NaN : Score(squared, count, w.Weights)).ToArray();
}

// Main pass: every gene of every needed bucket, every weighting.
var main = await PassAsync("select sequence_id, bucket, payload from gene_profile where bucket = any($1) and window_count > 0",
    () => new Tallies(queries.Count, W), (t, id, bucket, profile) =>
    {
        var squared = new double[64 * F];
        var candidate = taxa[id];
        var relevant = byBucket[bucket];
        for (int g = 0; g < profile.Genes.Count; g++)
        {
            var gene = profile.Genes[g];
            double[]? windows = null;
            foreach (var query in relevant)
            {
                if (!Comparable(query, gene) || IsSelf(query, id, gene)) continue;
                windows ??= profile.WindowsOf(g);
                if (squared.Length < windows.Length) squared = new double[windows.Length];
                int count = Differences(query.Axes, windows, squared);
                if (count == 0) continue;
                bool outside = OutsideGenus(query, candidate);
                bool same = gene.Family == query.Truth.Family;
                for (int v = 0; v < W; v++)
                {
                    double score = Score(squared, count, weightings[v].Weights);
                    ref double best = ref t.BestAll[query.Index * W + v];
                    if (score > best) best = score;
                    if (!outside) continue;
                    var top = t.Top[query.Index * W + v];
                    int k = v == 0 ? diagnosticK : topK;
                    if (top.Count < k) top.Enqueue((gene.Family, gene.Residues), score);
                    else if (top.TryPeek(out _, out double lowest) && score > lowest) top.EnqueueDequeue((gene.Family, gene.Residues), score);
                    if (v == 0 && same)
                    {
                        t.SameFamily[query.Index]++;
                        if (score > t.BestSame[query.Index]) (t.BestSame[query.Index], t.BestSameResidues[query.Index]) = (score, gene.Residues);
                    }
                }
            }
        }
    });
Console.WriteLine($"Main pass done, {total.Elapsed.TotalMinutes:F1} min.");

// Merge the workers.
var aggregators = new (string Name, Func<List<(int Family, int Residues)>, int> Pick)[]
{
    ("best hit", hits => hits.Count > 0 ? hits[0].Family : ZCurveGeneProfile.NoFamily),
    ("best with family", hits => hits.Select(h => h.Family).FirstOrDefault(f => f != ZCurveGeneProfile.NoFamily, ZCurveGeneProfile.NoFamily)),
    ("vote 3", hits => Vote(hits, 3)),
    ("vote 5", hits => Vote(hits, 5)),
    ("vote 10", hits => Vote(hits, 10)),
};
static int Vote(List<(int Family, int Residues)> hits, int k)
{
    var known = hits.Where(h => h.Family != ZCurveGeneProfile.NoFamily).Take(k).Select((h, i) => (h.Family, i)).ToList();
    return known.Count == 0 ? ZCurveGeneProfile.NoFamily
        : known.GroupBy(h => h.Family).OrderByDescending(g => g.Count()).ThenBy(g => g.Min(h => h.i)).First().Key;
}
int A = aggregators.Length;
var correct = new bool[queries.Count, W, A];
var exact = new bool[queries.Count, W];
var reachable = new bool[queries.Count];
var diagnostics = new (int Rank, double Gap, double LengthRatio)[queries.Count];
foreach (var query in queries)
{
    int q = query.Index;
    reachable[q] = main.Sum(t => t.SameFamily[q]) > 0;
    for (int v = 0; v < W; v++)
    {
        double bestAll = main.Max(t => t.BestAll[q * W + v]);
        exact[q, v] = query.Own[v] >= bestAll;
        var hits = main.SelectMany(t => t.Top[q * W + v].UnorderedItems).OrderByDescending(e => e.Priority).Select(e => e.Element).ToList();
        for (int a = 0; a < A; a++) correct[q, v, a] = aggregators[a].Pick(hits) == query.Truth.Family;
        if (v == 0)
        {
            int rank = hits.FindIndex(h => h.Family == query.Truth.Family);
            double topScore = main.SelectMany(t => t.Top[q * W].UnorderedItems).Select(e => e.Priority).DefaultIfEmpty(double.NaN).Max();
            int best = Enumerable.Range(0, main.Length).MaxBy(w => main[w].BestSame[q]);
            diagnostics[q] = (rank < 0 ? -1 : rank + 1, topScore - main[best].BestSame[q], (double)main[best].BestSameResidues[q] / query.Residues);
        }
    }
}

// Summary per half and band, weighting and aggregator.
string Percent(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
var lines = new List<string> { "Half,Band,Weighting,Choice,Reachable,Family top-1,Proteins,Exact top-1" };
foreach (var half in new[] { "tuning", "test", "all" })
    foreach (int band in Enumerable.Range(0, bandNames.Length).Prepend(-1))
    {
        var cell = queries.Where(q => (half == "all" || q.Tuning == (half == "tuning")) && (band < 0 || q.Band == band)).ToList();
        var reach = cell.Where(q => reachable[q.Index]).ToList();
        for (int v = 0; v < W; v++)
            for (int a = 0; a < A; a++)
                lines.Add(string.Join(",", half, band < 0 ? "all" : bandNames[band], weightings[v].Name, aggregators[a].Name, reach.Count,
                    Percent(reach.Count(q => correct[q.Index, v, a]), reach.Count), cell.Count, Percent(cell.Count(q => exact[q.Index, v]), cell.Count)));
    }

// Rule: best on the tuning half, judged on the test half against (unweighted, best hit), genome bootstrap. In --confirm, the fixed
// weights (best hit) against unweighted (best hit) on every query.
double Share(IEnumerable<Query> set, Func<Query, bool> hit) { var l = set.ToList(); return l.Count == 0 ? double.NaN : l.Count(hit) / (double)l.Count; }
var tuningReach = queries.Where(q => q.Tuning && reachable[q.Index]).ToList();
var (bestV, bestA) = confirm ? (1, 0) : (from v in Enumerable.Range(0, W) from a in Enumerable.Range(0, A) select (v, a))
    .MaxBy(c => Share(tuningReach, q => correct[q.Index, c.v, c.a]));
var testGenomes = queries.Where(q => !q.Tuning).GroupBy(q => q.Genome).Select(g => g.ToList()).ToList();
double Delta(IEnumerable<Query> set, Func<Query, bool> after, Func<Query, bool> before) =>
    Share(set, after) - Share(set, before);
var testReach = queries.Where(q => !q.Tuning && reachable[q.Index]).ToList();
var testAll = queries.Where(q => !q.Tuning).ToList();
double familyDelta = Delta(testReach, q => correct[q.Index, bestV, bestA], q => correct[q.Index, 0, 0]);
double exactDelta = Delta(testAll, q => exact[q.Index, bestV], q => exact[q.Index, 0]);
var rng = new Random(seed);
var deltas = new List<double>();
for (int r = 0; r < bootstrapRounds; r++)
{
    // Materialized: a lazy sample would be drawn again for each side of the difference, which unpairs it (fixed 2026-11-20; G4's
    // first CI, -8.7 to +15.3, came from that bug).
    var resampled = Enumerable.Range(0, testGenomes.Count).SelectMany(_ => testGenomes[rng.Next(testGenomes.Count)]).Where(q => reachable[q.Index]).ToList();
    double d = Delta(resampled, q => correct[q.Index, bestV, bestA], q => correct[q.Index, 0, 0]);
    if (!double.IsNaN(d)) deltas.Add(d);
}
deltas.Sort();
double low = deltas[(int)(0.025 * deltas.Count)], high = deltas[(int)(0.975 * deltas.Count) - 1];
string Points(double x) => (100 * x).ToString("+0.0;-0.0", CultureInfo.InvariantCulture);
lines.Add("");
if (confirm)
{
    // Exact top-1 over every query, with the same genome bootstrap.
    var exactDeltas = new List<double>();
    var exactRng = new Random(confirmSeed);
    for (int r = 0; r < bootstrapRounds; r++)
    {
        var resampled = Enumerable.Range(0, testGenomes.Count).SelectMany(_ => testGenomes[exactRng.Next(testGenomes.Count)]).ToList();
        exactDeltas.Add(Delta(resampled, q => exact[q.Index, 1], q => exact[q.Index, 0]));
    }
    exactDeltas.Sort();
    double exactLow = exactDeltas[(int)(0.025 * exactDeltas.Count)], exactHigh = exactDeltas[(int)(0.975 * exactDeltas.Count) - 1];
    bool confirmed = exactDelta >= 0.02 && exactLow > 0 && familyDelta >= -0.01;
    lines.Add($"Weights,{Vector(weightings[1].Weights)} (fixed from G4)");
    lines.Add($"Exact top-1 unweighted -> weighted,{100 * Share(testAll, q => exact[q.Index, 0]):F1}% -> {100 * Share(testAll, q => exact[q.Index, 1]):F1}%: " +
        $"{Points(exactDelta)} points (95% CI {Points(exactLow)} to {Points(exactHigh)}; {testAll.Count} proteins in {testGenomes.Count} genomes)");
    lines.Add($"Family top-1 (own genus removed) unweighted -> weighted,{100 * Share(testReach, q => correct[q.Index, 0, 0]):F1}% -> " +
        $"{100 * Share(testReach, q => correct[q.Index, 1, 0]):F1}%: {Points(familyDelta)} points (95% CI {Points(low)} to {Points(high)}; {testReach.Count} reachable)");
    lines.Add("Rule (exact >= +2 points with CI above 0, family >= -1 point)," + (confirmed ? "go: weights confirmed" : "stop: keep unweighted"));
}
else
{
    bool go = familyDelta >= 0.05 && low > 0 && exactDelta >= -0.01;
    lines.Add($"Weights 1/variance,{Vector(weightings[3].Weights)}");
    lines.Add($"Weights 1/within,{Vector(weightings[4].Weights)}");
    lines.Add($"Weights between/within,{Vector(weightings[5].Weights)}");
    lines.Add($"Best on tuning half,{weightings[bestV].Name} / {aggregators[bestA].Name}");
    lines.Add($"Test half vs unweighted / best hit,family top-1 {Points(familyDelta)} points (95% CI {Points(low)} to {Points(high)}; {testReach.Count} proteins in " +
        $"{testGenomes.Count} genomes); exact top-1 {Points(exactDelta)} points");
    lines.Add("Rule (>= +5 points, CI excludes 0, exact >= -1 point)," + (go ? "go" : "stop: report today's number"));
}

// Diagnostic of today's misses (unweighted, best hit), all queries.
var misses = queries.Where(q => reachable[q.Index] && !correct[q.Index, 0, 0]).ToList();
lines.Add("");
lines.Add($"Today's misses,{misses.Count} of {queries.Count(q => reachable[q.Index])} reachable");
lines.Add($"  first same-family gene at rank 2-10,{misses.Count(q => diagnostics[q.Index].Rank is >= 2 and <= 10)}");
lines.Add($"  at rank 11-{diagnosticK},{misses.Count(q => diagnostics[q.Index].Rank > 10)}");
lines.Add($"  below rank {diagnosticK},{misses.Count(q => diagnostics[q.Index].Rank < 0)}");
lines.Add($"  top hit has no family,{misses.Count(q => main.SelectMany(t => t.Top[q.Index * W].UnorderedItems).MaxBy(e => e.Priority).Element.Family == ZCurveGeneProfile.NoFamily)}");
foreach (int band in Enumerable.Range(0, bandNames.Length))
{
    var bandMisses = misses.Where(q => q.Band == band).Select(q => diagnostics[q.Index]).ToList();
    var bandHits = queries.Where(q => q.Band == band && reachable[q.Index] && correct[q.Index, 0, 0]).Select(q => diagnostics[q.Index]).ToList();
    string Median(IEnumerable<double> xs) { var l = xs.Order().ToList(); return l.Count == 0 ? "" : l[l.Count / 2].ToString("F3", CultureInfo.InvariantCulture); }
    lines.Add($"  {bandNames[band]},misses {bandMisses.Count}: median length ratio of best same-family gene {Median(bandMisses.Select(d => d.LengthRatio))} " +
        $"(hits {Median(bandHits.Select(d => d.LengthRatio))}); off by > 10% {bandMisses.Count(d => Math.Abs(d.LengthRatio - 1) > 0.1)}; median score gap {Median(bandMisses.Select(d => d.Gap))}");
}
lines.Add($"Total,{total.Elapsed.TotalMinutes:F1} min");
File.WriteAllLines(summaryPath, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines.Where(l => !l.StartsWith("tuning,") && !l.StartsWith("test,") && !(l.StartsWith("all,") && !l.StartsWith("all,all,")))));

var perQuery = new List<string> { "Genome,Half,Band,Residues,Family,Reachable,Rank first same family,Score gap,Length ratio," +
    string.Join(",", from v in Enumerable.Range(0, W) from a in Enumerable.Range(0, A) select $"{weightings[v].Name} / {aggregators[a].Name}") };
foreach (var q in queries)
    perQuery.Add(string.Join(",", q.Genome, q.Tuning ? "tuning" : "test", bandNames[q.Band], q.Residues, q.Truth.Family, reachable[q.Index] ? 1 : 0,
        diagnostics[q.Index].Rank, diagnostics[q.Index].Gap.ToString("F4", CultureInfo.InvariantCulture),
        diagnostics[q.Index].LengthRatio.ToString("F3", CultureInfo.InvariantCulture),
        string.Join(",", from v in Enumerable.Range(0, W) from a in Enumerable.Range(0, A) select correct[q.Index, v, a] ? 1 : 0)));
File.WriteAllLines(queriesPath, perQuery);

record Taxon(long Source, string? Genus);

sealed class Query(string genome, bool tuning, int band, int residues, long sequenceId, Taxon taxon, ZCurveGeneProfile.Gene truth, double[] axes,
    double[] truthWindows, HashSet<int> buckets)
{
    private static int next;
    public int Index { get; } = next++;
    public string Genome => genome;
    public bool Tuning => tuning;
    public int Band => band;
    public int Residues => residues;
    public long SequenceId => sequenceId;
    public Taxon Taxon => taxon;
    public ZCurveGeneProfile.Gene Truth => truth;
    public double[] Axes => axes;
    public double[] TruthWindows => truthWindows;
    public HashSet<int> Buckets => buckets;
    public double[] Own { get; set; } = [];
}

sealed class Moments
{
    public readonly double[] Sum = new double[ZCurveGeneProfile.FeatureCount], SumSquares = new double[ZCurveGeneProfile.FeatureCount];
    public readonly double[] Within = new double[ZCurveGeneProfile.FeatureCount], Between = new double[ZCurveGeneProfile.FeatureCount];
    public readonly long[] Counts = new long[2];
    public long Windows;

    public void Add(Moments other)
    {
        for (int f = 0; f < Sum.Length; f++)
        {
            Sum[f] += other.Sum[f]; SumSquares[f] += other.SumSquares[f]; Within[f] += other.Within[f]; Between[f] += other.Between[f];
        }
        Counts[0] += other.Counts[0]; Counts[1] += other.Counts[1]; Windows += other.Windows;
    }
}

sealed class Tallies
{
    public readonly double[] BestAll, BestSame;
    public readonly int[] SameFamily, BestSameResidues;
    public readonly PriorityQueue<(int Family, int Residues), double>[] Top;

    public Tallies(int queries, int weightings)
    {
        BestAll = Enumerable.Repeat(double.NegativeInfinity, queries * weightings).ToArray();
        BestSame = Enumerable.Repeat(double.NegativeInfinity, queries).ToArray();
        SameFamily = new int[queries];
        BestSameResidues = new int[queries];
        Top = Enumerable.Range(0, queries * weightings).Select(_ => new PriorityQueue<(int, int), double>()).ToArray();
    }
}

sealed class MomentTallies
{
    public readonly double[] Best;
    public readonly long[] Tied, TiedRight;
    public readonly int[] SameFamily;
    public readonly long[] Histogram, Undefined;

    public MomentTallies(int queries, int subScores, int bands, int bins)
    {
        Best = Enumerable.Repeat(double.NegativeInfinity, queries * subScores).ToArray();
        Tied = new long[queries * subScores];
        TiedRight = new long[queries * subScores];
        SameFamily = new int[queries];
        Histogram = new long[bands * subScores * 2 * bins];
        Undefined = new long[subScores];
    }
}
