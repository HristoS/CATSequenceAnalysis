#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Text;
using System.Text.Json;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Exact-gene miss check (2026-11-20). From the protein alone, the query's own gene is first for 85-87% of proteins and in the top 10 for
// 95-100% (ProteinSearch_CatalogueV2.csv; ProteinWeights_Confirm.csv). Ties don't explain it: the checks count a tie for the own gene, so
// every miss is a gene that scores strictly higher. Which genes are they? For G4's 1,143 proteins (protein_family_diagnostic.cs's rule
// and seed), every comparable gene of the catalogue is scored as ProteinSearch scores it (unweighted; from 500 aa against genes of 500 aa
// and up over their common windows, shorter within 10% of the length). Per query: the own gene's score, how many genes beat it, exact
// ties, and the top 10 winners. For each miss, the top winner's protein is fetched from BV-BRC (proteins, never stored) and aligned with
// the query (global, free end gaps, BLOSUM62), and the miss is classed, in this order:
//   identical protein   identity >= 99% over >= 95% of the longer protein;
//   near-identical      identity >= 90% over >= 90% of the longer;
//   shorter partial     500+ aa only: the winner is under 90% of the query's length and compared on fewer windows than the query has;
//   homolog             identity >= Rost's HSSP curve + 5 points over >= 60% of the shorter;
//   unrelated           otherwise ("no protein" when the winner isn't a BV-BRC gene).
// and by where the winner lies: same genome (a paralog), same species, same genus, other.
// Rule fixed before the run: identical + near-identical >= half of the misses -> the paper reports "right protein first" next to "own
// gene first"; shorter partial >= a third of the 500+ misses -> a card for penalising unmatched windows, measured on these queries.
// Args: output prefix (default ProteinExactMiss), writing <prefix>_Misses.csv and <prefix>_Summary.csv.
const int seed = 20261105, perBand = 5, maxAlignResidues = 6000, topWinners = 10;
string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 48, ["archaeon"] = 8, ["plasmid"] = 8 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
const int F = ZCurveGeneProfile.FeatureCount;
var prefix = args.Length > 0 ? args[0] : "ProteinExactMiss";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
using var http = new HttpClient { Timeout = TimeSpan.FromMinutes(5) };
http.DefaultRequestHeaders.UserAgent.ParseAdd("CATSequenceAnalysis-platform-check/1.0");
http.DefaultRequestHeaders.Accept.ParseAdd("application/json");
var total = Stopwatch.StartNew();
CheckBlosum();

async Task<JsonElement> GetAsync(string url)
{
    await Task.Delay(1000);
    using var response = await http.GetAsync(url);
    response.EnsureSuccessStatusCode();
    using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
    return document.RootElement.Clone();
}

var records = new Dictionary<long, RecordInfo>();
await using (var command = dataSource.CreateCommand(
    "select r.id, r.source_id, r.tax_genus, s.provider, s.accession, r.record_id, r.tax_order from sequence_record r join sequence_source s on s.id = r.source_id"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
        records[reader.GetInt64(0)] = new RecordInfo(reader.GetInt64(1), reader.IsDBNull(2) ? null : reader.GetString(2), reader.GetString(3),
            reader.IsDBNull(4) ? null : reader.GetString(4), reader.GetString(5), reader.IsDBNull(6) ? null : reader.GetString(6));

// Query selection: protein_family_diagnostic.cs unchanged (same genomes, same proteins, same order).
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
        var byRecord = genome.ToDictionary(s => s.Record);
        JsonElement features;
        try
        {
            features = await GetAsync($"https://www.bv-brc.org/api/genome_feature/?and(eq(genome_id,{genome.Key}),eq(feature_type,CDS)," +
                $"eq(annotation,PATRIC),ge(aa_length,{ZCurveGeneProfile.MinResidues}))&select(accession,start,end,strand,aa_length,aa_sequence_md5)&sort(+feature_id)&limit(25000)");
        }
        catch (HttpRequestException e) { Console.WriteLine($"{kind} {genome.Key}: skipped ({e.Message})"); continue; }
        var picker = new Random(seed);
        var picked = features.EnumerateArray()
            .Where(f => f.TryGetProperty("aa_sequence_md5", out _) && byRecord.ContainsKey(f.GetProperty("accession").GetString()!))
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
            var record = byRecord[feature.GetProperty("accession").GetString()!];
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
            if (profile.Genes[gene].Family == ZCurveGeneProfile.NoFamily || records[record.Id].Genus is null) continue;
            var code = record.Order is not null && mycoplasmaOrders.Contains(record.Order) ? ZCurveGeneProfile.mycoplasmaCode : ZCurveGeneProfile.standardCode;
            queries.Add(new Query(queries.Count, genome.Key, band, protein, code, record.Id, records[record.Id], profile.Genes[gene],
                ZCurveGeneProfile.QueryFromProtein(protein, null, code), GeneProfiles.BucketsFor(protein.Length, ProteinSearch.LengthBand).ToHashSet()));
        }
    }
}
Console.WriteLine($"{queries.Count} proteins with a PGfam and a genus (G4 had 1,143). {total.Elapsed.TotalMinutes:F1} min.");

// ---- the exact-miss check ----
// Two scores: today's (unweighted) and G4's fixed 1/variance weights (added on request before the run; the rule is judged unweighted).
var weightings = new (string Name, double[] Weights)[]
{
    ("unweighted", Enumerable.Repeat(1.0, F).ToArray()),
    ("1/variance", [0.741, 0.838, 0.333, 1.17, 1.16, 1.15, 1.19, 1.21, 1.2]),
};
int W = weightings.Length;
static double Score(double[] query, double[] gene, double[] weights)
{
    int count = Math.Min(query.Length, gene.Length) / F;
    if (count == 0) return double.NaN;
    double sum = 0;
    for (int w = 0; w < count; w++)
    {
        double s = 0;
        for (int f = 0; f < F; f++) { double d = query[w * F + f] - gene[w * F + f]; s += weights[f] * d * d; }
        sum += Math.Sqrt(s);
    }
    return 1 - sum / count;
}

var speciesOf = new Dictionary<long, string?>();
await using (var command = dataSource.CreateCommand("select id, tax_species from sequence_record"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync()) speciesOf[reader.GetInt64(0)] = reader.IsDBNull(1) ? null : reader.GetString(1);

// The own gene's score per weighting, as the pass scores every other gene; unweighted it must equal ZCurveGeneProfile.Score.
var own = new double[queries.Count, W];
foreach (var group in queries.GroupBy(q => q.SequenceId))
{
    var profile = (await GeneProfiles.ReadAsync(dataSource, group.Key))!;
    foreach (var query in group)
    {
        int gene = Enumerable.Range(0, profile.Genes.Count).First(g => profile.Genes[g].Start == query.Truth.Start && profile.Genes[g].End == query.Truth.End
            && profile.Genes[g].Reverse == query.Truth.Reverse);
        var windows = profile.WindowsOf(gene);
        for (int v = 0; v < W; v++) own[query.Index, v] = Score(query.Axes, windows, weightings[v].Weights);
        if (Math.Abs(own[query.Index, 0] - profile.Score(query.Axes, gene)) > 1e-12) throw new InvalidOperationException("Own score differs from ZCurveGeneProfile.Score.");
    }
}

bool Comparable(Query query, ZCurveGeneProfile.Gene gene) => query.Protein.Length >= ZCurveGeneProfile.LongResidues
    ? gene.Residues >= ZCurveGeneProfile.LongResidues
    : gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - query.Protein.Length) <= ProteinSearch.LengthBand * query.Protein.Length;
bool IsSelf(Query query, long id, ZCurveGeneProfile.Gene gene) =>
    id == query.SequenceId && gene.Start == query.Truth.Start && gene.End == query.Truth.End && gene.Reverse == query.Truth.Reverse;
var byBucket = Enumerable.Range(0, 64).Select(b => queries.Where(q => q.Buckets.Contains(b)).ToArray()).ToArray();
var needed = queries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();

var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
{
    var tally = new Tally(queries.Count, W);
    foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable())
    {
        var profile = ZCurveGeneProfile.FromBinary(payload);
        for (int g = 0; g < profile.Genes.Count; g++)
        {
            var gene = profile.Genes[g];
            double[]? windows = null;
            foreach (var query in byBucket[bucket])
            {
                if (!Comparable(query, gene) || IsSelf(query, id, gene)) continue;
                windows ??= profile.WindowsOf(g);
                for (int v = 0; v < W; v++)
                {
                    double score = Score(query.Axes, windows, weightings[v].Weights);
                    if (double.IsNaN(score)) continue;
                    double mine = own[query.Index, v];
                    if (score == mine) tally.Ties[query.Index, v]++;
                    if (score <= mine) continue;
                    tally.Beaten[query.Index, v]++;
                    var top = tally.Winners[query.Index, v];
                    var hit = new Hit(id, gene, score, []);
                    if (top.Count < topWinners) top.Enqueue(hit, score);
                    else if (top.TryPeek(out var peeked, out double lowest) && score > lowest) top.EnqueueDequeue(hit, score);
                }
            }
        }
    }
    return tally;
})).ToArray();
await using (var command = dataSource.CreateCommand("select sequence_id, bucket, payload from gene_profile where bucket = any($1) and window_count > 0"))
{
    command.Parameters.AddWithValue(needed);
    command.CommandTimeout = 0;
    await using var reader = await command.ExecuteReaderAsync();
    long read = 0;
    while (await reader.ReadAsync())
    {
        rows.Add((reader.GetInt64(0), reader.GetInt16(1), reader.GetFieldValue<byte[]>(2)));
        if (++read % 20_000 == 0) Console.WriteLine($"  {read:N0} rows, {total.Elapsed.TotalMinutes:F1} min");
    }
}
rows.CompleteAdding();
var perWorker = await Task.WhenAll(workers);
var beatenBy = new int[queries.Count, W];
var tiesOf = new int[queries.Count, W];
var winnersOf = new List<Hit>[queries.Count, W];
foreach (var query in queries)
    for (int v = 0; v < W; v++)
    {
        beatenBy[query.Index, v] = perWorker.Sum(t => t.Beaten[query.Index, v]);
        tiesOf[query.Index, v] = perWorker.Sum(t => t.Ties[query.Index, v]);
        winnersOf[query.Index, v] = perWorker.SelectMany(t => t.Winners[query.Index, v].UnorderedItems).OrderByDescending(e => e.Priority)
            .Select(e => e.Element).Take(topWinners).ToList();
    }
Console.WriteLine($"Pass done, {total.Elapsed.TotalMinutes:F1} min. Exact misses: " +
    string.Join(", ", Enumerable.Range(0, W).Select(v => $"{weightings[v].Name} {queries.Count(q => beatenBy[q.Index, v] > 0)}")));

var geneProteins = new Dictionary<(long, int, int, bool), string?>();
async Task<string?> GeneProteinAsync(Hit hit)
{
    var key = (hit.SequenceId, hit.Gene.Start, hit.Gene.End, hit.Gene.Reverse);
    if (geneProteins.TryGetValue(key, out var known)) return known;
    var record = records[hit.SequenceId];
    string? protein = null;
    if (record.Provider == "BV-BRC" && record.Accession is not null && hit.Gene.Chromosome == 0)
    {
        try
        {
            var found = await GetAsync($"https://www.bv-brc.org/api/genome_feature/?and(eq(genome_id,{record.Accession}),eq(accession,{Uri.EscapeDataString(record.RecordId)})," +
                $"eq(start,{hit.Gene.Start + 1}),eq(end,{hit.Gene.End}),eq(feature_type,CDS),eq(annotation,PATRIC))&select(strand,aa_sequence_md5)&limit(5)");
            var md5 = found.EnumerateArray().Where(f => (f.GetProperty("strand").GetString() == "-") == hit.Gene.Reverse && f.TryGetProperty("aa_sequence_md5", out _))
                .Select(f => f.GetProperty("aa_sequence_md5").GetString()).FirstOrDefault();
            if (md5 is not null)
                protein = (await GetAsync($"https://www.bv-brc.org/api/feature_sequence/?eq(md5,{md5})&select(sequence)&limit(1)"))
                    .EnumerateArray().Select(p => p.GetProperty("sequence").GetString()!.TrimEnd('*')).FirstOrDefault();
        }
        catch (HttpRequestException e) { Console.WriteLine($"  protein of record {hit.SequenceId} at {hit.Gene.Start}: {e.Message}"); }
    }
    return geneProteins[key] = protein;
}

string Where(Query query, long id) => records[id].Source == query.Record.Source ? "same genome"
    : speciesOf[id] is { } s && s == speciesOf[query.SequenceId] ? "same species"
    : records[id].Genus is { } g && g == query.Record.Genus ? "same genus" : "other";
string F4(double x) => double.IsNaN(x) ? "" : x.ToString("F4", CultureInfo.InvariantCulture);
var missLines = new List<string> { "Weighting,Query genome,Band,Query aa,Own score,Own rank,Exact ties,Winner score,Gap,Winner record,Where,Winner aa," +
    "Windows compared,Query windows,Identity,Cover of longer,Class,Winner family,Query family" };
var classOf = new Dictionary<(int, int), string>();
var whereOf = new Dictionary<(int, int), string>();
for (int v = 0; v < W; v++)
    foreach (var query in queries.Where(q => beatenBy[q.Index, v] > 0))
    {
        var winner = winnersOf[query.Index, v][0];
        var protein = await GeneProteinAsync(winner);
        int queryWindows = query.Axes.Length / F, compared = Math.Min(queryWindows, winner.Gene.WindowCount);
        double identity = double.NaN, coverLonger = double.NaN;
        string cls = "no protein";
        if (protein is not null && query.Protein.Length <= maxAlignResidues && protein.Length <= maxAlignResidues)
        {
            var align = Alignment(query.Protein, protein);
            int shorter = Math.Min(query.Protein.Length, protein.Length), longer = Math.Max(query.Protein.Length, protein.Length);
            double pairs = double.IsNaN(align.AlignedShare) ? 0 : align.AlignedShare * shorter;
            identity = align.Identity; coverLonger = pairs / longer;
            cls = identity >= 0.99 && coverLonger >= 0.95 ? "identical protein"
                : identity >= 0.90 && coverLonger >= 0.90 ? "near-identical"
                : query.Band == 4 && protein.Length < 0.9 * query.Protein.Length && compared < queryWindows ? "shorter partial"
                : pairs >= 0.6 * shorter && 100 * identity >= Hssp((int)pairs) + 5 ? "homolog"
                : "unrelated";
        }
        else if (protein is not null) cls = "too long to align";
        classOf[(query.Index, v)] = cls;
        whereOf[(query.Index, v)] = Where(query, winner.SequenceId);
        missLines.Add(string.Join(",", weightings[v].Name, query.Genome, bandNames[query.Band], query.Protein.Length, F4(own[query.Index, v]), 1 + beatenBy[query.Index, v],
            tiesOf[query.Index, v], F4(winner.Score), F4(winner.Score - own[query.Index, v]), winner.SequenceId, whereOf[(query.Index, v)], winner.Gene.Residues,
            compared, queryWindows, F4(identity), F4(coverLonger), cls, winner.Gene.Family, query.Truth.Family));
        if (classOf.Count % 25 == 0) Console.WriteLine($"  {classOf.Count} misses classed, {total.Elapsed.TotalMinutes:F1} min");
    }
File.WriteAllLines($"{prefix}_Misses.csv", missLines);

// Summary per weighting and band.
string Pct(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
string[] classes = ["identical protein", "near-identical", "shorter partial", "homolog", "unrelated", "no protein", "too long to align"];
string[] places = ["same genome", "same species", "same genus", "other"];
var summary = new List<string> { "Weighting,Band,Queries,Own gene first,Own gene in top 10,Right protein first,With exact ties,Misses," +
    string.Join(",", classes) + "," + string.Join(",", places) };
for (int v = 0; v < W; v++)
    foreach (int band in Enumerable.Range(0, bandNames.Length).Append(-1))
    {
        int vv = v;
        var cell = queries.Where(q => band < 0 || q.Band == band).ToList();
        var missed = cell.Where(q => beatenBy[q.Index, vv] > 0).ToList();
        int rightProtein = cell.Count(q => beatenBy[q.Index, vv] == 0 || classOf.GetValueOrDefault((q.Index, vv)) is "identical protein" or "near-identical");
        summary.Add(string.Join(",", weightings[v].Name, band < 0 ? "all" : bandNames[band], cell.Count, Pct(cell.Count(q => beatenBy[q.Index, vv] == 0), cell.Count),
            Pct(cell.Count(q => beatenBy[q.Index, vv] < 10), cell.Count), Pct(rightProtein, cell.Count), cell.Count(q => tiesOf[q.Index, vv] > 0), missed.Count,
            string.Join(",", classes.Select(c => missed.Count(q => classOf[(q.Index, vv)] == c))),
            string.Join(",", places.Select(p => missed.Count(q => whereOf[(q.Index, vv)] == p)))));
    }
// What the weights change: misses fixed and misses made, by today's class.
var fixedByWeights = queries.Where(q => beatenBy[q.Index, 0] > 0 && beatenBy[q.Index, 1] == 0).ToList();
var madeByWeights = queries.Where(q => beatenBy[q.Index, 0] == 0 && beatenBy[q.Index, 1] > 0).ToList();
summary.Add("");
summary.Add($"Weights fix (missed unweighted; own gene first weighted),{fixedByWeights.Count}: " +
    string.Join(" / ", classes.Select(c => $"{c} {fixedByWeights.Count(q => classOf[(q.Index, 0)] == c)}").Where(s => !s.EndsWith(" 0"))));
summary.Add($"Weights break (first unweighted; missed weighted),{madeByWeights.Count}: " +
    string.Join(" / ", classes.Select(c => $"{c} {madeByWeights.Count(q => classOf[(q.Index, 1)] == c)}").Where(s => !s.EndsWith(" 0"))));
var misses = queries.Where(q => beatenBy[q.Index, 0] > 0).ToList();
var longMisses = misses.Where(q => q.Band == 4).ToList();
int sameProtein = misses.Count(q => classOf[(q.Index, 0)] is "identical protein" or "near-identical");
int partial = longMisses.Count(q => classOf[(q.Index, 0)] == "shorter partial");
summary.Add($"Rule (unweighted): identical + near-identical >= half of the misses,{sameProtein} of {misses.Count} ({Pct(sameProtein, misses.Count)}): " +
    (2 * sameProtein >= misses.Count ? "report 'right protein first' beside 'own gene first'" : "no"));
summary.Add($"Rule (unweighted): shorter partial >= a third of the 500+ misses,{partial} of {longMisses.Count} ({Pct(partial, longMisses.Count)}): " +
    (3 * partial >= longMisses.Count && longMisses.Count > 0 ? "card for penalising unmatched windows" : "no"));
summary.Add($"Queries whose own gene ties with another gene (unweighted),{queries.Count(q => tiesOf[q.Index, 0] > 0)} (ties count for the own gene here; the API breaks them by record id)");
summary.Add($"Total,{total.Elapsed.TotalMinutes:F1} min");
File.WriteAllLines($"{prefix}_Summary.csv", summary);
Console.WriteLine(string.Join(Environment.NewLine, summary));

// Rost (1999) HSSP curve: the identity (%) above which an alignment of L residues is in the safe zone of homology.
static double Hssp(int length) => length <= 11 ? 100 : length > 450 ? 19.5 : 480 * Math.Pow(length, -0.32 * (1 + Math.Exp(-length / 1000.0)));

// Global alignment with free end gaps (Gotoh, BLOSUM62, a gap of k residues costs 11 + k).
static Align Alignment(string query, string gene)
{
    const int open = 12, extend = 1, negative = int.MinValue / 4;
    int n = query.Length, m = gene.Length;
    var trace = new byte[(n + 1) * (m + 1)];
    var mPrev = new int[m + 1]; var xPrev = new int[m + 1]; var yPrev = new int[m + 1];
    var mCur = new int[m + 1]; var xCur = new int[m + 1]; var yCur = new int[m + 1];
    for (int j = 0; j <= m; j++) { mPrev[j] = j == 0 ? 0 : negative; xPrev[j] = negative; yPrev[j] = j == 0 ? negative : 0; }
    int bestScore = negative, bestI = 0, bestJ = 0, bestState = 0;
    void Consider(int score, int i, int j, int state) { if (score > bestScore) (bestScore, bestI, bestJ, bestState) = (score, i, j, state); }
    for (int i = 1; i <= n; i++)
    {
        mCur[0] = negative; xCur[0] = 0; yCur[0] = negative;
        int qi = Blosum.Index(query[i - 1]);
        for (int j = 1; j <= m; j++)
        {
            byte t = 0;
            int fromM = mPrev[j - 1], fromX = xPrev[j - 1], fromY = yPrev[j - 1];
            int diagonal = fromM; if (fromX > diagonal) { diagonal = fromX; t = 1; } if (fromY > diagonal) { diagonal = fromY; t = 2; }
            mCur[j] = diagonal + Blosum.Score(qi, Blosum.Index(gene[j - 1]));
            int xOpen = mPrev[j] - open, xExtend = xPrev[j] - extend;
            if (xExtend > xOpen) { xCur[j] = xExtend; t |= 4; } else xCur[j] = xOpen;
            int yOpen = mCur[j - 1] - open, yExtend = yCur[j - 1] - extend;
            if (yExtend > yOpen) { yCur[j] = yExtend; t |= 8; } else yCur[j] = yOpen;
            trace[i * (m + 1) + j] = t;
            if (i == n || j == m) { Consider(mCur[j], i, j, 0); Consider(xCur[j], i, j, 1); Consider(yCur[j], i, j, 2); }
        }
        (mPrev, mCur) = (mCur, mPrev); (xPrev, xCur) = (xCur, xPrev); (yPrev, yCur) = (yCur, yPrev);
    }

    // Trace back: the aligned residue pairs and the gaps between the first and last of them.
    var pairs = new List<(int I, int J)>();
    int ci = bestI, cj = bestJ, state = bestState;
    var gapRuns = new List<int>();
    int run = 0;
    while (ci > 0 && cj > 0)
    {
        byte t = trace[ci * (m + 1) + cj];
        if (state == 0)
        {
            // Traced from the end, so a run met before any pair is a trailing end gap, not an inner indel.
            if (run > 0) { if (pairs.Count > 0) gapRuns.Add(run); run = 0; }
            pairs.Add((ci - 1, cj - 1));
            state = t & 3; ci--; cj--;
        }
        else if (state == 1) { run++; state = (t & 4) != 0 ? 1 : 0; ci--; }
        else { run++; state = (t & 8) != 0 ? 2 : 0; cj--; }
    }
    if (pairs.Count == 0) return Align.None;
    // A run still open when the trace stops is a leading end gap and isn't counted either.
    pairs.Reverse();
    int identical = pairs.Count(p => char.ToUpperInvariant(query[p.I]) == char.ToUpperInvariant(gene[p.J]));
    return new Align(identical / (double)pairs.Count, pairs.Count / (double)Math.Min(n, m), pairs[0].J - pairs[0].I, gapRuns.Count, gapRuns.Sum());
}

static void CheckBlosum()
{
    for (int a = 0; a < Blosum.Letters.Length; a++)
        for (int b = 0; b < Blosum.Letters.Length; b++)
            if (Blosum.Score(a, b) != Blosum.Score(b, a)) throw new InvalidOperationException($"BLOSUM62 isn't symmetric at {Blosum.Letters[a]}{Blosum.Letters[b]}.");
    if (string.Concat("ARNDCQEGHILKMFPSTWYV".Select(c => Blosum.Score(Blosum.Index(c), Blosum.Index(c)).ToString(CultureInfo.InvariantCulture) + ",")) != "4,5,6,6,9,5,5,6,8,4,4,5,5,6,7,4,5,11,7,4,")
        throw new InvalidOperationException("BLOSUM62's diagonal is wrong.");
}

record RecordInfo(long Source, string? Genus, string Provider, string? Accession, string RecordId, string? Order);
record Query(int Index, string Genome, int Band, string Protein, string Code, long SequenceId, RecordInfo Record, ZCurveGeneProfile.Gene Truth, double[] Axes, HashSet<int> Buckets);
record Hit(long SequenceId, ZCurveGeneProfile.Gene Gene, double Score, double[] Windows);
record Scan(double AtZero, int BestOffset, double Best);
record Align(double Identity, double AlignedShare, int? Offset, int InnerIndels, int InnerIndelResidues)
{
    public static readonly Align None = new(double.NaN, double.NaN, null, 0, 0);
}

static class ExpectedAxesCache
{
    private static readonly ConcurrentDictionary<string, (double X, double Y, double Z)[]> tables = new();
    public static (double X, double Y, double Z)[] For(string code) => tables.GetOrAdd(code, c => ZCurveGeneProfile.ExpectedAxes(null, c));
}

static class Blosum
{
    public const string Letters = "ARNDCQEGHILKMFPSTWYVBZX*";
    private static readonly int[] matrix = Parse(
        " 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n" +
        "-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n" +
        "-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n" +
        "-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n" +
        " 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n" +
        "-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n" +
        "-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n" +
        " 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n" +
        "-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n" +
        "-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n" +
        "-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n" +
        "-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n" +
        "-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n" +
        "-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n" +
        "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n" +
        " 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n" +
        " 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n" +
        "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n" +
        "-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n" +
        " 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n" +
        "-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n" +
        "-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n" +
        " 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n" +
        "-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1");

    private static int[] Parse(string text) =>
        text.Split((char[])[' ', '\n'], StringSplitOptions.RemoveEmptyEntries).Select(s => int.Parse(s, CultureInfo.InvariantCulture)).ToArray();

    public static int Index(char c)
    {
        int i = Letters.IndexOf(char.ToUpperInvariant(c));
        return i >= 0 ? i : Letters.IndexOf('X');
    }

    public static int Score(int a, int b) => matrix[a * Letters.Length + b];
}

sealed class Tally(int queries, int weightings)
{
    public readonly int[,] Beaten = new int[queries, weightings], Ties = new int[queries, weightings];
    public readonly PriorityQueue<Hit, double>[,] Winners = Make(queries, weightings);

    private static PriorityQueue<Hit, double>[,] Make(int queries, int weightings)
    {
        var winners = new PriorityQueue<Hit, double>[queries, weightings];
        for (int q = 0; q < queries; q++) for (int w = 0; w < weightings; w++) winners[q, w] = new PriorityQueue<Hit, double>();
        return winners;
    }
}
