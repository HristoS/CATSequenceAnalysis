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

// Symmetric homolog key (2026-11-20). Protein family numbers so far use BV-BRC's PGfam as the key, a function label built from
// similarity, and G8's alignments showed it cuts both ways: 23% of the hits that "beat" the right family were real homologs, and 16% of
// same-PGfam pairs were not detectably related. Those were misses only. Here every query is judged the same way, hits and misses.
// Queries: G4's 1,143 proteins (protein_family_diagnostic.cs's rule and seed). Score: ZCurveGeneProfile.Score as adopted 2026-11-20
// (FeatureWeights). Over the genes outside the query's genus, as ProteinSearch compares them: the top hit, and the best gene of the
// query's PGfam. Both proteins are fetched from BV-BRC (never stored) and aligned with the query (global, free end gaps, BLOSUM62); a
// gene is a homolog at identity >= Rost's HSSP curve + 5 points over >= 60% of the shorter protein.
// Reported, per band and overall, with 95% genome-bootstrap CIs: PGfam top-1 (the top hit's PGfam is the query's) and homolog top-1
// (the top hit is a homolog), over queries with a same-PGfam gene outside the genus; the identity curve (homolog top-1 by the identity of
// the closest homolog found, the same-PGfam gene or the top hit); and the key's quality on hits and misses alike (same-PGfam genes that
// are not detectably homologous).
// Rule fixed before the run: both keys are reported whatever they show. "Misses are partly label noise" is claimed only if the share of
// undetectable same-PGfam genes among PGfam misses exceeds that among PGfam hits by >= 10 points with the CI above 0.
// Args: output prefix (default ProteinHomologKey), writing <prefix>_Queries.csv and <prefix>_Summary.csv.
const int seed = 20261105, perBand = 5, maxAlignResidues = 6000, bootstrapRounds = 2000, batchSize = 20;
string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 48, ["archaeon"] = 8, ["plasmid"] = 8 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
const int F = ZCurveGeneProfile.FeatureCount;
var prefix = args.Length > 0 ? args[0] : "ProteinHomologKey";
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

// ---- the symmetric key ----
bool Comparable(Query query, ZCurveGeneProfile.Gene gene) => query.Protein.Length >= ZCurveGeneProfile.LongResidues
    ? gene.Residues >= ZCurveGeneProfile.LongResidues
    : gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - query.Protein.Length) <= ProteinSearch.LengthBand * query.Protein.Length;
bool OutsideGenus(Query query, RecordInfo candidate) => candidate.Source != query.Record.Source && (candidate.Genus is null || candidate.Genus != query.Record.Genus);
var byBucket = Enumerable.Range(0, 64).Select(b => queries.Where(q => q.Buckets.Contains(b)).ToArray()).ToArray();
var needed = queries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();

var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
{
    var top = new Hit?[queries.Count];
    var same = new Hit?[queries.Count];
    foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable())
    {
        var candidate = records[id];
        var profile = ZCurveGeneProfile.FromBinary(payload);
        for (int g = 0; g < profile.Genes.Count; g++)
        {
            var gene = profile.Genes[g];
            foreach (var query in byBucket[bucket])
            {
                if (!Comparable(query, gene) || !OutsideGenus(query, candidate)) continue;
                double score = profile.Score(query.Axes, g);
                if (double.IsNaN(score)) continue;
                if (top[query.Index] is not { } t || score > t.Score) top[query.Index] = new Hit(id, gene, score, []);
                if (gene.Family == query.Truth.Family && (same[query.Index] is not { } s || score > s.Score)) same[query.Index] = new Hit(id, gene, score, []);
            }
        }
    }
    return (top, same);
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
var topOf = queries.Select(q => perWorker.Select(w => w.top[q.Index]).Where(h => h is not null).MaxBy(h => h!.Score)).ToArray();
var sameOf = queries.Select(q => perWorker.Select(w => w.same[q.Index]).Where(h => h is not null).MaxBy(h => h!.Score)).ToArray();
Console.WriteLine($"Pass done, {total.Elapsed.TotalMinutes:F1} min.");

// Proteins of the genes to judge, fetched in batches: one genome_feature request for up to 20 genes (an or() of their coordinates),
// then their sequences by md5, 50 at a time. A failed batch falls back to one request per gene.
var geneProteins = new Dictionary<(long, int, int, bool), string?>();
var toFetch = topOf.Concat(sameOf).Where(h => h is not null).Select(h => h!)
    .Where(h => records[h.SequenceId].Provider == "BV-BRC" && records[h.SequenceId].Accession is not null && h.Gene.Chromosome == 0)
    .DistinctBy(h => (h.SequenceId, h.Gene.Start, h.Gene.End, h.Gene.Reverse)).ToList();
string Clause(Hit h) => $"and(eq(genome_id,{records[h.SequenceId].Accession}),eq(accession,{Uri.EscapeDataString(records[h.SequenceId].RecordId)})," +
    $"eq(start,{h.Gene.Start + 1}),eq(end,{h.Gene.End}))";
var md5Of = new Dictionary<(long, int, int, bool), string>();
async Task FeaturesAsync(List<Hit> batch)
{
    var byCoordinates = batch.GroupBy(h => (records[h.SequenceId].RecordId, h.Gene.Start + 1, h.Gene.End)).ToDictionary(g => g.Key, g => g.ToList());
    var found = await GetAsync("https://www.bv-brc.org/api/genome_feature/?and(or(" + string.Join(",", batch.Select(Clause)) +
        "),eq(feature_type,CDS),eq(annotation,PATRIC))&select(accession,start,end,strand,aa_sequence_md5)&limit(200)");
    foreach (var f in found.EnumerateArray())
    {
        if (!f.TryGetProperty("aa_sequence_md5", out var md5)) continue;
        var key = (f.GetProperty("accession").GetString()!, f.GetProperty("start").GetInt32(), f.GetProperty("end").GetInt32());
        if (!byCoordinates.TryGetValue(key, out var hits)) continue;
        bool reverse = f.GetProperty("strand").GetString() == "-";
        foreach (var h in hits.Where(h => h.Gene.Reverse == reverse)) md5Of[(h.SequenceId, h.Gene.Start, h.Gene.End, h.Gene.Reverse)] = md5.GetString()!;
    }
}
foreach (var batch in toFetch.Chunk(batchSize))
{
    try { await FeaturesAsync(batch.ToList()); }
    catch (HttpRequestException)
    {
        foreach (var h in batch)
            try { await FeaturesAsync([h]); } catch (HttpRequestException e) { Console.WriteLine($"  record {h.SequenceId} at {h.Gene.Start}: {e.Message}"); }
    }
    if (md5Of.Count % 200 < batchSize) Console.WriteLine($"  {md5Of.Count} of {toFetch.Count} genes located, {total.Elapsed.TotalMinutes:F1} min");
}
var sequences = new Dictionary<string, string>();
foreach (var chunk in md5Of.Values.Distinct().Chunk(50))
    foreach (var p in (await GetAsync($"https://www.bv-brc.org/api/feature_sequence/?in(md5,({string.Join(",", chunk)}))&select(md5,sequence)&limit(1000)")).EnumerateArray())
        sequences[p.GetProperty("md5").GetString()!] = p.GetProperty("sequence").GetString()!.TrimEnd('*');
string? ProteinOf(Hit h) => md5Of.TryGetValue((h.SequenceId, h.Gene.Start, h.Gene.End, h.Gene.Reverse), out var md5) && sequences.TryGetValue(md5, out var s) ? s : null;
Console.WriteLine($"{sequences.Count} proteins fetched, {total.Elapsed.TotalMinutes:F1} min.");

// Judge each query: is the top hit a homolog, is the same-PGfam gene one.
(bool? Homolog, double Identity) Judge(Query query, Hit? hit)
{
    if (hit is null || ProteinOf(hit) is not { } protein || query.Protein.Length > maxAlignResidues || protein.Length > maxAlignResidues) return (null, double.NaN);
    var align = Alignment(query.Protein, protein);
    int shorter = Math.Min(query.Protein.Length, protein.Length);
    double pairs = double.IsNaN(align.AlignedShare) ? 0 : align.AlignedShare * shorter;
    return (pairs >= 0.6 * shorter && 100 * align.Identity >= Hssp((int)pairs) + 5, align.Identity);
}
var outcomes = new List<Outcome>();
var queryLines = new List<string> { "Query genome,Band,Query aa,Query family,Top family,Top score,PGfam top-1,Top homolog,Top identity,Same-PGfam score,Same-PGfam homolog,Same-PGfam identity" };
string F4(double x) => double.IsNaN(x) ? "" : x.ToString("F4", CultureInfo.InvariantCulture);
string B(bool? b) => b is null ? "" : b.Value ? "1" : "0";
foreach (var query in queries)
{
    var top = topOf[query.Index];
    var sameHit = sameOf[query.Index];
    if (top is null || sameHit is null) continue;   // not reachable: no gene of its PGfam outside the genus
    bool pgfam = top.Gene.Family == query.Truth.Family;
    var (topHomolog, topIdentity) = Judge(query, top);
    var (sameHomolog, sameIdentity) = pgfam && top == sameHit ? (topHomolog, topIdentity) : Judge(query, sameHit);
    outcomes.Add(new Outcome(query, pgfam, topHomolog, topIdentity, sameHomolog, sameIdentity));
    queryLines.Add(string.Join(",", query.Genome, bandNames[query.Band], query.Protein.Length, query.Truth.Family, top.Gene.Family, F4(top.Score), pgfam ? 1 : 0,
        B(topHomolog), F4(topIdentity), F4(sameHit.Score), B(sameHomolog), F4(sameIdentity)));
}
File.WriteAllLines($"{prefix}_Queries.csv", queryLines);

// Shares with genome-bootstrap CIs.
string Pct(double x) => double.IsNaN(x) ? "" : (100 * x).ToString("F1", CultureInfo.InvariantCulture) + "%";
(double Value, double Low, double High) Boot(List<Outcome> set, Func<List<Outcome>, double> statistic)
{
    var genomes = set.GroupBy(o => o.Query.Genome).Select(g => g.ToList()).ToList();
    var random = new Random(seed);
    var values = new List<double>();
    for (int r = 0; r < bootstrapRounds; r++)
    {
        var sample = Enumerable.Range(0, genomes.Count).SelectMany(_ => genomes[random.Next(genomes.Count)]).ToList();
        double v = statistic(sample);
        if (!double.IsNaN(v)) values.Add(v);
    }
    values.Sort();
    return (statistic(set), values.Count == 0 ? double.NaN : values[(int)(0.025 * values.Count)], values.Count == 0 ? double.NaN : values[(int)(0.975 * values.Count) - 1]);
}
string Show((double Value, double Low, double High) b) => $"{Pct(b.Value)} ({Pct(b.Low)} to {Pct(b.High)})";
double Share(List<Outcome> s, Func<Outcome, bool?> pick) { var known = s.Select(pick).Where(v => v is not null).ToList(); return known.Count == 0 ? double.NaN : known.Count(v => v!.Value) / (double)known.Count; }

var summary = new List<string> { "Band,Queries,PGfam top-1 (95% CI),Homolog top-1 (95% CI),Top hit judged,Same-PGfam not detectably homologous" };
foreach (int band in Enumerable.Range(0, bandNames.Length).Append(-1))
{
    var cell = outcomes.Where(o => band < 0 || o.Query.Band == band).ToList();
    if (cell.Count == 0) continue;
    summary.Add(string.Join(",", band < 0 ? "all" : bandNames[band], cell.Count, Show(Boot(cell, s => s.Count(o => o.PgfamTop1) / (double)s.Count)),
        Show(Boot(cell, s => Share(s, o => o.TopHomolog))), cell.Count(o => o.TopHomolog is not null),
        Pct(Share(cell, o => o.SameHomolog is null ? null : !o.SameHomolog.Value))));
}
summary.Add("");
summary.Add("Closest homolog's identity,Queries,Homolog top-1");
foreach (var (low, high) in new[] { (0.0, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 0.95), (0.95, 1.01) })
{
    var bin = outcomes.Where(o => o.TopHomolog is not null)
        .Select(o => (o, Closest: new[] { o.TopHomolog == true ? o.TopIdentity : double.NaN, o.SameHomolog == true ? o.SameIdentity : double.NaN }.Where(x => !double.IsNaN(x)).DefaultIfEmpty(double.NaN).Max()))
        .Where(x => x.Closest >= low && x.Closest < high).ToList();
    if (bin.Count > 0) summary.Add($"{low:F2}-{Math.Min(high, 1):F2},{bin.Count},{Pct(bin.Count(x => x.o.TopHomolog == true) / (double)bin.Count)}");
}
var hits = outcomes.Where(o => o.PgfamTop1 && o.SameHomolog is not null).ToList();
var misses = outcomes.Where(o => !o.PgfamTop1 && o.SameHomolog is not null).ToList();
double Doubtful(List<Outcome> s) => s.Count == 0 ? double.NaN : s.Count(o => o.SameHomolog == false) / (double)s.Count;
var difference = Boot(outcomes.Where(o => o.SameHomolog is not null).ToList(),
    s => Doubtful(s.Where(o => !o.PgfamTop1).ToList()) - Doubtful(s.Where(o => o.PgfamTop1).ToList()));
summary.Add("");
summary.Add($"Same-PGfam gene not detectably homologous: PGfam hits,{Pct(Doubtful(hits))} of {hits.Count}");
summary.Add($"Same-PGfam gene not detectably homologous: PGfam misses,{Pct(Doubtful(misses))} of {misses.Count}");
summary.Add($"PGfam misses whose top hit is a homolog,{Pct(Share(misses, o => o.TopHomolog))}");
summary.Add($"Difference misses - hits,{(100 * difference.Value):+0.0;-0.0} points (95% CI {(100 * difference.Low):+0.0;-0.0} to {(100 * difference.High):+0.0;-0.0})");
summary.Add("Rule (>= 10 points, CI above 0)," + (difference.Value >= 0.10 && difference.Low > 0 ? "claim: misses are partly label noise" : "no claim; report both keys"));
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

record Outcome(Query Query, bool PgfamTop1, bool? TopHomolog, double TopIdentity, bool? SameHomolog, double SameIdentity);
