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

// Protein offset and feature check (2026-11-20), after protein_family_diagnostic.cs (G4): with the query's own genus removed, the right
// PGfam was first for 64.1% of 989 proteins; 139 misses had it at rank 2-10 and 150 at 11-200. Two questions about the 9 numbers per
// window (Legendre moments P0-P2 of x, y, z over codon positions 1 and 2):
//
// A. Which numbers matter? For every query with a same-family gene outside its genus, the best such gene ("true") and the top hit of
//    another family or none ("wrong", what a miss loses to), from their stored windows: per feature, the mean squared difference to the
//    query for true and wrong, and how often true is the closer of the two on that feature alone. Today's score (unweighted) and the
//    1/variance weights G4 printed are both used to pick true and wrong.
// B. Are the misses offset copies? For each miss under today's score whose true and wrong genes are BV-BRC genes, both proteins and the
//    query's are fetched from BV-BRC (never stored; proteins, not DNA). The query and each gene are aligned globally with free end gaps
//    (BLOSUM62, gap 12 + 1 per residue): identity, the N-terminal offset (target residue facing the query's first aligned one) and inner
//    indels. Then both are profiled from the protein at every offset d of -60..60 aa (d > 0 drops the gene's first d residues, d < 0 the
//    query's), compared over their common length with the query's window rule, under both weightings. A miss "flips" when the true gene's
//    best-offset score beats the wrong gene's best-offset score.
//    Rule fixed before the run (unweighted): flips in >= 30% of the misses with both proteins, and in most of those the best offset within
//    10 aa of the alignment's offset -> an offset search is worth a ranking test (G6); under 15% -> offsets don't explain the misses and
//    G6's offset arm is dropped; in between, reported as unclear. Flipping against the top wrong hit only is necessary, not sufficient,
//    for a rank-1 gain: the ranking test is the real check.
// Also written: for up to 6 offset misses, window-by-window similarity at a fine step (60 aa windows every 10 aa), to see shifted stretches.
// Queries: protein_family_diagnostic.cs's rule and seed exactly (48 bacteria, 8 archaea, 8 plasmids, seed 20261105).
// Args: output prefix (default ProteinOffset), writing <prefix>_Features.csv, <prefix>_Pairs.csv, <prefix>_Dotplots.csv, <prefix>_Summary.csv.
const int seed = 20261105, perBand = 5, maxOffset = 60, maxAlignResidues = 4000, dotWindow = 60, dotStep = 10, dotPairs = 6;
string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
var genomesPerKind = new Dictionary<string, int> { ["bacterium"] = 48, ["archaeon"] = 8, ["plasmid"] = 8 };
string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
const int F = ZCurveGeneProfile.FeatureCount;
string[] featureNames = ["P0 x", "P0 y", "P0 z", "P1 x", "P1 y", "P1 z", "P2 x", "P2 y", "P2 z"];
// G4's 1/variance weights (catalogue-wide feature variance from a 2% row sample, normalized to sum 9).
var weightings = new (string Name, double[] Weights)[]
{
    ("unweighted", Enumerable.Repeat(1.0, F).ToArray()),
    ("1/variance", [0.741, 0.838, 0.333, 1.17, 1.16, 1.15, 1.19, 1.21, 1.2]),
};
int W = weightings.Length;
var prefix = args.Length > 0 ? args[0] : "ProteinOffset";
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

// The script rebuilds windows itself (any window size, any offset); it must equal the library's for the whole protein.
foreach (var query in queries)
{
    var mine = Profile(query.Protein, query.Code, query.Protein.Length >= ZCurveGeneProfile.LongResidues);
    if (mine.Length != query.Axes.Length || mine.Zip(query.Axes).Any(p => Math.Abs(p.First - p.Second) > 1e-9))
        throw new InvalidOperationException($"Profile rebuild differs from QueryFromProtein for a {query.Protein.Length} aa query.");
}

bool Comparable(Query query, ZCurveGeneProfile.Gene gene) => query.Protein.Length >= ZCurveGeneProfile.LongResidues
    ? gene.Residues >= ZCurveGeneProfile.LongResidues
    : gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - query.Protein.Length) <= ProteinSearch.LengthBand * query.Protein.Length;
bool OutsideGenus(Query query, RecordInfo candidate) => candidate.Source != query.Record.Source && (candidate.Genus is null || candidate.Genus != query.Record.Genus);
var byBucket = Enumerable.Range(0, 64).Select(b => queries.Where(q => q.Buckets.Contains(b)).ToArray()).ToArray();
var needed = queries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();

// One pass over the needed buckets: per query and weighting, the best same-family gene and the top other gene outside the genus.
var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
{
    var best = new Hit?[queries.Count, W, 2];   // [query, weighting, 0 = same family, 1 = other]
    var squared = new double[64 * F];
    foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable())
    {
        var candidate = records[id];
        var profile = ZCurveGeneProfile.FromBinary(payload);
        for (int g = 0; g < profile.Genes.Count; g++)
        {
            var gene = profile.Genes[g];
            double[]? windows = null;
            foreach (var query in byBucket[bucket])
            {
                if (!Comparable(query, gene) || !OutsideGenus(query, candidate)) continue;
                windows ??= profile.WindowsOf(g);
                if (squared.Length < windows.Length) squared = new double[windows.Length];
                int count = Differences(query.Axes, windows, squared);
                if (count == 0) continue;
                int side = gene.Family == query.Truth.Family ? 0 : 1;
                for (int v = 0; v < W; v++)
                {
                    double score = Score(squared, count, weightings[v].Weights);
                    if (best[query.Index, v, side] is not { } current || score > current.Score)
                        best[query.Index, v, side] = new Hit(id, gene, score, windows);
                }
            }
        }
    }
    return best;
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
var bestHits = new Hit?[queries.Count, W, 2];
for (int q = 0; q < queries.Count; q++)
    for (int v = 0; v < W; v++)
        for (int side = 0; side < 2; side++)
            bestHits[q, v, side] = perWorker.Select(w => w[q, v, side]).Where(h => h is not null).MaxBy(h => h!.Score);
Console.WriteLine($"Pass done, {total.Elapsed.TotalMinutes:F1} min.");

string F3(double x) => double.IsNaN(x) ? "" : x.ToString("F4", CultureInfo.InvariantCulture);
string Pct(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";

// A. Which features separate the right family from what beats it.
var featureLines = new List<string> { "Picked by,Group,Queries,Feature,Mean d2 true,Mean d2 wrong,Wrong/true,True closer" };
for (int v = 0; v < W; v++)
{
    var reachable = queries.Where(q => bestHits[q.Index, v, 0] is not null && bestHits[q.Index, v, 1] is not null).ToList();
    foreach (var (group, members) in new[] { ("all", reachable),
        ("hits", reachable.Where(q => bestHits[q.Index, v, 0]!.Score > bestHits[q.Index, v, 1]!.Score).ToList()),
        ("misses", reachable.Where(q => bestHits[q.Index, v, 0]!.Score <= bestHits[q.Index, v, 1]!.Score).ToList()) })
    {
        if (members.Count == 0) continue;
        var trueD = members.Select(q => PerFeature(q.Axes, bestHits[q.Index, v, 0]!.Windows)).ToList();
        var wrongD = members.Select(q => PerFeature(q.Axes, bestHits[q.Index, v, 1]!.Windows)).ToList();
        void Line(string feature, int[] indices)
        {
            double t = trueD.Average(d => indices.Sum(f => d[f])), w = wrongD.Average(d => indices.Sum(f => d[f]));
            int closer = Enumerable.Range(0, members.Count).Count(i => indices.Sum(f => trueD[i][f]) < indices.Sum(f => wrongD[i][f]));
            featureLines.Add(string.Join(",", weightings[v].Name, group, members.Count, feature, F3(t), F3(w), F3(w / t), Pct(closer, members.Count)));
        }
        for (int f = 0; f < F; f++) Line(featureNames[f], [f]);
        Line("P0 (means)", [0, 1, 2]); Line("P1 (trend)", [3, 4, 5]); Line("P2 (curvature)", [6, 7, 8]);
        Line("x (purine)", [0, 3, 6]); Line("y (amino)", [1, 4, 7]); Line("z (strong)", [2, 5, 8]);
        Line("all 9", Enumerable.Range(0, F).ToArray());
    }
}
File.WriteAllLines($"{prefix}_Features.csv", featureLines);
Console.WriteLine($"Features written, {total.Elapsed.TotalMinutes:F1} min.");

// B. The misses under today's score: both genes' proteins from BV-BRC, alignment, offset scan.
var misses = queries.Where(q => bestHits[q.Index, 0, 0] is { } t && bestHits[q.Index, 0, 1] is { } w && t.Score <= w.Score).ToList();
Console.WriteLine($"{misses.Count} misses under today's score (G4: 355).");
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

var pairLines = new List<string> { "Query genome,Band,Query aa,Family,Pair,Gene record,Gene aa,Gene family,Identity,Aligned share,Alignment offset,Inner indels," +
    "Inner indel aa,Stored score,Overlap score at 0 unweighted,Best offset unweighted,Best score unweighted,Overlap score at 0 weighted,Best offset weighted,Best score weighted" };
var offsetFeatureShift = new double[2, F];   // true pairs of usable misses: mean d2 per feature at offset 0 and at the best offset (unweighted)
var outcomes = new List<(Query Query, Align TrueAlign, Scan[] TrueScan, Scan[] WrongScan)>();
foreach (var query in misses)
{
    var trueHit = bestHits[query.Index, 0, 0]!;
    var wrongHit = bestHits[query.Index, 0, 1]!;
    var trueProtein = await GeneProteinAsync(trueHit);
    var wrongProtein = await GeneProteinAsync(wrongHit);
    if (trueProtein is null || wrongProtein is null) continue;
    var pairs = new[] { ("true", trueHit, trueProtein), ("wrong", wrongHit, wrongProtein) };
    var aligns = new Align[2];
    var scans = new Scan[2][];
    for (int p = 0; p < 2; p++)
    {
        var (label, hit, protein) = pairs[p];
        aligns[p] = query.Protein.Length <= maxAlignResidues && protein.Length <= maxAlignResidues ? Alignment(query.Protein, protein) : Align.None;
        scans[p] = weightings.Select(w => OffsetScan(query, protein, w.Weights)).ToArray();
        pairLines.Add(string.Join(",", query.Genome, bandNames[query.Band], query.Protein.Length, query.Truth.Family, label, hit.SequenceId, protein.Length,
            hit.Gene.Family, F3(aligns[p].Identity), F3(aligns[p].AlignedShare), aligns[p].Offset?.ToString(CultureInfo.InvariantCulture) ?? "",
            aligns[p].InnerIndels, aligns[p].InnerIndelResidues, F3(hit.Score),
            F3(scans[p][0].AtZero), scans[p][0].BestOffset, F3(scans[p][0].Best), F3(scans[p][1].AtZero), scans[p][1].BestOffset, F3(scans[p][1].Best)));
    }
    var atZero = PerFeature(Profile(Overlap(query.Protein, trueProtein, 0).Query, query.Code, query.Protein.Length >= ZCurveGeneProfile.LongResidues),
        Profile(Overlap(query.Protein, trueProtein, 0).Gene, query.Code, query.Protein.Length >= ZCurveGeneProfile.LongResidues));
    var (bq, bg) = Overlap(query.Protein, trueProtein, scans[0][0].BestOffset);
    var atBest = PerFeature(Profile(bq, query.Code, query.Protein.Length >= ZCurveGeneProfile.LongResidues), Profile(bg, query.Code, query.Protein.Length >= ZCurveGeneProfile.LongResidues));
    for (int f = 0; f < F; f++) { offsetFeatureShift[0, f] += atZero[f]; offsetFeatureShift[1, f] += atBest[f]; }
    outcomes.Add((query, aligns[0], scans[0], scans[1]));
    if (outcomes.Count % 25 == 0) Console.WriteLine($"  {outcomes.Count} misses with both proteins, {total.Elapsed.TotalMinutes:F1} min");
}
File.WriteAllLines($"{prefix}_Pairs.csv", pairLines);

// Dot plots: the clearest offset misses (alignment offset >= 10 aa, identity >= 30%), window by window at a fine step.
var dotLines = new List<string> { "Query genome,Query aa,Gene aa,Alignment offset,Query window start aa,Gene window start aa,Similarity" };
foreach (var (query, align, _, _) in outcomes.Where(o => o.TrueAlign.Offset is { } d && Math.Abs(d) >= 10 && o.TrueAlign.Identity >= 0.3)
    .OrderByDescending(o => o.TrueAlign.Identity).Take(dotPairs))
{
    var gene = geneProteins[(bestHits[query.Index, 0, 0]!.SequenceId, bestHits[query.Index, 0, 0]!.Gene.Start, bestHits[query.Index, 0, 0]!.Gene.End, bestHits[query.Index, 0, 0]!.Gene.Reverse)]!;
    var qWindows = FineWindows(query.Protein, query.Code);
    var gWindows = FineWindows(gene, query.Code);
    for (int i = 0; i < qWindows.Count; i++)
        for (int j = 0; j < gWindows.Count; j++)
        {
            double d = Math.Sqrt(Enumerable.Range(0, F).Sum(f => Math.Pow(qWindows[i][f] - gWindows[j][f], 2)));
            dotLines.Add(string.Join(",", query.Genome, query.Protein.Length, gene.Length, align.Offset, i * dotStep, j * dotStep, F3(1 - d)));
        }
}
File.WriteAllLines($"{prefix}_Dotplots.csv", dotLines);

// Summary and the rule.
var summary = new List<string> { "Item,Value" };
summary.Add($"Queries,{queries.Count}");
summary.Add($"Misses under today's score (own genus removed),{misses.Count}");
summary.Add($"Misses with both proteins from BV-BRC,{outcomes.Count}");
for (int v = 0; v < W; v++)
{
    int flips = outcomes.Count(o => o.TrueScan[v].Best > o.WrongScan[v].Best);
    int atZeroAlready = outcomes.Count(o => o.TrueScan[v].AtZero > o.WrongScan[v].AtZero);
    summary.Add($"{weightings[v].Name}: true beats wrong over the common length at offset 0,{Pct(atZeroAlready, outcomes.Count)} ({atZeroAlready})");
    summary.Add($"{weightings[v].Name}: true beats wrong at their best offsets (flip),{Pct(flips, outcomes.Count)} ({flips})");
}
var flipped = outcomes.Where(o => o.TrueScan[0].Best > o.WrongScan[0].Best).ToList();
var withOffset = flipped.Where(o => o.TrueAlign.Offset is not null).ToList();
int agree = withOffset.Count(o => Math.Abs(o.TrueScan[0].BestOffset - o.TrueAlign.Offset!.Value) <= 10);
summary.Add($"Flips whose best offset is within 10 aa of the alignment's,{Pct(agree, withOffset.Count)} ({agree} of {withOffset.Count})");
var identities = outcomes.Where(o => !double.IsNaN(o.TrueAlign.Identity)).Select(o => o.TrueAlign.Identity).Order().ToList();
if (identities.Count > 0)
    summary.Add($"True pairs' identity: quartiles,{F3(identities[identities.Count / 4])} / {F3(identities[identities.Count / 2])} / {F3(identities[3 * identities.Count / 4])}");
summary.Add($"True pairs with alignment offset >= 10 aa,{Pct(outcomes.Count(o => o.TrueAlign.Offset is { } d && Math.Abs(d) >= 10), outcomes.Count)}");
summary.Add($"True pairs with inner indels,{Pct(outcomes.Count(o => o.TrueAlign.InnerIndels > 0), outcomes.Count)}");
int used = Math.Max(outcomes.Count, 1);
foreach (var (name, indices) in new[] { ("P0", new[] { 0, 1, 2 }), ("P1", new[] { 3, 4, 5 }), ("P2", new[] { 6, 7, 8 }) })
    summary.Add($"True pairs' mean d2 of {name}: offset 0 -> best offset,{F3(indices.Sum(f => offsetFeatureShift[0, f]) / used)} -> {F3(indices.Sum(f => offsetFeatureShift[1, f]) / used)}");
double flipShare = outcomes.Count == 0 ? 0 : flipped.Count / (double)outcomes.Count;
string verdict = flipShare >= 0.30 && withOffset.Count > 0 && agree * 2 > withOffset.Count ? "go: an offset search is worth a ranking test (G6)"
    : flipShare < 0.15 ? "stop: offsets don't explain the misses; G6's offset arm is dropped" : "unclear: report, no ranking test yet";
summary.Add($"Rule (unweighted: flips >= 30% and most within 10 aa -> ranking test; < 15% -> drop),{verdict}");
summary.Add($"Total,{total.Elapsed.TotalMinutes:F1} min");
File.WriteAllLines($"{prefix}_Summary.csv", summary);
Console.WriteLine(string.Join(Environment.NewLine, summary));
Console.WriteLine(string.Join(Environment.NewLine, featureLines.Where(l => l.Contains(",all,") && (l.Contains("(") || l.Contains("all 9")))));

// ---- helpers ----

// The query and gene from the offset on, cut to their common length: d > 0 drops the gene's first d residues, d < 0 the query's.
static (string Query, string Gene) Overlap(string query, string gene, int d)
{
    var q = d < 0 ? query[Math.Min(-d, query.Length)..] : query;
    var g = d > 0 ? gene[Math.Min(d, gene.Length)..] : gene;
    int length = Math.Min(q.Length, g.Length);
    return (q[..length], g[..length]);
}

// The query's window rule kept whatever the overlap's length: 333-base windows for a long query, 5 equal windows otherwise.
Scan OffsetScan(Query query, string gene, double[] weights)
{
    bool longRule = query.Protein.Length >= ZCurveGeneProfile.LongResidues;
    double atZero = double.NaN, best = double.NegativeInfinity;
    int bestOffset = 0;
    var squared = new double[64 * F];
    for (int d = -maxOffset; d <= maxOffset; d++)
    {
        var (q, g) = Overlap(query.Protein, gene, d);
        if (q.Length < ZCurveGeneProfile.MinResidues) continue;
        var a = Profile(q, query.Code, longRule);
        var b = Profile(g, query.Code, longRule);
        if (squared.Length < a.Length) squared = new double[a.Length];
        int count = Differences(a, b, squared);
        if (count == 0) continue;
        double score = Score(squared, count, weights);
        if (d == 0) atZero = score;
        if (score > best || (score == best && Math.Abs(d) < Math.Abs(bestOffset))) (best, bestOffset) = (score, d);
    }
    return new Scan(atZero, bestOffset, best);
}

// Windows rebuilt from a protein exactly as ZCurveGeneProfile does (checked against QueryFromProtein above), with the window rule chosen.
static double[] Profile(string protein, string code, bool longRule)
{
    int window = longRule ? ZCurveGeneProfile.WindowBases : protein.Length >= ZCurveGeneProfile.MinResidues ? 2 * protein.Length / ZCurveGeneProfile.SplitWindows : 0;
    return Moments(protein, code, 0, 2 * protein.Length, window);
}

static double[] Moments(string protein, string code, int startBase, int bases, int window)
{
    var table = ExpectedAxesCache.For(code);
    int count = window == 0 ? 0 : bases / window;
    var result = new double[count * F];
    for (int w = 0; w < count; w++)
    {
        int slot = w * F;
        for (int i = 0; i < window; i++)
        {
            int b = startBase + w * window + i;
            int letter = char.ToUpperInvariant(protein[b / 2]) - 'A';
            var (x, y, z) = letter is >= 0 and < 26 ? table[2 * letter + b % 2] : (0, 0, 0);
            double t = -1 + 2 * (i + 0.5) / window, p2 = (3 * t * t - 1) / 2;
            result[slot] += x; result[slot + 1] += y; result[slot + 2] += z;
            result[slot + 3] += x * t; result[slot + 4] += y * t; result[slot + 5] += z * t;
            result[slot + 6] += x * p2; result[slot + 7] += y * p2; result[slot + 8] += z * p2;
        }
        for (int f = 0; f < F; f++) result[slot + f] *= Math.Sqrt(2 * (f / 3) + 1) / window;
    }
    return result;
}

// 60 aa windows every 10 aa, for the dot plots.
List<double[]> FineWindows(string protein, string code)
{
    var list = new List<double[]>();
    for (int start = 0; start + dotWindow <= protein.Length; start += dotStep)
        list.Add(Moments(protein, code, 2 * start, 2 * dotWindow, 2 * dotWindow));
    return list;
}

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

// Mean squared difference per feature over the common windows.
static double[] PerFeature(double[] query, double[] gene)
{
    var result = new double[F];
    int count = Math.Min(query.Length, gene.Length) / F;
    for (int w = 0; w < count; w++)
        for (int f = 0; f < F; f++) result[f] += Math.Pow(query[w * F + f] - gene[w * F + f], 2) / count;
    return result;
}

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
