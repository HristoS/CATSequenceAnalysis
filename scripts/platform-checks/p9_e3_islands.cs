#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Formats.Tar;
using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using System.Text.RegularExpressions;
using SequenceAnalyses;
using SequencePlatform.Core;

// P9 S9.3 (E3): foreign islands from the stored 500 bp windows, the thesis method's own question. A detector scores every window by
// how far its local composition sits from the genome's own; islands from IslandViewer 4 are the answer key (D-P9c, 2026-11-04):
// - Main key: IslandPick (comparative genomics, so independent of composition), matched to the corpus's chromosomes by name (RefSeq
//   NZ_X to GenBank X), every region inside our length, same version where ours is known. Windows it doesn't flag count as negatives,
//   which understates precision for every detector alike.
// - Second key: the curated 104-genome reference (Sup_tables 1 positives, 2 negatives; windows in neither are left out), RefSeq NC_
//   mapped to GenBank through NCBI's nuccore summary (cached in <prefix>_RefSeqMap.csv) and kept when the lengths are equal.
// Fixed before the run:
// - Detector: per axis, the window value minus a baseline (detrended: a running median over 201 windows, about 100 kbp, so the
//   replication skew's flip at origin and terminus isn't read as an island; genomeMean: the genome's mean, the classic form, as a
//   check); the residuals averaged over 8 windows (4 kbp); score = Mahalanobis distance squared under the genome's own covariance
//   of those averages. zcurve uses x, y and z; gc uses z alone (z is the GC axis), so the difference is what x and y add; xy is context.
// - A window is positive when at least 250 bp of it lies in a region.
// - Genomes split in half (seeded); each detector's F1 threshold is chosen on the first half and tested on the second.
// - Pass: zcurve's test F1 minus gc's, detrended, has a bootstrap 95 % lower bound (over test genomes) above 0. Average precision
//   over all genomes is reported with its own bootstrap.
// Writes <prefix>_Islands.csv and <prefix>_RefSeqMap.csv. The reference is downloaded into --reference <dir> (default: the temp
// folder) and never committed: the site states no licence.
const int window = RecordStore.FineWindowSize, medianWidth = 201, medianStride = 20, smoothWidth = 8, minOverlapBp = 250;
const int bins = 2048, bootstraps = 1000, seed = 20261104;
const string sourceUrl = "https://www.pathogenomics.sfu.ca/islandviewer/download/datasets/";
int referenceAt = Array.IndexOf(args, "--reference");
var prefix = args.Where((a, i) => !a.StartsWith("--", StringComparison.Ordinal) && i != referenceAt + 1).FirstOrDefault() ?? "P9_E3";
var referenceDir = referenceAt >= 0 ? args[referenceAt + 1] : Path.Combine(Path.GetTempPath(), "p9_island_reference");
string[] detectors = ["zcurve", "gc", "xy"];
string[] variants = ["detrended", "genomeMean"];

await DownloadAsync(referenceDir, ["all_gis_islandpick_iv4.txt.tar.gz", "Sup_tables.tar.gz"]);
var islandPick = Regions(Path.Combine(referenceDir, "all_gis_islandpick_iv4.txt"), startColumn: 1);
var curatedPositive = Regions(Path.Combine(referenceDir, "SupTable1.txt"), startColumn: 3);
var curatedNegative = Regions(Path.Combine(referenceDir, "SupTable2.txt"), startColumn: 3);

await using var dataSource = ExperimentCorpus.OpenSnapshot();
var corpusIds = (await ExperimentCorpus.LoadAsync(dataSource)).Where(g => !g.Excluded && g.Group == SignatureIndex.Cellular).Select(g => g.Id).ToArray();
var ours = new Dictionary<string, Chromosome>(StringComparer.Ordinal);
await using (var command = dataSource.CreateCommand(
    "select id, record_id, length_bp, coalesce(r.tax_genus, ''), substring(annotations->'contig'->>'description' from '^[A-Z_]+[0-9]+\\.([0-9]+)$') " +
    "from sequence_record r where id = any($1)"))
{
    command.Parameters.AddWithValue(corpusIds);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync())
        ours[reader.GetString(1).Split('.')[0]] = new Chromosome(reader.GetInt64(0), reader.GetString(1), reader.GetInt64(2), reader.GetString(3),
            reader.IsDBNull(4) ? null : reader.GetString(4));
}

var pickGenomes = islandPick.Where(e => ours.TryGetValue(e.Key, out var c) && e.Value.All(r => r.End <= c.LengthBp)
        && (c.Version is null || c.Version == e.Value[0].Version))
    .Select(e => (Chromosome: ours[e.Key], Positive: e.Value, Negative: (List<Region>?)null)).OrderBy(g => g.Chromosome.Id).ToList();
var map = await RefSeqMapAsync($"{prefix}_RefSeqMap.csv", curatedPositive.Keys.Union(curatedNegative.Keys).Order().ToArray());
// Stored under the NC_ name itself, else under the GenBank record it copies; either way the length must be NCBI's.
Chromosome? CuratedMatch(string k) => map.TryGetValue(k, out var m)
    ? (ours.TryGetValue(k, out var c) || ours.TryGetValue(m.GenBank, out c)) && c.LengthBp == m.LengthBp ? c : null
    : null;
var curatedGenomes = curatedPositive.Keys.Union(curatedNegative.Keys)
    .Where(k => CuratedMatch(k) is not null)
    .Select(k => (Chromosome: CuratedMatch(k)!, Positive: curatedPositive.GetValueOrDefault(k) ?? [], Negative: (List<Region>?)(curatedNegative.GetValueOrDefault(k) ?? [])))
    .OrderBy(g => g.Chromosome.Id).ToList();
Console.WriteLine($"IslandPick: {pickGenomes.Count} chromosomes, {pickGenomes.Sum(g => g.Positive.Count)} regions. " +
    $"Curated: {curatedGenomes.Count} of {map.Count} mapped chromosomes.");

var ids = pickGenomes.Concat(curatedGenomes).Select(g => g.Chromosome.Id).Distinct().ToArray();
var profiles = await ProfileLevels.ReadManyAsync(dataSource, ids, RecordStore.Method, RecordStore.MethodVersion, window);
var pick = pickGenomes.Where(g => profiles.ContainsKey(g.Chromosome.Id)).Select(g => Label(g.Chromosome, profiles[g.Chromosome.Id], g.Positive, null)).ToList();
var curated = curatedGenomes.Where(g => profiles.ContainsKey(g.Chromosome.Id)).Select(g => Label(g.Chromosome, profiles[g.Chromosome.Id], g.Positive, g.Negative)).ToList();
Console.WriteLine($"Profiles: {pick.Count} IslandPick, {curated.Count} curated; {pick.Sum(g => g.Labels.Length):N0} windows, " +
    $"{pick.Sum(g => g.Labels.Count(l => l == 1)):N0} positive.");

var order = Enumerable.Range(0, pick.Count).ToArray();
new Random(seed).Shuffle(order);
var train = order.Take(pick.Count / 2).ToHashSet();
var test = Enumerable.Range(0, pick.Count).Where(i => !train.Contains(i)).ToArray();
var draws = Draws(pick.Count, bootstraps, seed);
var testDraws = Draws(test.Length, bootstraps, seed + 1);

var rows = new List<string>();
bool pass = false;
foreach (var variant in variants)
{
    var evaluated = new Dictionary<string, Evaluation>();
    foreach (var detector in detectors)
    {
        var scores = pick.Select(g => Scores(g.Profile, detector, variant == "detrended")).ToArray();
        var edges = Edges(scores, bins, seed);
        var histograms = pick.Select((g, i) => Histogram(scores[i], g.Labels, edges)).ToArray();
        var ap = histograms.Aggregate(new Histogram(new long[bins], new long[bins]), Add);
        int threshold = BestThreshold(train.Select(i => histograms[i]).Aggregate(new Histogram(new long[bins], new long[bins]), Add));
        var testCounts = test.Select(i => Counts(histograms[i], threshold)).ToArray();
        var rounds = draws.AsParallel().AsOrdered().Select(d => AveragePrecision(Weighted(histograms, d))).ToArray();
        evaluated[detector] = new Evaluation(histograms, testCounts, threshold, rounds);

        var apRounds = rounds.Order().ToArray();
        var (tp, fp, fn) = Sum(testCounts, null);
        var curatedResult = Curated(curated, detector, variant == "detrended", edges, threshold);
        rows.Add(string.Create(CultureInfo.InvariantCulture,
            $"islandPick, {detector}, {variant}, {pick.Count}, {ap.Positive.Sum() + ap.Negative.Sum()}, {ap.Positive.Sum()}, " +
            $"{(double)ap.Positive.Sum() / (ap.Positive.Sum() + ap.Negative.Sum()):F4}, {AveragePrecision(ap):F4}, {Low(apRounds):F4}, {High(apRounds):F4}, " +
            $"{test.Length}, {Precision(tp, fp):F4}, {Recall(tp, fn):F4}, {F1(tp, fp, fn):F4}, , , "));
        Console.WriteLine(rows[^1]);
        rows.Add(curatedResult);
        Console.WriteLine(rows[^1]);
    }

    foreach (var other in new[] { "gc", "xy" })
    {
        var a = evaluated["zcurve"];
        var b = evaluated[other];
        var apDifference = a.ApRounds.Zip(b.ApRounds, (x, y) => x - y).Order().ToArray();
        var f1Difference = testDraws.Select(d => F1Of(Sum(a.TestCounts, d)) - F1Of(Sum(b.TestCounts, d))).Order().ToArray();
        double apPoint = AveragePrecision(Weighted(a.Histograms, null)) - AveragePrecision(Weighted(b.Histograms, null));
        double f1Point = F1Of(Sum(a.TestCounts, null)) - F1Of(Sum(b.TestCounts, null));
        bool rowPass = other == "gc" && Low(f1Difference) > 0;
        if (variant == "detrended" && other == "gc") pass = rowPass;
        rows.Add(string.Create(CultureInfo.InvariantCulture,
            $"islandPick, zcurveMinus{char.ToUpperInvariant(other[0])}{other[1..]}, {variant}, {pick.Count}, , , , {apPoint:F4}, {Low(apDifference):F4}, {High(apDifference):F4}, " +
            $"{test.Length}, , , {f1Point:F4}, {Low(f1Difference):F4}, {High(f1Difference):F4}, {(other == "gc" ? rowPass.ToString() : "")}"));
        Console.WriteLine(rows[^1]);
    }
}
Console.WriteLine($"Verdict: {(pass ? "pass" : "fail")}.");

using (var csv = File.CreateText($"{prefix}_Islands.csv"))
{
    csv.WriteLine("Key, Detector, Variant, Genomes, Windows, Positives, BaseRate, AveragePrecision, ApLow, ApHigh, TestGenomes, Precision, Recall, F1, F1Low, F1High, Pass");
    foreach (var row in rows) csv.WriteLine(row);
}

// Each window's label: 1 when >= minOverlapBp of it lies in a positive region; with negatives given, 0 only inside a negative region
// and -1 (left out) elsewhere or in both; without them, 0 elsewhere.
static Genome Label(Chromosome chromosome, ZCurveWindowedProfile profile, List<Region> positive, List<Region>? negative)
{
    int n = profile.WindowCount;
    var inPositive = Covered(positive, n);
    var inNegative = negative is null ? null : Covered(negative, n);
    var labels = new int[n];
    for (int i = 0; i < n; i++)
        labels[i] = inNegative is null ? (inPositive[i] ? 1 : 0)
            : inPositive[i] == inNegative[i] ? -1 : inPositive[i] ? 1 : 0;
    return new Genome(chromosome, profile, labels);
}

// Regions are 1-based and inclusive; window i covers bases [i * window, (i + 1) * window) counted from 0.
static bool[] Covered(List<Region> regions, int windows)
{
    var overlap = new int[windows];
    foreach (var r in regions)
    {
        long start = r.Start - 1, end = r.End;
        for (long w = start / window; w < windows && w * window < end; w++)
            overlap[w] += (int)(Math.Min(end, (w + 1) * window) - Math.Max(start, w * window));
    }
    return overlap.Select(o => o >= minOverlapBp).ToArray();
}

static double[] Scores(ZCurveWindowedProfile profile, string detector, bool detrend)
{
    var windows = profile.Windows;
    int n = windows.Count;
    var axes = new List<double[]>();
    if (detector != "gc") { axes.Add(windows.Select(w => w.X).ToArray()); axes.Add(windows.Select(w => w.Y).ToArray()); }
    if (detector != "xy") axes.Add(windows.Select(w => w.Z).ToArray());
    var smoothed = axes.Select(values => Smooth(Subtract(values, detrend ? RunningMedian(values) : null))).ToArray();
    int d = smoothed.Length;
    var means = smoothed.Select(s => s.Average()).ToArray();
    var covariance = new double[d, d];
    for (int a = 0; a < d; a++)
        for (int b = 0; b < d; b++)
        {
            double sum = 0;
            for (int i = 0; i < n; i++) sum += (smoothed[a][i] - means[a]) * (smoothed[b][i] - means[b]);
            covariance[a, b] = sum / n + (a == b ? 1e-12 : 0);
        }
    var inverse = Invert(covariance);
    var scores = new double[n];
    var r = new double[d];
    for (int i = 0; i < n; i++)
    {
        for (int a = 0; a < d; a++) r[a] = smoothed[a][i] - means[a];
        double q = 0;
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++) q += r[a] * inverse[a, b] * r[b];
        scores[i] = q;
    }
    return scores;
}

static double[] Subtract(double[] values, double[]? baseline)
{
    double mean = values.Average();
    return values.Select((v, i) => v - (baseline?[i] ?? mean)).ToArray();
}

// The median of each medianWidth-window stretch at every medianStride-th window, linearly interpolated between them.
static double[] RunningMedian(double[] values)
{
    int n = values.Length;
    var centers = new List<int>();
    for (int c = 0; c < n; c += medianStride) centers.Add(c);
    if (centers[^1] != n - 1) centers.Add(n - 1);
    var medians = centers.Select(c =>
    {
        int from = Math.Max(0, c - medianWidth / 2), to = Math.Min(n, c + medianWidth / 2 + 1);
        var stretch = values[from..to];
        Array.Sort(stretch);
        return stretch[stretch.Length / 2];
    }).ToArray();
    var baseline = new double[n];
    for (int k = 0; k + 1 < centers.Count; k++)
        for (int i = centers[k]; i <= centers[k + 1]; i++)
        {
            double t = centers[k + 1] == centers[k] ? 0 : (double)(i - centers[k]) / (centers[k + 1] - centers[k]);
            baseline[i] = medians[k] + t * (medians[k + 1] - medians[k]);
        }
    if (centers.Count == 1) Array.Fill(baseline, medians[0]);
    return baseline;
}

// Mean over the smoothWidth windows starting smoothWidth / 2 before i, clipped at the ends.
static double[] Smooth(double[] values)
{
    int n = values.Length;
    var prefix = new double[n + 1];
    for (int i = 0; i < n; i++) prefix[i + 1] = prefix[i] + values[i];
    var result = new double[n];
    for (int i = 0; i < n; i++)
    {
        int from = Math.Max(0, i - smoothWidth / 2), to = Math.Min(n, i - smoothWidth / 2 + smoothWidth);
        result[i] = (prefix[to] - prefix[from]) / (to - from);
    }
    return result;
}

static double[,] Invert(double[,] m)
{
    int d = m.GetLength(0);
    var a = new double[d, 2 * d];
    for (int i = 0; i < d; i++) { for (int j = 0; j < d; j++) a[i, j] = m[i, j]; a[i, d + i] = 1; }
    for (int col = 0; col < d; col++)
    {
        int pivot = Enumerable.Range(col, d - col).MaxBy(r => Math.Abs(a[r, col]));
        for (int j = 0; j < 2 * d; j++) (a[col, j], a[pivot, j]) = (a[pivot, j], a[col, j]);
        double p = a[col, col];
        for (int j = 0; j < 2 * d; j++) a[col, j] /= p;
        for (int r = 0; r < d; r++)
            if (r != col) { double f = a[r, col]; for (int j = 0; j < 2 * d; j++) a[r, j] -= f * a[col, j]; }
    }
    var inverse = new double[d, d];
    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) inverse[i, j] = a[i, d + j];
    return inverse;
}

// Quantile bin edges over every window's score (from a seeded sample); labels are not used, so the test half doesn't leak.
static double[] Edges(double[][] scores, int count, int seed)
{
    var rng = new Random(seed);
    var all = scores.SelectMany(s => s).ToArray();
    var sample = all.Length <= 2_000_000 ? all : Enumerable.Range(0, 2_000_000).Select(_ => all[rng.Next(all.Length)]).ToArray();
    Array.Sort(sample);
    return Enumerable.Range(1, count - 1).Select(k => sample[(long)k * sample.Length / count]).ToArray();
}

static Histogram Histogram(double[] scores, int[] labels, double[] edges)
{
    var positive = new long[edges.Length + 1];
    var negative = new long[edges.Length + 1];
    for (int i = 0; i < scores.Length; i++)
    {
        if (labels[i] < 0) continue;
        int at = Array.BinarySearch(edges, scores[i]);
        int bin = at >= 0 ? at + 1 : ~at;
        if (labels[i] == 1) positive[bin]++; else negative[bin]++;
    }
    return new Histogram(positive, negative);
}

static Histogram Add(Histogram a, Histogram b) =>
    new(a.Positive.Zip(b.Positive, (x, y) => x + y).ToArray(), a.Negative.Zip(b.Negative, (x, y) => x + y).ToArray());

static Histogram Weighted(Histogram[] histograms, int[]? weights)
{
    int width = histograms[0].Positive.Length;
    var positive = new long[width];
    var negative = new long[width];
    for (int g = 0; g < histograms.Length; g++)
    {
        int w = weights?[g] ?? 1;
        if (w == 0) continue;
        for (int b = 0; b < width; b++) { positive[b] += w * histograms[g].Positive[b]; negative[b] += w * histograms[g].Negative[b]; }
    }
    return new Histogram(positive, negative);
}

// Windows in one bin share a score, so they enter the ranking together.
static double AveragePrecision(Histogram h)
{
    long total = h.Positive.Sum(), tp = 0, fp = 0;
    if (total == 0) return double.NaN;
    double ap = 0;
    for (int b = h.Positive.Length - 1; b >= 0; b--)
    {
        tp += h.Positive[b];
        fp += h.Negative[b];
        if (h.Positive[b] > 0) ap += (double)h.Positive[b] / total * tp / (tp + fp);
    }
    return ap;
}

// The bin from which windows are called islands that maximizes F1.
static int BestThreshold(Histogram h)
{
    long total = h.Positive.Sum(), tp = 0, fp = 0;
    double best = -1;
    int at = h.Positive.Length - 1;
    for (int b = h.Positive.Length - 1; b >= 0; b--)
    {
        tp += h.Positive[b];
        fp += h.Negative[b];
        double f1 = F1(tp, fp, total - tp);
        if (f1 > best) { best = f1; at = b; }
    }
    return at;
}

static (long Tp, long Fp, long Fn) Counts(Histogram h, int threshold)
{
    long tp = 0, fp = 0, fn = 0;
    for (int b = 0; b < h.Positive.Length; b++)
        if (b >= threshold) { tp += h.Positive[b]; fp += h.Negative[b]; } else fn += h.Positive[b];
    return (tp, fp, fn);
}

static (long Tp, long Fp, long Fn) Sum((long Tp, long Fp, long Fn)[] counts, int[]? weights)
{
    long tp = 0, fp = 0, fn = 0;
    for (int g = 0; g < counts.Length; g++)
    {
        int w = weights?[g] ?? 1;
        tp += w * counts[g].Tp; fp += w * counts[g].Fp; fn += w * counts[g].Fn;
    }
    return (tp, fp, fn);
}

static int[][] Draws(int genomes, int rounds, int seed)
{
    var rng = new Random(seed);
    return Enumerable.Range(0, rounds).Select(_ =>
    {
        var weights = new int[genomes];
        for (int j = 0; j < genomes; j++) weights[rng.Next(genomes)]++;
        return weights;
    }).ToArray();
}

static double F1Of((long Tp, long Fp, long Fn) c) => F1(c.Tp, c.Fp, c.Fn);
static double F1(long tp, long fp, long fn) => tp == 0 ? 0 : 2.0 * tp / (2.0 * tp + fp + fn);
static double Precision(long tp, long fp) => tp + fp == 0 ? double.NaN : (double)tp / (tp + fp);
static double Recall(long tp, long fn) => tp + fn == 0 ? double.NaN : (double)tp / (tp + fn);
static double Low(double[] sorted) => sorted[(int)(0.025 * sorted.Length)];
static double High(double[] sorted) => sorted[(int)(0.975 * sorted.Length) - 1];

// The curated key, scored with the threshold IslandPick's training half chose (no tuning on it).
static string Curated(List<Genome> genomes, string detector, bool detrend, double[] edges, int threshold)
{
    var histograms = genomes.Select(g => Histogram(Scores(g.Profile, detector, detrend), g.Labels, edges)).ToArray();
    var total = Weighted(histograms, null);
    var (tp, fp, fn) = Counts(total, threshold);
    long positives = total.Positive.Sum(), labelled = positives + total.Negative.Sum();
    return string.Create(CultureInfo.InvariantCulture,
        $"curated, {detector}, {(detrend ? "detrended" : "genomeMean")}, {genomes.Count}, {labelled}, {positives}, {(double)positives / labelled:F4}, " +
        $"{AveragePrecision(total):F4}, , , {genomes.Count}, {Precision(tp, fp):F4}, {Recall(tp, fn):F4}, {F1(tp, fp, fn):F4}, , , ");
}

// Regions per sequence, keyed by accession without version and without RefSeq's NZ_ prefix (an NZ_ record copies the GenBank one).
// The curated tables end their lines with a bare CR.
static Dictionary<string, List<Region>> Regions(string path, int startColumn)
{
    var regions = new Dictionary<string, List<Region>>(StringComparer.Ordinal);
    var lines = File.ReadAllText(path).Replace("\r\n", "\n").Replace('\r', '\n').Split('\n').Skip(1);
    foreach (var line in lines.Where(l => l.Trim() != ""))
    {
        var fields = line.Split('\t');
        var accession = fields[0].Trim().Split('.');
        var key = accession[0].StartsWith("NZ_", StringComparison.Ordinal) ? accession[0][3..] : accession[0];
        if (!regions.TryGetValue(key, out var list)) regions[key] = list = [];
        list.Add(new Region(long.Parse(fields[startColumn], CultureInfo.InvariantCulture), long.Parse(fields[startColumn + 1], CultureInfo.InvariantCulture),
            accession.Length > 1 ? accession[1] : null));
    }
    return regions;
}

static async Task DownloadAsync(string directory, string[] archives)
{
    Directory.CreateDirectory(directory);
    using var http = new HttpClient();
    foreach (var archive in archives)
    {
        var path = Path.Combine(directory, archive);
        if (File.Exists(path)) continue;
        Console.WriteLine($"Downloading {sourceUrl}{archive}");
        await File.WriteAllBytesAsync(path, await http.GetByteArrayAsync(sourceUrl + archive));
        await using var gzip = new GZipStream(File.OpenRead(path), CompressionMode.Decompress);
        await using var tar = new TarReader(gzip);
        while (await tar.GetNextEntryAsync() is { } entry)
            if (entry.EntryType is TarEntryType.RegularFile or TarEntryType.V7RegularFile && !Path.GetFileName(entry.Name).StartsWith("._", StringComparison.Ordinal))
                await entry.ExtractToFileAsync(Path.Combine(directory, Path.GetFileName(entry.Name)), overwrite: true);
    }
}

// RefSeq NC_ to the GenBank record it copies (nuccore summary's assemblyacc) and its length, cached so a rerun doesn't call NCBI.
static async Task<Dictionary<string, RefSeqLink>> RefSeqMapAsync(string cache, string[] refSeq)
{
    var map = new Dictionary<string, RefSeqLink>(StringComparer.Ordinal);
    if (File.Exists(cache))
    {
        foreach (var line in File.ReadLines(cache).Skip(1))
        {
            var f = line.Split(", ");
            map[f[0]] = new RefSeqLink(f[1], long.Parse(f[2], CultureInfo.InvariantCulture));
        }
        return map;
    }
    using var http = new HttpClient();
    foreach (var batch in refSeq.Chunk(100))
    {
        var url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&tool=cat-platform-p9&id=" + string.Join(",", batch);
        using var document = JsonDocument.Parse(await http.GetStringAsync(url));
        var result = document.RootElement.GetProperty("result");
        foreach (var uid in result.GetProperty("uids").EnumerateArray())
        {
            var summary = result.GetProperty(uid.GetString()!);
            var accession = summary.GetProperty("accessionversion").GetString()!.Split('.')[0];
            var genBank = summary.TryGetProperty("assemblyacc", out var a) ? a.GetString() : null;
            if (!string.IsNullOrEmpty(genBank)) map[accession] = new RefSeqLink(genBank.Split('.')[0], summary.GetProperty("slen").GetInt64());
        }
        await Task.Delay(400);   // NCBI asks for at most 3 requests a second without a key
    }
    using var csv = File.CreateText(cache);
    csv.WriteLine("RefSeq, GenBank, LengthBp");
    foreach (var (refSeqAccession, link) in map.OrderBy(m => m.Key, StringComparer.Ordinal))
        csv.WriteLine($"{refSeqAccession}, {link.GenBank}, {link.LengthBp}");
    return map;
}

record Region(long Start, long End, string? Version);

record Chromosome(long Id, string RecordId, long LengthBp, string Genus, string? Version);

record Genome(Chromosome Chromosome, ZCurveWindowedProfile Profile, int[] Labels);

record Histogram(long[] Positive, long[] Negative);

record Evaluation(Histogram[] Histograms, (long Tp, long Fp, long Fn)[] TestCounts, int Threshold, double[] ApRounds);

record RefSeqLink(string GenBank, long LengthBp);
