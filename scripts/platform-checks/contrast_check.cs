#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.Text.RegularExpressions;
using SequencePlatform.Core;

// Check for a group-contrast runner ("can the signal tell A from B?") that catalogue visitors would run. Reads the live database,
// writes nothing to it. Fixed before the run (2026-11-12):
// - Runner: k nearest within A and B pooled (log rho*, Euclidean; baseline GC + log10 length, standardized), self excluded. A genome's
//   score is the share of its neighbours in A; AUC (Mann-Whitney) of that score, A against B; 1,000 bootstraps of genomes (neighbours
//   fixed) for the 95 % interval; 200 label shuffles for the null. "Separates": interval low above the shuffles' 95th percentile.
//   "Beyond GC": the signature's interval low above the baseline's AUC.
// - Replay of E5 (pass): per provider genus, the two largest habitats by E5's rules, >= 5 genomes each, k = 5. Must pass as E5 did:
//   excess (AUC minus shuffle mean) > 0 in more genera than < 0, one-sided sign test p < 0.05.
// - Negative control (stands in for E1, which scores single genomes, not two groups): each genus with >= 20 genomes cut into random
//   halves, k = 10. "Separates" in at most 10 % of genera, else the runner misleads.
// - Timing: 2,000 against 2,000, both spaces, at 4 and 2 threads (the server has 2 vCPU). At most 2 s, else it must be a worker job.
// - Coverage: within-genus habitat pairs and within-viral-family host-genus pairs with >= 50 genomes per side. At least 20 contrasts.
// Writes <prefix>_Replay.csv, <prefix>_Control.csv and <prefix>_Coverage.csv.
const int seed = 20261112, minSide = 50, timingSide = 2000;
var prefix = args.FirstOrDefault() ?? "ContrastCheck";
var rules = new (string Habitat, Regex Pattern)[]
{
    ("none", new(@"^(missing|not available.*|not applicable|not collected|not provided|unknown|n/?a|-|none|environmental|environment|water|r2a|envo:\d+|other|sample)$")),
    ("clinical", new(@"blood|sputum|urine|wound|patient|hospital|clinical|abscess|cerebrospinal|\bcsf\b|bronch|trache|catheter|\bpus\b|swab|peritoneal|pleural|aspirate|biopsy")),
    ("food", new(@"\bfood|milk|cheese|dairy|meat|sausage|salami|kimchi|fermented|yog(h)?urt|kefir|\bwine|\bbeer|vinegar|sourdough|pickle|sauerkraut")),
    ("host-associated", new(@"^[a-z ]+:|fec(es|al)|faec|stool|\bgut\b|intestin|ca?ecal|cecum|caecum|rumen|oral|skin|nasal|saliva|nodule|\broot|\bleaf|phyllosphere|endophyt|insect|beetle|\btick|mosquito|larva|\bhost\b")),
    ("marine", new(@"marine|\bsea\b|seawater|ocean|coastal|hydrothermal|cold seep|estuar|mangrove|coral|deep-sea")),
    ("freshwater", new(@"freshwater|fresh water|\blake|\briver|\bpond|stream|groundwater|drinking water")),
    ("soil", new(@"soil|rhizosphere|compost|sediment|\bmud\b|permafrost|desert")),
};
string? Habitat(CorpusGenome g) => g.Field("isolation_source") is { } source
    ? rules.FirstOrDefault(r => r.Pattern.IsMatch(source.Trim().ToLowerInvariant())).Habitat is { } h and not ("none" or "other") ? h : null
    : null;

await using var dataSource = ExperimentCorpus.OpenLive();
var corpus = (await ExperimentCorpus.LoadAsync(dataSource)).Where(g => !g.Excluded).ToList();
var cellular = corpus.Where(g => g.Group == SignatureIndex.Cellular).ToList();
var viral = corpus.Where(g => g.Group == SignatureIndex.Viral).ToList();
Console.WriteLine($"{cellular.Count} cellular, {viral.Count} viral genomes.");

// E5 replay.
var replay = new List<(string Genus, string A, string B, int NA, int NB, Contrast Result)>();
foreach (var genus in cellular.Where(g => g.Provider.Genus is not null && Habitat(g) is not null).GroupBy(g => g.Provider.Genus!)
    .OrderBy(g => g.Key, StringComparer.Ordinal))
{
    var habitats = genus.GroupBy(g => Habitat(g)!).Where(h => h.Count() >= 5).OrderByDescending(h => h.Count()).ThenBy(h => h.Key, StringComparer.Ordinal).Take(2).ToList();
    if (habitats.Count < 2) continue;
    replay.Add((genus.Key, habitats[0].Key, habitats[1].Key, habitats[0].Count(), habitats[1].Count(),
        Run(habitats[0].ToList(), habitats[1].ToList(), 5, seed + replay.Count, -1)));
}
int positive = replay.Count(r => r.Result.Signature.Auc - r.Result.Signature.NullMean > 1e-9);
int negative = replay.Count(r => r.Result.Signature.Auc - r.Result.Signature.NullMean < -1e-9);
double replayP = SignTest(positive, negative);
bool replayPass = positive > negative && replayP < 0.05;
Console.WriteLine($"E5 replay: {replay.Count} genera; excess > 0 in {positive}, < 0 in {negative} (sign test p {replayP:G3}); " +
    $"separates {replay.Count(r => r.Result.Signature.Separates)}, beyond GC {replay.Count(r => r.Result.BeyondGc)}, GC separates {replay.Count(r => r.Result.Gc.Separates)}. " +
    $"Reproduces E5's pass: {replayPass}.");
using (var csv = File.CreateText($"{prefix}_Replay.csv"))
{
    csv.WriteLine("Genus, HabitatA, HabitatB, GenomesA, GenomesB, " + Contrast.Header);
    foreach (var r in replay) csv.WriteLine($"\"{r.Genus}\", {r.A}, {r.B}, {r.NA}, {r.NB}, {r.Result.Row}");
}

// Negative control: random halves of one genus.
var control = new List<(string Genus, int N, Contrast Result)>();
foreach (var genus in cellular.Where(g => g.Provider.Genus is not null).GroupBy(g => g.Provider.Genus!).Where(g => g.Count() >= 20)
    .OrderBy(g => g.Key, StringComparer.Ordinal))
{
    var members = genus.OrderBy(g => g.Id).ToList();
    new Random(seed + 7919 * control.Count).Shuffle(System.Runtime.InteropServices.CollectionsMarshal.AsSpan(members));
    int half = members.Count / 2;
    control.Add((genus.Key, members.Count, Run(members.GetRange(0, half), members.GetRange(half, members.Count - half), 10, seed + control.Count, -1)));
}
double falseSignature = control.Count(c => c.Result.Signature.Separates) / (double)control.Count;
double falseGc = control.Count(c => c.Result.Gc.Separates) / (double)control.Count;
double falseBeyond = control.Count(c => c.Result.BeyondGc) / (double)control.Count;
bool controlPass = falseSignature <= 0.10;
Console.WriteLine($"Control: {control.Count} genera split at random; signature separates {falseSignature:P1}, GC {falseGc:P1}, beyond GC {falseBeyond:P1}. " +
    $"Within 10 %: {controlPass}.");
using (var csv = File.CreateText($"{prefix}_Control.csv"))
{
    csv.WriteLine("Genus, Genomes, " + Contrast.Header);
    foreach (var c in control) csv.WriteLine($"\"{c.Genus}\", {c.N}, {c.Result.Row}");
}

// Timing: the two largest cellular phyla, 2,000 each.
var phyla = cellular.Where(g => g.Provider.Phylum is not null).GroupBy(g => g.Provider.Phylum!).OrderByDescending(p => p.Count()).Take(2).ToList();
var sideA = phyla[0].OrderBy(g => g.Id).Take(timingSide).ToList();
var sideB = phyla[1].OrderBy(g => g.Id).Take(timingSide).ToList();
Run(sideA, sideB, 10, seed, -1);   // warm-up (JIT)
var timings = new List<string>();
foreach (int threads in new[] { 4, 2 })
{
    var seconds = new List<double>();
    for (int repeat = 0; repeat < 5; repeat++)
    {
        var watch = Stopwatch.StartNew();
        Run(sideA, sideB, 10, seed, threads);
        seconds.Add(watch.Elapsed.TotalSeconds);
    }
    double median = seconds.Order().ElementAt(2);
    timings.Add(string.Create(CultureInfo.InvariantCulture, $"{threads} threads {median:F2} s"));
}
var timed = Run(sideA, sideB, 10, seed, -1);
bool timingPass = timings.All(t => double.Parse(t.Split(' ')[2], CultureInfo.InvariantCulture) <= 2.0);
Console.WriteLine($"Timing, {phyla[0].Key} against {phyla[1].Key}, {sideA.Count} + {sideB.Count}: {string.Join(", ", timings)} (median of 5). " +
    $"Signature AUC {timed.Signature.Auc:F3}, GC {timed.Gc.Auc:F3}. Within 2 s: {timingPass}.");

// Coverage: contrasts a visitor could ask with >= 50 genomes per side.
var coverage = new List<string>();
int habitatPairs = 0, hostPairs = 0, habitatGenera = 0, hostFamilies = 0;
foreach (var genus in cellular.Where(g => g.Provider.Genus is not null && Habitat(g) is not null).GroupBy(g => g.Provider.Genus!))
{
    var sides = genus.GroupBy(g => Habitat(g)!).Where(h => h.Count() >= minSide).Select(h => $"{h.Key} {h.Count()}").ToList();
    int pairs = sides.Count * (sides.Count - 1) / 2;
    if (pairs == 0) continue;
    habitatPairs += pairs;
    habitatGenera++;
    coverage.Add($"habitat, \"{genus.Key}\", {pairs}, \"{string.Join("; ", sides)}\"");
}
foreach (var family in viral.Where(g => g.Provider.Family is not null && HostGenus(g) is not null).GroupBy(g => g.Provider.Family!))
{
    var sides = family.GroupBy(g => HostGenus(g)!).Where(h => h.Count() >= minSide).Select(h => $"{h.Key} {h.Count()}").ToList();
    int pairs = sides.Count * (sides.Count - 1) / 2;
    if (pairs == 0) continue;
    hostPairs += pairs;
    hostFamilies++;
    coverage.Add($"host, \"{family.Key}\", {pairs}, \"{string.Join("; ", sides)}\"");
}
bool coveragePass = habitatPairs + hostPairs >= 20;
Console.WriteLine($"Coverage (>= {minSide} per side): {habitatPairs} habitat pairs in {habitatGenera} genera, {hostPairs} host-genus pairs in {hostFamilies} viral families. " +
    $"At least 20: {coveragePass}.");
using (var csv = File.CreateText($"{prefix}_Coverage.csv"))
{
    csv.WriteLine("Kind, Taxon, Pairs, Sides");
    foreach (var line in coverage) csv.WriteLine(line);
}
Console.WriteLine($"Verdict: replay {replayPass}, control {controlPass}, timing {timingPass}, coverage {coveragePass}.");

static string? HostGenus(CorpusGenome g) =>
    g.Field("host_name") is { } host && host.Split(' ', StringSplitOptions.RemoveEmptyEntries) is [var first, ..] && first.Length > 2 && char.IsUpper(first[0])
        && first.All(char.IsLetter) && first is not ("Homo" or "Unknown" or "Environmental") ? first : null;

static Contrast Run(List<CorpusGenome> a, List<CorpusGenome> b, int k, int seed, int threads)
{
    var members = a.Concat(b).ToArray();
    var inA = members.Select((_, i) => i < a.Count).ToArray();
    var signature = Space.Evaluate(Neighbours(members.Select(g => g.Signature).ToArray(), k, threads), inA, seed);
    var gc = Space.Evaluate(Neighbours(Standardize(members.Select(g => new[] { g.Gc, Math.Log10(g.LengthBp) }).ToArray()), k, threads), inA, seed + 1);
    return new Contrast(signature, gc);
}

static double[][] Standardize(double[][] points)
{
    int dimensions = points[0].Length;
    var result = points.Select(p => (double[])p.Clone()).ToArray();
    for (int d = 0; d < dimensions; d++)
    {
        double mean = points.Average(p => p[d]);
        double sd = Math.Sqrt(points.Average(p => (p[d] - mean) * (p[d] - mean)));
        foreach (var p in result) p[d] = sd > 0 ? (p[d] - mean) / sd : 0;
    }
    return result;
}

// The k nearest of each point among the others, brute force; ties by index.
static int[][] Neighbours(double[][] points, int k, int threads)
{
    int n = points.Length;
    k = Math.Min(k, n - 1);
    var result = new int[n][];
    var options = new ParallelOptions { MaxDegreeOfParallelism = threads };
    Parallel.For(0, n, options, i =>
    {
        var best = new int[k];
        var bestDistance = new double[k];
        Array.Fill(bestDistance, double.PositiveInfinity);
        var p = points[i];
        for (int j = 0; j < n; j++)
        {
            if (j == i) continue;
            var q = points[j];
            double sum = 0;
            for (int d = 0; d < p.Length; d++) { double diff = p[d] - q[d]; sum += diff * diff; }
            if (sum >= bestDistance[k - 1]) continue;
            int slot = k - 1;
            while (slot > 0 && bestDistance[slot - 1] > sum) { bestDistance[slot] = bestDistance[slot - 1]; best[slot] = best[slot - 1]; slot--; }
            bestDistance[slot] = sum;
            best[slot] = j;
        }
        result[i] = best;
    });
    return result;
}

static double SignTest(int positive, int negative)
{
    int n = positive + negative;
    if (n == 0) return 1;
    double p = 0;
    for (int x = positive; x <= n; x++) p += Math.Exp(LogChoose(n, x) - n * Math.Log(2));
    return Math.Min(1, p);
}

static double LogChoose(int n, int x) => LogFactorial(n) - LogFactorial(x) - LogFactorial(n - x);

static double LogFactorial(int n)
{
    double sum = 0;
    for (int i = 2; i <= n; i++) sum += Math.Log(i);
    return sum;
}

internal sealed record Space(double Auc, double Low, double High, double NullMean, double Null95)
{
    private const int shuffles = 200, bootstraps = 1000;

    public bool Separates => Low > Null95;

    public static Space Evaluate(int[][] neighbours, bool[] inA, int seed)
    {
        double auc = MannWhitney(Scores(neighbours, inA), inA);
        var rng = new Random(seed);
        var shuffled = (bool[])inA.Clone();
        var nulls = new double[shuffles];
        for (int s = 0; s < shuffles; s++)
        {
            rng.Shuffle(shuffled);
            nulls[s] = MannWhitney(Scores(neighbours, shuffled), shuffled);
        }
        Array.Sort(nulls);
        var scores = Scores(neighbours, inA);
        var aIndex = Enumerable.Range(0, inA.Length).Where(i => inA[i]).ToArray();
        var bIndex = Enumerable.Range(0, inA.Length).Where(i => !inA[i]).ToArray();
        var boot = new double[bootstraps];
        for (int r = 0; r < bootstraps; r++)
        {
            var sampleA = aIndex.Select(_ => scores[aIndex[rng.Next(aIndex.Length)]]).ToArray();
            var sampleB = bIndex.Select(_ => scores[bIndex[rng.Next(bIndex.Length)]]).ToArray();
            boot[r] = MannWhitney(sampleA, sampleB);
        }
        Array.Sort(boot);
        return new Space(auc, boot[(int)(0.025 * bootstraps)], boot[(int)(0.975 * bootstraps) - 1], nulls.Average(), nulls[(int)(0.95 * shuffles)]);
    }

    private static double[] Scores(int[][] neighbours, bool[] inA) =>
        neighbours.Select(n => n.Count(j => inA[j]) / (double)n.Length).ToArray();

    private static double MannWhitney(double[] scores, bool[] inA) =>
        MannWhitney(scores.Where((_, i) => inA[i]).ToArray(), scores.Where((_, i) => !inA[i]).ToArray());

    // Mann-Whitney: P(score of an A genome > score of a B genome) + half the ties. Scores take at most k + 1 values, so count them.
    private static double MannWhitney(double[] a, double[] b)
    {
        var bSorted = (double[])b.Clone();
        Array.Sort(bSorted);
        double sum = 0;
        foreach (var x in a)
        {
            int below = LowerBound(bSorted, x), notAbove = UpperBound(bSorted, x);
            sum += below + 0.5 * (notAbove - below);
        }
        return sum / ((double)a.Length * b.Length);
    }

    private static int LowerBound(double[] sorted, double x)
    {
        int low = 0, high = sorted.Length;
        while (low < high) { int mid = (low + high) / 2; if (sorted[mid] < x) low = mid + 1; else high = mid; }
        return low;
    }

    private static int UpperBound(double[] sorted, double x)
    {
        int low = 0, high = sorted.Length;
        while (low < high) { int mid = (low + high) / 2; if (sorted[mid] <= x) low = mid + 1; else high = mid; }
        return low;
    }
}

internal sealed record Contrast(Space Signature, Space Gc)
{
    public const string Header = "SignatureAuc, SignatureLow, SignatureHigh, SignatureNullMean, SignatureNull95, SignatureSeparates, " +
        "GcAuc, GcLow, GcHigh, GcNullMean, GcNull95, GcSeparates, BeyondGc";

    public bool BeyondGc => Signature.Low > Gc.Auc;

    public string Row => string.Create(CultureInfo.InvariantCulture,
        $"{Signature.Auc:F4}, {Signature.Low:F4}, {Signature.High:F4}, {Signature.NullMean:F4}, {Signature.Null95:F4}, {Signature.Separates}, " +
        $"{Gc.Auc:F4}, {Gc.Low:F4}, {Gc.High:F4}, {Gc.NullMean:F4}, {Gc.Null95:F4}, {Gc.Separates}, {BeyondGc}");
}
