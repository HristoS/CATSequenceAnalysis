#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Text.Json;
using Npgsql;
using SequencePlatform.Core;

// P9 S9.11 (E2b), arm A6: a phage's host genus from stored 6-mer frequencies. E2 (round 1) passed weakly on log rho* (top-1 0.057 vs
// GC 0.024, but always guessing Escherichia scored 0.110); P9_KmerCheck put 6-mer frequencies at 0.256 against 0.150 on a balanced
// sample. Fixed before the run (2026-11-12):
// - Corpus: round 2 (records <= RoundTwoMaxId, fingerprint verified). Phages as E2: viruses whose host_name's first word is a stored
//   cellular genus (provider lineage). Hosts: every cellular corpus genome with a genus (its largest chromosome).
// - A6: both-strand 6-mer frequencies (each 6-mer counted with its reverse complement, 2,080 values; stored counts are one strand),
//   Euclidean. Baselines recomputed on the same set: log rho* (E2's arm, Euclidean), |delta GC|, the most common host, and A6 against
//   host labels shuffled across phages.
// - Prediction as E2: genera ranked by their nearest genome; top-1 and top-5 over every stored genus.
// - Pass: A6 minus log rho* in top-1 genus accuracy, bootstrap 95 % lower bound above 0, both overall and macro (mean over host genera
//   with >= 5 phages).
// Writes <prefix>_PhageHost.csv (each row with the corpus fingerprint).
// U3 (2026-11-12), after the run: --publish writes the run to ExperimentStore as e2b-phage-host with A6's five nearest host genera for
// every phage with a recorded host (kind "recorded") and for viruses without a host_name whose provider class has >= 50 recorded hosts,
// >= 90 % of them a stored cellular genus (kind "unrecorded"; the qualifying classes are printed first, and listed in <prefix>_Classes.csv).
// Each finding carries how often A6's first guess of that genus was right among recorded phages, the phage itself left out.
const int minHostPhages = 5, shuffles = 200, bootstraps = 1000, seed = 20261112, canonicalCount = 2080, guesses = 5, minClassRecorded = 50;
const double minClassBacterial = 0.9;
var prefix = args.FirstOrDefault(a => !a.StartsWith("--", StringComparison.Ordinal)) ?? "P9_E2b";
bool publish = args.Contains("--publish");

await using var dataSource = ExperimentCorpus.OpenLive();
var verified = ExperimentCorpus.VerifyRoundTwoAsync(dataSource);
var corpus = (await ExperimentCorpus.LoadAsync(dataSource, ExperimentCorpus.RoundTwoMaxId)).Where(g => !g.Excluded).ToList();
var cellular = corpus.Where(g => g.Group == SignatureIndex.Cellular && g.Provider.Genus is not null).ToList();
var generaIndex = cellular.Select(g => g.Provider.Genus!).Distinct().Order(StringComparer.Ordinal).Select((g, i) => (g, i))
    .ToDictionary(x => x.g, x => x.i, StringComparer.Ordinal);
var phages = corpus.Where(g => g.Group == SignatureIndex.Viral && HostGenus(g) is { } h && generaIndex.ContainsKey(h)).ToList();
var truth = phages.Select(p => generaIndex[HostGenus(p)!]).ToArray();
var genusOf = cellular.Select(g => generaIndex[g.Provider.Genus!]).ToArray();
Console.WriteLine($"{phages.Count} phages over {truth.Distinct().Count()} host genera; {cellular.Count} cellular genomes in {generaIndex.Count} genera.");

// Viruses without a recorded host, in classes whose recorded hosts are nearly all stored cellular genera (U3's rule).
var viral = corpus.Where(g => g.Group == SignatureIndex.Viral).ToList();
var classes = viral.Where(g => g.Provider.Class is not null && g.Field("host_name") is not null).GroupBy(g => g.Provider.Class!)
    .Select(c => (Class: c.Key, Recorded: c.Count(), Bacterial: c.Count(g => HostGenus(g) is { } h && generaIndex.ContainsKey(h)),
        Unrecorded: viral.Count(g => g.Provider.Class == c.Key && g.Field("host_name") is null)))
    .OrderByDescending(c => c.Recorded).ToList();
var qualifying = classes.Where(c => c.Recorded >= minClassRecorded && c.Bacterial >= minClassBacterial * c.Recorded).Select(c => c.Class).ToHashSet(StringComparer.Ordinal);
var unrecorded = viral.Where(g => g.Field("host_name") is null && g.Provider.Class is { } c && qualifying.Contains(c)).ToList();
using (var csv = File.CreateText($"{prefix}_Classes.csv"))
{
    csv.WriteLine("Class, Recorded, BacterialShare, Unrecorded, Qualifies");
    foreach (var c in classes.Where(c => c.Recorded >= minClassRecorded))
    {
        var line = string.Create(CultureInfo.InvariantCulture, $"{c.Class}, {c.Recorded}, {c.Bacterial / (double)c.Recorded:F3}, {c.Unrecorded}, {qualifying.Contains(c.Class)}");
        csv.WriteLine(line);
        Console.WriteLine("  " + line);
    }
}
Console.WriteLine($"Classes qualifying (>= {minClassRecorded} recorded hosts, >= {minClassBacterial:P0} a stored cellular genus): {string.Join(", ", qualifying.Order())}; " +
    $"{unrecorded.Count} viruses without a recorded host in them.");

var watch = Stopwatch.StartNew();
var frequencies = await FrequenciesAsync(dataSource, cellular.Select(g => g.Id).Concat(phages.Select(p => p.Id)).Concat(unrecorded.Select(u => u.Id)).ToArray());
var hostVectors = cellular.Select(g => frequencies[g.Id]).ToArray();
var phageVectors = phages.Select(p => frequencies[p.Id]).ToArray();
var unrecordedVectors = unrecorded.Select(u => frequencies[u.Id]).ToArray();
Console.WriteLine($"6-mer frequencies read in {watch.Elapsed.TotalSeconds:F0} s.");

watch.Restart();
var a6 = GenusOrders(phages.Count, generaIndex.Count, genusOf, (i, g) => SquaredFloat(phageVectors[i], hostVectors[g]));
var arms = new (string Name, int[][] Order)[]
{
    ("a6", a6.Orders),
    ("signature", GenusOrders(phages.Count, generaIndex.Count, genusOf, (i, g) => Squared(phages[i].Signature, cellular[g].Signature)).Orders),
    ("gc", GenusOrders(phages.Count, generaIndex.Count, genusOf, (i, g) => Math.Abs(phages[i].Gc - cellular[g].Gc)).Orders),
};
var a6Unrecorded = GenusOrders(unrecorded.Count, generaIndex.Count, genusOf, (i, g) => SquaredFloat(unrecordedVectors[i], hostVectors[g]));
Console.WriteLine($"Rankings in {watch.Elapsed.TotalSeconds:F0} s.");
var fingerprint = await verified;

var hostName = phages.Select(p => HostGenus(p)!).ToArray();
var macroGenera = hostName.GroupBy(h => h).Where(h => h.Count() >= minHostPhages).Select(h => h.Key).ToHashSet(StringComparer.Ordinal);
var rankOfTruth = arms.ToDictionary(a => a.Name, a => a.Order.Select((order, i) => Array.IndexOf(order, truth[i]) + 1).ToArray());
var common = truth.GroupBy(t => t).OrderByDescending(t => t.Count()).ThenBy(t => t.Key).Select(t => t.Key).ToArray();
rankOfTruth["mostCommon"] = truth.Select(t => Array.IndexOf(common, t) + 1).ToArray();

var rows = new List<string>();
foreach (var (name, ranks) in rankOfTruth)
    rows.Add(Row(name, ranks.Average(r => r <= 1 ? 1.0 : 0), ranks.Average(r => r <= 5 ? 1.0 : 0), Macro(ranks, 1), Macro(ranks, 5), null, null, null, null, null));

// A6 against host labels shuffled across phages.
var rng = new Random(seed);
var labels = truth.ToArray();
var shuffled1 = new double[shuffles];
for (int s = 0; s < shuffles; s++)
{
    rng.Shuffle(labels);
    shuffled1[s] = Enumerable.Range(0, phages.Count).Average(i => arms[0].Order[i][0] == labels[i] ? 1.0 : 0);
}
Array.Sort(shuffled1);
rows.Add(Row("a6Shuffled97.5", shuffled1[(int)(0.975 * shuffles) - 1], null, null, null, null, null, null, null, null));

bool pass = false;
var differences = new List<object>();
foreach (var baseline in new[] { "signature", "gc", "mostCommon" })
{
    var (overall, overallLow, overallHigh, macro, macroLow, macroHigh) = BootstrapDifference(
        rankOfTruth["a6"].Select(r => r == 1).ToArray(), rankOfTruth[baseline].Select(r => r == 1).ToArray(), hostName, macroGenera);
    bool beats = overallLow > 0 && macroLow > 0;
    if (baseline == "signature") pass = beats;
    differences.Add(new { baseline, overall, overallLow, overallHigh, macro, macroLow, macroHigh });
    rows.Add(Row($"a6Minus_{baseline}", overall, null, macro, null, overallLow, overallHigh, macroLow, macroHigh, baseline == "signature" ? beats : null));
}
foreach (var row in rows) Console.WriteLine(row);
Console.WriteLine($"Verdict (A6 beats log rho*, overall and macro): {(pass ? "pass" : "fail")}. {macroGenera.Count} host genera in the macro mean.");

using (var csv = File.CreateText($"{prefix}_PhageHost.csv"))
{
    csv.WriteLine("Predictor, Phages, Top1, Top5, MacroTop1, MacroTop5, MacroGenera, DiffLow, DiffHigh, MacroDiffLow, MacroDiffHigh, Pass, Fingerprint");
    foreach (var row in rows) csv.WriteLine($"{row}, {fingerprint}");
}

// How often A6's first guess of a genus was right among the recorded phages, and the nearest genera with their family, for findings.
var genusName = generaIndex.OrderBy(g => g.Value).Select(g => g.Key).ToArray();
var familyOfGenus = cellular.GroupBy(g => g.Provider.Genus!).ToDictionary(g => g.Key, g => g.Select(x => x.Provider.Family).Where(f => f is not null)
    .GroupBy(f => f!).OrderByDescending(f => f.Count()).ThenBy(f => f.Key, StringComparer.Ordinal).FirstOrDefault()?.Key, StringComparer.Ordinal);
var guessed = new int[genusName.Length];
var right = new int[genusName.Length];
for (int i = 0; i < phages.Count; i++)
{
    guessed[a6.Orders[i][0]]++;
    if (a6.Orders[i][0] == truth[i]) right[truth[i]]++;
}
object Guesses(int[] order, double[] distances) =>
    order.Take(guesses).Select((g, j) => new { genus = genusName[g], family = familyOfGenus[genusName[g]], distance = distances[j] }).ToArray();
Console.WriteLine($"First guesses: {guessed.Count(c => c > 0)} genera; the most guessed {string.Join(", ", Enumerable.Range(0, genusName.Length)
    .OrderByDescending(g => guessed[g]).Take(5).Select(g => $"{genusName[g]} {right[g]}/{guessed[g]}"))}.");

if (publish)
{
    var json = new JsonSerializerOptions { PropertyNamingPolicy = JsonNamingPolicy.CamelCase };
    var predictors = rankOfTruth.Select(r => new
    {
        predictor = r.Key,
        top1 = r.Value.Average(x => x <= 1 ? 1.0 : 0),
        top5 = r.Value.Average(x => x <= 5 ? 1.0 : 0),
        macroTop1 = Macro(r.Value, 1),
        macroTop5 = Macro(r.Value, 5),
    }).ToList();
    var summary = JsonSerializer.Serialize(new
    {
        rule = "Pass when 6-mers beat log ρ* in top-1 host genus, the bootstrap 95 % lower bound of the difference above 0, both over all phages " +
            $"and averaged over host genera with ≥ {minHostPhages} phages.",
        analysed = phages.Count,
        analysedNote = $"phages with a recorded host among {generaIndex.Count:N0} stored genera; {unrecorded.Count:N0} more without a recorded host get a guess",
        keyNumbers = new object[]
        {
            new { label = "top-1 host genus, 6-mers", value = predictors.First(p => p.predictor == "a6").top1 },
            new { label = "top-1, log ρ* (round 1)", value = predictors.First(p => p.predictor == "signature").top1 },
            new { label = "top-1, always the commonest host", value = predictors.First(p => p.predictor == "mostCommon").top1 },
        },
        predictors,
        shuffled = shuffled1[(int)(0.975 * shuffles) - 1],
        differences,
        macroGenera = macroGenera.Count,
        hostGenera = truth.Distinct().Count(),
        storedGenera = generaIndex.Count,
        classes = classes.Where(c => c.Recorded >= minClassRecorded).Select(c => new { name = c.Class, recorded = c.Recorded, bacterialShare = c.Bacterial / (double)c.Recorded,
            unrecorded = c.Unrecorded, qualifies = qualifying.Contains(c.Class) }),
        classRule = $"A virus without a recorded host gets a guess when its class has ≥ {minClassRecorded} recorded hosts, ≥ {minClassBacterial:P0} of them a stored bacterial or archaeal genus.",
        reading = "A guess is the nearest stored genera by 6-mer composition, not evidence of infection. The right genus comes first for about one " +
            "phage in five and is among the five for one in three.",
        fingerprint,
    }, json);
    var findings = new List<ExperimentFinding>();
    for (int i = 0; i < phages.Count; i++)
    {
        int first = a6.Orders[i][0];
        bool self = first == truth[i];
        findings.Add(new ExperimentFinding(phages[i].Id, "recorded", a6.Distances[i][0], JsonSerializer.Serialize(new
        {
            host = hostName[i], hostFamily = familyOfGenus[hostName[i]], rank = rankOfTruth["a6"][i], @class = phages[i].Provider.Class,
            guesses = Guesses(a6.Orders[i], a6.Distances[i]), firstGuess = new { guessed = guessed[first] - 1, right = right[first] - (self ? 1 : 0) },
        }, json)));
    }
    for (int i = 0; i < unrecorded.Count; i++)
    {
        int first = a6Unrecorded.Orders[i][0];
        findings.Add(new ExperimentFinding(unrecorded[i].Id, "unrecorded", a6Unrecorded.Distances[i][0], JsonSerializer.Serialize(new
        {
            @class = unrecorded[i].Provider.Class, guesses = Guesses(a6Unrecorded.Orders[i], a6Unrecorded.Distances[i]),
            firstGuess = new { guessed = guessed[first], right = right[first] },
        }, json)));
    }
    var (stored, skipped) = await ExperimentStore.PublishAsync(dataSource, "e2b-phage-host", "Phage host from 6-mers",
        "Which bacterial or archaeal genus does a phage's 6-mer composition point to, and how often is that its recorded host?",
        pass ? ExperimentStore.Pass : ExperimentStore.Fail, $"live, records <= {ExperimentCorpus.RoundTwoMaxId}, fingerprint {fingerprint}",
        "scripts/platform-checks/p9_e2b_phage_host.cs", summary, findings);
    Console.WriteLine($"Published e2b-phage-host: {stored} findings stored, {skipped} skipped.");
}

double? Macro(int[] ranks, int within) =>
    Enumerable.Range(0, ranks.Length).Where(i => macroGenera.Contains(hostName[i])).GroupBy(i => hostName[i]).Average(g => g.Average(i => ranks[i] <= within ? 1.0 : 0));

string Row(string predictor, double top1, double? top5, double? macro1, double? macro5, double? low, double? high, double? macroLow, double? macroHigh, bool? beats) =>
    string.Create(CultureInfo.InvariantCulture,
        $"{predictor}, {phages.Count}, {top1:F4}, {top5:F4}, {macro1:F4}, {macro5:F4}, {macroGenera.Count}, {low:F4}, {high:F4}, {macroLow:F4}, {macroHigh:F4}, {beats}");

static string? HostGenus(CorpusGenome g) => g.Field("host_name")?.Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault();

// Both-strand 6-mer frequencies per record: each one-strand count added to its canonical word (the smaller of a word and its reverse
// complement), then divided by the total.
static async Task<Dictionary<long, float[]>> FrequenciesAsync(NpgsqlDataSource dataSource, long[] ids)
{
    var canonical = new int[4096];
    var slot = new Dictionary<int, int>();
    for (int word = 0; word < 4096; word++)
    {
        int reverse = 0;
        for (int p = 0, w = word; p < 6; p++, w >>= 2) reverse = (reverse << 2) | (3 - (w & 3));
        int key = Math.Min(word, reverse);
        if (!slot.TryGetValue(key, out var s)) slot[key] = s = slot.Count;
        canonical[word] = s;
    }
    if (slot.Count != canonicalCount) throw new InvalidOperationException($"{slot.Count} canonical 6-mers, expected {canonicalCount}.");
    var result = new Dictionary<long, float[]>();
    await using var command = dataSource.CreateCommand("select sequence_id, counts from record_kmer where k = 6 and status = 'done' and sequence_id = any($1)");
    command.Parameters.AddWithValue(ids);
    command.CommandTimeout = 3600;
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync())
    {
        var counts = reader.GetFieldValue<int[]>(1);
        var folded = new double[canonicalCount];
        double total = 0;
        for (int w = 0; w < counts.Length; w++) { folded[canonical[w]] += counts[w]; total += counts[w]; }
        result[reader.GetInt64(0)] = folded.Select(c => (float)(c / total)).ToArray();
    }
    var missing = ids.Where(id => !result.ContainsKey(id)).ToList();
    if (missing.Count > 0) throw new InvalidOperationException($"{missing.Count} records without 6-mers, e.g. {missing[0]}.");
    return result;
}

// Per phage: genera ordered by their nearest genome (ties by genus index), and the first few genera's distances (square root of the
// distance function's value).
static (int[][] Orders, double[][] Distances) GenusOrders(int phageCount, int genusCount, int[] genusOf, Func<int, int, double> distance)
{
    var orders = new int[phageCount][];
    var distances = new double[phageCount][];
    Parallel.For(0, phageCount, i =>
    {
        var nearest = new double[genusCount];
        Array.Fill(nearest, double.PositiveInfinity);
        for (int g = 0; g < genusOf.Length; g++)
        {
            double d = distance(i, g);
            if (d < nearest[genusOf[g]]) nearest[genusOf[g]] = d;
        }
        orders[i] = Enumerable.Range(0, genusCount).OrderBy(t => nearest[t]).ThenBy(t => t).ToArray();
        distances[i] = orders[i].Take(guesses).Select(t => Math.Sqrt(nearest[t])).ToArray();
    });
    return (orders, distances);
}

static double Squared(double[] a, double[] b)
{
    double sum = 0;
    for (int f = 0; f < a.Length; f++) { double d = a[f] - b[f]; sum += d * d; }
    return sum;
}

static double SquaredFloat(float[] a, float[] b)
{
    var sum = Vector<float>.Zero;
    int i = 0;
    for (; i <= a.Length - Vector<float>.Count; i += Vector<float>.Count)
    {
        var d = new Vector<float>(a, i) - new Vector<float>(b, i);
        sum += d * d;
    }
    float total = Vector.Sum(sum);
    for (; i < a.Length; i++) { float d = a[i] - b[i]; total += d * d; }
    return total;
}

// A minus B in top-1 accuracy, overall and macro, with a percentile bootstrap over phages. As E2.
static (double Overall, double OverallLow, double OverallHigh, double Macro, double MacroLow, double MacroHigh) BootstrapDifference(
    bool[] a, bool[] b, string[] host, HashSet<string> macroGenera)
{
    int n = a.Length;
    (double Overall, double Macro) Difference(int[] picks)
    {
        double overall = picks.Average(j => (a[j] ? 1.0 : 0) - (b[j] ? 1.0 : 0));
        double macro = picks.Where(j => macroGenera.Contains(host[j])).GroupBy(j => host[j]).Average(g => g.Average(j => (a[j] ? 1.0 : 0) - (b[j] ? 1.0 : 0)));
        return (overall, macro);
    }
    var (overallAll, macroAll) = Difference(Enumerable.Range(0, n).ToArray());
    var rng = new Random(seed);
    var overalls = new double[bootstraps];
    var macros = new double[bootstraps];
    var draw = new int[n];
    for (int s = 0; s < bootstraps; s++)
    {
        for (int j = 0; j < n; j++) draw[j] = rng.Next(n);
        (overalls[s], macros[s]) = Difference(draw);
    }
    Array.Sort(overalls);
    Array.Sort(macros);
    int low = (int)(0.025 * bootstraps), high = (int)(0.975 * bootstraps) - 1;
    return (overallAll, overalls[low], overalls[high], macroAll, macros[low], macros[high]);
}
