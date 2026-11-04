#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;
using System.Text.RegularExpressions;
using SequencePlatform.Core;

// P9 S9.5 (E5): traits beyond ancestry. Within a genus, do a genome's nearest neighbours by signature share its habitat more than
// a shuffle of habitat labels inside that genus predicts? Shuffling only inside the genus keeps the genus's ancestry in the null.
// Fixed before the run:
// - Habitat from BV-BRC isolation_source by the keyword rules below, first match wins: clinical, food, host-associated, marine,
//   freshwater, soil, other. Vague values ("water", "environmental", "missing", ...) get no label. "other" is matched but not scored.
// - Per genus (provider lineage) with >= 10 labelled genomes and >= 2 habitats: each genome's 5 nearest labelled genomes of the genus
//   (log rho*, Euclidean), the mean share with the same habitat; against 1,000 shuffles of the genus's labels. Excess = observed
//   minus the shuffle mean.
// - Stop before scoring if fewer than 20 genera qualify.
// - Pass: more genera with positive excess than negative, one-sided sign test p < 0.05.
// - Context, not the rule: the same with |delta GC| as the distance, and the same per species (ancestry within a species remains).
// Writes <prefix>_Habitat.csv (per genus and species) and <prefix>_Rules.csv (coverage of each rule).
const int k = 5, minGenomes = 10, minGenera = 20, shuffles = 1000, seed = 20261104;
var prefix = args.FirstOrDefault() ?? "P9_E5";
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
await using var dataSource = ExperimentCorpus.OpenSnapshot();
var cellular = (await ExperimentCorpus.LoadAsync(dataSource)).Where(g => !g.Excluded && g.Group == SignatureIndex.Cellular).ToList();
var labelled = new List<(CorpusGenome Genome, string Habitat)>();
var coverage = new Dictionary<string, int>();
int withSource = 0;
foreach (var g in cellular)
{
    if (g.Field("isolation_source") is not { } source) continue;
    withSource++;
    var text = source.Trim().ToLowerInvariant();
    var habitat = rules.FirstOrDefault(r => r.Pattern.IsMatch(text)).Habitat ?? "other";
    coverage[habitat] = coverage.GetValueOrDefault(habitat) + 1;
    if (habitat is not ("none" or "other")) labelled.Add((g, habitat));
}
using (var csv = File.CreateText($"{prefix}_Rules.csv"))
{
    csv.WriteLine("Habitat, Genomes, ShareOfWithSource");
    foreach (var (habitat, count) in coverage.OrderByDescending(c => c.Value))
        csv.WriteLine(string.Create(CultureInfo.InvariantCulture, $"{habitat}, {count}, {(double)count / withSource:F4}"));
}
Console.WriteLine($"{cellular.Count} cellular genomes, {withSource} with an isolation source, {labelled.Count} labelled: " +
    string.Join(", ", coverage.OrderByDescending(c => c.Value).Select(c => $"{c.Key} {c.Value}")));

var rows = new List<string>();
bool pass = false;
foreach (var rank in new[] { "genus", "species" })
{
    var groups = labelled.Where(l => (rank == "genus" ? l.Genome.Provider.Genus : l.Genome.Provider.Species) is not null)
        .GroupBy(l => rank == "genus" ? l.Genome.Provider.Genus! : l.Genome.Provider.Species!)
        .Where(g => g.Count() >= minGenomes && g.Select(l => l.Habitat).Distinct().Count() >= 2)
        .OrderBy(g => g.Key, StringComparer.Ordinal).ToList();
    Console.WriteLine($"{rank}: {groups.Count} qualify (>= {minGenomes} labelled genomes, >= 2 habitats).");
    if (rank == "genus" && groups.Count < minGenera)
    {
        Console.WriteLine($"Stop: fewer than {minGenera} genera qualify.");
        break;
    }
    if (groups.Count == 0) continue;   // species: the catalogue holds about one reference genome per bacterial species

    var results = groups.Select((group, index) => (group, index)).AsParallel().Select(item =>
    {
        var (group, index) = item;
        var members = group.ToArray();
        var labels = members.Select(m => m.Habitat).ToArray();
        var (signature, signatureNull, signatureP) = Score(Neighbours(members, (a, b) => Squared(a.Signature, b.Signature)), labels, seed + 2 * index);
        var (gc, gcNull, _) = Score(Neighbours(members, (a, b) => Math.Abs(a.Gc - b.Gc)), labels, seed + 2 * index + 1);
        return (Name: group.Key, Genomes: members.Length, Habitats: labels.Distinct().Count(),
            Majority: labels.GroupBy(l => l).Max(l => l.Count()) / (double)members.Length,
            Signature: signature, SignatureNull: signatureNull, SignatureP: signatureP, Gc: gc, GcNull: gcNull);
    }).OrderBy(r => r.Name, StringComparer.Ordinal).ToList();

    foreach (var r in results)
        rows.Add(string.Create(CultureInfo.InvariantCulture,
            $"{rank}, \"{r.Name}\", {r.Genomes}, {r.Habitats}, {r.Majority:F3}, {r.Signature:F4}, {r.SignatureNull:F4}, {r.Signature - r.SignatureNull:F4}, {r.SignatureP:F4}, " +
            $"{r.Gc:F4}, {r.GcNull:F4}, {r.Gc - r.GcNull:F4}, , "));
    int positive = results.Count(r => r.Signature - r.SignatureNull > 1e-9), negative = results.Count(r => r.Signature - r.SignatureNull < -1e-9);
    int gcPositive = results.Count(r => r.Gc - r.GcNull > 1e-9), gcNegative = results.Count(r => r.Gc - r.GcNull < -1e-9);
    int better = results.Count(r => r.Signature - r.SignatureNull > r.Gc - r.GcNull + 1e-9), worse = results.Count(r => r.Signature - r.SignatureNull < r.Gc - r.GcNull - 1e-9);
    double p = SignTest(positive, negative), gcP = SignTest(gcPositive, gcNegative), versusGc = SignTest(better, worse);
    bool rankPass = positive > negative && p < 0.05;
    if (rank == "genus") pass = rankPass;
    var summary = string.Create(CultureInfo.InvariantCulture,
        $"{rank}, all, {results.Sum(r => r.Genomes)}, {results.Count}, , {results.Average(r => r.Signature):F4}, {results.Average(r => r.SignatureNull):F4}, " +
        $"{results.Average(r => r.Signature - r.SignatureNull):F4}, {p:G3}, {results.Average(r => r.Gc):F4}, {results.Average(r => r.GcNull):F4}, " +
        $"{results.Average(r => r.Gc - r.GcNull):F4}, {versusGc:G3}, {(rank == "genus" ? rankPass.ToString() : "")}");
    rows.Add(summary);
    Console.WriteLine($"{rank}: signature excess > 0 in {positive}, < 0 in {negative} (sign test p {p:G3}; {results.Count(r => r.SignatureP < 0.05)} with p < 0.05 alone); " +
        $"GC {gcPositive} / {gcNegative} (p {gcP:G3}); signature beats GC in {better}, loses in {worse} (p {versusGc:G3}).");
    Console.WriteLine(summary);
}
Console.WriteLine($"Verdict: {(pass ? "pass" : "fail")}.");

using (var csv = File.CreateText($"{prefix}_Habitat.csv"))
{
    csv.WriteLine("Rank, Taxon, Genomes, Habitats, MajorityShare, SignatureShare, SignatureShuffled, SignatureExcess, SignatureP, GcShare, GcShuffled, GcExcess, SignVsGcP, Pass");
    foreach (var row in rows) csv.WriteLine(row);
}

static double Squared(double[] a, double[] b)
{
    double sum = 0;
    for (int f = 0; f < a.Length; f++) { double d = a[f] - b[f]; sum += d * d; }
    return sum;
}

// The k nearest of each member within its group, nearest first, itself excluded; ties by index.
static int[][] Neighbours((CorpusGenome Genome, string Habitat)[] members, Func<CorpusGenome, CorpusGenome, double> distance) =>
    members.Select((m, i) => Enumerable.Range(0, members.Length).Where(j => j != i)
        .OrderBy(j => distance(m.Genome, members[j].Genome)).ThenBy(j => j).Take(k).ToArray()).ToArray();

// Mean share of same-habitat neighbours; its mean over within-group shuffles of the labels; the one-sided p of the observed share.
static (double Observed, double Shuffled, double P) Score(int[][] neighbours, string[] labels, int seed)
{
    double Share(string[] l) => neighbours.Select((n, i) => n.Count(j => l[j] == l[i]) / (double)n.Length).Average();
    double observed = Share(labels);
    var rng = new Random(seed);
    var shuffled = labels.ToArray();
    double sum = 0;
    int atLeast = 0;
    for (int s = 0; s < shuffles; s++)
    {
        rng.Shuffle(shuffled);
        double value = Share(shuffled);
        sum += value;
        if (value >= observed - 1e-12) atLeast++;
    }
    return (observed, sum / shuffles, (atLeast + 1.0) / (shuffles + 1));
}

// One-sided binomial sign test: P(X >= positive) with X ~ Binomial(positive + negative, 1/2).
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
