#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;
using System.Numerics;
using System.Text.RegularExpressions;
using Npgsql;
using SequencePlatform.Core;

// P9 S9.12 (E5b): traits beyond ancestry with richer signatures. E5 (round 1, p9_e5_habitat.cs) passed on log rho*: within a genus,
// a genome's nearest neighbours share its habitat more than a shuffle of the genus's labels predicts. Which signature carries more?
// Fixed before the run (2026-11-22):
// - Corpus: round 2 (records <= RoundTwoMaxId, fingerprint verified), bacteria and archaea, partial left out. Habitat labels, the
//   >= 10 genomes / >= 2 habitats rule, k = 5 and 1,000 shuffles exactly as E5.
// - Signatures on the same genomes: log rho* (Euclidean; the round-1 arm, recomputed), 6-mer frequency (both strands folded, Euclidean;
//   E2b's A6), amino-acid composition (1-mers derived from the stored 4-mers, frequencies, Euclidean), amino-acid CVTree (Core's
//   CvTreeSignature, the relatives signal); GC (|delta GC|) as the baseline. Genomes without amino-acid 4-mers are left out of all arms.
// - Per genus, every signature sees the same 1,000 shuffles, so their excesses (observed share minus shuffle mean) pair.
// - Pass, per challenger (6-mers, amino-acid composition, CVTree) on provider genus: more genera where its excess beats log rho*'s than
//   where it loses, one-sided sign test p < 0.05 / 3. None passes: E5's result stands.
// - Context, not the rule: the same on GTDB genus (as GTDB writes it, so Bacillus_A is its own genus), which splits many provider genera.
//   A challenger that passes on provider genus but not here is reported as finer ancestry, not lifestyle.
// Writes <prefix>_Habitat.csv (per grouping and genus, each row with the corpus fingerprint).
const int k = 5, minGenomes = 10, shuffles = 1000, batch = 1000, canonicalCount = 2080, seed = 20261122;
const double alpha = 0.05 / 3;
var prefix = args.FirstOrDefault() ?? "P9_E5b";
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
var challengers = new[] { "kmer6", "aminoAcid", "cvtree" };

await using var dataSource = ExperimentCorpus.OpenRoundTwo();
var verified = ExperimentCorpus.VerifyRoundTwoAsync(dataSource);
var cellular = (await ExperimentCorpus.LoadAsync(dataSource, ExperimentCorpus.RoundTwoMaxId))
    .Where(g => !g.Excluded && g.Group == SignatureIndex.Cellular).ToList();
var labelled = new List<(CorpusGenome Genome, string Habitat)>();
foreach (var g in cellular)
{
    if (g.Field("isolation_source") is not { } source) continue;
    var habitat = rules.FirstOrDefault(r => r.Pattern.IsMatch(source.Trim().ToLowerInvariant())).Habitat ?? "other";
    if (habitat is not ("none" or "other")) labelled.Add((g, habitat));
}
var ids = labelled.Select(l => l.Genome.Id).ToArray();
var proteins = await ProteinSignaturesAsync(dataSource, ids);
int withoutProteins = labelled.Count(l => !proteins.ContainsKey(l.Genome.Id));
labelled = labelled.Where(l => proteins.ContainsKey(l.Genome.Id)).ToList();
var kmers = await FrequenciesAsync(dataSource, labelled.Select(l => l.Genome.Id).ToArray());
Console.WriteLine($"{cellular.Count} cellular genomes, {labelled.Count} labelled with amino-acid 4-mers ({withoutProteins} labelled without, left out).");

var distances = new (string Name, Func<CorpusGenome, CorpusGenome, double> Distance)[]
{
    ("logRho", (a, b) => Squared(a.Signature, b.Signature)),
    ("kmer6", (a, b) => SquaredFloat(kmers[a.Id], kmers[b.Id])),
    ("aminoAcid", (a, b) => Squared(proteins[a.Id].Composition, proteins[b.Id].Composition)),
    ("cvtree", (a, b) => CvTreeSignature.Distance(proteins[a.Id].CvTree, proteins[b.Id].CvTree)),
    ("gc", (a, b) => Math.Abs(a.Gc - b.Gc)),
};

var rows = new List<string>();
var verdicts = new List<string>();
foreach (var grouping in new[] { "provider", "gtdb" })
{
    var groups = labelled.Where(l => Genus(l.Genome, grouping) is not null).GroupBy(l => Genus(l.Genome, grouping)!)
        .Where(g => g.Count() >= minGenomes && g.Select(l => l.Habitat).Distinct().Count() >= 2)
        .OrderBy(g => g.Key, StringComparer.Ordinal).ToList();
    var results = groups.Select((group, index) => (group, index)).AsParallel().Select(item =>
    {
        var members = item.group.ToArray();
        var labels = members.Select(m => m.Habitat).ToArray();
        var excess = distances.ToDictionary(d => d.Name, d =>
        {
            var (observed, shuffled, _) = Score(Neighbours(members, d.Distance), labels, seed + item.index);
            return (Observed: observed, Excess: observed - shuffled);
        });
        return (Name: item.group.Key, Genomes: members.Length, Habitats: labels.Distinct().Count(), Excess: excess);
    }).OrderBy(r => r.Name, StringComparer.Ordinal).ToList();

    foreach (var r in results)
        rows.Add(string.Create(CultureInfo.InvariantCulture, $"{grouping}, {Quote(r.Name)}, {r.Genomes}, {r.Habitats}, ") +
            string.Join(", ", distances.Select(d => string.Create(CultureInfo.InvariantCulture, $"{r.Excess[d.Name].Observed:F4}, {r.Excess[d.Name].Excess:F4}"))));

    Console.WriteLine($"{grouping} genus: {results.Count} qualify, {results.Sum(r => r.Genomes)} genomes.");
    foreach (var d in distances)
    {
        int positive = results.Count(r => r.Excess[d.Name].Excess > 1e-9), negative = results.Count(r => r.Excess[d.Name].Excess < -1e-9);
        Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
            $"  {d.Name,-10} excess > 0 in {positive}, < 0 in {negative} (p {SignTest(positive, negative):G3}); mean excess {results.Average(r => r.Excess[d.Name].Excess):F4}"));
    }
    foreach (var name in challengers.Append("gc"))
    {
        int better = results.Count(r => r.Excess[name].Excess > r.Excess["logRho"].Excess + 1e-9);
        int worse = results.Count(r => r.Excess[name].Excess < r.Excess["logRho"].Excess - 1e-9);
        double p = SignTest(better, worse);
        bool passes = name != "gc" && better > worse && p < alpha;
        double meanDifference = results.Average(r => r.Excess[name].Excess - r.Excess["logRho"].Excess);
        var line = string.Create(CultureInfo.InvariantCulture,
            $"{grouping}, {name} vs logRho, better {better}, worse {worse}, p {p:G3}, mean difference {meanDifference:F4}, {(name == "gc" ? "baseline" : passes ? "pass" : "no")}");
        verdicts.Add(line);
        Console.WriteLine("  " + line);
    }
}

var fingerprint = await verified;
Console.WriteLine($"Corpus fingerprint {fingerprint} verified.");
using (var csv = File.CreateText($"{prefix}_Habitat.csv"))
{
    csv.WriteLine("Grouping, Genus, Genomes, Habitats, " + string.Join(", ", distances.Select(d => $"{d.Name}Share, {d.Name}Excess")) + ", Fingerprint");
    foreach (var row in rows) csv.WriteLine(row + ", " + fingerprint);
    foreach (var line in verdicts) csv.WriteLine("# " + line);
}

static string? Genus(CorpusGenome g, string grouping) => grouping == "provider" ? g.Provider.Genus : g.Gtdb?.Genus;

static string Quote(string text) => "\"" + text.Replace("\"", "\"\"") + "\"";

// Per record, from its stored amino-acid 4-mers: the 1-mer frequencies (each gene's last three residues aren't counted) and the CVTree
// vector; records with no 4-mers, or none with a word, are left out. Read in batches so the payloads never sit in memory at once.
static async Task<Dictionary<long, (double[] Composition, short[] CvTree)>> ProteinSignaturesAsync(NpgsqlDataSource dataSource, long[] ids)
{
    var result = new Dictionary<long, (double[] Composition, short[] CvTree)>();
    foreach (var chunk in ids.Chunk(batch))
    {
        var payloads = new List<(long Id, byte[] Payload)>();
        await using (var command = dataSource.CreateCommand("select sequence_id, payload from protein_kmer where k = $1 and status = $2 and sequence_id = any($3)"))
        {
            command.Parameters.AddWithValue((short)CvTreeSignature.K);
            command.Parameters.AddWithValue(ProteinKmers.Done);
            command.Parameters.AddWithValue(chunk);
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) payloads.Add((reader.GetInt64(0), reader.GetFieldValue<byte[]>(1)));
        }
        foreach (var (id, composition, cvtree) in payloads.AsParallel().Select(p =>
        {
            var counts = ProteinKmers.Decode(p.Payload, CvTreeSignature.K);
            var singles = ProteinKmers.Shorter(ProteinKmers.Shorter(ProteinKmers.Shorter(counts)));
            double total = singles.Sum(c => (double)c);
            return (p.Id, Composition: singles.Select(c => total > 0 ? c / total : 0).ToArray(), CvTree: CvTreeSignature.FromFourMers(counts));
        }).ToList())
            if (cvtree is not null) result[id] = (composition, cvtree);
    }
    return result;
}

// Both-strand 6-mer frequencies per record, as E2b: each one-strand count added to its canonical word (the smaller of a word and its
// reverse complement), then divided by the total.
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

// The k nearest of each member within its group, nearest first, itself excluded; ties by index. As E5.
static int[][] Neighbours((CorpusGenome Genome, string Habitat)[] members, Func<CorpusGenome, CorpusGenome, double> distance) =>
    members.Select((m, i) => Enumerable.Range(0, members.Length).Where(j => j != i)
        .OrderBy(j => distance(m.Genome, members[j].Genome)).ThenBy(j => j).Take(k).ToArray()).ToArray();

// Mean share of same-habitat neighbours; its mean over within-group shuffles of the labels; the one-sided p of the observed share. As E5.
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
