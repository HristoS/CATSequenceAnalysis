#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;
using System.Numerics;
using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequencePlatform.Core;

// P9 S9.10 (E1b): misplaced genomes on amino-acid CVTree neighbours. E1 (round 1) failed on log rho*; its GTDB key turned out to be
// 87 % wholesale renames (p9_e1b_check.cs). Fixed before the run (2026-11-12):
// - Corpus: round 2 (records <= RoundTwoMaxId, fingerprint verified), bacteria and archaea with amino-acid 4-mers, partial left out.
//   Viruses left out: their key keeps 202 minority positives in only 11 families.
// - Arms on the same genomes: hashed CVTree (4-mers, 3- and 2-mers derived, 8,192 dimensions, cosine), log rho* (Euclidean; the
//   round-1 arm, recomputed), GC + log10 length (standardized, Euclidean).
// - Score as E1: k = 10 nearest in the group; a genome's family (provider lineage, >= 5 genomes) median same-family share minus its own.
// - Keys: checkmPoor as E1 (contamination >= 5 % or completeness < 90 %). gtdbFamilyMinority: GTDB family (suffix stripped) differs
//   from the provider's, counted only in provider families where positives are under half of the genomes with the key; every genome of
//   other families has no key value. Context, not the rule: gtdbIsolated, the same with positives also under half of their genus.
// - Pass: for at least one of checkmPoor and gtdbFamilyMinority, CVTree's AUC bootstrap 95 % lower bound is above 0.6 and above the
//   GC + length AUC. Beats the round-1 arm (round-2 stop rule): CVTree's lower bound above log rho*'s AUC on that key.
// - If it fails, A1 (the anomaly agent) stays off for good.
// Writes <prefix>_Outliers.csv and <prefix>_Top.csv, each row with the corpus fingerprint; --publish writes the run to the live database.
const int k = 10, minFamily = 5, bootstraps = 1000, topCount = 200, dimensions = 8192, batch = 1000, seed = 20261112;
const double requiredAuc = 0.6;
var prefix = args.FirstOrDefault(a => !a.StartsWith("--", StringComparison.Ordinal)) ?? "P9_E1b";
bool publish = args.Contains("--publish");
var suffix = new Regex(@"_[A-Z]+$", RegexOptions.Compiled);
string? Strip(string? name) => name is null ? null : suffix.Replace(name, "");

await using var dataSource = ExperimentCorpus.OpenLive();
var verified = ExperimentCorpus.VerifyRoundTwoAsync(dataSource);
var corpus = await ExperimentCorpus.LoadAsync(dataSource, ExperimentCorpus.RoundTwoMaxId);
var withKmers = new HashSet<long>();
await using (var command = dataSource.CreateCommand("select sequence_id from protein_kmer where k = 4 and status = 'done' and sequence_id <= $1"))
{
    command.Parameters.AddWithValue(ExperimentCorpus.RoundTwoMaxId);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) withKmers.Add(reader.GetInt64(0));
}
var members = corpus.Where(g => !g.Excluded && g.Group == SignatureIndex.Cellular && withKmers.Contains(g.Id)).ToList();
Console.WriteLine($"{members.Count} bacteria and archaea with amino-acid 4-mers (of {corpus.Count(g => !g.Excluded && g.Group == SignatureIndex.Cellular)}).");

var cvtree = await CvTreeVectorsAsync(dataSource, members.Select(g => g.Id).ToArray());
var fingerprint = await verified;
Console.WriteLine($"CVTree vectors built; corpus fingerprint {fingerprint} verified.");

var family = members.Select(g => g.Provider.Family).ToArray();
var genus = members.Select(g => g.Provider.Genus).ToArray();
var spaces = new (string Name, int[][] Neighbours)[]
{
    ("cvtree", CosineNeighbours(cvtree)),
    ("signature", EuclideanNeighbours(members.Select(g => g.Signature).ToArray())),
    ("gcLength", EuclideanNeighbours(Standardize(members.Select(g => new[] { g.Gc, Math.Log10(g.LengthBp) }).ToArray()))),
};
var scores = spaces.ToDictionary(s => s.Name, s => Scores(s.Neighbours, family));
int scored = scores["cvtree"].Score.Count(s => s is not null);
Console.WriteLine($"{scored} genomes scored (provider family with >= {minFamily} genomes).");

var gtdbDiffers = members.Select((g, i) => scores["cvtree"].Score[i] is not null && g.Gtdb?.Family is { } gf ? Strip(gf) != g.Provider.Family : (bool?)null).ToArray();
var familyRate = Rates(family, gtdbDiffers);
var genusRate = Rates(genus, gtdbDiffers);
var keys = new (string Name, bool InRule, bool?[] Positive)[]
{
    ("checkmPoor", true, members.Select(g => Number(g.Field("checkm_contamination")) is { } c && Number(g.Field("checkm_completeness")) is { } m ? c >= 5 || m < 90 : (bool?)null).ToArray()),
    ("gtdbFamilyMinority", true, gtdbDiffers.Select((d, i) => d is not null && familyRate[family[i]!] < 0.5 ? d : null).ToArray()),
    ("gtdbIsolated", false, gtdbDiffers.Select((d, i) => d is not null && familyRate[family[i]!] < 0.5
        && (d == false || genus[i] is not { } gn || genusRate[gn] < 0.5) ? d : null).ToArray()),
};

var results = new List<ResultRow>();
bool pass = false, beatsRoundOne = false;
using (var csv = File.CreateText($"{prefix}_Outliers.csv"))
{
    csv.WriteLine("Key, InRule, Space, Scored, WithKey, Positives, BaseRate, Auc, AucLow, AucHigh, Top1PctRate, Enrichment, Pass, BeatsSignature, Fingerprint");
    foreach (var (keyName, inRule, positive) in keys)
    {
        var aucs = new Dictionary<string, ResultRow>();
        foreach (var (space, _) in spaces)
        {
            var score = scores[space].Score;
            var pairs = Enumerable.Range(0, members.Count).Where(i => score[i] is not null && positive[i] is not null)
                .Select(i => (Score: score[i]!.Value, Positive: positive[i]!.Value)).ToArray();
            int positives = pairs.Count(p => p.Positive);
            double baseRate = (double)positives / pairs.Length;
            var (auc, low, high) = BootstrapAuc(pairs, bootstraps, seed);
            var topPct = pairs.OrderByDescending(p => p.Score).Take(Math.Max(1, pairs.Length / 100)).ToArray();
            double topRate = (double)topPct.Count(p => p.Positive) / topPct.Length;
            aucs[space] = new ResultRow(keyName, inRule, space, scored, pairs.Length, positives, baseRate, auc, low, high, topRate, topRate / baseRate);
        }
        var cv = aucs["cvtree"];
        bool keyPass = cv.AucLow > requiredAuc && cv.AucLow > aucs["gcLength"].Auc;
        bool keyBeats = cv.AucLow > aucs["signature"].Auc;
        if (inRule) { pass |= keyPass; beatsRoundOne |= keyBeats; }
        foreach (var row in aucs.Values)
        {
            results.Add(row);
            var line = string.Create(CultureInfo.InvariantCulture,
                $"{row.Key}, {row.InRule}, {row.Space}, {row.Scored}, {row.WithKey}, {row.Positives}, {row.BaseRate:F4}, {row.Auc:F4}, {row.AucLow:F4}, {row.AucHigh:F4}, " +
                $"{row.Top1PctRate:F4}, {row.Enrichment:F2}, {(row.Space == "cvtree" ? keyPass.ToString() : "")}, {(row.Space == "cvtree" ? keyBeats.ToString() : "")}, {fingerprint}");
            csv.WriteLine(line);
            Console.WriteLine(line);
        }
    }
}
Console.WriteLine($"Verdict: {(pass ? "pass" : "fail")}; CVTree beats the round-1 arm (log rho*) on a rule key: {beatsRoundOne}.");

// The top of CVTree's list with its evidence, and (after the run, not the rule) whom GTDB sides with when CVTree flags a genome.
var cvScores = scores["cvtree"];
var cvNeighbours = spaces[0].Neighbours;
var top = Enumerable.Range(0, members.Count).Where(i => cvScores.Score[i] is not null).OrderByDescending(i => cvScores.Score[i]).ThenBy(i => members[i].Id)
    .Take(topCount).Select(i =>
    {
        var majority = cvNeighbours[i].Select(j => family[j]).Where(f => f is not null).GroupBy(f => f!)
            .OrderByDescending(f => f.Count()).ThenBy(f => f.Key, StringComparer.Ordinal).FirstOrDefault();
        return new TopRow(members[i], cvScores.Score[i]!.Value, cvScores.Share[i]!.Value, cvScores.Median[i]!.Value, majority?.Key, majority?.Count() ?? 0,
            cvNeighbours[i].Select(j => members[j].Id).ToArray(), Strip(members[i].Gtdb?.Family), keys[0].Positive[i], keys[1].Positive[i]);
    }).ToList();
var withGtdb = top.Where(t => t.GtdbFamily is not null).ToList();
var diagnostic = new
{
    topWithGtdb = withGtdb.Count,
    gtdbAgreesWithProvider = withGtdb.Count(t => t.GtdbFamily == t.Genome.Provider.Family),
    gtdbAgreesWithNeighbours = withGtdb.Count(t => t.GtdbFamily == t.Majority),
    topGcBelow035 = top.Count(t => t.Genome.Gc < 0.35),
    topGcAbove065 = top.Count(t => t.Genome.Gc > 0.65),
};
Console.WriteLine($"After the run: of {diagnostic.topWithGtdb} top genomes with a GTDB family, GTDB agrees with the provider in " +
    $"{diagnostic.gtdbAgreesWithProvider} and with the neighbours' majority in {diagnostic.gtdbAgreesWithNeighbours}; GC < 0.35: {diagnostic.topGcBelow035}, > 0.65: {diagnostic.topGcAbove065}.");
using (var csv = File.CreateText($"{prefix}_Top.csv"))
{
    csv.WriteLine("Id, Name, ProviderFamily, Score, Share, FamilyMedianShare, NeighbourMajorityFamily, MajorityCount, GtdbFamily, CheckmPoor, GtdbFamilyMinority, Gc, LengthBp, Fingerprint");
    foreach (var t in top)
        csv.WriteLine(string.Create(CultureInfo.InvariantCulture,
            $"{t.Genome.Id}, {Quote(t.Genome.Name)}, {t.Genome.Provider.Family}, {t.Score:F2}, {t.Share:F2}, {t.FamilyMedian:F2}, {t.Majority}, {t.MajorityCount}, " +
            $"{t.GtdbFamily}, {t.CheckmPoor}, {t.GtdbMinority}, {t.Genome.Gc:F3}, {t.Genome.LengthBp}, {fingerprint}"));
}

if (publish)
{
    var json = new JsonSerializerOptions { PropertyNamingPolicy = JsonNamingPolicy.CamelCase };
    var gtdbRow = results.First(r => r.Key == "gtdbFamilyMinority" && r.Space == "cvtree");
    var summary = JsonSerializer.Serialize(new
    {
        rule = $"Pass when, for CheckM poor or GTDB family differs (minority families), CVTree's AUC 95 % lower bound is above {requiredAuc} and above GC + length.",
        analysed = members.Count,
        analysedNote = $"bacteria and archaea with amino-acid 4-mers; the {topCount} highest CVTree scores are listed",
        keyNumbers = new object[]
        {
            new { label = "AUC, CVTree, GTDB family differs (minority families)", value = gtdbRow.Auc },
            new { label = "its 95 % lower bound", value = gtdbRow.AucLow },
            new { label = "enrichment in the top 1 %", value = gtdbRow.Enrichment },
        },
        results,
        afterTheRun = diagnostic,
        fingerprint,
    }, json);
    var findings = top.Select(t => new ExperimentFinding(t.Genome.Id, "outlier", t.Score, JsonSerializer.Serialize(new
    {
        family = t.Genome.Provider.Family, share = t.Share, familyMedianShare = t.FamilyMedian, majorityFamily = t.Majority, majorityCount = t.MajorityCount,
        neighbours = t.NeighbourIds, gtdbFamily = t.GtdbFamily,
    }, json))).ToList();
    var (stored, skipped) = await ExperimentStore.PublishAsync(dataSource, "e1b-outliers-cvtree", "Misplaced genomes (amino-acid CVTree)",
        "Do genomes whose amino-acid CVTree neighbours belong to another family point to mislabels or contamination more than GC and length do?",
        pass ? ExperimentStore.Pass : ExperimentStore.Fail, $"live, records <= {ExperimentCorpus.RoundTwoMaxId}, fingerprint {fingerprint}",
        "scripts/platform-checks/p9_e1b_cvtree.cs", summary, findings);
    Console.WriteLine($"Published e1b-outliers-cvtree: {stored} findings stored, {skipped} skipped.");
}

static double? Number(string? text) => double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out var value) ? value : null;
static string Quote(string text) => "\"" + text.Replace("\"", "\"\"") + "\"";

// Share of key positives per taxon among the genomes with a key value.
static Dictionary<string, double> Rates(string?[] taxon, bool?[] positive) =>
    Enumerable.Range(0, taxon.Length).Where(i => taxon[i] is not null && positive[i] is not null).GroupBy(i => taxon[i]!)
        .ToDictionary(g => g.Key, g => g.Count(i => positive[i] == true) / (double)g.Count());

// Hashed CVTree vectors (unit length) in the order of ids, read in batches so the 4-mer payloads (about 1.6 GB) never sit in memory at once.
static async Task<float[][]> CvTreeVectorsAsync(NpgsqlDataSource dataSource, long[] ids)
{
    var vectors = new Dictionary<long, float[]>();
    foreach (var chunk in ids.Chunk(batch))
    {
        var payloads = new List<(long Id, byte[] Payload)>();
        await using (var command = dataSource.CreateCommand("select sequence_id, payload from protein_kmer where k = 4 and status = 'done' and sequence_id = any($1)"))
        {
            command.Parameters.AddWithValue(chunk);
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) payloads.Add((reader.GetInt64(0), reader.GetFieldValue<byte[]>(1)));
        }
        foreach (var (id, vector) in payloads.AsParallel().Select(p => (p.Id, Normalize(Hash(Composition(ProteinKmers.Decode(p.Payload, 4)))))).ToList())
            vectors[id] = vector;
    }
    return ids.Select(id => vectors[id]).ToArray();
}

// Hao's composition vector of the 4-mers, with the 3-mers and 2-mers derived from the 4-mer counts (ProteinKmers.Shorter).
static float[] Composition(int[] words)
{
    var shorterCounts = ProteinKmers.Shorter(words);
    var shortestCounts = ProteinKmers.Shorter(shorterCounts);
    double total = words.Sum(w => (double)w), shorterTotal = shorterCounts.Sum(w => (double)w), shortestTotal = shortestCounts.Sum(w => (double)w);
    var vector = new float[words.Length];
    int suffixSpace = shorterCounts.Length, middleSpace = shortestCounts.Length;
    for (int w = 0; w < words.Length; w++)
    {
        int prefix = w / 20, suffix = w % suffixSpace, middle = prefix % middleSpace;
        double shortest = shortestCounts[middle] / shortestTotal;
        double expected = shortest > 0 ? (shorterCounts[prefix] / shorterTotal) * (shorterCounts[suffix] / shorterTotal) / shortest : 0;
        vector[w] = expected > 0 ? (float)((words[w] / total - expected) / expected) : 0;
    }
    return vector;
}

// Feature hashing exactly as cvtree_compression_check.cs (commit f5cb655), so its 8,192-dimension result applies.
static float[] Hash(float[] vector)
{
    var hashed = new float[dimensions];
    for (int w = 0; w < vector.Length; w++)
    {
        ulong h = (ulong)w * 0x9E3779B97F4A7C15UL;
        h ^= h >> 29;
        h *= 0xBF58476D1CE4E5B9UL;
        h ^= h >> 32;
        hashed[(int)(h % dimensions)] += (h & (1UL << 63)) == 0 ? vector[w] : -vector[w];
    }
    return hashed;
}

static float[] Normalize(float[] vector)
{
    double norm = Math.Sqrt(vector.Sum(x => (double)x * x));
    if (norm > 0) for (int i = 0; i < vector.Length; i++) vector[i] = (float)(vector[i] / norm);
    return vector;
}

// The k nearest by cosine (unit vectors: largest dot product), nearest first, itself excluded; ties by index.
static int[][] CosineNeighbours(float[][] vectors)
{
    int n = vectors.Length;
    var result = new int[n][];
    Parallel.For(0, n, i =>
    {
        var best = new int[k];
        var bestDot = new float[k];
        Array.Fill(bestDot, float.NegativeInfinity);
        var a = vectors[i];
        for (int j = 0; j < n; j++)
        {
            if (j == i) continue;
            float dot = Dot(a, vectors[j]);
            if (dot <= bestDot[k - 1]) continue;
            int slot = k - 1;
            while (slot > 0 && bestDot[slot - 1] < dot) { bestDot[slot] = bestDot[slot - 1]; best[slot] = best[slot - 1]; slot--; }
            bestDot[slot] = dot;
            best[slot] = j;
        }
        result[i] = best;
    });
    return result;
}

static float Dot(float[] a, float[] b)
{
    var sum = Vector<float>.Zero;
    int i = 0;
    for (; i <= a.Length - Vector<float>.Count; i += Vector<float>.Count)
        sum += new Vector<float>(a, i) * new Vector<float>(b, i);
    float total = Vector.Sum(sum);
    for (; i < a.Length; i++) total += a[i] * b[i];
    return total;
}

// The k nearest of every point (Euclidean), nearest first, itself excluded; as E1.
static int[][] EuclideanNeighbours(double[][] points)
{
    int n = points.Length, d = points[0].Length;
    var flat = new double[n * d];
    for (int i = 0; i < n; i++) Array.Copy(points[i], 0, flat, i * d, d);
    var result = new int[n][];
    Parallel.For(0, n, i =>
    {
        var bestDistance = new double[k];
        var bestIndex = new int[k];
        Array.Fill(bestDistance, double.PositiveInfinity);
        for (int j = 0; j < n; j++)
        {
            if (j == i) continue;
            double sum = 0;
            for (int f = 0; f < d; f++) { double diff = flat[i * d + f] - flat[j * d + f]; sum += diff * diff; }
            if (sum >= bestDistance[k - 1]) continue;
            int at = k - 1;
            while (at > 0 && bestDistance[at - 1] > sum) { bestDistance[at] = bestDistance[at - 1]; bestIndex[at] = bestIndex[at - 1]; at--; }
            bestDistance[at] = sum;
            bestIndex[at] = j;
        }
        result[i] = bestIndex;
    });
    return result;
}

// Per genome: its family's median same-family share minus its own (null outside families with >= minFamily genomes), its share, the median. As E1.
static (double?[] Score, double?[] Share, double?[] Median) Scores(int[][] neighbours, string?[] taxon)
{
    int n = taxon.Length;
    var size = taxon.Where(t => t is not null).GroupBy(t => t!).ToDictionary(g => g.Key, g => g.Count());
    var share = new double?[n];
    for (int i = 0; i < n; i++)
        if (taxon[i] is { } t && size[t] >= minFamily) share[i] = neighbours[i].Count(j => taxon[j] == t) / (double)k;
    var median = Enumerable.Range(0, n).Where(i => share[i] is not null).GroupBy(i => taxon[i]!)
        .ToDictionary(g => g.Key, g => { var s = g.Select(i => share[i]!.Value).Order().ToArray(); return s[s.Length / 2]; });
    var score = new double?[n];
    var medianOf = new double?[n];
    for (int i = 0; i < n; i++)
        if (share[i] is { } s) { medianOf[i] = median[taxon[i]!]; score[i] = medianOf[i] - s; }
    return (score, share, medianOf);
}

// AUC (ties count half) with a percentile bootstrap over genomes. As E1.
static (double Auc, double Low, double High) BootstrapAuc((double Score, bool Positive)[] pairs, int rounds, int seed)
{
    double auc = Auc(pairs);
    var rng = new Random(seed);
    var samples = new double[rounds];
    var draw = new (double, bool)[pairs.Length];
    for (int r = 0; r < rounds; r++)
    {
        for (int i = 0; i < draw.Length; i++) draw[i] = pairs[rng.Next(pairs.Length)];
        samples[r] = Auc(draw);
    }
    Array.Sort(samples);
    return (auc, samples[(int)(0.025 * rounds)], samples[(int)(0.975 * rounds) - 1]);
}

static double Auc((double Score, bool Positive)[] pairs)
{
    var sorted = pairs.OrderBy(p => p.Score).ToArray();
    double rankSum = 0;
    long positives = 0;
    for (int i = 0; i < sorted.Length;)
    {
        int j = i;
        while (j < sorted.Length && sorted[j].Score == sorted[i].Score) j++;
        double rank = (i + j + 1) / 2.0;
        for (int m = i; m < j; m++) if (sorted[m].Positive) { rankSum += rank; positives++; }
        i = j;
    }
    long negatives = sorted.Length - positives;
    return positives == 0 || negatives == 0 ? double.NaN : (rankSum - positives * (positives + 1) / 2.0) / (positives * (double)negatives);
}

static double[][] Standardize(double[][] rows)
{
    int width = rows[0].Length;
    var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
    var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(rows.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
    return rows.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray();
}

record ResultRow(string Key, bool InRule, string Space, int Scored, int WithKey, int Positives, double BaseRate, double Auc, double AucLow, double AucHigh,
    double Top1PctRate, double Enrichment);

record TopRow(CorpusGenome Genome, double Score, double Share, double FamilyMedian, string? Majority, int MajorityCount, long[] NeighbourIds,
    string? GtdbFamily, bool? CheckmPoor, bool? GtdbMinority);
