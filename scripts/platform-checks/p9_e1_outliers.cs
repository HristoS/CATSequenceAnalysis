#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;
using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequencePlatform.Core;

// P9 S9.1 (E1): misplaced genomes. A genome whose nearest neighbours by signature don't share its family, although its family's
// genomes usually do, is a candidate mislabel, contamination or odd composition. Fixed before the run:
// - Score: within its group, the median share of same-family genomes among the 10 nearest over the genome's family (provider lineage,
//   families with >= 5 genomes), minus the genome's own share. Genus the same way, reported alongside.
// - Keys, never used by the score: cellular, GTDB family differs from the provider's (GTDB's _A-style suffixes stripped); cellular,
//   CheckM contamination >= 5 % or completeness < 90 %; viral, NCBI's current family differs from the provider's (reclassified since).
// - Pass: for at least one key, the AUC's bootstrap 95 % lower bound is above 0.6 and above the AUC of the same score computed in
//   GC + log10 length space.
// Eukaryotes are left out: one genome per family, so there is no family share to fall short of.
// Writes <prefix>_Outliers.csv (per group, key and space) and <prefix>_Top.csv; --publish writes the run into the live database.
const int k = 10, minFamily = 5, bootstraps = 1000, topPerGroup = 200;
const double requiredAuc = 0.6;
var prefix = args.FirstOrDefault(a => !a.StartsWith("--", StringComparison.Ordinal)) ?? "P9_E1";
bool publish = args.Contains("--publish");
var suffix = new Regex(@"_[A-Z]+$", RegexOptions.Compiled);
string? Strip(string? name) => name is null ? null : suffix.Replace(name, "");

await using var dataSource = ExperimentCorpus.OpenSnapshot();
var corpus = (await ExperimentCorpus.LoadAsync(dataSource)).Where(g => !g.Excluded && g.Group != SignatureIndex.Eukaryotic).ToList();
Console.WriteLine($"{corpus.Count} genomes (eukaryotes and partial records left out).");

var results = new List<ResultRow>();
var top = new List<TopRow>();
using var csv = File.CreateText($"{prefix}_Outliers.csv");
csv.WriteLine("Group, Key, Space, Scored, WithKey, Positives, BaseRate, Auc, AucLow, AucHigh, Top1PctRate, Enrichment, Pass");
bool anyPass = false;

foreach (var group in new[] { SignatureIndex.Cellular, SignatureIndex.Viral })
{
    var members = corpus.Where(g => g.Group == group).ToList();
    var signature = Neighbours(members.Select(g => g.Signature).ToArray());
    var gcLength = Neighbours(Standardize(members.Select(g => new[] { g.Gc, Math.Log10(g.LengthBp) }).ToArray()));
    var family = members.Select(g => g.Provider.Family).ToArray();
    var genus = members.Select(g => g.Provider.Genus).ToArray();
    var (familyScore, familyShare, familyMedian) = Scores(signature, family);
    var (genusScore, _, _) = Scores(signature, genus);
    var (gcScore, _, _) = Scores(gcLength, family);
    int scored = familyScore.Count(s => s is not null);
    Console.WriteLine($"{group}: {members.Count} genomes, {scored} scored at family.");

    var keys = new List<(string Name, bool?[] Positive)>();
    if (group == SignatureIndex.Cellular)
    {
        keys.Add(("gtdbFamilyDiffers", members.Select(g => g.Gtdb?.Family is { } gf && g.Provider.Family is { } pf ? Strip(gf) != pf : (bool?)null).ToArray()));
        keys.Add(("checkmPoor", members.Select(g => Number(g.Field("checkm_contamination")) is { } c && Number(g.Field("checkm_completeness")) is { } m ? c >= 5 || m < 90 : (bool?)null).ToArray()));
    }
    else
        keys.Add(("ncbiFamilyDiffers", members.Select(g => g.Ncbi?.Family is { } nf && g.Provider.Family is { } pf ? nf != pf : (bool?)null).ToArray()));

    foreach (var (keyName, positive) in keys)
    {
        var aucs = new Dictionary<string, (double Auc, double Low, double High)>();
        foreach (var (space, score) in new[] { ("signature", familyScore), ("gcLength", gcScore) })
        {
            var pairs = Enumerable.Range(0, members.Count).Where(i => score[i] is not null && positive[i] is not null)
                .Select(i => (Score: score[i]!.Value, Positive: positive[i]!.Value)).ToArray();
            int positives = pairs.Count(p => p.Positive);
            double baseRate = (double)positives / pairs.Length;
            var (auc, low, high) = BootstrapAuc(pairs, bootstraps, seed: 20261103);
            var topPct = pairs.OrderByDescending(p => p.Score).Take(Math.Max(1, pairs.Length / 100)).ToArray();
            double topRate = (double)topPct.Count(p => p.Positive) / topPct.Length;
            aucs[space] = (auc, low, high);
            var line = string.Create(CultureInfo.InvariantCulture,
                $"{group}, {keyName}, {space}, {scored}, {pairs.Length}, {positives}, {baseRate:F4}, {auc:F4}, {low:F4}, {high:F4}, {topRate:F4}, {topRate / baseRate:F2}");
            results.Add(new ResultRow(group, keyName, space, scored, pairs.Length, positives, baseRate, auc, low, high, topRate, topRate / baseRate));
            if (space == "gcLength")
            {
                bool keyPass = aucs["signature"].Low > requiredAuc && aucs["signature"].Low > aucs["gcLength"].Auc;
                anyPass |= keyPass;
                csv.WriteLine(line + $", {keyPass}");
                Console.WriteLine(line + $", pass {keyPass}");
            }
            else
            {
                csv.WriteLine(line + ", ");
                Console.WriteLine(line);
            }
        }
    }

    // The top of the list with its evidence: what the neighbours are, and what the keys say.
    foreach (var i in Enumerable.Range(0, members.Count).Where(i => familyScore[i] is not null).OrderByDescending(i => familyScore[i]).Take(topPerGroup))
    {
        var g = members[i];
        var neighbourFamilies = signature[i].Select(j => family[j]).ToArray();
        var majority = neighbourFamilies.Where(f => f is not null).GroupBy(f => f).OrderByDescending(f => f.Count()).ThenBy(f => f.Key, StringComparer.Ordinal).FirstOrDefault();
        top.Add(new TopRow(group, g, familyScore[i]!.Value, familyShare[i]!.Value, familyMedian[i]!.Value, genusScore[i], majority?.Key, majority?.Count() ?? 0,
            signature[i].Select(j => members[j].Id).ToArray(), signature[i].Select(j => family[j]).ToArray(), Strip(g.Gtdb?.Family), g.Ncbi?.Family,
            Number(g.Field("checkm_completeness")), Number(g.Field("checkm_contamination"))));
    }
}
Console.WriteLine($"Verdict: {(anyPass ? "pass" : "fail")}.");

// Added after the run (not part of the rule): when the score flags a cellular genome, whom does GTDB side with, the provider or the
// neighbours? Names are compared as they stand, so a nomenclature difference counts as neither.
var flagged = top.Where(t => t.Group == SignatureIndex.Cellular && t.GtdbFamily is not null).ToList();
var diagnostic = new
{
    cellularTopWithGtdb = flagged.Count,
    gtdbAgreesWithProvider = flagged.Count(t => t.GtdbFamily == t.Genome.Provider.Family),
    gtdbAgreesWithNeighbours = flagged.Count(t => t.GtdbFamily == t.Majority),
    cellularTopGcBelow035 = top.Count(t => t.Group == SignatureIndex.Cellular && t.Genome.Gc < 0.35),
    cellularTopGcAbove065 = top.Count(t => t.Group == SignatureIndex.Cellular && t.Genome.Gc > 0.65),
};
Console.WriteLine($"After the run: of {diagnostic.cellularTopWithGtdb} top cellular genomes with a GTDB family, GTDB agrees with the provider in " +
    $"{diagnostic.gtdbAgreesWithProvider} and with the neighbours' majority in {diagnostic.gtdbAgreesWithNeighbours}.");

using (var topCsv = File.CreateText($"{prefix}_Top.csv"))
{
    topCsv.WriteLine("Group, Id, Name, ProviderFamily, Score, Share, FamilyMedianShare, GenusScore, NeighbourMajorityFamily, MajorityCount, GtdbFamily, NcbiFamily, CheckmCompleteness, CheckmContamination, Gc, LengthBp");
    foreach (var t in top)
        topCsv.WriteLine(string.Create(CultureInfo.InvariantCulture,
            $"{t.Group}, {t.Genome.Id}, {Quote(t.Genome.Name)}, {t.Genome.Provider.Family}, {t.Score:F2}, {t.Share:F2}, {t.FamilyMedian:F2}, {t.GenusScore:F2}, {t.Majority}, {t.MajorityCount}, " +
            $"{t.GtdbFamily}, {t.NcbiFamily}, {t.Completeness}, {t.Contamination}, {t.Genome.Gc:F3}, {t.Genome.LengthBp}"));
}

if (publish)
{
    var json = new JsonSerializerOptions { PropertyNamingPolicy = JsonNamingPolicy.CamelCase };
    var cellularGtdb = results.First(r => r.Group == SignatureIndex.Cellular && r.Key == "gtdbFamilyDiffers" && r.Space == "signature");
    var summary = JsonSerializer.Serialize(new
    {
        rule = $"Pass when, for at least one key, the AUC's 95 % lower bound is above {requiredAuc} and above the same score in GC + length space.",
        analysed = corpus.Count,
        analysedNote = $"every bacterium, archaeon and virus; the {topPerGroup} highest scores per group are listed",
        keyNumbers = new object[]
        {
            new { label = "AUC, cellular, GTDB family differs", value = cellularGtdb.Auc },
            new { label = "its 95 % lower bound", value = cellularGtdb.AucLow },
            new { label = "enrichment in the top 1 %", value = cellularGtdb.Enrichment },
        },
        results,
        afterTheRun = diagnostic,
        reading = "The flagged genomes are mostly composition outliers (endosymbionts, extreme GC, small viral genomes), not mislabels: " +
            "where GTDB has a family for them, it agrees with the provider, not with the neighbours.",
    }, json);
    var findings = top.Select(t => new ExperimentFinding(t.Genome.Id, "outlier", t.Score, JsonSerializer.Serialize(new
    {
        group = t.Group, family = t.Genome.Provider.Family, share = t.Share, familyMedianShare = t.FamilyMedian, genusScore = t.GenusScore,
        majorityFamily = t.Majority, majorityCount = t.MajorityCount, neighbours = t.NeighbourIds, neighbourFamilies = t.NeighbourFamilies,
        gtdbFamily = t.GtdbFamily, ncbiFamily = t.NcbiFamily, checkmCompleteness = t.Completeness, checkmContamination = t.Contamination,
    }, json))).ToList();
    await using var live = NpgsqlDataSource.Create(ConnectionString.Resolve());
    var (stored, skipped) = await ExperimentStore.PublishAsync(live, "e1-outliers", "Misplaced genomes",
        "Which genomes sit among neighbours of another family, and does that point to mislabels or contamination more than GC and length do?",
        anyPass ? ExperimentStore.Pass : ExperimentStore.Fail, ExperimentCorpus.SnapshotDatabase, "scripts/platform-checks/p9_e1_outliers.cs", summary, findings);
    Console.WriteLine($"Published e1-outliers: {stored} findings stored, {skipped} skipped.");
}

static double? Number(string? text) => double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out var value) ? value : null;
static string Quote(string text) => "\"" + text.Replace("\"", "\"\"") + "\"";

// The k nearest of every point (Euclidean), nearest first, itself excluded.
static int[][] Neighbours(double[][] points)
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

// Per genome: its taxon's median same-taxon share minus its own (null outside taxa with >= minFamily genomes), its own share, the median.
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

// AUC (ties count half) with a percentile bootstrap over genomes.
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
        double rank = (i + j + 1) / 2.0;   // mean of 1-based ranks i+1..j
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

record ResultRow(string Group, string Key, string Space, int Scored, int WithKey, int Positives, double BaseRate, double Auc, double AucLow, double AucHigh,
    double Top1PctRate, double Enrichment);

record TopRow(string Group, CorpusGenome Genome, double Score, double Share, double FamilyMedian, double? GenusScore, string? Majority, int MajorityCount,
    long[] NeighbourIds, string?[] NeighbourFamilies, string? GtdbFamily, string? NcbiFamily, double? Completeness, double? Contamination);
