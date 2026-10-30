#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S5.1 (P5, R2): candidate tree distances on gc_stratum_check.cs's sample (same rule and seed: 437 GTDB-labelled genomes, one per
// species, up to 4 species in each of 150 genera). Every level is downsampled from the stored 500 bp profile. Candidates:
//  - compare:     Compare at the level, as stored (the baseline, d = 1 - s).
//  - strandOrigin: the best Compare over both genomes' versions: forward and flipped, each as stored and, for circular chromosomes
//                 OriginRotation accepts, rotated to the origin (as find-similar scores, but every version against every version).
//  - coverage:    Compare times shorter / longer window count, so the longer genome's uncovered windows count as unmatched.
//  - summary:     Euclidean distance of standardized per-genome features at the level: per axis |mean| (x and y change sign with
//                 the strand), SD and lag-1 autocorrelation across windows.
// Baselines: gc (-|ΔGC| from dinucleotide counts), lengthRatio (shorter / longer, the length term coverage adds on its own) and
// gcLength (Euclidean distance of standardized GC and log10 length: R4's GC-plus-length baseline).
// AUC per rank and |ΔGC| stratum with 95% bootstrap CIs over genomes (200 resamples), paired against gc, gcLength and compare;
// triangle violations over all triples per distance. Args: output path (default S51_Candidates.csv), coarse factor over 500 bp
// (default 64 = 32 kbp). The triangle counts go to the same path with a _Triangle suffix.
const int maxGenera = 150, maxSpeciesPerGenus = 4, resamples = 200, minLength = 256_000;
var output = args.Length > 0 ? args[0] : "S51_Candidates.csv";
int factor = args.Length > 1 ? int.Parse(args[1]) : 64;
int windowSize = RecordStore.FineWindowSize * factor;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var all = new List<Genome>();
await using (var cmd = dataSource.CreateCommand(
    "select distinct on (r.source_id) r.id, r.length_bp, r.dinucleotide_counts, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species, " +
    "r.annotations->'contig'->>'topology' " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null " +
    "order by r.source_id, r.length_bp desc, r.id"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        if (r.GetInt64(1) < minLength || r.IsDBNull(6) || r.IsDBNull(5) || r.IsDBNull(4) || r.IsDBNull(3)) continue;
        var counts = r.GetFieldValue<long[]>(2);
        double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
        all.Add(new Genome(r.GetInt64(0), r.GetInt64(1), gc, r.GetString(3), r.GetString(4), r.GetString(5), r.GetString(6),
            !r.IsDBNull(7) && r.GetString(7) == "circular"));
    }

var rng = new Random(20261029);
var perSpecies = all.GroupBy(g => g.Species).Select(s => s.OrderBy(_ => rng.Next()).First()).ToList();
var sample = perSpecies.GroupBy(g => g.Genus).Where(s => s.Count() >= 2).OrderBy(_ => rng.Next()).Take(maxGenera)
    .SelectMany(s => s.OrderBy(_ => rng.Next()).Take(maxSpeciesPerGenus)).ToList();

var fine = await ProfileLevels.ReadManyAsync(dataSource, sample.Select(g => g.Id).ToArray(),
    RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
sample = sample.Where(g => fine.ContainsKey(g.Id)).ToList();
int n = sample.Count;

// Each genome's versions at the level: [0] as stored, then flipped, then (when rotated) rotated and rotated-flipped.
var versions = sample.Select(g =>
{
    var stored = fine[g.Id];
    var list = new List<ZCurveWindowedProfile> { stored, Flip(stored) };
    if (g.Circular && OriginRotation.ToOrigin(stored, g.Length) is { } turned)
    {
        var reverse = Flip(turned);
        list.Add(turned);
        list.Add(reverse.RotateTo(reverse.GcSkewExtremes().OriginWindow));
    }
    return list.Select(v => factor == 1 ? v : v.Downsample(factor)).ToArray();
}).ToArray();
Console.Error.WriteLine($"sample {n} genomes in {sample.Select(g => g.Genus).Distinct().Count()} genera, " +
    $"{sample.Select(g => g.Family).Distinct().Count()} families, {sample.Select(g => g.Phylum).Distinct().Count()} phyla; " +
    $"{versions.Count(v => v.Length > 2)} rotated; level {windowSize} bp");

var features = Standardize(versions.Select(v => Features(v[0])).ToArray());
var gcLength = Standardize(sample.Select(g => new[] { g.Gc, Math.Log10(g.Length) }).ToArray());

string[] candidates = { "compare", "strandOrigin", "coverage", "summary", "gc", "gcLength", "lengthRatio" };
var score = candidates.Select(_ => new double[n, n]).ToArray();
var watch = System.Diagnostics.Stopwatch.StartNew();
Parallel.For(0, n, a =>
{
    for (int b = a + 1; b < n; b++)
    {
        double plain = ZCurveWindowedProfile.Compare(versions[a][0], versions[b][0]);
        double best = plain;
        foreach (var va in versions[a])
            foreach (var vb in versions[b])
                best = Math.Max(best, ZCurveWindowedProfile.Compare(va, vb));
        int wa = versions[a][0].WindowCount, wb = versions[b][0].WindowCount;
        double ratio = (double)Math.Min(wa, wb) / Math.Max(wa, wb);
        double[] values = { plain, best, plain * ratio, -Euclidean(features[a], features[b]), -Math.Abs(sample[a].Gc - sample[b].Gc),
            -Euclidean(gcLength[a], gcLength[b]), ratio };
        for (int c = 0; c < candidates.Length; c++) score[c][a, b] = score[c][b, a] = values[c];
    }
});
Console.Error.WriteLine($"{n * (n - 1) / 2:N0} pairs scored in {watch.Elapsed.TotalSeconds:F1} s");

// Triangle violations: d = 1 - s for the Compare-based candidates, d = -s otherwise (Euclidean, |ΔGC|; lengthRatio isn't a distance,
// and gcLength is Euclidean like summary).
using (var triangle = File.CreateText(Path.ChangeExtension(output, null) + "_Triangle.csv"))
{
    triangle.WriteLine("WindowSize, Candidate, Genomes, InequalitiesChecked, InequalitiesViolated, MaxViolation");
    for (int c = 0; c < 6; c++)
    {
        var s = score[c];
        Func<int, int, double> d = c < 3 ? (a, b) => 1 - s[a, b] : (a, b) => -s[a, b];
        long violated = 0;
        double maxViolation = 0;
        var gate = new object();
        Parallel.For(0, n, a =>
        {
            long localViolated = 0;
            double localMax = 0;
            for (int b = a + 1; b < n; b++)
                for (int e = b + 1; e < n; e++)
                {
                    double ab = d(a, b), be = d(b, e), ae = d(a, e);
                    double excess = Math.Max(ab - be - ae, Math.Max(be - ab - ae, ae - ab - be));   // at most one side can exceed
                    if (excess > 1e-9) { localViolated++; localMax = Math.Max(localMax, excess); }
                }
            lock (gate) { violated += localViolated; maxViolation = Math.Max(maxViolation, localMax); }
        });
        long checkedCount = 3L * n * (n - 1) * (n - 2) / 6;
        string line = $"{windowSize}, {candidates[c]}, {n}, {checkedCount}, {violated}, {maxViolation:G4}";
        triangle.WriteLine(line);
        Console.Error.WriteLine(line);
    }
}

var ranks = new (string Name, Func<Genome, string> Of)[] { ("genus", g => g.Genus), ("family", g => g.Family), ("phylum", g => g.Phylum) };
var strata = new (string Name, double Min, double Max)[] { ("all", 0, 1), ("dGC<0.02", 0, 0.02), ("dGC 0.02-0.10", 0.02, 0.10), ("dGC>=0.10", 0.10, 1) };

// AUCs of one draw of genomes (a genome drawn twice is not paired with itself): [rank, stratum, candidate].
(double Auc, int Positives, int Negatives)[,,] Aucs(int[] drawn)
{
    var result = new (double, int, int)[ranks.Length, strata.Length, candidates.Length];
    for (int k = 0; k < ranks.Length; k++)
        for (int st = 0; st < strata.Length; st++)
        {
            var pairs = new List<(int A, int B, bool Positive)>();
            for (int i = 0; i < drawn.Length; i++)
                for (int j = i + 1; j < drawn.Length; j++)
                {
                    int a = drawn[i], b = drawn[j];
                    if (a == b) continue;
                    double dGc = Math.Abs(sample[a].Gc - sample[b].Gc);
                    if (dGc < strata[st].Min || dGc >= strata[st].Max) continue;
                    pairs.Add((a, b, ranks[k].Of(sample[a]) == ranks[k].Of(sample[b])));
                }
            for (int c = 0; c < candidates.Length; c++)
                result[k, st, c] = Auc(pairs.Select(p => (score[c][p.A, p.B], p.Positive)).ToList());
        }
    return result;
}

var observed = Aucs(Enumerable.Range(0, n).ToArray());
var draws = new (double Auc, int Positives, int Negatives)[resamples][,,];
Parallel.For(0, resamples, r =>
{
    var local = new Random(20261029 + r);
    draws[r] = Aucs(Enumerable.Range(0, n).Select(_ => local.Next(n)).ToArray());
});

string Interval(IEnumerable<double> values)
{
    var sorted = values.Where(v => !double.IsNaN(v)).Order().ToArray();
    if (sorted.Length == 0) return "NaN, NaN";
    double At(double p) => sorted[(int)Math.Round(p * (sorted.Length - 1))];
    return $"{At(0.025):F4}, {At(0.975):F4}";
}

var references = new[] { "gc", "gcLength", "compare" }.Select(name => Array.IndexOf(candidates, name)).ToArray();
using var report = File.CreateText(output);
report.WriteLine("WindowSize, Candidate, Rank, Stratum, Positives, Negatives, Auc, CiLow, CiHigh, " +
    "MinusGc, MinusGcCiLow, MinusGcCiHigh, MinusGcLength, MinusGcLengthCiLow, MinusGcLengthCiHigh, " +
    "MinusCompare, MinusCompareCiLow, MinusCompareCiHigh");
for (int k = 0; k < ranks.Length; k++)
    for (int st = 0; st < strata.Length; st++)
        for (int c = 0; c < candidates.Length; c++)
        {
            var (auc, positives, negatives) = observed[k, st, c];
            string line = $"{windowSize}, {candidates[c]}, {ranks[k].Name}, {strata[st].Name}, {positives}, {negatives}, {auc:F4}, " +
                $"{Interval(draws.Select(d => d[k, st, c].Auc))}, " + string.Join(", ", references.Select(x =>
                    $"{auc - observed[k, st, x].Auc:F4}, {Interval(draws.Select(d => d[k, st, c].Auc - d[k, st, x].Auc))}"));
            report.WriteLine(line);
            Console.WriteLine(line);
        }

static ZCurveWindowedProfile Flip(ZCurveWindowedProfile profile) =>
    ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize, profile.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

static double Euclidean(double[] p, double[] q) => Math.Sqrt(p.Zip(q, (u, v) => (u - v) * (u - v)).Sum());

static double[] Features(ZCurveWindowedProfile level)
{
    var result = new List<double>();
    foreach (var (values, signed) in new[] { (level.Windows.Select(w => w.X).ToArray(), false),
        (level.Windows.Select(w => w.Y).ToArray(), false), (level.Windows.Select(w => w.Z).ToArray(), true) })
    {
        double mean = values.Average();
        double variance = values.Sum(v => (v - mean) * (v - mean)) / values.Length;
        double lag1 = 0;
        for (int i = 1; i < values.Length; i++) lag1 += (values[i] - mean) * (values[i - 1] - mean);
        result.Add(signed ? mean : Math.Abs(mean));   // z is GC and keeps its sign on both strands
        result.Add(Math.Sqrt(variance));
        result.Add(variance == 0 ? 0 : lag1 / ((values.Length - 1) * variance));
    }
    return result.ToArray();
}

static double[][] Standardize(double[][] rows)
{
    int width = rows[0].Length;
    var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
    var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(rows.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
    return rows.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray();
}

// Rank-sum AUC with average ranks for ties, as RetrievalClassificationBenchmark.ComputeAuc; NaN when a class is empty.
static (double Auc, int Positives, int Negatives) Auc(List<(double Score, bool Positive)> pairs)
{
    int positives = pairs.Count(p => p.Positive), negatives = pairs.Count - positives;
    if (positives == 0 || negatives == 0) return (double.NaN, positives, negatives);
    var sorted = pairs.OrderBy(p => p.Score).ToArray();
    double rankSum = 0;
    for (int i = 0; i < sorted.Length;)
    {
        int j = i;
        while (j + 1 < sorted.Length && sorted[j + 1].Score == sorted[i].Score) j++;
        double averageRank = (i + j + 2) / 2.0;
        for (int t = i; t <= j; t++) if (sorted[t].Positive) rankSum += averageRank;
        i = j + 1;
    }
    return ((rankSum - positives * (positives + 1) / 2.0) / ((double)positives * negatives), positives, negatives);
}

record Genome(long Id, long Length, double Gc, string Phylum, string Family, string Genus, string Species, bool Circular);
