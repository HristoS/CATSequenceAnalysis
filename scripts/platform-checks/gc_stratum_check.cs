#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S5.1 gate (P5, R2): does today's distance, 1 - Compare, separate taxa beyond GC content on pairs whose GC is nearly equal?
// R1 predicts the added value has to show in the low-|ΔGC| stratum, where GC-only is near chance. Genomes per D9: the largest
// bacterial or archaeal record of each source, labelled by GTDB (gtdb_taxon), at least 256 kbp so even the 32 kbp level keeps
// 8 windows; the floor doesn't depend on the level, so every level scores the same sample.
// One genome per GTDB species; up to 150 genera with at least 2 species, up to 4 species each, so same-genus pairs are
// different species. All pairs compared at one fixed level (default 32 kbp, the level find-similar uses for chromosomes),
// as stored (no strand flip or origin rotation: those are R2 candidates, not the baseline). GC from the stored dinucleotide
// counts (first base of each pair). AUC by rank sum per rank (genus, family, phylum: positive = same taxon at that rank) and
// stratum of |ΔGC|; 95% CIs from 200 bootstrap resamples of genomes, and the paired CI of AUC(Z-curve) - AUC(GC).
// Args: output path (default S51_GcStratumCheck.csv), window size (default 32000).
const int maxGenera = 150, maxSpeciesPerGenus = 4, resamples = 200, minLength = 256_000;
var output = args.Length > 0 ? args[0] : "S51_GcStratumCheck.csv";
int windowSize = args.Length > 1 ? int.Parse(args[1]) : 32_000;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var all = new List<Genome>();
await using (var cmd = dataSource.CreateCommand(
    "select distinct on (r.source_id) r.id, r.length_bp, r.dinucleotide_counts, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species " +
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
        all.Add(new Genome(r.GetInt64(0), r.GetInt64(1), gc, r.GetString(3), r.GetString(4), r.GetString(5), r.GetString(6)));
    }

var rng = new Random(20261029);
var perSpecies = all.GroupBy(g => g.Species).Select(s => s.OrderBy(_ => rng.Next()).First()).ToList();
var sample = perSpecies.GroupBy(g => g.Genus).Where(s => s.Count() >= 2).OrderBy(_ => rng.Next()).Take(maxGenera)
    .SelectMany(s => s.OrderBy(_ => rng.Next()).Take(maxSpeciesPerGenus)).ToList();

var levels = await ProfileLevels.ReadManyAsync(dataSource, sample.Select(g => g.Id).ToArray(),
    RecordStore.Method, RecordStore.MethodVersion, windowSize);
sample = sample.Where(g => levels.ContainsKey(g.Id)).ToList();
Console.Error.WriteLine($"{all.Count:N0} labelled genomes, {perSpecies.Count:N0} species; sample {sample.Count} genomes in " +
    $"{sample.Select(g => g.Genus).Distinct().Count()} genera, {sample.Select(g => g.Family).Distinct().Count()} families, " +
    $"{sample.Select(g => g.Phylum).Distinct().Count()} phyla");

int n = sample.Count;
var zCurve = new double[n, n];
var watch = System.Diagnostics.Stopwatch.StartNew();
Parallel.For(0, n, a =>
{
    for (int b = a + 1; b < n; b++)
        zCurve[a, b] = zCurve[b, a] = ZCurveWindowedProfile.Compare(levels[sample[a].Id], levels[sample[b].Id]);
});
Console.Error.WriteLine($"{n * (n - 1) / 2:N0} pairs compared in {watch.Elapsed.TotalSeconds:F1} s");

var ranks = new (string Name, Func<Genome, string> Of)[] { ("genus", g => g.Genus), ("family", g => g.Family), ("phylum", g => g.Phylum) };
var strata = new (string Name, double Min, double Max)[] { ("all", 0, 1), ("dGC<0.02", 0, 0.02), ("dGC 0.02-0.10", 0.02, 0.10), ("dGC>=0.10", 0.10, 1) };

// AUCs of one draw of genomes (indexes into sample; a genome drawn twice is not paired with itself): [rank, stratum, method].
(double Auc, int Positives, int Negatives)[,,] Aucs(int[] drawn)
{
    var result = new (double, int, int)[ranks.Length, strata.Length, 2];
    for (int k = 0; k < ranks.Length; k++)
        for (int s = 0; s < strata.Length; s++)
        {
            var zPairs = new List<(double Score, bool Positive)>();
            var gcPairs = new List<(double Score, bool Positive)>();
            for (int i = 0; i < drawn.Length; i++)
                for (int j = i + 1; j < drawn.Length; j++)
                {
                    int a = drawn[i], b = drawn[j];
                    if (a == b) continue;
                    double dGc = Math.Abs(sample[a].Gc - sample[b].Gc);
                    if (dGc < strata[s].Min || dGc >= strata[s].Max) continue;
                    bool positive = ranks[k].Of(sample[a]) == ranks[k].Of(sample[b]);
                    zPairs.Add((zCurve[a, b], positive));
                    gcPairs.Add((-dGc, positive));
                }
            result[k, s, 0] = Auc(zPairs);
            result[k, s, 1] = Auc(gcPairs);
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

double Percentile(IEnumerable<double> values, double p)
{
    var sorted = values.Where(v => !double.IsNaN(v)).Order().ToArray();
    return sorted.Length == 0 ? double.NaN : sorted[(int)Math.Clamp(Math.Round(p * (sorted.Length - 1)), 0, sorted.Length - 1)];
}

using var report = File.CreateText(output);
report.WriteLine("WindowSize, Rank, Stratum, Positives, Negatives, AucZCurve, ZCurveCiLow, ZCurveCiHigh, AucGc, GcCiLow, GcCiHigh, " +
    "AucDifference, DifferenceCiLow, DifferenceCiHigh");
for (int k = 0; k < ranks.Length; k++)
    for (int s = 0; s < strata.Length; s++)
    {
        var (zAuc, positives, negatives) = observed[k, s, 0];
        double gcAuc = observed[k, s, 1].Auc;
        var zDraws = draws.Select(d => d[k, s, 0].Auc).ToList();
        var gcDraws = draws.Select(d => d[k, s, 1].Auc).ToList();
        var differences = draws.Select(d => d[k, s, 0].Auc - d[k, s, 1].Auc).ToList();
        string line = $"{windowSize}, {ranks[k].Name}, {strata[s].Name}, {positives}, {negatives}, " +
            $"{zAuc:F4}, {Percentile(zDraws, 0.025):F4}, {Percentile(zDraws, 0.975):F4}, " +
            $"{gcAuc:F4}, {Percentile(gcDraws, 0.025):F4}, {Percentile(gcDraws, 0.975):F4}, " +
            $"{zAuc - gcAuc:F4}, {Percentile(differences, 0.025):F4}, {Percentile(differences, 0.975):F4}";
        report.WriteLine(line);
        Console.WriteLine(line);
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

record Genome(long Id, long Length, double Gc, string Phylum, string Family, string Genus, string Species);
