#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// P5, E3 (2026-10-30): does the Z-curve add anything to Karlin's delta*? r1_baselines.cs's sample and signatures (437 genomes, stored
// data only). Combined distances delta* + w * d_summary, each distance divided by its SD over all pairs, for w in {0.1, 0.25, 0.5, 1}
// (w = 0 is Karlin alone); and, as the control, delta* + w * d_gcLength on the same grid: if the window summary adds no more than GC
// and length do, the Z-curve contributes nothing of its own. The grid and the control were fixed before the run.
// AUC per GTDB rank, all pairs and |ΔGC| < 0.02, with 95% bootstrap CIs (200 resamples of genomes) and the paired difference against
// Karlin alone; NJ and UPGMA clades recovered per GTDB rank. Args: output path (default S52b_E3_KarlinPlusZCurve.csv); trees to *_Trees.csv.
const int maxGenera = 150, maxSpeciesPerGenus = 4, resamples = 200, minLength = 256_000, factor = 64;
double[] weights = { 0.1, 0.25, 0.5, 1 };
var output = args.Length > 0 ? args[0] : "S52b_E3_KarlinPlusZCurve.csv";
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
        all.Add(new Genome(r.GetInt64(0), r.GetInt64(1), gc, counts, r.GetString(3), r.GetString(4), r.GetString(5), r.GetString(6)));
    }

var rng = new Random(20261029);
var perSpecies = all.GroupBy(g => g.Species).Select(s => s.OrderBy(_ => rng.Next()).First()).ToList();
var sample = perSpecies.GroupBy(g => g.Genus).Where(s => s.Count() >= 2).OrderBy(_ => rng.Next()).Take(maxGenera)
    .SelectMany(s => s.OrderBy(_ => rng.Next()).Take(maxSpeciesPerGenus)).ToList();
var fine = await ProfileLevels.ReadManyAsync(dataSource, sample.Select(g => g.Id).ToArray(),
    RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
sample = sample.Where(g => fine.ContainsKey(g.Id)).ToList();
int n = sample.Count;
var coarse = sample.Select(g => fine[g.Id].Downsample(factor)).ToArray();
var rho = sample.Select(g => GenomeSignatures.KarlinRho(g.Dinucleotides)).ToArray();
Console.Error.WriteLine($"sample {n} genomes");

var karlin = new double[n, n];
Parallel.For(0, n, a =>
{
    for (int b = a + 1; b < n; b++)
        karlin[a, b] = karlin[b, a] = GenomeSignatures.KarlinDelta(rho[a], rho[b]);
});
var summary = StandardizedEuclidean(coarse.Select(GenomeSignatures.WindowSummary));
var gcLength = StandardizedEuclidean(sample.Select(g => new[] { g.Gc, Math.Log10(g.Length) }));
var (karlinUnit, summaryUnit, gcLengthUnit) = (UnitSd(karlin), UnitSd(summary), UnitSd(gcLength));
var distances = new List<(string Name, double[,] Matrix)> { ("karlin", karlin), ("summary", summary), ("gcLength", gcLength) };
foreach (double w in weights)
{
    distances.Add(($"karlin+{w}*summary", Sum(karlinUnit, summaryUnit, w)));
    distances.Add(($"karlin+{w}*gcLength", Sum(karlinUnit, gcLengthUnit, w)));
}
string[] names = distances.Select(d => d.Name).ToArray();

var ranks = new (string Name, Func<Genome, string> Of)[] { ("genus", g => g.Genus), ("family", g => g.Family), ("phylum", g => g.Phylum) };
var strata = new (string Name, double Min, double Max)[] { ("all", 0, 1), ("dGC<0.02", 0, 0.02) };

// AUCs of one draw of genomes (a genome drawn twice is not paired with itself): [rank, stratum, distance]; score = -distance.
(double Auc, int Positives, int Negatives)[,,] Aucs(int[] drawn)
{
    var result = new (double, int, int)[ranks.Length, strata.Length, distances.Count];
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
            for (int c = 0; c < distances.Count; c++)
                result[k, st, c] = Auc(pairs.Select(p => (-distances[c].Matrix[p.A, p.B], p.Positive)).ToList());
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
    return string.Create(CultureInfo.InvariantCulture, $"{At(0.025):F4}, {At(0.975):F4}");
}

int reference = Array.IndexOf(names, "karlin");
using (var report = File.CreateText(output))
{
    report.WriteLine("Distance, Rank, Stratum, Positives, Negatives, Auc, CiLow, CiHigh, MinusKarlin, MinusKarlinCiLow, MinusKarlinCiHigh");
    for (int c = 0; c < distances.Count; c++)
        for (int k = 0; k < ranks.Length; k++)
            for (int st = 0; st < strata.Length; st++)
            {
                var (auc, positives, negatives) = observed[k, st, c];
                string line = string.Create(CultureInfo.InvariantCulture,
                    $"{names[c]}, {ranks[k].Name}, {strata[st].Name}, {positives}, {negatives}, {auc:F4}, " +
                    $"{Interval(draws.Select(d => d[k, st, c].Auc))}, {auc - observed[k, st, reference].Auc:F4}, ") +
                    Interval(draws.Select(d => d[k, st, c].Auc - d[k, st, reference].Auc));
                report.WriteLine(line);
                Console.WriteLine(line);
            }
}

using (var trees = File.CreateText(Path.ChangeExtension(output, null) + "_Trees.csv"))
{
    trees.WriteLine("Distance, Method, Rank, Taxa, Clades, MeanConsistency");
    var labels = sample.Select(g => $"{g.Id}|{g.Genus}").ToArray();
    foreach (var (name, matrix) in distances)
        foreach (var (method, tree) in new[] { ("NJ", Phylogeny.NeighbourJoining(matrix, labels)), ("UPGMA", Phylogeny.Upgma(matrix, labels)) })
            foreach (var (rank, of) in ranks)
            {
                var taxa = sample.Select(of).ToArray();
                var recovery = tree.Recovery(taxa);
                string line = string.Create(CultureInfo.InvariantCulture,
                    $"{name}, {method}, {rank}, {recovery.Count}, {recovery.Count(x => x.Monophyletic)}, {recovery.Average(x => x.Consistency):F3}");
                trees.WriteLine(line);
                Console.WriteLine(line);
            }
}

// A distance divided by its SD over all pairs a < b, so two distances can be added on one scale.
static double[,] UnitSd(double[,] matrix)
{
    int count = matrix.GetLength(0);
    var values = new List<double>();
    for (int a = 0; a < count; a++)
        for (int b = a + 1; b < count; b++)
            values.Add(matrix[a, b]);
    double mean = values.Average(), sd = Math.Sqrt(values.Average(v => (v - mean) * (v - mean)));
    var result = new double[count, count];
    for (int a = 0; a < count; a++)
        for (int b = 0; b < count; b++)
            result[a, b] = matrix[a, b] / sd;
    return result;
}

static double[,] Sum(double[,] first, double[,] second, double weight)
{
    int count = first.GetLength(0);
    var result = new double[count, count];
    for (int a = 0; a < count; a++)
        for (int b = 0; b < count; b++)
            result[a, b] = first[a, b] + weight * second[a, b];
    return result;
}

static double[,] StandardizedEuclidean(IEnumerable<double[]> rows)
{
    var points = rows.ToArray();
    int width = points[0].Length;
    var means = Enumerable.Range(0, width).Select(f => points.Average(r => r[f])).ToArray();
    var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(points.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
    return Euclidean(points.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray());
}

static double[,] Euclidean(double[][] points)
{
    int count = points.Length;
    var result = new double[count, count];
    for (int a = 0; a < count; a++)
        for (int b = a + 1; b < count; b++)
            result[a, b] = result[b, a] = Math.Sqrt(points[a].Zip(points[b], (p, q) => (p - q) * (p - q)).Sum());
    return result;
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

record Genome(long Id, long Length, double Gc, long[] Dinucleotides, string Phylum, string Family, string Genus, string Species);
