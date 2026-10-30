#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S5.2 (P5, R3) on the catalogue: neighbour-joining and UPGMA trees from Z-curve distances only, on r2_candidates.cs's sample (same
// rule and seed: 437 GTDB-labelled genomes, one per species, up to 4 species in each of 150 genera), at 32 kbp downsampled from the
// stored 500 bp profile. Distances: compare (1 - Compare as stored), strandOrigin (1 - the best Compare over both strands and, for
// circular chromosomes OriginRotation accepts, origin rotation) and summary (standardized window statistics, Euclidean); gc and
// gcLength trees as context. Per tree and rank (genus, family, phylum; GTDB and NCBI's current lineage from taxon): clades
// recovered among taxa with at least 2 leaves, mean consistency, and a 1,000-permutation label null on the same tree. Also how
// much each distance follows GC: Spearman correlation of the distance with |ΔGC| over all pairs.
// Args: output path (default S52_Trees_Catalogue.csv); Newick files go to trees/ beside it, leaf orders to *_LeafOrder.csv.
const int maxGenera = 150, maxSpeciesPerGenus = 4, minLength = 256_000, factor = 64, permutations = 1000;
var output = args.Length > 0 ? args[0] : "S52_Trees_Catalogue.csv";
var treeDirectory = Path.Combine(Path.GetDirectoryName(Path.GetFullPath(output))!, "trees");
Directory.CreateDirectory(treeDirectory);
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var all = new List<Genome>();
await using (var cmd = dataSource.CreateCommand(
    "select distinct on (r.source_id) r.id, r.length_bp, r.dinucleotide_counts, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species, " +
    "r.annotations->'contig'->>'topology', t.tax_phylum, t.tax_family, t.tax_genus, coalesce(r.name, '') " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "left join taxon t on t.taxid = r.ncbi_taxid " +
    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null " +
    "order by r.source_id, r.length_bp desc, r.id"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        if (r.GetInt64(1) < minLength || r.IsDBNull(6) || r.IsDBNull(5) || r.IsDBNull(4) || r.IsDBNull(3)) continue;
        var counts = r.GetFieldValue<long[]>(2);
        double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
        string? Ncbi(int column) => r.IsDBNull(column) ? null : r.GetString(column);
        all.Add(new Genome(r.GetInt64(0), r.GetInt64(1), gc, r.GetString(3), r.GetString(4), r.GetString(5), r.GetString(6),
            !r.IsDBNull(7) && r.GetString(7) == "circular", Ncbi(8), Ncbi(9), Ncbi(10), r.GetString(11)));
    }

var rng = new Random(20261029);
var perSpecies = all.GroupBy(g => g.Species).Select(s => s.OrderBy(_ => rng.Next()).First()).ToList();
var sample = perSpecies.GroupBy(g => g.Genus).Where(s => s.Count() >= 2).OrderBy(_ => rng.Next()).Take(maxGenera)
    .SelectMany(s => s.OrderBy(_ => rng.Next()).Take(maxSpeciesPerGenus)).ToList();

var fine = await ProfileLevels.ReadManyAsync(dataSource, sample.Select(g => g.Id).ToArray(),
    RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
sample = sample.Where(g => fine.ContainsKey(g.Id)).ToList();
int n = sample.Count;

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
    return list.Select(v => v.Downsample(factor)).ToArray();
}).ToArray();
Console.Error.WriteLine($"sample {n} genomes, {versions.Count(v => v.Length > 2)} rotated");

var plain = new double[n, n];
var best = new double[n, n];
Parallel.For(0, n, a =>
{
    for (int b = a + 1; b < n; b++)
    {
        plain[a, b] = plain[b, a] = 1 - ZCurveWindowedProfile.Compare(versions[a][0], versions[b][0]);
        double top = double.NegativeInfinity;
        foreach (var x in versions[a])
            foreach (var y in versions[b])
                top = Math.Max(top, ZCurveWindowedProfile.Compare(x, y));
        best[a, b] = best[b, a] = 1 - top;
    }
});
var distances = new (string Name, double[,] Matrix)[]
{
    ("compare", plain),
    ("strandOrigin", best),
    ("summary", Euclidean(Standardize(versions.Select(v => Features(v[0])).ToArray()))),
    ("gc", Euclidean(sample.Select(g => new[] { g.Gc }).ToArray())),
    ("gcLength", Euclidean(Standardize(sample.Select(g => new[] { g.Gc, Math.Log10(g.Length) }).ToArray()))),
};

// Missing NCBI ranks get a label of their own, so they count as singletons and drop out.
var ranks = new (string Name, string[] Taxa)[]
{
    ("gtdb genus", sample.Select(g => g.Genus).ToArray()),
    ("gtdb family", sample.Select(g => g.Family).ToArray()),
    ("gtdb phylum", sample.Select(g => g.Phylum).ToArray()),
    ("ncbi genus", sample.Select(g => g.NcbiGenus ?? $"none {g.Id}").ToArray()),
    ("ncbi family", sample.Select(g => g.NcbiFamily ?? $"none {g.Id}").ToArray()),
    ("ncbi phylum", sample.Select(g => g.NcbiPhylum ?? $"none {g.Id}").ToArray()),
};
var labels = sample.Select(g => $"{g.Name} {g.Id}|{g.Genus}|{g.Phylum}").ToArray();   // the record id keeps leaf labels unique

using var report = File.CreateText(output);
report.WriteLine("Distance, SpearmanWithAbsGcDifference, Method, Rank, Taxa, Monophyletic, MeanConsistency, NullMeanMonophyletic, NullP, NjVsUpgmaNrf");
using var order = File.CreateText(Path.ChangeExtension(output, null) + "_LeafOrder.csv");
order.WriteLine("Distance, Method, Position, Id, Name, GtdbGenus, GtdbFamily, GtdbPhylum, Gc, LengthBp");
foreach (var (name, matrix) in distances)
{
    double gcFollow = Spearman(matrix, (a, b) => Math.Abs(sample[a].Gc - sample[b].Gc));
    var trees = new[] { ("NJ", Phylogeny.NeighbourJoining(matrix, labels)), ("UPGMA", Phylogeny.Upgma(matrix, labels)) };
    double disagreement = trees[0].Item2.NormalizedRobinsonFoulds(trees[1].Item2);
    foreach (var (method, tree) in trees)
    {
        File.WriteAllText(Path.Combine(treeDirectory, $"S52_catalogue_32kbp_{name}_{method}.nwk"), tree.ToNewick() + Environment.NewLine);
        int position = 0;
        foreach (int leaf in tree.LeafOrder())
        {
            var g = sample[leaf];
            order.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"{name}, {method}, {position++}, {g.Id}, \"{g.Name}\", {g.Genus}, {g.Family}, {g.Phylum}, {g.Gc:F4}, {g.Length}"));
        }
        foreach (var (rank, taxa) in ranks)
        {
            var recovery = tree.Recovery(taxa);
            int observed = recovery.Count(x => x.Monophyletic);
            var nullCounts = new int[permutations];
            Parallel.For(0, permutations, p =>
            {
                var local = new Random(20261029 + p);
                nullCounts[p] = tree.Recovery(taxa.OrderBy(_ => local.Next()).ToArray()).Count(x => x.Monophyletic);
            });
            double pValue = (1 + nullCounts.Count(c => c >= observed)) / (permutations + 1.0);
            string line = string.Create(CultureInfo.InvariantCulture,
                $"{name}, {gcFollow:F3}, {method}, {rank}, {recovery.Count}, {observed}, {recovery.Average(x => x.Consistency):F3}, " +
                $"{nullCounts.Average():F2}, {pValue:F4}, {disagreement:F3}");
            report.WriteLine(line);
            Console.WriteLine(line);
        }
    }
}

// Spearman's rank correlation over all pairs a < b of the matrix against another pair value (average ranks for ties).
double Spearman(double[,] matrix, Func<int, int, double> other)
{
    var xs = new List<double>();
    var ys = new List<double>();
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
        {
            xs.Add(matrix[a, b]);
            ys.Add(other(a, b));
        }
    var rx = Ranks(xs);
    var ry = Ranks(ys);
    double mx = rx.Average(), my = ry.Average();
    double sxy = 0, sxx = 0, syy = 0;
    for (int i = 0; i < rx.Length; i++)
    {
        sxy += (rx[i] - mx) * (ry[i] - my);
        sxx += (rx[i] - mx) * (rx[i] - mx);
        syy += (ry[i] - my) * (ry[i] - my);
    }
    return sxy / Math.Sqrt(sxx * syy);
}

static double[] Ranks(List<double> values)
{
    var index = Enumerable.Range(0, values.Count).OrderBy(i => values[i]).ToArray();
    var ranks = new double[values.Count];
    for (int i = 0; i < index.Length;)
    {
        int j = i;
        while (j + 1 < index.Length && values[index[j + 1]] == values[index[i]]) j++;
        for (int t = i; t <= j; t++) ranks[index[t]] = (i + j) / 2.0;
        i = j + 1;
    }
    return ranks;
}

static ZCurveWindowedProfile Flip(ZCurveWindowedProfile profile) =>
    ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize, profile.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

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

static double[,] Euclidean(double[][] points)
{
    int count = points.Length;
    var result = new double[count, count];
    for (int a = 0; a < count; a++)
        for (int b = a + 1; b < count; b++)
            result[a, b] = result[b, a] = Math.Sqrt(points[a].Zip(points[b], (p, q) => (p - q) * (p - q)).Sum());
    return result;
}

record Genome(long Id, long Length, double Gc, string Phylum, string Family, string Genus, string Species, bool Circular,
    string? NcbiPhylum, string? NcbiFamily, string? NcbiGenus, string Name);
