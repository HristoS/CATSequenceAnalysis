#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S5.3 (P5, R4): a map of the catalogue from the second-order Z-curve (GenomeSignatures.SecondOrder from the stored dinucleotide
// counts; no DNA). Two corpora, one genome per species (seed 20261030): cellular (bacteria and archaea, the largest record per source
// as D9, GTDB v232 labels) and viruses (NCBI's current lineage from taxon). Spaces: secondOrder (9-D) and its PCA to 2-D (centred, not
// rescaled: all nine are correlations); karlin (the 16 rho*, Euclidean) and its PCA to 2-D; gcLength (standardized GC and log10 length,
// 2-D), as context. Per space and rank: kNN purity (share of the 10 nearest genomes in the same taxon, over genomes with a label), its
// label-permutation null (50 shuffles of whole lineages); mean silhouette over taxa with at least 2 genomes; and the gap test: for every
// genome of a genus with at least 2 genomes, is its nearest genome outside that genus in the same family, and the same phylum? (null:
// 20 shuffles). Writes S53_Map_Metrics.csv, S53_Map_Loadings.csv and S53_Map_Coordinates_{cellular,viruses}.csv beside the given path.
const int k = 10, neighbourList = 300, purityShuffles = 50, gapShuffles = 20;
var output = args.Length > 0 ? args[0] : "S53_Map_Metrics.csv";
var directory = Path.GetDirectoryName(Path.GetFullPath(output))!;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
string[] secondOrderNames = { "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz" };
string[] karlinNames = Enumerable.Range(0, 16).Select(i => $"{"ACGT"[i / 4]}{"ACGT"[i % 4]}").ToArray();

var cellular = await LoadAsync(
    "select distinct on (r.source_id) r.id, r.name, r.length_bp, r.dinucleotide_counts, g.tax_domain, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null order by r.source_id, r.length_bp desc, r.id");
var viruses = await LoadAsync(
    "select r.id, r.name, r.length_bp, r.dinucleotide_counts, t.tax_realm, t.tax_phylum, t.tax_family, t.tax_genus, t.tax_species " +
    "from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid " +
    "where r.kind = 'virus' and r.dinucleotide_counts is not null and r.length_bp >= 1000 order by r.id");

using var metrics = File.CreateText(output);
metrics.WriteLine("Corpus, Space, Metric, Rank, Genomes, Value, Null");
using var loadings = File.CreateText(Path.Combine(directory, "S53_Map_Loadings.csv"));
loadings.WriteLine("Corpus, Space, Component, ExplainedVariance, " + string.Join(", ", karlinNames.Select((name, i) => $"F{i}")) + ", FeatureNames");

foreach (var (corpus, genomes) in new[] { ("cellular", cellular), ("viruses", viruses) })
{
    int n = genomes.Count;
    var secondOrder = genomes.Select(g => GenomeSignatures.SecondOrder(g.Dinucleotides)).ToArray();
    var karlin = genomes.Select(g => GenomeSignatures.KarlinRho(g.Dinucleotides)).ToArray();
    var secondOrderPca = Ordination.Pca(secondOrder, 2);
    var karlinPca = Ordination.Pca(karlin, 2);
    var gcLength = Standardize(genomes.Select(g => new[] { g.Gc, Math.Log10(g.Length) }).ToArray());
    Console.Error.WriteLine($"{corpus}: {n} genomes (one per species); second-order PC1/PC2 explain " +
        $"{secondOrderPca.ExplainedVariance[0]:P1} / {secondOrderPca.ExplainedVariance[1]:P1}, Karlin {karlinPca.ExplainedVariance[0]:P1} / {karlinPca.ExplainedVariance[1]:P1}");

    foreach (var (space, pca, names) in new[] { ("secondOrder2D", secondOrderPca, secondOrderNames), ("karlin2D", karlinPca, karlinNames) })
        for (int c = 0; c < 2; c++)
            loadings.WriteLine(string.Create(CultureInfo.InvariantCulture, $"{corpus}, {space}, PC{c + 1}, {pca.ExplainedVariance[c]:F4}, ") +
                string.Join(", ", Enumerable.Range(0, 16).Select(f => f < names.Length ? pca.Loadings[c][f].ToString("F4", CultureInfo.InvariantCulture) : "")) +
                $", {string.Join(" ", names)}");

    using (var coordinates = File.CreateText(Path.Combine(directory, $"S53_Map_Coordinates_{corpus}.csv")))
    {
        coordinates.WriteLine("Id, Name, Top, Phylum, Family, Genus, Species, Gc, LengthBp, SecondOrderPc1, SecondOrderPc2, KarlinPc1, KarlinPc2");
        for (int i = 0; i < n; i++)
        {
            var g = genomes[i];
            coordinates.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"{g.Id}, {Quote(g.Name)}, {Quote(g.Top)}, {Quote(g.Phylum)}, {Quote(g.Family)}, {Quote(g.Genus)}, {Quote(g.Species)}, {g.Gc:F4}, {g.Length}, " +
                $"{secondOrderPca.Scores[i][0]:F5}, {secondOrderPca.Scores[i][1]:F5}, {karlinPca.Scores[i][0]:F5}, {karlinPca.Scores[i][1]:F5}"));
        }
    }

    var ranks = new (string Name, Func<Genome, string?> Of)[] { ("genus", g => g.Genus), ("family", g => g.Family), ("phylum", g => g.Phylum) };
    foreach (var (space, points) in new[] { ("secondOrder2D", secondOrderPca.Scores), ("secondOrder9D", secondOrder),
        ("karlin2D", karlinPca.Scores), ("karlin16D", karlin), ("gcLength2D", gcLength) })
    {
        var neighbours = Neighbours(points);
        foreach (var (rank, of) in ranks)
        {
            var labels = genomes.Select(of).ToArray();
            double purity = Purity(neighbours, labels);
            var rng = new Random(20261030);
            double purityNull = Enumerable.Range(0, purityShuffles).Average(_ => Purity(neighbours, Shuffle(labels, rng)));
            double silhouette = Silhouette(points, labels);
            int queries = labels.Count(l => l is not null);
            Write(corpus, space, "knnPurity10", rank, queries, purity, purityNull);
            Write(corpus, space, "silhouette", rank, labels.GroupBy(l => l).Where(gr => gr.Key is not null && gr.Count() >= 2).Sum(gr => gr.Count()), silhouette, double.NaN);
        }

        // Gap test on whole lineages (genus, family, phylum together), so the null keeps them consistent.
        var lineages = genomes.Select(g => (g.Genus, g.Family, g.Phylum)).ToArray();
        var (tested, sameFamily, samePhylum) = Gap(neighbours, lineages);
        var gapRng = new Random(20261031);
        var shuffled = Enumerable.Range(0, gapShuffles).Select(_ => Gap(neighbours, Shuffle(lineages, gapRng))).ToList();
        Write(corpus, space, "gapSameFamily", "genus", tested, sameFamily, shuffled.Average(s => s.SameFamily));
        Write(corpus, space, "gapSamePhylum", "genus", tested, samePhylum, shuffled.Average(s => s.SamePhylum));
    }
}

void Write(string corpus, string space, string metric, string rank, int count, double value, double nullValue)
{
    string line = string.Create(CultureInfo.InvariantCulture, $"{corpus}, {space}, {metric}, {rank}, {count}, {value:F4}, {nullValue:F4}");
    metrics.WriteLine(line);
    Console.WriteLine(line);
}

async Task<List<Genome>> LoadAsync(string sql)
{
    var rows = new List<Genome>();
    await using (var cmd = dataSource.CreateCommand(sql))
    await using (var r = await cmd.ExecuteReaderAsync())
        while (await r.ReadAsync())
        {
            string? Text(int column) => r.IsDBNull(column) || r.GetString(column) == "" ? null : r.GetString(column);
            var counts = r.GetFieldValue<long[]>(3);
            if (counts.Sum() == 0) continue;
            double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
            rows.Add(new Genome(r.GetInt64(0), Text(1) ?? "", r.GetInt64(2), gc, counts, Text(4), Text(5), Text(6), Text(7), Text(8)));
        }
    // One genome per species; a genome without a species stands for itself.
    var rng = new Random(20261030);
    return rows.GroupBy(g => g.Species ?? $"id {g.Id}").Select(s => s.OrderBy(_ => rng.Next()).First()).OrderBy(g => g.Id).ToList();
}

// The nearest genomes of each genome, nearest first, up to neighbourList of them.
static int[][] Neighbours(double[][] points)
{
    int n = points.Length;
    var result = new int[n][];
    Parallel.For(0, n, i =>
    {
        var distances = new double[n];
        for (int j = 0; j < n; j++)
        {
            double sum = 0;
            for (int f = 0; f < points[i].Length; f++)
            {
                double d = points[i][f] - points[j][f];
                sum += d * d;
            }
            distances[j] = j == i ? double.PositiveInfinity : sum;
        }
        var order = Enumerable.Range(0, n).ToArray();
        Array.Sort(distances, order);
        result[i] = order.Take(Math.Min(neighbourList, n - 1)).ToArray();
    });
    return result;
}

static double Purity(int[][] neighbours, string?[] labels)
{
    double sum = 0;
    int count = 0;
    for (int i = 0; i < labels.Length; i++)
    {
        if (labels[i] is null) continue;
        sum += neighbours[i].Take(k).Count(j => labels[j] == labels[i]) / (double)k;
        count++;
    }
    return sum / count;
}

static (int Tested, double SameFamily, double SamePhylum) Gap(int[][] neighbours, (string? Genus, string? Family, string? Phylum)[] lineages)
{
    var genusSize = lineages.Where(l => l.Genus is not null).GroupBy(l => l.Genus).ToDictionary(g => g.Key!, g => g.Count());
    int tested = 0, family = 0, phylum = 0;
    for (int i = 0; i < lineages.Length; i++)
    {
        var l = lineages[i];
        if (l.Genus is null || genusSize[l.Genus] < 2) continue;
        int position = Array.FindIndex(neighbours[i], j => lineages[j].Genus != l.Genus);
        if (position < 0) continue;
        int outside = neighbours[i][position];
        tested++;
        if (l.Family is not null && lineages[outside].Family == l.Family) family++;
        if (l.Phylum is not null && lineages[outside].Phylum == l.Phylum) phylum++;
    }
    return (tested, (double)family / tested, (double)phylum / tested);
}

// Mean silhouette over genomes whose taxon has at least 2 genomes: (b - a) / max(a, b), a the mean distance to its own taxon, b the
// smallest mean distance to another such taxon.
static double Silhouette(double[][] points, string?[] labels)
{
    var taxa = labels.Where(l => l is not null).GroupBy(l => l).Where(g => g.Count() >= 2).Select(g => g.Key!).ToList();
    var index = taxa.Select((t, i) => (t, i)).ToDictionary(x => x.t, x => x.i);
    var member = labels.Select(l => l is not null && index.TryGetValue(l, out int t) ? t : -1).ToArray();
    var size = new int[taxa.Count];
    foreach (int t in member) if (t >= 0) size[t]++;
    int n = points.Length;
    var values = new double[n];
    Parallel.For(0, n, i =>
    {
        values[i] = double.NaN;
        if (member[i] < 0) return;
        var sums = new double[taxa.Count];
        for (int j = 0; j < n; j++)
        {
            if (j == i || member[j] < 0) continue;
            double sum = 0;
            for (int f = 0; f < points[i].Length; f++)
            {
                double d = points[i][f] - points[j][f];
                sum += d * d;
            }
            sums[member[j]] += Math.Sqrt(sum);
        }
        double a = sums[member[i]] / (size[member[i]] - 1), b = double.PositiveInfinity;
        for (int t = 0; t < taxa.Count; t++)
            if (t != member[i]) b = Math.Min(b, sums[t] / size[t]);
        values[i] = Math.Max(a, b) == 0 ? 0 : (b - a) / Math.Max(a, b);
    });
    return values.Where(v => !double.IsNaN(v)).Average();
}

static T[] Shuffle<T>(T[] values, Random rng) => values.OrderBy(_ => rng.Next()).ToArray();

static double[][] Standardize(double[][] rows)
{
    int width = rows[0].Length;
    var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
    var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(rows.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
    return rows.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray();
}

static string Quote(string? text) => text is null ? "" : "\"" + text.Replace("\"", "'") + "\"";

record Genome(long Id, string Name, long Length, double Gc, long[] Dinucleotides, string? Top, string? Phylum, string? Family, string? Genus, string? Species);
