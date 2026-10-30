#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using System.IO.Compression;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// P5, measure B (2026-10-30): does each whole-genome signature track evolutionary distance? GTDB's reference trees (bac120 and ar53,
// built from 120 / 53 concatenated marker proteins; the release the catalogue's labels come from) give a patristic distance between
// species representatives. Each catalogue genome (one per GTDB species, as r4_map.cs) is placed at its species' representative leaf.
// Per domain, Spearman correlation of signature distance with patristic distance over genome pairs, overall and within depth bands
// (same family, same phylum but another family, other phyla), for the second-order Z-curve (9, Euclidean), Karlin's rho* (16,
// Euclidean), |ΔGC| and GC + length; with a 95% bootstrap CI of the Z-curve minus Karlin difference (50 resamples of 1,000 genomes).
// Point estimates on a random 2,500 genomes per domain when there are more. Args: the folder holding bac120.tree.gz, ar53.tree.gz,
// bac120_taxonomy.tsv.gz and ar53_taxonomy.tsv.gz (GTDB's release files; not kept in the repository), and the output path (default
// S54_GtdbTreeDistance.csv).
const int pointGenomes = 2500, resampleGenomes = 1000, resamples = 50;
var gtdb = args[0];
var output = args.Length > 1 ? args[1] : "S54_GtdbTreeDistance.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var genomes = new List<Genome>();
await using (var cmd = dataSource.CreateCommand(
    "select distinct on (r.source_id) r.id, r.length_bp, r.dinucleotide_counts, g.tax_domain, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null order by r.source_id, r.length_bp desc, r.id"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        if (r.IsDBNull(7) || r.IsDBNull(3)) continue;
        var counts = r.GetFieldValue<long[]>(2);
        if (counts.Sum() == 0) continue;
        double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
        string? Text(int column) => r.IsDBNull(column) ? null : r.GetString(column);
        genomes.Add(new Genome(r.GetInt64(0), r.GetInt64(1), gc, counts, r.GetString(3), Text(4), Text(5), Text(6), r.GetString(7)));
    }
var rng = new Random(20261030);
genomes = genomes.GroupBy(g => g.Species).Select(s => s.OrderBy(_ => rng.Next()).First()).OrderBy(g => g.Id).ToList();

using var report = File.CreateText(output);
report.WriteLine("Domain, Band, Pairs, Distance, Spearman, MinusKarlin, MinusKarlinCiLow, MinusKarlinCiHigh");
foreach (var (domain, prefix) in new[] { ("Bacteria", "bac120"), ("Archaea", "ar53") })
{
    // Species -> representative leaf: the tree's leaves are the representatives, the taxonomy file names each one's species.
    var speciesOf = new Dictionary<string, string>();
    using (var reader = new StreamReader(new GZipStream(File.OpenRead(Path.Combine(gtdb, $"{prefix}_taxonomy.tsv.gz")), CompressionMode.Decompress)))
        for (string? line; (line = reader.ReadLine()) != null;)
        {
            var parts = line.Split('\t');
            speciesOf[parts[0]] = parts[1][(parts[1].LastIndexOf("s__", StringComparison.Ordinal) + 3)..];
        }
    var tree = Tree.Parse(new StreamReader(new GZipStream(File.OpenRead(Path.Combine(gtdb, $"{prefix}.tree.gz")), CompressionMode.Decompress)).ReadToEnd());
    var leafOfSpecies = new Dictionary<string, int>();
    foreach (var (label, node) in tree.Leaves)
        if (speciesOf.TryGetValue(label, out var species)) leafOfSpecies[species] = node;

    var placed = genomes.Where(g => g.Domain == $"d__{domain}" || g.Domain == domain)
        .Select(g => (Genome: g, Leaf: leafOfSpecies.TryGetValue(StripPrefix(g.Species), out int leaf) ? leaf : -1)).Where(x => x.Leaf >= 0).ToList();
    Console.Error.WriteLine($"{domain}: {tree.Leaves.Count} tree leaves, {placed.Count} catalogue species placed");

    var all = placed.Select(p => p.Genome).ToArray();
    var leaves = placed.Select(p => p.Leaf).ToArray();
    var features = new Func<Genome, double[]>[]
    {
        g => GenomeSignatures.SecondOrder(g.Dinucleotides),
        g => GenomeSignatures.KarlinRho(g.Dinucleotides),
        g => new[] { g.Gc },
    };
    var names = new[] { "secondOrder", "karlin", "gc", "gcLength" };
    var vectors = features.Select(f => all.Select(f).ToArray()).ToList();
    vectors.Add(Standardize(all.Select(g => new[] { g.Gc, Math.Log10(g.Length) }).ToArray()));

    string Band(int a, int b) =>
        all[a].Family is not null && all[a].Family == all[b].Family ? "same family"
        : all[a].Phylum is not null && all[a].Phylum == all[b].Phylum ? "same phylum, other family" : "other phyla";

    // Spearman per band and distance over the pairs of one set of genomes (indexes into all; a genome drawn twice is not paired with itself).
    Dictionary<string, (int Pairs, double[] Rho)> Correlations(int[] chosen)
    {
        var byBand = new Dictionary<string, (List<double> Tree, List<double>[] Signature)>();
        for (int i = 0; i < chosen.Length; i++)
            for (int j = i + 1; j < chosen.Length; j++)
            {
                int a = chosen[i], b = chosen[j];
                if (a == b) continue;
                double patristic = tree.Distance(leaves[a], leaves[b]);
                foreach (var band in new[] { "all", Band(a, b) })
                {
                    if (!byBand.TryGetValue(band, out var lists)) byBand[band] = lists = (new List<double>(), names.Select(_ => new List<double>()).ToArray());
                    lists.Tree.Add(patristic);
                    for (int d = 0; d < names.Length; d++) lists.Signature[d].Add(Euclidean(vectors[d][a], vectors[d][b]));
                }
            }
        return byBand.ToDictionary(kv => kv.Key, kv => (kv.Value.Tree.Count, kv.Value.Signature.Select(s => Spearman(kv.Value.Tree, s)).ToArray()));
    }

    var local = new Random(20261030);
    var point = Correlations(Enumerable.Range(0, all.Length).OrderBy(_ => local.Next()).Take(pointGenomes).ToArray());
    var draws = new Dictionary<string, (int, double[])>[resamples];
    Parallel.For(0, resamples, new ParallelOptions { MaxDegreeOfParallelism = 4 }, r =>
    {
        var draw = new Random(20261100 + r);
        draws[r] = Correlations(Enumerable.Range(0, resampleGenomes).Select(_ => draw.Next(all.Length)).ToArray());
    });
    foreach (var band in new[] { "all", "same family", "same phylum, other family", "other phyla" })
    {
        if (!point.TryGetValue(band, out var values)) continue;
        for (int d = 0; d < names.Length; d++)
        {
            var differences = draws.Where(x => x.ContainsKey(band)).Select(x => x[band].Item2[d] - x[band].Item2[1]).Order().ToArray();
            string line = string.Create(CultureInfo.InvariantCulture,
                $"{domain}, {band}, {values.Pairs}, {names[d]}, {values.Rho[d]:F4}, {values.Rho[d] - values.Rho[1]:F4}, " +
                $"{differences[(int)(0.025 * (differences.Length - 1))]:F4}, {differences[(int)(0.975 * (differences.Length - 1))]:F4}");
            report.WriteLine(line);
            Console.WriteLine(line);
        }
    }
}

static string StripPrefix(string species) => species.StartsWith("s__", StringComparison.Ordinal) ? species[3..] : species;

static double Euclidean(double[] a, double[] b)
{
    double sum = 0;
    for (int f = 0; f < a.Length; f++) sum += (a[f] - b[f]) * (a[f] - b[f]);
    return Math.Sqrt(sum);
}

static double[][] Standardize(double[][] rows)
{
    int width = rows[0].Length;
    var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
    var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(rows.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
    return rows.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray();
}

static double Spearman(List<double> xs, List<double> ys)
{
    var (rx, ry) = (Ranks(xs), Ranks(ys));
    double mean = (xs.Count - 1) / 2.0, sxy = 0, sxx = 0, syy = 0;
    for (int i = 0; i < rx.Length; i++)
    {
        sxy += (rx[i] - mean) * (ry[i] - mean);
        sxx += (rx[i] - mean) * (rx[i] - mean);
        syy += (ry[i] - mean) * (ry[i] - mean);
    }
    return sxy / Math.Sqrt(sxx * syy);
}

static double[] Ranks(List<double> values)
{
    var index = Enumerable.Range(0, values.Count).ToArray();
    var keys = values.ToArray();
    Array.Sort(keys, index);
    var ranks = new double[values.Count];
    for (int i = 0; i < index.Length;)
    {
        int j = i;
        while (j + 1 < index.Length && keys[j + 1] == keys[i]) j++;
        for (int t = i; t <= j; t++) ranks[index[t]] = (i + j) / 2.0;
        i = j + 1;
    }
    return ranks;
}

record Genome(long Id, long Length, double Gc, long[] Dinucleotides, string Domain, string? Phylum, string? Family, string? Genus, string Species);

// A rooted tree read from Newick, with patristic distances by lowest common ancestor (Euler tour and a sparse table of minima).
sealed class Tree
{
    private readonly List<int> parent = new();
    private readonly List<double> length = new();
    private readonly List<List<int>> children = new();
    private double[] depth = Array.Empty<double>();
    private int[] level = Array.Empty<int>();
    private int[] first = Array.Empty<int>();
    private int[][] table = Array.Empty<int[]>();
    private int[] euler = Array.Empty<int>();

    public List<(string Label, int Node)> Leaves { get; } = new();

    private int Add(int up)
    {
        parent.Add(up);
        length.Add(0);
        children.Add(new List<int>());
        if (up >= 0) children[up].Add(parent.Count - 1);
        return parent.Count - 1;
    }

    // Iterative, so the depth of GTDB's trees can't overflow the stack. Labels may be quoted ('100.0:p__Bacillota').
    public static Tree Parse(string text)
    {
        var tree = new Tree();
        int current = tree.Add(-1);
        for (int i = 0; i < text.Length;)
        {
            char c = text[i];
            if (c == '(') { current = tree.Add(current); i++; }
            else if (c == ',') { current = tree.Add(tree.parent[current]); i++; }
            else if (c == ')') { current = tree.parent[current]; i++; }
            else if (c == ';' || char.IsWhiteSpace(c)) i++;
            else if (c == ':')
            {
                int start = ++i;
                while (i < text.Length && ",);".IndexOf(text[i]) < 0) i++;
                tree.length[current] = double.Parse(text.AsSpan(start, i - start), CultureInfo.InvariantCulture);
            }
            else
            {
                string label;
                if (c == '\'')
                {
                    int end = text.IndexOf('\'', i + 1);
                    label = text.Substring(i + 1, end - i - 1);
                    i = end + 1;
                }
                else
                {
                    int start = i;
                    while (i < text.Length && ":,);".IndexOf(text[i]) < 0) i++;
                    label = text[start..i];
                }
                if (tree.children[current].Count == 0) tree.Leaves.Add((label, current));
            }
        }
        tree.Index();
        return tree;
    }

    private void Index()
    {
        int n = parent.Count;
        depth = new double[n];
        level = new int[n];
        first = new int[n];
        var tour = new List<int>(2 * n);
        var stack = new Stack<(int Node, int Child)>();
        stack.Push((0, 0));
        first[0] = 0;
        tour.Add(0);
        while (stack.Count > 0)
        {
            var (node, child) = stack.Pop();
            if (child < children[node].Count)
            {
                stack.Push((node, child + 1));
                int next = children[node][child];
                depth[next] = depth[node] + length[next];
                level[next] = level[node] + 1;
                first[next] = tour.Count;
                tour.Add(next);
                stack.Push((next, 0));
            }
            else if (stack.Count > 0) tour.Add(stack.Peek().Node);
        }
        euler = tour.ToArray();
        int logs = 1;
        while ((1 << logs) <= euler.Length) logs++;
        table = new int[logs][];
        table[0] = Enumerable.Range(0, euler.Length).ToArray();
        for (int k = 1; k < logs; k++)
        {
            table[k] = new int[euler.Length];
            for (int i = 0; i + (1 << k) <= euler.Length; i++)
            {
                int a = table[k - 1][i], b = table[k - 1][i + (1 << (k - 1))];
                table[k][i] = level[euler[a]] <= level[euler[b]] ? a : b;
            }
        }
    }

    public double Distance(int a, int b)
    {
        int lo = Math.Min(first[a], first[b]), hi = Math.Max(first[a], first[b]);
        int k = 31 - System.Numerics.BitOperations.LeadingZeroCount((uint)(hi - lo + 1));
        int x = table[k][lo], y = table[k][hi - (1 << k) + 1];
        int ancestor = euler[level[euler[x]] <= level[euler[y]] ? x : y];
        return depth[a] + depth[b] - 2 * depth[ancestor];
    }
}
