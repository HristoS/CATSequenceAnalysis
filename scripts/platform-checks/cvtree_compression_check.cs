#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Globalization;
using System.IO.Compression;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// CVTree compression check (2026-11-04): relatives from amino-acid CVTree 4-mers need a signature per genome, and the full composition
// vector has 20^4 = 160,000 dimensions (about 22 GB of floats for the catalogue, computed). Does a hashed vector keep T2's accuracy?
// Input: genomes the S7.R pass has finished (their 3- and 4-mer counts are in protein_kmer), bacteria with a GTDB v232 placement, the
// largest record >= 1 Mbp, one genome per species, as T2. The 4-mer composition vector (Hao: a(w) = (f(w) - f0(w)) / f0(w), f0 from the
// two 3-mers and the shared 2-mer; 2-mers are the 3-mer counts summed over the last letter) is computed from the stored counts.
// Compressed: feature hashing into d buckets (each word to one bucket with a fixed pseudo-random sign), d = 256 ... 8192, and d = 2048
// rounded to 16-bit halves. Context: CVTree 3-mers (8,000 dimensions) and Karlin's rho* from the stored dinucleotides. Distance (1 - cosine)
// / 2 (Euclidean for Karlin, as S54). Measures as T2: Spearman with GTDB patristic distance per band, and 5-nearest family purity.
// Rule fixed before the run: the smallest d whose within-family Spearman and purity are both within 0.01 of the full vector's.
// Second run (2026-11-05): can the 3-mer table go? Arms full4derived and hashed8192derived take the 3-mers from the 4-mer counts summed
// over the last letter; rule fixed before that run: drop the table if hashed8192derived is within 0.005 of hashed8192 on both measures.
// Args: the folder with GTDB's bac120.tree.gz and bac120_taxonomy.tsv.gz, the output path (default CvTree_Compression.csv), and
// optionally how many species to draw (the finished catalogue has too many for full 4-mer vectors in memory).
const int seed = 20261104;
var gtdb = args[0];
var output = args.Length > 1 ? args[1] : "CvTree_Compression.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var speciesOf = new Dictionary<string, string>();
using (var reader = new StreamReader(new GZipStream(File.OpenRead(Path.Combine(gtdb, "bac120_taxonomy.tsv.gz")), CompressionMode.Decompress)))
    for (string? line; (line = reader.ReadLine()) != null;)
    {
        var parts = line.Split('\t');
        speciesOf[parts[0]] = parts[1][(parts[1].LastIndexOf("s__", StringComparison.Ordinal) + 3)..];
    }
var tree = Tree.Parse(new StreamReader(new GZipStream(File.OpenRead(Path.Combine(gtdb, "bac120.tree.gz")), CompressionMode.Decompress)).ReadToEnd());
var leafOfSpecies = new Dictionary<string, int>();
foreach (var (label, node) in tree.Leaves)
    if (speciesOf.TryGetValue(label, out var species)) leafOfSpecies[species] = node;

var genomes = new List<Genome>();
await using (var cmd = dataSource.CreateCommand(
    "select distinct on (r.source_id) r.id, r.dinucleotide_counts, g.tax_phylum, g.tax_family, g.tax_species, k3.payload, k4.payload " +
    "from sequence_record r join protein_kmer k3 on k3.sequence_id = r.id and k3.k = 3 and k3.status = 'done' " +
    "join protein_kmer k4 on k4.sequence_id = r.id and k4.k = 4 and k4.status = 'done' " +
    "join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind = 'bacterium' and r.length_bp >= 1000000 and r.dinucleotide_counts is not null order by r.source_id, r.length_bp desc"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        if (r.IsDBNull(2) || r.IsDBNull(3) || r.IsDBNull(4)) continue;
        if (!leafOfSpecies.TryGetValue(StripPrefix(r.GetString(4)), out int leaf)) continue;
        genomes.Add(new Genome(r.GetInt64(0), r.GetFieldValue<long[]>(1), r.GetString(2), r.GetString(3), r.GetString(4), leaf,
            ProteinKmers.Decode(r.GetFieldValue<byte[]>(5), 3), ProteinKmers.Decode(r.GetFieldValue<byte[]>(6), 4)));
    }
var rng = new Random(seed);
genomes = genomes.GroupBy(g => g.Species).OrderBy(s => s.Key, StringComparer.Ordinal).Select(s => s.OrderBy(_ => rng.Next()).First()).ToList();
// The finished catalogue has thousands of species; full 4-mer vectors for all of them don't fit in memory.
if (args.Length > 2) genomes = genomes.OrderBy(_ => rng.Next()).Take(int.Parse(args[2], CultureInfo.InvariantCulture)).ToList();
int n = genomes.Count;
Console.WriteLine($"{n} genomes, {genomes.Select(g => g.Family).Distinct().Count()} families, {genomes.Select(g => g.Phylum).Distinct().Count()} phyla.");

var full = genomes.AsParallel().AsOrdered().Select(g => Composition(g.Counts4, g.Counts3, 4)).ToArray();
var signatures = new List<(string Name, int Dimensions, Func<int, int, double> Distance)>
{
    ("full4", full[0].Length, (a, b) => CosineDistance(full[a], full[b])),
};
var three = genomes.AsParallel().AsOrdered().Select(g => Composition(g.Counts3, Marginal(g.Counts3), 3)).ToArray();
signatures.Add(("full3", three[0].Length, (a, b) => CosineDistance(three[a], three[b])));
foreach (var d in new[] { 256, 512, 1024, 2048, 4096, 8192 })
{
    var hashed = full.Select(v => Hash(v, d)).ToArray();
    signatures.Add(($"hashed{d}", d, (a, b) => CosineDistance(hashed[a], hashed[b])));
    if (d == 2048)
    {
        var halves = hashed.Select(v => v.Select(x => (float)(Half)x).ToArray()).ToArray();
        signatures.Add(("hashed2048half", d, (a, b) => CosineDistance(halves[a], halves[b])));
    }
}
// Without the 3-mer table (2026-11-05): 3-mers as the 4-mer counts summed over the last letter, short only of each gene's last 3-mer.
var derived = genomes.AsParallel().AsOrdered().Select(g => Composition(g.Counts4, Marginal(g.Counts4), 4)).ToArray();
signatures.Add(("full4derived", derived[0].Length, (a, b) => CosineDistance(derived[a], derived[b])));
var derivedHashed = derived.Select(v => Hash(v, 8192)).ToArray();
signatures.Add(("hashed8192derived", 8192, (a, b) => CosineDistance(derivedHashed[a], derivedHashed[b])));
var karlin = genomes.Select(g => GenomeSignatures.KarlinRho(g.Dinucleotides)).ToArray();
signatures.Add(("karlin", 16, (a, b) => Math.Sqrt(karlin[a].Zip(karlin[b], (x, y) => (x - y) * (x - y)).Sum())));

var pairs = new List<(int A, int B)>();
for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
        pairs.Add((a, b));
var patristic = pairs.Select(p => tree.Distance(genomes[p.A].Leaf, genomes[p.B].Leaf)).ToArray();
string Band((int A, int B) p) => genomes[p.A].Family == genomes[p.B].Family ? "same family"
    : genomes[p.A].Phylum == genomes[p.B].Phylum ? "same phylum, other family" : "other phyla";
var bands = pairs.Select(Band).ToArray();

using var report = File.CreateText(output);
report.WriteLine("Signature, Dimensions, Band, Pairs, Spearman, Purity5Family");
var summary = new Dictionary<string, (double WithinFamily, double Purity, int Dimensions)>();
foreach (var (name, dimensions, distance) in signatures)
{
    var values = new double[pairs.Count];
    Parallel.For(0, pairs.Count, p => values[p] = distance(pairs[p].A, pairs[p].B));
    var matrix = new double[n, n];
    for (int p = 0; p < pairs.Count; p++)
        matrix[pairs[p].A, pairs[p].B] = matrix[pairs[p].B, pairs[p].A] = values[p];
    double purity = Enumerable.Range(0, n).Average(a => Enumerable.Range(0, n).Where(b => b != a).OrderBy(b => matrix[a, b]).Take(5)
        .Average(b => genomes[b].Family == genomes[a].Family ? 1.0 : 0));
    foreach (var band in new[] { "all", "same family", "same phylum, other family", "other phyla" })
    {
        var chosen = Enumerable.Range(0, pairs.Count).Where(p => band == "all" || bands[p] == band).ToArray();
        double rho = Spearman(chosen.Select(p => patristic[p]).ToList(), chosen.Select(p => values[p]).ToList());
        if (band == "same family") summary[name] = (rho, purity, dimensions);
        var line = string.Create(CultureInfo.InvariantCulture, $"{name}, {dimensions}, {band}, {chosen.Length}, {rho:F4}, {purity:F4}");
        report.WriteLine(line);
        Console.WriteLine(line);
    }
}
var reference = summary["full4"];
var chosenSignature = summary.Where(s => s.Key.StartsWith("hashed", StringComparison.Ordinal) && !s.Key.EndsWith("half", StringComparison.Ordinal)
        && !s.Key.EndsWith("derived", StringComparison.Ordinal)
        && s.Value.WithinFamily >= reference.WithinFamily - 0.01 && s.Value.Purity >= reference.Purity - 0.01)
    .OrderBy(s => s.Value.Dimensions).Select(s => s.Key).FirstOrDefault();
Console.WriteLine($"Rule (full4 within-family {reference.WithinFamily:F4}, purity {reference.Purity:F4}): {chosenSignature ?? "no hashed size qualifies"}");
var (withStored, withDerived) = (summary["hashed8192"], summary["hashed8192derived"]);
bool dropThree = Math.Abs(withDerived.WithinFamily - withStored.WithinFamily) <= 0.005 && Math.Abs(withDerived.Purity - withStored.Purity) <= 0.005;
Console.WriteLine($"3-mer rule (hashed8192derived within 0.005 of hashed8192 on both): {withDerived.WithinFamily:F4} vs {withStored.WithinFamily:F4}, " +
    $"{withDerived.Purity:F4} vs {withStored.Purity:F4}: {(dropThree ? "drop the 3-mer table" : "keep it")}");

static string StripPrefix(string species) => species.StartsWith("s__", StringComparison.Ordinal) ? species[3..] : species;

// Counts of the (k-1)-mers summed over the last letter: the shorter words, short only of each protein's last word.
static int[] Marginal(int[] counts) => Enumerable.Range(0, counts.Length / 20).Select(w => Enumerable.Range(0, 20).Sum(l => counts[20 * w + l])).ToArray();

// Hao's composition vector from k-mer and (k-1)-mer counts; the (k-2)-mers are the (k-1)-mer counts summed over the last letter.
static float[] Composition(int[] words, int[] shorterCounts, int k)
{
    var shortestCounts = Marginal(shorterCounts);
    double total = words.Sum(), shorterTotal = shorterCounts.Sum(), shortestTotal = shortestCounts.Sum();
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

// Feature hashing: word w goes to bucket h(w) mod d with sign from another bit of h(w); the same h for every genome.
static float[] Hash(float[] vector, int d)
{
    var hashed = new float[d];
    for (int w = 0; w < vector.Length; w++)
    {
        ulong h = (ulong)w * 0x9E3779B97F4A7C15UL;
        h ^= h >> 29;
        h *= 0xBF58476D1CE4E5B9UL;
        h ^= h >> 32;
        hashed[(int)(h % (ulong)d)] += (h & (1UL << 63)) == 0 ? vector[w] : -vector[w];
    }
    return hashed;
}

static double CosineDistance(float[] a, float[] b)
{
    double dot = 0, na = 0, nb = 0;
    for (int i = 0; i < a.Length; i++)
    {
        dot += a[i] * (double)b[i];
        na += a[i] * (double)a[i];
        nb += b[i] * (double)b[i];
    }
    return (1 - dot / Math.Sqrt(na * nb)) / 2;
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

record Genome(long Id, long[] Dinucleotides, string Phylum, string Family, string Species, int Leaf, int[] Counts3, int[] Counts4);

// A rooted tree read from Newick, with patristic distances by lowest common ancestor (copied from gtdb_tree_distance.cs).
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
