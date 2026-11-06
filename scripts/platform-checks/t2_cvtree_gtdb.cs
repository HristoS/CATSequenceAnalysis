#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.IO.Compression;
using System.Text;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// P7 T2 (2026-11-04): does the amino-acid CVTree signature, which beat Karlin on the bundled trees (R6e_Trees_Bundled.csv), also track
// GTDB's tree better than Karlin on catalogue genomes? Asked before the re-download, because it would need gene calls and amino-acid
// k-mer counts at ingest. Sample: bacterial GTDB families with >= perFamily species whose largest stored record is >= 1 Mbp (one
// genome per species, as gtdb_tree_distance.cs), `families` of them drawn at random, perFamily species each; Mycoplasma-like orders
// (genetic code 4) are left out. Each record is streamed from NCBI efetch and held in memory only while its signatures are computed:
// no file, nothing kept but the signatures. The stored dinucleotide counts are compared with the streamed ones, so a wrong download shows.
// Genes as in GeneLevelSignatureBenchmark: ATG-to-stop ORFs >= 300 bp on both strands, longest first, dropping overlaps > 60 bp.
// Distances: karlin (rho*, Euclidean, as S54), gc, aminoAcids (standardized Euclidean), cvTree3 and cvTree4 (Hao's Markov-corrected
// amino-acid k-mers, (1 - cosine) / 2), dnaCv6 (the same on both-strand DNA 6-mers, which the platform stores). Measures: Spearman with
// GTDB patristic distance per band (all, same family, same phylum but another family, other phyla); 95% CI of each minus karlin within
// a family from resampling families; 5-nearest family purity as context.
// Rule fixed before the run: pass if cvTree3 or cvTree4 beats karlin's within-family Spearman with the CI's lower bound above 0; else stop.
// Args: the folder with GTDB's bac120.tree.gz and bac120_taxonomy.tsv.gz (v232, not kept here), and the output path (default T2_CvTreeGtdb.csv).
const int families = 40, perFamily = 8, minLengthBp = 1_000_000, minGeneLength = 300, maxOverlap = 60, resamples = 1000, seed = 20261104;
const string aminoAcidLetters = "ACDEFGHIKLMNPQRSTVWY";
// NCBI table 11 in TCAG order (index = 16 * first + 4 * second + third).
const string code11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
var gtdb = args[0];
var output = args.Length > 1 ? args[1] : "T2_CvTreeGtdb.csv";
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

var candidates = new List<Candidate>();
await using (var cmd = dataSource.CreateCommand(
    "select distinct on (r.source_id) r.id, r.record_id, r.length_bp, r.dinucleotide_counts, g.tax_phylum, g.tax_order, g.tax_family, g.tax_species " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind = 'bacterium' and r.dinucleotide_counts is not null order by r.source_id, r.length_bp desc, r.id"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        if (r.IsDBNull(4) || r.IsDBNull(5) || r.IsDBNull(6) || r.IsDBNull(7) || r.GetInt64(2) < minLengthBp) continue;
        var order = r.GetString(5);
        if (order.Contains("Mycoplasm", StringComparison.Ordinal) || order.Contains("Entomoplasm", StringComparison.Ordinal)) continue;
        if (!leafOfSpecies.TryGetValue(StripPrefix(r.GetString(7)), out int leaf)) continue;
        candidates.Add(new Candidate(r.GetInt64(0), r.GetString(1), r.GetInt64(2), r.GetFieldValue<long[]>(3), r.GetString(4), r.GetString(6), r.GetString(7), leaf));
    }
var rng = new Random(seed);
var oneEach = candidates.GroupBy(c => c.Species).OrderBy(s => s.Key, StringComparer.Ordinal).Select(s => s.OrderBy(_ => rng.Next()).First()).ToList();
var eligible = oneEach.GroupBy(c => c.Family).Where(f => f.Count() >= perFamily).OrderBy(f => f.Key, StringComparer.Ordinal).ToList();
var sample = eligible.OrderBy(_ => rng.Next()).Take(families)
    .SelectMany(f => f.OrderBy(_ => rng.Next()).Take(perFamily)).ToList();
Console.WriteLine($"{candidates.Count} placed records, {oneEach.Count} species, {eligible.Count} families with >= {perFamily}; " +
    $"sampling {sample.Count} genomes from {sample.Select(s => s.Family).Distinct().Count()} families, {sample.Select(s => s.Phylum).Distinct().Count()} phyla.");

using var http = new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
var watch = Stopwatch.StartNew();
long bases = 0;
int countMismatches = 0;
var profiles = new List<Profile>();
foreach (var (genome, index) in sample.Select((g, i) => (g, i)))
{
    var sequence = await FetchAsync(genome.Accession);
    bases += sequence.Length;
    if (!GenomeSignatures.DinucleotideCounts(sequence).SequenceEqual(genome.StoredDinucleotides)) countMismatches++;
    var reverse = ReverseComplement(sequence);
    var proteins = CallGenes(sequence, reverse).Select(Translate).ToList();
    profiles.Add(new Profile(genome,
        GenomeSignatures.KarlinRho(GenomeSignatures.DinucleotideCounts(sequence)),
        sequence.Count(c => c is 'G' or 'C') / (double)sequence.Count(c => c is 'A' or 'C' or 'G' or 'T'),
        AminoAcidFrequencies(proteins),
        CvTreeVector(proteins, 3, aminoAcidLetters),
        CvTreeVector(proteins, 4, aminoAcidLetters),
        CvTreeVector(new[] { sequence, reverse }, 6, "ACGT"),
        proteins.Count));
    if ((index + 1) % 20 == 0)
        Console.WriteLine($"  {index + 1}/{sample.Count}, {watch.Elapsed.TotalMinutes:F1} min, {bases / 1e6:F0} Mbp, median genes {profiles.Select(p => p.Genes).Order().ElementAt(profiles.Count / 2)}");
}
Console.WriteLine($"Streamed {bases / 1e6:F0} Mbp in {watch.Elapsed.TotalMinutes:F1} min; streamed dinucleotides differ from the stored ones for {countMismatches} of {sample.Count}.");

var aminoAcids = Standardize(profiles.Select(p => p.AminoAcids).ToArray());
var names = new[] { "karlin", "gc", "aminoAcids", "cvTree3", "cvTree4", "dnaCv6" };
Func<int, int, double>[] distance =
{
    (a, b) => Euclidean(profiles[a].Karlin, profiles[b].Karlin),
    (a, b) => Math.Abs(profiles[a].Gc - profiles[b].Gc),
    (a, b) => Euclidean(aminoAcids[a], aminoAcids[b]),
    (a, b) => CosineDistance(profiles[a].CvTree3, profiles[b].CvTree3),
    (a, b) => CosineDistance(profiles[a].CvTree4, profiles[b].CvTree4),
    (a, b) => CosineDistance(profiles[a].DnaCv6, profiles[b].DnaCv6),
};
int n = profiles.Count;
var pairs = new List<(int A, int B)>();
for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
        pairs.Add((a, b));
var matrix = new double[names.Length][];
for (int d = 0; d < names.Length; d++)
{
    var values = new double[pairs.Count];
    int column = d;
    Parallel.For(0, pairs.Count, p => values[p] = distance[column](pairs[p].A, pairs[p].B));
    matrix[d] = values;
}
var patristic = pairs.Select(p => tree.Distance(profiles[p.A].Genome.Leaf, profiles[p.B].Genome.Leaf)).ToArray();
string Band((int A, int B) p) =>
    profiles[p.A].Genome.Family == profiles[p.B].Genome.Family ? "same family"
    : profiles[p.A].Genome.Phylum == profiles[p.B].Genome.Phylum ? "same phylum, other family" : "other phyla";
var bands = pairs.Select(Band).ToArray();

using var report = File.CreateText(output);
report.WriteLine("Measure, Band, Pairs, Distance, Value, MinusKarlin, MinusKarlinCiLow, MinusKarlinCiHigh");
void Write(string measure, string band, int count, int d, double value, double minusKarlin, string ciLow, string ciHigh)
{
    var line = string.Create(CultureInfo.InvariantCulture, $"{measure}, {band}, {count}, {names[d]}, {value:F4}, {minusKarlin:F4}, {ciLow}, {ciHigh}");
    report.WriteLine(line);
    Console.WriteLine(line);
}

var familyNames = profiles.Select(p => p.Genome.Family).Distinct().OrderBy(f => f, StringComparer.Ordinal).ToArray();
var pairsOfFamily = familyNames.ToDictionary(f => f, f => Enumerable.Range(0, pairs.Count)
    .Where(p => bands[p] == "same family" && profiles[pairs[p].A].Genome.Family == f).ToArray());
var draws = new double[resamples][];
var drawRng = new Random(seed + 1);
var drawFamilies = Enumerable.Range(0, resamples).Select(_ => Enumerable.Range(0, familyNames.Length).Select(_ => familyNames[drawRng.Next(familyNames.Length)]).ToArray()).ToArray();
Parallel.For(0, resamples, r =>
{
    var chosen = drawFamilies[r].SelectMany(f => pairsOfFamily[f]).ToArray();
    var treeDistances = chosen.Select(p => patristic[p]).ToList();
    draws[r] = Enumerable.Range(0, names.Length).Select(d => Spearman(treeDistances, chosen.Select(p => matrix[d][p]).ToList())).ToArray();
});

var withinFamily = new double[names.Length];
foreach (var band in new[] { "all", "same family", "same phylum, other family", "other phyla" })
{
    var chosen = Enumerable.Range(0, pairs.Count).Where(p => band == "all" || bands[p] == band).ToArray();
    var treeDistances = chosen.Select(p => patristic[p]).ToList();
    var rho = Enumerable.Range(0, names.Length).Select(d => Spearman(treeDistances, chosen.Select(p => matrix[d][p]).ToList())).ToArray();
    for (int d = 0; d < names.Length; d++)
    {
        string low = "", high = "";
        if (band == "same family")
        {
            withinFamily[d] = rho[d];
            var differences = draws.Select(x => x[d] - x[0]).Order().ToArray();
            low = differences[(int)(0.025 * (resamples - 1))].ToString("F4", CultureInfo.InvariantCulture);
            high = differences[(int)(0.975 * (resamples - 1))].ToString("F4", CultureInfo.InvariantCulture);
        }
        Write("spearman", band, chosen.Length, d, rho[d], rho[d] - rho[0], low, high);
    }
}

// 5-nearest family purity: the share of each genome's 5 nearest (by the distance) that share its family; chance is (perFamily - 1) / (n - 1).
var purity = new double[names.Length];
for (int d = 0; d < names.Length; d++)
{
    var full = new double[n, n];
    for (int p = 0; p < pairs.Count; p++)
        full[pairs[p].A, pairs[p].B] = full[pairs[p].B, pairs[p].A] = matrix[d][p];
    purity[d] = Enumerable.Range(0, n).Average(a => Enumerable.Range(0, n).Where(b => b != a).OrderBy(b => full[a, b]).Take(5)
        .Average(b => profiles[b].Genome.Family == profiles[a].Genome.Family ? 1.0 : 0));
    Write("purity5", "family", n, d, purity[d], purity[d] - purity[0], "", "");
}

var passing = new[] { 3, 4 }.Where(d => draws.Select(x => x[d] - x[0]).Order().ElementAt((int)(0.025 * (resamples - 1))) > 0).Select(d => names[d]).ToList();
Console.WriteLine($"Rule (karlin within-family Spearman {withinFamily[0]:F4}): " +
    (passing.Count > 0 ? $"pass for {string.Join(", ", passing)}" : "stop, no cvTree beats karlin within a family with the CI above 0"));

async Task<string> FetchAsync(string accession)
{
    for (int attempt = 1; ; attempt++)
    {
        try
        {
            var url = $"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&tool=cat-platform-t2&id={accession}";
            await using var stream = await http.GetStreamAsync(url);
            using var reader = new StreamReader(stream);
            var sequence = new StringBuilder();
            for (string? line; (line = await reader.ReadLineAsync()) != null;)
                if (!line.StartsWith('>')) sequence.Append(line.Trim().ToUpperInvariant());
            await Task.Delay(400);   // NCBI asks for at most 3 requests a second without a key
            return sequence.ToString();
        }
        catch (Exception e) when (attempt < 4)
        {
            Console.WriteLine($"  retry {accession}: {e.Message}");
            await Task.Delay(2000 * attempt);
        }
    }
}

static string StripPrefix(string species) => species.StartsWith("s__", StringComparison.Ordinal) ? species[3..] : species;

// Longest first, so at a locus the real gene usually wins over shorter antisense or nested ORFs. Returned in coding orientation.
static List<string> CallGenes(string forward, string reverse)
{
    int length = forward.Length;
    var candidates = AtgGenes(forward).Select(g => (g.Start, g.End, Reverse: false))
        .Concat(AtgGenes(reverse).Select(g => (Start: length - g.End, End: length - g.Start, Reverse: true)))
        .OrderByDescending(g => g.End - g.Start)
        .ToList();
    var covered = new bool[length];
    var genes = new List<string>();
    foreach (var (start, end, isReverse) in candidates)
    {
        int overlap = 0;
        for (int i = start; i < end; i++)
            if (covered[i]) overlap++;
        if (overlap > maxOverlap) continue;
        for (int i = start; i < end; i++)
            covered[i] = true;
        genes.Add(isReverse ? reverse.Substring(length - end, end - start) : forward.Substring(start, end - start));
    }
    return genes;
}

// A stop-free stretch per frame (a non-ACGT codon ends one like a stop), trimmed to its first in-frame ATG.
static IEnumerable<(int Start, int End)> AtgGenes(string strand)
{
    for (int frame = 0; frame < 3; frame++)
    {
        int start = frame, i = frame;
        for (; i + 3 <= strand.Length; i += 3)
        {
            int codon = CodonIndex(strand, i);
            if (codon >= 0 && code11[codon] != '*') continue;
            if (Trimmed(strand, start, i) is { } gene) yield return gene;
            start = i + 3;
        }
        if (Trimmed(strand, start, i) is { } last) yield return last;
    }
}

static (int Start, int End)? Trimmed(string strand, int start, int end)
{
    for (int i = start; i + minGeneLength <= end; i += 3)
        if (strand[i] == 'A' && strand[i + 1] == 'T' && strand[i + 2] == 'G') return (i, end);
    return null;
}

static int CodonIndex(string sequence, int i)
{
    int first = "TCAG".IndexOf(sequence[i]), second = "TCAG".IndexOf(sequence[i + 1]), third = "TCAG".IndexOf(sequence[i + 2]);
    return first < 0 || second < 0 || third < 0 ? -1 : 16 * first + 4 * second + third;
}

static string Translate(string gene)
{
    var protein = new char[gene.Length / 3];
    for (int i = 0; i < protein.Length; i++)
        protein[i] = code11[CodonIndex(gene, 3 * i)];
    return new string(protein);
}

static string ReverseComplement(string sequence)
{
    var result = new char[sequence.Length];
    for (int i = 0; i < sequence.Length; i++)
        result[sequence.Length - 1 - i] = sequence[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
    return new string(result);
}

static double[] AminoAcidFrequencies(List<string> proteins)
{
    var counts = new double[aminoAcidLetters.Length];
    foreach (var protein in proteins)
        foreach (var aminoAcid in protein)
            counts[aminoAcidLetters.IndexOf(aminoAcid)]++;
    var total = counts.Sum();
    return counts.Select(c => c / total).ToArray();
}

// Hao's composition vector: a(w) = (f(w) - f0(w)) / f0(w), f0 predicting w from its two (k-1)-mers and their shared (k-2)-mer.
static float[] CvTreeVector(IReadOnlyList<string> texts, int k, string alphabet)
{
    var words = KmerFrequencies(texts, k, alphabet);
    var shorter = KmerFrequencies(texts, k - 1, alphabet);
    var shortest = KmerFrequencies(texts, k - 2, alphabet);
    var vector = new float[words.Length];
    for (int w = 0; w < words.Length; w++)
    {
        int prefix = w / alphabet.Length, suffix = w % shorter.Length, middle = prefix % shortest.Length;
        double expected = shortest[middle] > 0 ? shorter[prefix] * shorter[suffix] / shortest[middle] : 0;
        vector[w] = expected > 0 ? (float)((words[w] - expected) / expected) : 0;
    }
    return vector;
}

static double[] KmerFrequencies(IReadOnlyList<string> texts, int k, string alphabet)
{
    int space = 1;
    for (int i = 0; i < k; i++)
        space *= alphabet.Length;
    var lookup = new int[128];
    Array.Fill(lookup, -1);
    for (int i = 0; i < alphabet.Length; i++)
        lookup[alphabet[i]] = i;
    var counts = new double[space];
    long total = 0;
    foreach (var text in texts)
        for (int start = 0; start + k <= text.Length; start++)
        {
            int word = 0;
            for (int i = 0; i < k && word >= 0; i++)
            {
                int letter = text[start + i] < 128 ? lookup[text[start + i]] : -1;
                word = letter < 0 ? -1 : word * alphabet.Length + letter;
            }
            if (word < 0) continue;
            counts[word]++;
            total++;
        }
    for (int w = 0; w < space; w++)
        counts[w] /= total;
    return counts;
}

static double CosineDistance(float[] a, float[] b)
{
    double dot = 0, normA = 0, normB = 0;
    for (int i = 0; i < a.Length; i++)
    {
        dot += a[i] * (double)b[i];
        normA += a[i] * (double)a[i];
        normB += b[i] * (double)b[i];
    }
    return (1 - dot / Math.Sqrt(normA * normB)) / 2;
}

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

record Candidate(long Id, string Accession, long Length, long[] StoredDinucleotides, string Phylum, string Family, string Species, int Leaf);

record Profile(Candidate Genome, double[] Karlin, double Gc, double[] AminoAcids, float[] CvTree3, float[] CvTree4, float[] DnaCv6, int Genes);

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
