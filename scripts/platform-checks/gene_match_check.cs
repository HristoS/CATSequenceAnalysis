#:project ../../SequenceAnalise/SequenceAnalise.csproj
#:property PublishAot=false
using System.Collections.Concurrent;
using System.IO.Compression;
using System.Text.RegularExpressions;
using SequenceAnalyses;

// Run from the repo root after fetch_eukaryotes.cs (the genomes in data/genomes/eukaryotes/); no database needed. Release: dotnet run -c Release.
// Gene matching between two related eukaryotes by the Z-curve gene score alone (what the platform stores), as reciprocal best hits
// (RBH) under the platform's candidate rule (>= 500 aa: long genes; else within +-10% of length). Judged by NCBI's protein names:
// pass if >= 70% of RBH pairs with informative names on both sides share the normalized name, in both pairs. Reported with the chance
// rate (a random candidate under the same rule sharing the name) and an amino-acid 4-mer cosine RBH as a sequence-based reference.
var folder = args.Length > 0 ? args[0] : "data/genomes/eukaryotes";
var output = args.Length > 1 ? args[1] : "results/GeneMatch_Check.csv";
var pairs = new[]
{
    ("Podospora", "GCF_035222275.1_CBS112042p_1.0", "GCF_035222375.1_CBS415.72m_1.0"),
    ("Nelumbo", "GCF_055504915.1_China_Antique-T2T", "GCF_056320005.1_American_lotus-T2T"),
};
var rows = new List<string> { "pair,method,length_class,genes_a,genes_b,rbh,rbh_named,rbh_same_name,share_same_name,chance_same_name,rbh_pfam,rbh_same_pfam,share_same_pfam" };
foreach (var (label, a, b) in pairs)
{
    var genesA = Read(Path.Combine(folder, a + "_cds_from_genomic.fna.gz"));
    var genesB = Read(Path.Combine(folder, b + "_cds_from_genomic.fna.gz"));
    Console.WriteLine($"{label}: {genesA.Count} / {genesB.Count} genes >= 50 aa (longest isoform)");
    var profileA = Profile(genesA); var profileB = Profile(genesB);
    var queryA = genesA.Select(g => ZCurveGeneProfile.QueryFromGene(g.Dna)).ToArray();
    var queryB = genesB.Select(g => ZCurveGeneProfile.QueryFromGene(g.Dna)).ToArray();

    var stopwatch = System.Diagnostics.Stopwatch.StartNew();
    var zAB = Best(genesA, genesB, (i, j) => profileB.Score(queryA[i], j));
    var zBA = Best(genesB, genesA, (j, i) => profileA.Score(queryB[j], i));
    Console.WriteLine($"  z-curve scores {stopwatch.Elapsed.TotalSeconds:F0} s");
    stopwatch.Restart();
    var kA = genesA.Select(g => Kmers(g.Protein)).ToArray(); var kB = genesB.Select(g => Kmers(g.Protein)).ToArray();
    var kAB = Best(genesA, genesB, (i, j) => Cosine(kA[i], kB[j]));
    var kBA = Best(genesB, genesA, (j, i) => Cosine(kB[j], kA[i]));
    Console.WriteLine($"  4-mer scores {stopwatch.Elapsed.TotalSeconds:F0} s");

    foreach (var (method, ab, ba) in new[] { ("zcurve", zAB, zBA), ("aa4mer", kAB, kBA) })
        foreach (var lengthClass in new[] { "all", "short", "long" })
        {
            bool In(Gene g) => lengthClass == "all" || (lengthClass == "long") == (g.Residues >= ZCurveGeneProfile.LongResidues);
            var rbh = Enumerable.Range(0, genesA.Count).Where(i => In(genesA[i]) && ab[i] >= 0 && ba[ab[i]] == i).Select(i => (i, j: ab[i])).ToList();
            var named = rbh.Where(p => genesA[p.i].Name is not null && genesB[p.j].Name is not null).ToList();
            int same = named.Count(p => genesA[p.i].Name == genesB[p.j].Name);
            var pfam = rbh.Where(p => genesA[p.i].Pfam is not null && genesB[p.j].Pfam is not null).ToList();
            int samePfam = pfam.Count(p => genesA[p.i].Pfam == genesB[p.j].Pfam);
            double chance = Chance(genesA, genesB, In);
            var row = $"{label},{method},{lengthClass},{genesA.Count(In)},{genesB.Count(In)},{rbh.Count},{named.Count},{same},{Share(same, named.Count)},{chance:F4},{pfam.Count},{samePfam},{Share(samePfam, pfam.Count)}";
            rows.Add(row);
            Console.WriteLine("  " + row);
        }
    // A few Z-curve RBH pairs whose names differ, to see what the misses are.
    var misses = Enumerable.Range(0, genesA.Count).Where(i => zAB[i] >= 0 && zBA[zAB[i]] == i && genesA[i].Name is not null
        && genesB[zAB[i]].Name is not null && genesA[i].Name != genesB[zAB[i]].Name).Take(8);
    foreach (var i in misses) Console.WriteLine($"    miss {genesA[i].Residues} aa '{genesA[i].Name}' ~ {genesB[zAB[i]].Residues} aa '{genesB[zAB[i]].Name}'");
}
File.WriteAllLines(output, rows);
Console.WriteLine($"wrote {output}");

static string Share(int part, int whole) => whole == 0 ? "" : (part / (double)whole).ToString("F4");

static ZCurveGeneProfile Profile(List<Gene> genes) =>
    ZCurveGeneProfile.FromCodingSequences(genes.Select((g, i) => (0, i * 10, i * 10 + g.Dna.Length, false, i, g.Dna)));

// For each gene of `from`, the index of the best-scoring candidate in `to` under the platform's rule, or -1.
static int[] Best(List<Gene> from, List<Gene> to, Func<int, int, double> score)
{
    var order = Enumerable.Range(0, to.Count).OrderBy(j => to[j].Residues).ToArray();
    var lengths = order.Select(j => to[j].Residues).ToArray();
    var best = new int[from.Count];
    Parallel.For(0, from.Count, i =>
    {
        int r = from[i].Residues;
        (int lo, int hi) = r >= ZCurveGeneProfile.LongResidues ? (ZCurveGeneProfile.LongResidues, int.MaxValue)
            : ((int)Math.Ceiling(r * 0.9), Math.Min(ZCurveGeneProfile.LongResidues - 1, (int)Math.Floor(r * 1.1)));
        int start = LowerBound(lengths, lo);
        double top = double.NegativeInfinity; int at = -1;
        for (int k = start; k < order.Length && lengths[k] <= hi; k++)
        {
            double s = score(i, order[k]);
            if (s > top) { top = s; at = order[k]; }
        }
        best[i] = at;
    });
    return best;
}

static int LowerBound(int[] sorted, int value)
{
    int lo = 0, hi = sorted.Length;
    while (lo < hi) { int mid = (lo + hi) / 2; if (sorted[mid] < value) lo = mid + 1; else hi = mid; }
    return lo;
}

// Mean over named genes of A of the share of its candidates in B with the same name.
static double Chance(List<Gene> a, List<Gene> b, Func<Gene, bool> include)
{
    var order = b.OrderBy(g => g.Residues).ToArray();
    var lengths = order.Select(g => g.Residues).ToArray();
    var shares = new ConcurrentBag<double>();
    Parallel.ForEach(a.Where(g => g.Name is not null && include(g)), g =>
    {
        int r = g.Residues;
        (int lo, int hi) = r >= ZCurveGeneProfile.LongResidues ? (ZCurveGeneProfile.LongResidues, int.MaxValue)
            : ((int)Math.Ceiling(r * 0.9), Math.Min(ZCurveGeneProfile.LongResidues - 1, (int)Math.Floor(r * 1.1)));
        int total = 0, same = 0;
        for (int k = LowerBound(lengths, lo); k < order.Length && lengths[k] <= hi; k++)
        {
            if (order[k].Name is null) continue;
            total++;
            if (order[k].Name == g.Name) same++;
        }
        if (total > 0) shares.Add(same / (double)total);
    });
    return shares.Count == 0 ? double.NaN : shares.Average();
}

static (int[] Codes, double[] Weights) Kmers(string protein)
{
    const string letters = "ACDEFGHIKLMNPQRSTVWY";
    var counts = new Dictionary<int, int>();
    for (int p = 0; p + 4 <= protein.Length; p++)
    {
        int code = 0; bool ok = true;
        for (int q = 0; q < 4 && ok; q++) { int l = letters.IndexOf(protein[p + q]); if (l < 0) ok = false; else code = code * 20 + l; }
        if (ok) counts[code] = counts.GetValueOrDefault(code) + 1;
    }
    var codes = counts.Keys.Order().ToArray();
    var weights = codes.Select(c => (double)counts[c]).ToArray();
    double norm = Math.Sqrt(weights.Sum(w => w * w));
    if (norm > 0) for (int k = 0; k < weights.Length; k++) weights[k] /= norm;
    return (codes, weights);
}

static double Cosine((int[] Codes, double[] Weights) a, (int[] Codes, double[] Weights) b)
{
    double dot = 0; int i = 0, j = 0;
    while (i < a.Codes.Length && j < b.Codes.Length)
    {
        if (a.Codes[i] == b.Codes[j]) dot += a.Weights[i++] * b.Weights[j++];
        else if (a.Codes[i] < b.Codes[j]) i++; else j++;
    }
    return dot;
}

static List<Gene> Read(string path)
{
    var attribute = new Regex(@"\[(\w+)=([^\]]*)\]");
    var best = new Dictionary<string, Gene>();
    using var reader = new StreamReader(new GZipStream(File.OpenRead(path), CompressionMode.Decompress), bufferSize: 1 << 20);
    string? header = null; var dna = new System.Text.StringBuilder();
    void Flush()
    {
        if (header is null) return;
        var fields = new Dictionary<string, string>();
        foreach (Match m in attribute.Matches(header)) fields.TryAdd(m.Groups[1].Value, m.Groups[2].Value);
        if (fields.ContainsKey("pseudo") || fields.ContainsKey("partial") || fields.ContainsKey("exception")) return;
        var sequence = dna.ToString().ToUpperInvariant();
        if (sequence.Length % 3 != 0) return;
        var protein = ZCurveGeneProfile.Translate(sequence);
        if (protein.Length < ZCurveGeneProfile.MinResidues) return;
        var xref = fields.GetValueOrDefault("db_xref", "");
        var key = Regex.Match(xref, @"GeneID:(\d+)") is { Success: true } g ? g.Value : fields.GetValueOrDefault("locus_tag", header);
        var pfam = string.Join("+", Regex.Matches(xref, @"PFAM:(PF\d+)").Select(m => m.Groups[1].Value).Distinct().Order());
        var gene = new Gene(sequence, protein, protein.Length, Normalize(fields.GetValueOrDefault("protein")), pfam.Length == 0 ? null : pfam);
        if (!best.TryGetValue(key, out var other) || other.Dna.Length < sequence.Length) best[key] = gene;
    }
    for (string? line; (line = reader.ReadLine()) is not null;)
    {
        if (line.StartsWith('>')) { Flush(); header = line; dna.Clear(); }
        else dna.Append(line.Trim());
    }
    Flush();
    return best.Values.ToList();
}

static string? Normalize(string? name)
{
    if (name is null) return null;
    var n = name.ToLowerInvariant();
    n = n.Replace("low quality protein: ", "");
    n = Regex.Replace(n, @"\s+isoform\s+\S+$", "");
    n = n.Replace("-like", "").Replace("putative ", "").Replace("probable ", "");
    n = Regex.Replace(n, @"\s+homolog$", "");
    n = Regex.Replace(n, @"\s+", " ").Trim();
    return n.Length == 0 || n.StartsWith("uncharacterized") || n.StartsWith("hypothetical") ? null : n;
}

record Gene(string Dna, string Protein, int Residues, string? Name, string? Pfam);
