#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.Text;
using Npgsql;
using SequencePlatform.Core;

// P7 T3 (2026-11-04): protein search with the chosen per-CDS profile on a pool 10x the bundled one. The bundled checks measured 30
// genomes and extrapolated 510x to the catalogue with a model that over-predicted by 4-5 points over 8 -> 30 genomes; this measures
// 300 and shortens the extrapolation to about 51x. Sample: `genomes` bacterial species drawn at random, one genome each (as T2), whose
// largest stored record is >= 1 Mbp, Mycoplasma-like orders left out (code 4). Each record is streamed from NCBI efetch and held in
// memory only while it is profiled: no file, nothing kept but gene profiles, proteins (for the homolog count) and codon usage.
// Genes and profile as R6f and the CDS sweep's chosen cell: stop-free stretches from the first in-frame ATG, >= 1,500 bp, both strands;
// the first 1,500 bp, codon positions 1 and 2, windows of 333 such bases, Legendre moments P0-P2 over x, y, z scaled by sqrt(2k + 1)
// (27 numbers per gene); distance = mean over windows of the Euclidean distance. Queries: `queriesPerGenome` genes per genome; from the
// protein with uniform codon weights (ExpectedUniform) or each candidate genome's own codon usage (ExpectedPerTarget); Dna as control.
// Measured: strict top-1 (the query's own gene first) and homolog top-1 (the first hit shares >= 0.3 of the query's amino-acid 5-mers)
// at 30, 100 and 300 genomes. Model (R6f's, corrected as in the CDS sweep): fitted on the query's genome plus 29 others, predict 300.
// Rule fixed before the run: the 30 -> 300 prediction is within 5 points of the measured strict top-1, and the extrapolation from 300 to
// catalogueGenomes is >= 50% for ExpectedUniform or ExpectedPerTarget. Args: output path (default T3_Pool300.csv).
const int genomes = 300, queriesPerGenome = 2, geneLength = 1500, window = 333, catalogueGenomes = 15_314, tailNeighbours = 100, seed = 20261104;
const int codonUsageMinLength = 900;
const double orthologMinShare = 0.3, modelTolerance = 0.05;
const string code11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
var output = args.FirstOrDefault() ?? "T3_Pool300.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var candidates = new List<(string Accession, string Species)>();
await using (var cmd = dataSource.CreateCommand(
    "select distinct on (r.source_id) r.record_id, r.length_bp, g.tax_order, g.tax_species " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind = 'bacterium' order by r.source_id, r.length_bp desc, r.id"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        if (r.IsDBNull(2) || r.IsDBNull(3) || r.GetInt64(1) < 1_000_000) continue;
        var order = r.GetString(2);
        if (order.Contains("Mycoplasm", StringComparison.Ordinal) || order.Contains("Entomoplasm", StringComparison.Ordinal)) continue;
        candidates.Add((r.GetString(0), r.GetString(3)));
    }
var rng = new Random(seed);
var sample = candidates.GroupBy(c => c.Species).OrderBy(s => s.Key, StringComparer.Ordinal).Select(s => s.OrderBy(_ => rng.Next()).First())
    .OrderBy(_ => rng.Next()).Take(genomes).ToList();
Console.WriteLine($"{candidates.Count} records, sampling {sample.Count} species.");

var synonyms = Enumerable.Range(0, 64).Where(i => code11[i] != '*').GroupBy(i => code11[i]).ToDictionary(g => g.Key, g => g.ToArray());
var uniformTable = ExpectedTable(Enumerable.Repeat(1L, 64).ToArray());
var genes = new List<(int Genome, string Protein, float[] Features)>();
var tables = new List<(double X, double Y, double Z)[]>();
var queries = new List<int>();
using var http = new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
var watch = Stopwatch.StartNew();
long bases = 0;
for (int genome = 0; genome < sample.Count; genome++)
{
    var forward = await FetchAsync(sample[genome].Accession);
    bases += forward.Length;
    var reverse = ReverseComplement(forward);
    var usage = new long[64];
    foreach (var strand in new[] { forward, reverse })
        foreach (var (start, end) in Stretches(strand, codonUsageMinLength))
            for (int i = start; i < end; i += 3)
                usage[CodonIndex(strand, i)]++;
    tables.Add(ExpectedTable(usage));
    int first = genes.Count;
    foreach (var strand in new[] { forward, reverse })
        foreach (var (start, end) in Stretches(strand, geneLength))
            for (int i = start; i + geneLength <= end; i += 3)
                if (strand[i] == 'A' && strand[i + 1] == 'T' && strand[i + 2] == 'G')
                {
                    var gene = strand.Substring(i, geneLength);
                    var protein = new string(Enumerable.Range(0, geneLength / 3).Select(c => code11[CodonIndex(gene, 3 * c)]).ToArray());
                    genes.Add((genome, protein, Features(p => DnaAxes(gene[3 * (p / 2) + p % 2]))));
                    break;
                }
    var picker = new Random(seed + genome);
    queries.AddRange(Enumerable.Range(first, genes.Count - first).OrderBy(_ => picker.Next()).Take(queriesPerGenome));
    if ((genome + 1) % 25 == 0)
        Console.WriteLine($"  {genome + 1}/{sample.Count}, {watch.Elapsed.TotalMinutes:F1} min, {bases / 1e6:F0} Mbp, {genes.Count} genes");
}
Console.WriteLine($"Streamed {bases / 1e6:F0} Mbp in {watch.Elapsed.TotalMinutes:F1} min: {genes.Count} genes, {queries.Count} queries.");

var poolSizes = new[] { 30, 100, 300 };
using var report = File.CreateText(output);
report.WriteLine("Variant,Queries,Pool," + string.Join(",", poolSizes.Select(s => $"StrictTop1At{s}")) + "," +
    string.Join(",", poolSizes.Select(s => $"HomologTop1At{s}")) + ",ModelFrom30AtFull,ModelValid,ModelTop1AtCatalogue");
var verdicts = new List<(string Variant, bool Valid, double Catalogue)>();
foreach (var variant in new[] { "Dna", "ExpectedUniform", "ExpectedPerTarget" })
{
    var strict = new double[queries.Count, poolSizes.Length];
    var homolog = new double[queries.Count, poolSizes.Length];
    var from30 = new double[queries.Count];
    var catalogue = new double[queries.Count];
    Parallel.For(0, queries.Count, j =>
    {
        int truth = queries[j];
        int own = genes[truth].Genome;
        var byGenome = new float[sample.Count][];
        float[] Query(int genome) => byGenome[variant == "ExpectedPerTarget" ? genome : 0] ??= variant switch
        {
            "Dna" => genes[truth].Features,
            "ExpectedUniform" => Features(p => Expected(genes[truth].Protein, uniformTable, p)),
            _ => Features(p => Expected(genes[truth].Protein, tables[genome], p)),
        };
        double trueDistance = Distance(Query(own), genes[truth].Features);
        var rngQuery = new Random(seed ^ truth);
        var order = Enumerable.Range(0, sample.Count).Where(g => g != own).OrderBy(_ => rngQuery.Next()).ToArray();
        var rankOfGenome = new int[sample.Count];
        for (int r = 0; r < order.Length; r++) rankOfGenome[order[r]] = r + 1;   // own genome is rank 0: always in the pool
        // For each pool size s, the genomes in the pool are own + the first s - 1 of `order`.
        int closerOwn = 0;
        var closest = Enumerable.Repeat((Distance: double.PositiveInfinity, Gene: -1), poolSizes.Length).ToArray();
        var others30 = new List<double>();
        var othersAll = new List<double>();
        for (int i = 0; i < genes.Count; i++)
        {
            if (i == truth) continue;
            int genome = genes[i].Genome;
            double distance = Distance(Query(genome), genes[i].Features);
            if (genome == own)
            {
                if (distance < trueDistance) closerOwn++;
                for (int s = 0; s < poolSizes.Length; s++)
                    if (distance < closest[s].Distance) closest[s] = (distance, i);
                continue;
            }
            othersAll.Add(distance);
            if (rankOfGenome[genome] < 30) others30.Add(distance);
            for (int s = 0; s < poolSizes.Length; s++)
                if (rankOfGenome[genome] < poolSizes[s] && distance < closest[s].Distance) closest[s] = (distance, i);
        }
        for (int s = 0; s < poolSizes.Length; s++)
        {
            bool top = closest[s].Distance >= trueDistance;
            strict[j, s] = top ? 1 : 0;
            homolog[j, s] = top || KmerShare(genes[truth].Protein, genes[closest[s].Gene].Protein) >= orthologMinShare ? 1 : 0;
        }
        var sorted30 = others30.Order().ToArray();
        var sortedAll = othersAll.Order().ToArray();
        from30[j] = Survive(closerOwn, sorted30, trueDistance, sortedAll.Length);
        catalogue[j] = Survive(closerOwn, sortedAll, trueDistance, catalogueGenomes / (double)sample.Count * genes.Count);
    });
    double Column(double[,] values, int s) => Enumerable.Range(0, queries.Count).Average(j => values[j, s]);
    double measured = Column(strict, poolSizes.Length - 1), predicted = from30.Average();
    bool valid = Math.Abs(predicted - measured) <= modelTolerance;
    var line = string.Join(",", new[] { variant, queries.Count.ToString(), genes.Count.ToString() }
        .Concat(poolSizes.Select((_, s) => Format(Column(strict, s))))
        .Concat(poolSizes.Select((_, s) => Format(Column(homolog, s))))
        .Concat(new[] { Format(predicted), valid.ToString().ToLowerInvariant(), Format(catalogue.Average()) }));
    report.WriteLine(line);
    Console.WriteLine(line);
    verdicts.Add((variant, valid, catalogue.Average()));
}
var passing = verdicts.Where(v => v.Variant != "Dna" && v.Valid && v.Catalogue >= 0.5).Select(v => v.Variant).ToList();
Console.WriteLine($"Rule: {(passing.Count > 0 ? "pass for " + string.Join(", ", passing) : "stop")}");

// The query's own genome first: if any of its genes is closer, a miss at every scale; else the measured share of competitors closer
// when there is one, and only otherwise the tail extrapolation below the nearest competitor (Levina-Bickel local dimension of the k nearest).
static double Survive(int closerOwn, double[] sorted, double trueDistance, double competitors)
{
    if (closerOwn > 0) return 0;
    int closer = 0;
    while (closer < sorted.Length && sorted[closer] < trueDistance) closer++;
    double share;
    if (closer > 0) share = closer / (double)sorted.Length;
    else
    {
        int k = Math.Min(tailNeighbours, sorted.Length);
        double kth = Math.Max(sorted[k - 1], 1e-12), logSum = 0;
        for (int i = 0; i < k - 1; i++) logSum += Math.Log(kth / Math.Max(sorted[i], 1e-12));
        share = k / (double)sorted.Length * Math.Pow(Math.Max(trueDistance, 0) / kth, (k - 1) / logSum);
    }
    return Math.Exp(-competitors * share);
}

// 3 windows x 9 numbers over the first 1,000 positions-1-and-2 bases, laid out window by window.
static float[] Features(Func<int, (double X, double Y, double Z)> axesAt)
{
    const int windows = geneLength / 3 * 2 / window;
    var features = new float[windows * 9];
    for (int w = 0; w < windows; w++)
    {
        var sums = new double[9];
        for (int i = 0; i < window; i++)
        {
            var (x, y, z) = axesAt(w * window + i);
            double t = -1 + 2 * (i + 0.5) / window;
            for (int k = 0; k < 3; k++)
            {
                double weight = k switch { 0 => 1, 1 => t, _ => (3 * t * t - 1) / 2 };
                sums[3 * k] += x * weight;
                sums[3 * k + 1] += y * weight;
                sums[3 * k + 2] += z * weight;
            }
        }
        for (int f = 0; f < 9; f++)
            features[9 * w + f] = (float)(sums[f] * Math.Sqrt(2 * (f / 3) + 1) / window);
    }
    return features;
}

static double Distance(float[] a, float[] b)
{
    double sum = 0;
    for (int w = 0; w < a.Length / 9; w++)
    {
        double squared = 0;
        for (int f = 0; f < 9; f++)
        {
            double d = a[9 * w + f] - b[9 * w + f];
            squared += d * d;
        }
        sum += Math.Sqrt(squared);
    }
    return sum / (a.Length / 9);
}

static (double X, double Y, double Z) DnaAxes(char b) => b switch
{
    'A' => (1, 1, -1),
    'C' => (-1, 1, 1),
    'G' => (1, -1, 1),
    _ => (-1, -1, -1),
};

static (double X, double Y, double Z) Expected(string protein, (double X, double Y, double Z)[] table, int p) =>
    table[2 * (protein[p / 2] - 'A') + p % 2];

// Per amino acid and codon position 1 or 2, the usage-weighted mean of its synonyms' axes (uniform when the table has none).
(double X, double Y, double Z)[] ExpectedTable(long[] usage)
{
    var table = new (double X, double Y, double Z)[2 * 26];
    foreach (var (aminoAcid, codons) in synonyms)
    {
        long total = codons.Sum(c => usage[c]);
        foreach (var codon in codons)
        {
            double weight = total > 0 ? usage[codon] / (double)total : 1.0 / codons.Length;
            for (int position = 0; position < 2; position++)
            {
                var (x, y, z) = DnaAxes("TCAG"[position == 0 ? codon / 16 : codon / 4 % 4]);
                ref var cell = ref table[2 * (aminoAcid - 'A') + position];
                cell = (cell.X + weight * x, cell.Y + weight * y, cell.Z + weight * z);
            }
        }
    }
    return table;
}

static IEnumerable<(int Start, int End)> Stretches(string strand, int minLength)
{
    for (int frame = 0; frame < 3; frame++)
    {
        int start = frame, i = frame;
        for (; i + 3 <= strand.Length; i += 3)
        {
            int codon = CodonIndex(strand, i);
            if (codon >= 0 && code11[codon] != '*') continue;
            if (i - start >= minLength) yield return (start, i);
            start = i + 3;
        }
        if (i - start >= minLength) yield return (start, i);
    }
}

static int CodonIndex(string s, int i)
{
    int a = "TCAG".IndexOf(s[i]), b = "TCAG".IndexOf(s[i + 1]), c = "TCAG".IndexOf(s[i + 2]);
    return a < 0 || b < 0 || c < 0 ? -1 : 16 * a + 4 * b + c;
}

static string ReverseComplement(string s)
{
    var r = new char[s.Length];
    for (int i = 0; i < s.Length; i++)
        r[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
    return new string(r);
}

static double KmerShare(string query, string other)
{
    var a = Kmers(query);
    if (a.Count == 0) return 0;
    var b = Kmers(other);
    return a.Count(b.Contains) / (double)a.Count;
}

static HashSet<string> Kmers(string protein)
{
    var set = new HashSet<string>();
    for (int i = 0; i + 5 <= protein.Length; i++) set.Add(protein.Substring(i, 5));
    return set;
}

static string Format(double v) => v.ToString("F4", CultureInfo.InvariantCulture);

async Task<string> FetchAsync(string accession)
{
    for (int attempt = 1; ; attempt++)
    {
        try
        {
            var url = $"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&tool=cat-platform-t3&id={accession}";
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
