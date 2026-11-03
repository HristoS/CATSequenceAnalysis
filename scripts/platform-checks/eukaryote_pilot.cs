#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Globalization;
using Npgsql;
using SequencePlatform.Core;

// The eukaryote pilot's quality check: does the dinucleotide signature (log rho*, nuclear counts summed) keep eukaryotic lineages
// together? Over all eukaryote records and within each kingdom: 5-nearest purity (share of a genome's 5 nearest in the same taxon, over
// genomes whose taxon has at least 2 members), for the signature and, as context, for GC + log10 length, against 500 shuffles of the
// labels. A rank passes when no shuffle reaches the observed purity and the purity is at least twice the shuffles' mean; the pilot's rule
// is that class and order must pass in the onePerFamily sample. That sample keeps the first stored genome of each family, because
// re-picks after NCBI's list changed left some families with two, and a same-family pair would count as a class and order hit for
// free. Writes the given CSV.
const int k = 5, shuffles = 500;
const double requiredRatio = 2;
var output = args.Length > 0 ? args[0] : "P8_EukaryotePilot.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var genomes = new List<Genome>();
await using (var command = dataSource.CreateCommand(
    "select r.name, r.length_bp, r.dinucleotide_counts, coalesce(t.tax_kingdom, r.tax_kingdom), coalesce(t.tax_phylum, r.tax_phylum), " +
    "coalesce(t.tax_class, r.tax_class), coalesce(t.tax_order, r.tax_order), coalesce(t.tax_family, r.tax_family) " +
    "from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid where r.kind = 'eukaryote' and r.dinucleotide_counts is not null order by r.id"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
    {
        string? Text(int column) => reader.IsDBNull(column) || reader.GetString(column) == "" ? null : reader.GetString(column);
        var counts = reader.GetFieldValue<long[]>(2);
        double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();
        genomes.Add(new Genome(reader.GetString(0), reader.GetInt64(1), gc, SignatureIndex.Signature(counts), Text(3), Text(4), Text(5), Text(6), Text(7)));
    }
var onePerFamily = genomes.GroupBy(g => g.Family ?? g.Name).Select(g => g.First()).ToList();
Console.WriteLine($"{genomes.Count} eukaryote genomes, {onePerFamily.Count} with one per family.");

using var csv = File.CreateText(output);
csv.WriteLine("Sample, Corpus, Space, Rank, Genomes, Purity, NullMean, NullMax, Pass");
var ranks = new (string Name, Func<Genome, string?> Of)[] { ("kingdom", g => g.Kingdom), ("phylum", g => g.Phylum), ("class", g => g.Class), ("order", g => g.Order), ("family", g => g.Family) };
var corpora = new[] { (Sample: "onePerFamily", Pool: onePerFamily), (Sample: "all", Pool: genomes) }.SelectMany(s =>
    new[] { (s.Sample, Corpus: "all", Members: s.Pool) }.Concat(s.Pool.Where(g => g.Kingdom is not null).GroupBy(g => g.Kingdom!)
        .Where(g => g.Count() > k).Select(g => (s.Sample, Corpus: g.Key, Members: g.ToList()))));
foreach (var (sample, corpus, members) in corpora)
{
    var spaces = new (string Name, double[][] Points)[]
    {
        ("signature", members.Select(g => g.Signature).ToArray()),
        ("gcLength", Standardize(members.Select(g => new[] { g.Gc, Math.Log10(g.Length) }).ToArray())),
    };
    foreach (var (space, points) in spaces)
    {
        var neighbours = Neighbours(points);
        foreach (var (rank, of) in ranks)
        {
            var labels = members.Select(of).ToArray();
            var size = labels.Where(l => l is not null).GroupBy(l => l).ToDictionary(g => g.Key!, g => g.Count());
            for (int i = 0; i < labels.Length; i++) if (labels[i] is { } label && size[label] < 2) labels[i] = null;
            int scored = labels.Count(l => l is not null);
            if (scored == 0 || size.Count(s => s.Value >= 2) < 2) continue;   // one taxon can't be separated from anything

            double purity = Purity(neighbours, labels);
            var rng = new Random(20261030);
            var nulls = Enumerable.Range(0, shuffles).Select(_ => Purity(neighbours, labels.OrderBy(_ => rng.Next()).ToArray())).ToArray();
            bool pass = purity > nulls.Max() && purity >= requiredRatio * nulls.Average();
            var line = string.Create(CultureInfo.InvariantCulture, $"{sample}, {corpus}, {space}, {rank}, {scored}, {purity:F4}, {nulls.Average():F4}, {nulls.Max():F4}, {pass}");
            csv.WriteLine(line);
            Console.WriteLine(line);
        }
    }
}

static int[][] Neighbours(double[][] points)
{
    int n = points.Length;
    var result = new int[n][];
    for (int i = 0; i < n; i++)
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
        result[i] = order.Take(Math.Min(k, n - 1)).ToArray();
    }
    return result;
}

static double Purity(int[][] neighbours, string?[] labels)
{
    double sum = 0;
    int count = 0;
    for (int i = 0; i < labels.Length; i++)
    {
        if (labels[i] is null) continue;
        sum += neighbours[i].Count(j => labels[j] == labels[i]) / (double)neighbours[i].Length;
        count++;
    }
    return sum / count;
}

static double[][] Standardize(double[][] rows)
{
    int width = rows[0].Length;
    var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
    var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(rows.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
    return rows.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray();
}

record Genome(string Name, long Length, double Gc, double[] Signature, string? Kingdom, string? Phylum, string? Class, string? Order, string? Family);
