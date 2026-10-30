#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// P5, where the second-order Z-curve's gap to Karlin comes from (2026-10-30): r4_map.cs's corpora, 10-nearest purity per rank with a
// paired bootstrap against karlin (as second_order_log_neighbours.cs). Spaces, fixed before the run, all Euclidean: karlin (rho*, ratio
// scale, 16), karlinLog (log rho*, 0.5 added per count, 16), zcurveLogFull (all 16 coefficients of that log rho* in Z-curve
// coordinates: a rotation, so it must equal karlinLog), zcurveLogInteraction (only its 9 neighbour-interaction coefficients, the part
// the second-order Z-curve keeps; drops the constant and the 6 main-effect coefficients) and secondOrder (the 9 correlations). karlin
// against karlinLog isolates the scale; karlinLog against zcurveLogInteraction the dropped coefficients. Args: output path (default
// S52e_GapDecomposition.csv).
const int k = 10, neighbourList = 300, resamples = 1000;
var output = args.Length > 0 ? args[0] : "S52e_GapDecomposition.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var cellular = await LoadAsync(
    "select distinct on (r.source_id) r.id, r.name, r.length_bp, r.dinucleotide_counts, g.tax_domain, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null order by r.source_id, r.length_bp desc, r.id");
var viruses = await LoadAsync(
    "select r.id, r.name, r.length_bp, r.dinucleotide_counts, t.tax_realm, t.tax_phylum, t.tax_family, t.tax_genus, t.tax_species " +
    "from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid " +
    "where r.kind = 'virus' and r.dinucleotide_counts is not null and r.length_bp >= 1000 order by r.id");

using var report = File.CreateText(output);
report.WriteLine("Corpus, Space, Rank, Queries, Purity, MinusKarlin, MinusKarlinCiLow, MinusKarlinCiHigh");
foreach (var (corpus, genomes) in new[] { ("cellular", cellular), ("viruses", viruses) })
{
    var spaces = new (string Name, double[][] Points)[]
    {
        ("karlin", genomes.Select(g => GenomeSignatures.KarlinRho(g.Dinucleotides)).ToArray()),
        ("karlinLog", genomes.Select(g => GenomeSignatures.KarlinLogRho(g.Dinucleotides)).ToArray()),
        ("zcurveLogFull", genomes.Select(g => GenomeSignatures.ToZCurveBasis(GenomeSignatures.KarlinLogRho(g.Dinucleotides))).ToArray()),
        ("zcurveLogInteraction", genomes.Select(g => Interactions(GenomeSignatures.ToZCurveBasis(GenomeSignatures.KarlinLogRho(g.Dinucleotides)))).ToArray()),
        ("secondOrder", genomes.Select(g => GenomeSignatures.SecondOrder(g.Dinucleotides)).ToArray()),
    };
    var neighbours = spaces.Select(s => Neighbours(s.Points)).ToArray();
    foreach (var (rank, of) in new (string, Func<Genome, string?>)[] { ("genus", g => g.Genus), ("family", g => g.Family), ("phylum", g => g.Phylum) })
    {
        var labels = genomes.Select(of).ToArray();
        var queries = Enumerable.Range(0, labels.Length).Where(i => labels[i] is not null).ToArray();
        var perGenome = neighbours.Select(nb => queries.Select(i => nb[i].Take(k).Count(j => labels[j] == labels[i]) / (double)k).ToArray()).ToArray();
        var rng = new Random(20261030);
        var draws = Enumerable.Range(0, resamples).Select(_ => Enumerable.Range(0, queries.Length).Select(_ => rng.Next(queries.Length)).ToArray()).ToArray();
        for (int s = 0; s < spaces.Length; s++)
        {
            var differences = draws.Select(d => d.Average(q => perGenome[s][q] - perGenome[0][q])).Order().ToArray();
            string line = string.Create(CultureInfo.InvariantCulture,
                $"{corpus}, {spaces[s].Name}, {rank}, {queries.Length}, {perGenome[s].Average():F4}, {perGenome[s].Average() - perGenome[0].Average():F4}, " +
                $"{differences[(int)(0.025 * (resamples - 1))]:F4}, {differences[(int)(0.975 * (resamples - 1))]:F4}");
            report.WriteLine(line);
            Console.WriteLine(line);
        }
    }
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

// The 9 coefficients on h_i(n) h_j(n + 1) with i, j in (x, y, z): no constant, no main effects (ToZCurveBasis puts them at 4 * i + j).
static double[] Interactions(double[] coefficients) =>
    new[] { 1, 2, 3 }.SelectMany(i => new[] { 1, 2, 3 }.Select(j => coefficients[4 * i + j])).ToArray();

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

record Genome(long Id, string Name, long Length, double Gc, long[] Dinucleotides, string? Top, string? Phylum, string? Family, string? Genus, string? Species);
