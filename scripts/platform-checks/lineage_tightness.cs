#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// P5, measure A (2026-10-30): which lineages does each whole-genome dinucleotide signature keep together better, the second-order
// Z-curve (9 correlations) or Karlin's rho* (16 values)? r4_map.cs's corpora (one genome per species, stored dinucleotide counts),
// each signature in its full space, Euclidean. For every genus, family and phylum with at least 5 genomes: purity (the mean share of
// its genomes' 10 nearest genomes in the same taxon) and spread (mean distance of its genomes to their centroid over the mean distance
// of all genomes to the global centroid, in the same space). A taxon is a Z-curve win when its purity is higher with the Z-curve
// (ties excluded); per stratum a two-sided sign test. Strata: cellular by domain and by the taxon's mean GC (< 40, 40-60, > 60 %);
// viruses by realm and by host (phage when the recorded host's first word is a catalogue bacterial or archaeal genus, eukaryote when a
// host is recorded otherwise). Writes S54_LineageTightness_Taxa.csv and S54_LineageTightness.csv (pass the summary's path).
const int k = 10, minMembers = 5;
var output = args.Length > 0 ? args[0] : "S54_LineageTightness.csv";
var directory = Path.GetDirectoryName(Path.GetFullPath(output))!;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var cellular = await LoadAsync(
    "select distinct on (r.source_id) r.id, r.dinucleotide_counts, g.tax_domain, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species, null " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null order by r.source_id, r.length_bp desc, r.id");
var viruses = await LoadAsync(
    "select r.id, r.dinucleotide_counts, t.tax_realm, t.tax_phylum, t.tax_family, t.tax_genus, t.tax_species, r.annotations->'bvbrc'->>'host_name' " +
    "from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid " +
    "where r.kind = 'virus' and r.dinucleotide_counts is not null and r.length_bp >= 1000 order by r.id");

// Genera of the catalogue's bacteria and archaea, GTDB and NCBI spellings both, for telling phages from eukaryotic viruses.
var prokaryoticGenera = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
await using (var cmd = dataSource.CreateCommand(
    "select distinct split_part(g.tax_genus, '_', 1) from gtdb_taxon g where g.tax_genus is not null union " +
    "select distinct t.tax_genus from sequence_record r join taxon t on t.taxid = r.ncbi_taxid where r.kind in ('bacterium', 'archaeon', 'plasmid') and t.tax_genus is not null"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync()) prokaryoticGenera.Add(r.GetString(0));

using var taxaReport = File.CreateText(Path.Combine(directory, "S54_LineageTightness_Taxa.csv"));
taxaReport.WriteLine("Corpus, Rank, Taxon, Genomes, Stratum1, Stratum2, PurityZCurve, PurityKarlin, SpreadZCurve, SpreadKarlin");
var rows = new List<(string Corpus, string Rank, string Stratum, string Value, double Difference)>();

foreach (var (corpus, genomes) in new[] { ("cellular", cellular), ("viruses", viruses) })
{
    var spaces = new[]
    {
        genomes.Select(g => GenomeSignatures.SecondOrder(g.Dinucleotides)).ToArray(),
        genomes.Select(g => GenomeSignatures.KarlinRho(g.Dinucleotides)).ToArray(),
    };
    var neighbours = spaces.Select(Neighbours).ToArray();
    var overall = spaces.Select(points => { var c = Centroid(points); return points.Average(p => Distance(p, c)); }).ToArray();
    Console.Error.WriteLine($"{corpus}: {genomes.Count} genomes");

    foreach (var (rank, of) in new (string, Func<Genome, string?>)[] { ("genus", g => g.Genus), ("family", g => g.Family), ("phylum", g => g.Phylum) })
    {
        var labels = genomes.Select(of).ToArray();
        foreach (var taxon in Enumerable.Range(0, genomes.Count).Where(i => labels[i] is not null).GroupBy(i => labels[i]!).Where(g => g.Count() >= minMembers))
        {
            var members = taxon.ToArray();
            var purity = neighbours.Select(nb => members.Average(i => nb[i].Count(j => labels[j] == taxon.Key) / (double)k)).ToArray();
            var spread = spaces.Select((points, s) =>
            {
                var c = Centroid(members.Select(i => points[i]).ToArray());
                return members.Average(i => Distance(points[i], c)) / overall[s];
            }).ToArray();
            var (first, second) = corpus == "cellular"
                ? (Majority(members.Select(i => genomes[i].Top)), GcClass(members.Average(i => genomes[i].Gc)))
                : (RealmGroup(Majority(members.Select(i => genomes[i].Top))), Majority(members.Select(i => HostType(genomes[i].Host))));
            taxaReport.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"{corpus}, {rank}, \"{taxon.Key}\", {members.Length}, {first}, {second}, {purity[0]:F4}, {purity[1]:F4}, {spread[0]:F4}, {spread[1]:F4}"));
            double difference = purity[0] - purity[1];
            rows.Add((corpus, rank, "all", "all", difference));
            rows.Add((corpus, rank, corpus == "cellular" ? "domain" : "realm", first, difference));
            rows.Add((corpus, rank, corpus == "cellular" ? "gc" : "host", second, difference));
        }
    }
}

using (var report = File.CreateText(output))
{
    report.WriteLine("Corpus, Rank, Stratum, Value, Taxa, ZCurveWins, KarlinWins, Ties, SignTestP, MedianPurityDifference");
    foreach (var group in rows.GroupBy(r => (r.Corpus, r.Rank, r.Stratum, r.Value)).OrderBy(g => g.Key.Corpus).ThenBy(g => g.Key.Rank).ThenBy(g => g.Key.Stratum).ThenBy(g => g.Key.Value))
    {
        var differences = group.Select(r => r.Difference).ToList();
        int wins = differences.Count(d => d > 1e-12), losses = differences.Count(d => d < -1e-12), ties = differences.Count - wins - losses;
        string line = string.Create(CultureInfo.InvariantCulture,
            $"{group.Key.Corpus}, {group.Key.Rank}, {group.Key.Stratum}, {group.Key.Value}, {differences.Count}, {wins}, {losses}, {ties}, " +
            $"{SignTest(wins, losses):F4}, {differences.Order().ElementAt(differences.Count / 2):F4}");
        report.WriteLine(line);
        Console.WriteLine(line);
    }
}

string HostType(string? host)
{
    if (string.IsNullOrWhiteSpace(host)) return "unknown";
    var first = host.Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries)[0];
    return prokaryoticGenera.Contains(first) ? "phage" : "eukaryote";
}

static string RealmGroup(string realm) => realm is "Riboviria" or "Duplodnaviria" or "Floreoviria" or "Varidnaviria" ? realm : "other";

static string GcClass(double gc) => gc < 0.40 ? "GC<40" : gc <= 0.60 ? "GC40-60" : "GC>60";

static string Majority(IEnumerable<string?> values) =>
    values.Select(v => v ?? "none").GroupBy(v => v).OrderByDescending(g => g.Count()).ThenBy(g => g.Key).First().Key;

// Two-sided exact binomial test of wins against losses at p = 1/2.
static double SignTest(int wins, int losses)
{
    int n = wins + losses, extreme = Math.Min(wins, losses);
    if (n == 0) return 1;
    double tail = 0;
    for (int i = 0; i <= extreme; i++) tail += Math.Exp(LogChoose(n, i) - n * Math.Log(2));
    return Math.Min(1, 2 * tail);
}

static double LogChoose(int n, int r)
{
    double sum = 0;
    for (int i = 1; i <= r; i++) sum += Math.Log(n - r + i) - Math.Log(i);
    return sum;
}

static double[] Centroid(double[][] points) => Enumerable.Range(0, points[0].Length).Select(f => points.Average(p => p[f])).ToArray();

static double Distance(double[] a, double[] b)
{
    double sum = 0;
    for (int f = 0; f < a.Length; f++) sum += (a[f] - b[f]) * (a[f] - b[f]);
    return Math.Sqrt(sum);
}

// The k nearest genomes of each genome.
static int[][] Neighbours(double[][] points)
{
    int n = points.Length;
    var result = new int[n][];
    Parallel.For(0, n, i =>
    {
        var distances = new double[n];
        for (int j = 0; j < n; j++) distances[j] = j == i ? double.PositiveInfinity : Distance(points[i], points[j]);
        var order = Enumerable.Range(0, n).ToArray();
        Array.Sort(distances, order);
        result[i] = order.Take(k).ToArray();
    });
    return result;
}

async Task<List<Genome>> LoadAsync(string sql)
{
    var loaded = new List<Genome>();
    await using (var cmd = dataSource.CreateCommand(sql))
    await using (var r = await cmd.ExecuteReaderAsync())
        while (await r.ReadAsync())
        {
            string? Text(int column) => r.IsDBNull(column) || r.GetString(column) == "" ? null : r.GetString(column);
            var counts = r.GetFieldValue<long[]>(1);
            if (counts.Sum() == 0) continue;
            double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
            loaded.Add(new Genome(r.GetInt64(0), gc, counts, Text(2), Text(3), Text(4), Text(5), Text(6), Text(7)));
        }
    // One genome per species, as r4_map.cs (same seed); a genome without a species stands for itself.
    var rng = new Random(20261030);
    return loaded.GroupBy(g => g.Species ?? $"id {g.Id}").Select(s => s.OrderBy(_ => rng.Next()).First()).OrderBy(g => g.Id).ToList();
}

record Genome(long Id, double Gc, long[] Dinucleotides, string? Top, string? Phylum, string? Family, string? Genus, string? Species, string? Host);
