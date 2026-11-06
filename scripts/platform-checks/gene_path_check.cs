#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Two checks before the compare page's gene-built path and gene highlighting (2026-11-06).
// 1. Origin offset: the stored "1-origin" levels don't record where the origin is. Recomputed from today's stored 500 bp level
//    (OriginRotation.OriginWindow), does rotating by it reproduce the stored 2 kbp origin level exactly? Rule: all match.
// 2. Gene-built path: each record's stored genes in chromosome order from that origin, each window adding window positions x its
//    means (codon positions 1 and 2), a reverse-strand gene flipped onto the forward strand (windows in reverse order, x and y
//    negated); intergenic DNA left out. The overlay gap of compare_overlay_check.cs (drift removed, 400 points over the path, mean
//    distance over the curves' mean radius, both strands of the second record tried), K-12 against the same three genomes.
//    Rule fixed before the run: build the Genes source if the two relatives' gaps are no worse than on the genome-window path from
//    the origin (0.30 and 0.17) and the control stays above both.
long[] ids = [9354, 37323, 38539, 4907];
string[] names = ["K-12", "same genus: Escherichia sp. TC-EC600-tetX4", "other genus: Salmonella sp. S13", "control: B. subtilis 168"];
const int points = 400;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var lengths = new Dictionary<long, long>();
await using (var command = dataSource.CreateCommand("select id, length_bp from sequence_record where id = any($1)"))
{
    command.Parameters.AddWithValue(ids);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) lengths[reader.GetInt64(0)] = reader.GetInt64(1);
}

var origins = new Dictionary<long, long>();
Console.WriteLine("1. Origin offset");
foreach (var id in ids)
{
    var fine = (await ProfileLevels.ReadAsync(dataSource, id, RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize))!;
    var stored = await ProfileLevels.ReadAsync(dataSource, id, RecordStore.Method, OriginRotation.RotatedVersion(RecordStore.MethodVersion), 2000);
    var window = OriginRotation.OriginWindow(fine, lengths[id]);
    bool same = stored is not null && window is int w && RecordStore.Ladder(fine.RotateTo(w)).Level(4).Windows.SequenceEqual(stored.Windows);
    origins[id] = (window ?? 0) * (long)RecordStore.FineWindowSize;
    Console.WriteLine($"  {id}: origin window {window?.ToString() ?? "none"} ({origins[id]:N0} bp), stored origin level {(stored is null ? "none" : "present")}, rotation reproduces it: {same}");
}

double[][] Stretch(List<double[]> sums)
{
    int m = sums.Count;
    double meanK = (m - 1) / 2.0, varK = Enumerable.Range(0, m).Sum(k => (k - meanK) * (k - meanK));
    var fit = Enumerable.Range(0, 3).Select(a =>
    {
        double mean = sums.Average(s => s[a]);
        return (mean, slope: Enumerable.Range(0, m).Sum(k => (k - meanK) * (sums[k][a] - mean)) / varK);
    }).ToArray();
    return Enumerable.Range(0, points).Select(i =>
    {
        int k = (int)Math.Round(i * (m - 1) / (double)(points - 1));
        return Enumerable.Range(0, 3).Select(a => sums[k][a] - fit[a].mean - fit[a].slope * (k - meanK)).ToArray();
    }).ToArray();
}

async Task<double[][]> GenePathAsync(long id)
{
    var profile = (await GeneProfiles.ReadAsync(dataSource, id))!;
    long origin = origins[id], length = lengths[id];
    var order = Enumerable.Range(0, profile.Genes.Count).OrderBy(g => ((profile.Genes[g].Start - origin) % length + length) % length).ToList();
    var sums = new List<double[]> { new double[3] };
    foreach (int g in order)
    {
        var gene = profile.Genes[g];
        int positions = ZCurveGeneProfile.WindowOf(gene.Residues);
        var windows = profile.WindowsOf(g);
        var steps = Enumerable.Range(0, gene.WindowCount).Select(w => new[] { windows[9 * w], windows[9 * w + 1], windows[9 * w + 2] }).ToList();
        if (gene.Reverse) steps = steps.AsEnumerable().Reverse().Select(s => new[] { -s[0], -s[1], s[2] }).ToList();
        foreach (var s in steps)
        {
            var last = sums[^1];
            sums.Add([last[0] + positions * s[0], last[1] + positions * s[1], last[2] + positions * s[2]]);
        }
    }
    return Stretch(sums);
}

static double[][] ReverseStrand(double[][] p) => p.Reverse().Select(q => new[] { -q[0], -q[1], q[2] }).ToArray();
static double Radius(double[][] p) => p.Average(q => Math.Sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]));
static double Gap(double[][] a, double[][] b) =>
    a.Zip(b).Average(z => Math.Sqrt(Enumerable.Range(0, 3).Sum(i => Math.Pow(z.First[i] - z.Second[i], 2)))) / ((Radius(a) + Radius(b)) / 2);

Console.WriteLine("2. Gene-built path from the origin (genome-window path from the origin: 0.30 / 0.17 / 0.98)");
var lines = new List<string> { "Pair,Gene path gap (forward),Gene path gap (reverse strand),Best,Genome-window path from the origin" };
double[] windowGaps = [0.30, 0.17, 0.98];
var reference = await GenePathAsync(ids[0]);
var best = new double[3];
for (int i = 1; i < ids.Length; i++)
{
    var other = await GenePathAsync(ids[i]);
    double forward = Gap(reference, other), reverse = Gap(reference, ReverseStrand(other));
    best[i - 1] = Math.Min(forward, reverse);
    lines.Add(string.Create(CultureInfo.InvariantCulture, $"K-12 vs {names[i]},{forward:F2},{reverse:F2},{best[i - 1]:F2},{windowGaps[i - 1]:F2}"));
}
bool pass = best[0] <= windowGaps[0] && best[1] <= windowGaps[1] && best[2] > Math.Max(best[0], best[1]);
lines.Add("");
lines.Add($"Rule (relatives no worse than on the genome-window path, control above both): {(pass ? "build" : "stop")}");
Console.WriteLine(string.Join(Environment.NewLine, lines));
File.WriteAllLines(args.Length > 0 ? args[0] : "GenePath_Check.csv", lines);
