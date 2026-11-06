#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Globalization;
using System.Text;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Compare-page check (2026-11-06): do two genomes' Z-curve paths line up when drawn on one canvas? Each path is the running sum of a
// stored 2 kbp level, drift removed (each axis's least-squares line), stretched to 400 points over its length. Gap = mean distance
// between the two curves' points over the mean radius of the two curves. Both strands of the second genome are tried (the reverse
// strand walks backwards with x and y negated). Levels: as deposited ("1") and laid out from the replication origin ("1-origin").
// Rule fixed before the run: if close relatives don't line up better than the unrelated control even after origin rotation, the page
// shows curves side by side instead of overlaid. Args: output SVG folder (default .).
var folder = args.Length > 0 ? args[0] : ".";
(long, string)[] pairs = [(37323, "same genus: Escherichia sp. TC-EC600-tetX4"), (38539, "other genus: Salmonella sp. S13"), (4907, "control: B. subtilis 168")];
const long k12 = 9354;
const int points = 400;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

async Task<double[][]> PathAsync(long id, string version)
{
    var level = await ProfileLevels.ReadAsync(dataSource, id, RecordStore.Method, version, 2000) ?? throw new InvalidOperationException($"{id} has no {version} level");
    var w = level.Windows;
    var sums = new double[w.Count + 1][];
    sums[0] = [0, 0, 0];
    for (int k = 0; k < w.Count; k++) sums[k + 1] = [sums[k][0] + 2000 * w[k].X, sums[k][1] + 2000 * w[k].Y, sums[k][2] + 2000 * w[k].Z];
    int m = sums.Length;
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

static double[][] ReverseStrand(double[][] p) => p.Reverse().Select(q => new[] { -q[0], -q[1], q[2] }).ToArray();
static double Radius(double[][] p) => p.Average(q => Math.Sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]));
static double Gap(double[][] a, double[][] b) =>
    a.Zip(b).Average(z => Math.Sqrt(Enumerable.Range(0, 3).Sum(i => Math.Pow(z.First[i] - z.Second[i], 2)))) / ((Radius(a) + Radius(b)) / 2);

void Svg(string name, double[][] a, double[][] b, string title)
{
    double r = Math.Max(a.Concat(b).Max(q => Math.Abs(q[0])), a.Concat(b).Max(q => Math.Abs(q[1])));
    string Line(double[][] p, string colour) => $"<polyline fill='none' stroke='{colour}' stroke-width='2' points='" +
        string.Join(" ", p.Select(q => string.Create(CultureInfo.InvariantCulture, $"{250 + 220 * q[0] / r:F1},{270 - 220 * q[1] / r:F1}"))) + "'/>";
    File.WriteAllText(Path.Combine(folder, name), $"<svg xmlns='http://www.w3.org/2000/svg' width='500' height='520'><rect width='500' height='520' fill='white'/>" +
        $"<text x='10' y='20' font-family='sans-serif' font-size='13'>{title}</text>{Line(a, "#2a78d6")}{Line(b, "#eb6834")}</svg>");
}

var lines = new List<string> { "Pair,Level,Gap (forward),Gap (reverse strand),Best" };
foreach (var version in new[] { "1", "1-origin" })
{
    var reference = await PathAsync(k12, version);
    foreach (var (id, label) in pairs)
    {
        var other = await PathAsync(id, version);
        double forward = Gap(reference, other), reverse = Gap(reference, ReverseStrand(other));
        lines.Add(string.Create(CultureInfo.InvariantCulture, $"K-12 vs {label},{version},{forward:F2},{reverse:F2},{Math.Min(forward, reverse):F2}"));
        Svg($"overlay_{id}_{version}.svg", reference, forward <= reverse ? other : ReverseStrand(other), $"K-12 (blue) vs {label}, level {version}, x-y");
    }
}
Console.WriteLine(string.Join(Environment.NewLine, lines));
