#:project ../../SequenceAnalise/SequenceAnalise.csproj
#:property PublishAot=false
using SequenceAnalyses;

// Can a gene's Z-curve over codon positions 1 and 2 be redrawn from its stored windows (mean, trend, curvature per axis, float16)?
// Truth: the running sums of the gene's own bases. Compared: the curve rebuilt from the three moments, and straight segments from the
// means alone. Genes: ATG-to-stop ORFs of at least 50 codons on both strands of E. coli K-12, longest first, non-overlapping.
var record = FastaReader.ReadPrimaryRecord(args[0]).Sequence.ToUpperInvariant();
var genes = new List<string>();
foreach (var strand in new[] { record, Reverse(record) })
    for (int frame = 0; frame < 3; frame++)
    {
        int start = -1;
        for (int i = frame; i + 3 <= strand.Length; i += 3)
        {
            var codon = strand.AsSpan(i, 3);
            if (start < 0 && codon.SequenceEqual("ATG")) start = i;
            else if (start >= 0 && (codon.SequenceEqual("TAA") || codon.SequenceEqual("TAG") || codon.SequenceEqual("TGA")))
            {
                if ((i - start) / 3 >= 50) genes.Add(strand.Substring(start, i + 3 - start));
                start = -1;
            }
        }
    }
var rng = new Random(20261106);
var sample = genes.OrderBy(_ => rng.Next()).Take(3000).ToList();

var rows = new List<(bool Long, double Moments, double Straight, double Extent, double Step)>();
foreach (var gene in sample)
{
    var w = ZCurveGeneProfile.QueryFromGene(gene).Select(v => (double)(Half)v).ToArray();
    int residues = gene.Length / 3 - 1, window = residues >= 500 ? 333 : 2 * residues / 5, count = w.Length / 9;
    double[] truth = new double[3], rebuilt = new double[3], straight = new double[3];
    double worstMoments = 0, worstStraight = 0, extent = 0, step = 0;
    var path = new List<double[]>();
    for (int k = 0; k < count; k++)
    {
        var start = (double[])truth.Clone();
        double stepLength = 0;
        for (int a = 0; a < 3; a++) stepLength += Math.Pow(window * w[9 * k + a], 2);
        step += Math.Sqrt(stepLength) / count;
        for (int s = 1; s <= window; s++)
        {
            int i = k * window + s - 1;
            var (x, y, z) = Axes(gene[3 * (i / 2) + i % 2]);
            truth[0] += x; truth[1] += y; truth[2] += z;
            double tau = -1 + 2.0 * s / window, d2m = 0, d2s = 0, e = 0;
            for (int a = 0; a < 3; a++)
            {
                double m0 = w[9 * k + a], m1 = w[9 * k + 3 + a], m2 = w[9 * k + 6 + a];
                double r = start[a] + window / 2.0 * (m0 * (tau + 1) + Math.Sqrt(3) * m1 * (tau * tau - 1) / 2 + Math.Sqrt(5) * m2 * (tau * tau * tau - tau) / 2);
                double l = start[a] + m0 * s;
                d2m += (r - truth[a]) * (r - truth[a]);
                d2s += (l - truth[a]) * (l - truth[a]);
                e += truth[a] * truth[a];
            }
            worstMoments = Math.Max(worstMoments, Math.Sqrt(d2m));
            path.Add((double[])truth.Clone());
            worstStraight = Math.Max(worstStraight, Math.Sqrt(d2s));
            extent = Math.Max(extent, Math.Sqrt(e));
        }
        // The rebuilt curve is anchored at the true window start, as a drawing chains windows from exact window sums.
        for (int a = 0; a < 3; a++) truth[a] = start[a] + Enumerable.Range(0, window).Sum(s =>
        {
            int i = k * window + s; var ax = Axes(gene[3 * (i / 2) + i % 2]); return a == 0 ? ax.X : a == 1 ? ax.Y : ax.Z;
        });
    }
    int m = path.Count; double meanK = (m - 1) / 2.0, varK = Enumerable.Range(0, m).Sum(k => (k - meanK) * (k - meanK));
    double residual = 0;
    var fit = Enumerable.Range(0, 3).Select(a => { double mean = path.Average(q => q[a]); double slope = Enumerable.Range(0, m).Sum(k => (k - meanK) * (path[k][a] - mean)) / varK; return (mean, slope); }).ToArray();
    for (int k = 0; k < m; k++) residual = Math.Max(residual, Math.Sqrt(Enumerable.Range(0, 3).Sum(a => Math.Pow(path[k][a] - fit[a].mean - fit[a].slope * (k - meanK), 2))));
    rows.Add((residues >= 500, worstMoments, worstStraight, residual, step));
}

Console.WriteLine($"{sample.Count} genes ({rows.Count(r => r.Long)} of 500 aa or more)");
foreach (var group in rows.GroupBy(r => r.Long).OrderBy(g => g.Key))
{
    double Q(IEnumerable<double> v, double q) { var s = v.Order().ToList(); return s[(int)(q * (s.Count - 1))]; }
    var g = group.ToList();
    Console.WriteLine($"{(group.Key ? "500+ aa (333-position windows)" : "50-499 aa (5 windows)")}: {g.Count} genes");
    Console.WriteLine($"  worst gap to the true curve, bases: moments median {Q(g.Select(r => r.Moments), 0.5):F2}, 95% {Q(g.Select(r => r.Moments), 0.95):F2} | " +
        $"straight median {Q(g.Select(r => r.Straight), 0.5):F2}, 95% {Q(g.Select(r => r.Straight), 0.95):F2}");
    Console.WriteLine($"  as a share of the drift-removed extent: moments median {Q(g.Select(r => r.Moments / r.Extent), 0.5):P1}, 95% {Q(g.Select(r => r.Moments / r.Extent), 0.95):P1} | " +
        $"straight median {Q(g.Select(r => r.Straight / r.Extent), 0.5):P1}, 95% {Q(g.Select(r => r.Straight / r.Extent), 0.95):P1}");
    Console.WriteLine($"  as a share of the mean window step: moments median {Q(g.Select(r => r.Moments / r.Step), 0.5):P1} | straight median {Q(g.Select(r => r.Straight / r.Step), 0.5):P1}");
}

static string Reverse(string s)
{
    var r = new char[s.Length];
    for (int i = 0; i < s.Length; i++) r[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
    return new string(r);
}

static (int X, int Y, int Z) Axes(char b) => b switch
{
    'A' => (1, 1, -1), 'C' => (-1, 1, 1), 'G' => (1, -1, 1), 'T' => (-1, -1, -1), _ => (0, 0, 0),
};
