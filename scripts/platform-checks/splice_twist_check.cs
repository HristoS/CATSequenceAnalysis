#:project ../../SequenceAnalise/SequenceAnalise.csproj
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.IO.Compression;
using System.Text.RegularExpressions;
using SequenceAnalyses;

// Run from the repo root (genomes in data/genomes/eukaryotes/, the DNAshape table in its shape/ folder); no database. Release:
//   dotnet run -c Release scripts/platform-checks/splice_twist_check.cs [folder] [output.csv]
// Check (g), 2026-11-12: does the DNA's 3D geometry, summed step by step from the gene's ATG, mark where genes are spliced and stop?
// Each base step gets its six step parameters (tilt, roll, twist; shift, slide, rise) from the DNAshape Monte Carlo pentamer table
// (Zhou et al. 2013, via DNAshapeR's TableCompiler.cpp: Tilt, Roll_mc, HelT_mc, Shift, Slide, Rise, the step from the pentamer's centre to
// the next base). Steps are composed as in 3DNA (Lu & Olson 2003): rotation Rz(twist/2 - phi) Ry(gamma) Rz(twist/2 + phi), gamma = the
// bend sqrt(tilt^2 + roll^2), phi = atan2(tilt, roll); the displacement in the mid-step frame Rz(twist/2 - phi) Ry(gamma/2) Rz(phi). That
// gives every base a frame (orientation) and an origin (position, in angstrom). Between two points p -> q the features are the relative
// rotation (9 entries), the direction to q in p's frame (3) and log distance; the 1D variant is the summed twist alone (cos, sin).
// Pairs, real against decoy, on the genes of check (f) (longest isoform, canonical GT..AG, ATG..stop, no internal stop):
//   site_donor / site_acceptor  ATG -> each real donor / acceptor; decoys: other GT / AG cuts in the gene span (up to 5 per real site)
//   intron_span                 donor -> its acceptor; decoys: the same donor -> another AG at least the 1st-percentile intron on
//   exon_span                   exon start (ATG or acceptor) -> its donor; decoys: the same start -> another GT downstream
//   stop                        ATG -> stop codon: no decoys, only the concentration of the summed twist (R)
// The summed geometry is mostly length, so every test has a length-only control: the same features on an ideal helix with constant step
// parameters (the genome's means) plus log length. Logistic models (trained on one genome, AUC on its relative): control alone,
// control + 1D twist, control + 3D. Rule fixed before the run: the angles are worth using for cutting only if, in both Nelumbo
// directions, control + twist or control + 3D reaches AUC >= 0.55 on some test and beats the control alone by >= 0.02.
var folder = args.Length > 0 ? args[0] : "data/genomes/eukaryotes";
var output = args.Length > 1 ? args[1] : "results/SpliceTwist_Check.csv";
var pairs = new[]
{
    ("Podospora", ("Podospora bellae-mahoneyi", "GCF_035222275.1_CBS112042p_1.0"), ("Podospora pseudocomata", "GCF_035222375.1_CBS415.72m_1.0")),
    ("Nelumbo", ("Nelumbo nucifera", "GCF_055504915.1_China_Antique-T2T"), ("Nelumbo lutea", "GCF_056320005.1_American_lotus-T2T")),
};
const int pad = 8, decoysPerSite = 5;
string[] kinds = ["site_donor", "site_acceptor", "intron_span", "exon_span"];
var clock = Stopwatch.StartNew();
void Log(string message) => Console.WriteLine($"[{clock.Elapsed:hh\\:mm\\:ss}] {message}");

// ---- step parameters per pentamer (the step from its centre base to the next) ---------------------------------------------------------------
static int Code(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 };
var step = new double[1024, 6]; var known = new bool[1024]; // tilt, roll, twist (degrees); shift, slide, rise (angstrom)
foreach (var line in File.ReadLines(Path.Combine(folder, "shape", "pentamers.txt")))
{
    var f = line.Split(' '); int k = 0; foreach (char c in f[0]) k = k * 4 + Code(c);
    double V(int i) => double.Parse(f[1 + i], CultureInfo.InvariantCulture);
    step[k, 0] = V(33); step[k, 1] = V(37); step[k, 2] = V(50); step[k, 3] = V(43); step[k, 4] = V(48); step[k, 5] = V(41); known[k] = true;
    int rc = 0; for (int i = 4; i >= 0; i--) rc = rc * 4 + (3 - Code(f[0][i]));
    // the reverse complement's step 2 is this pentamer's step 1, read from the other strand: tilt and shift change sign
    if (!known[rc]) { step[rc, 0] = -V(32); step[rc, 1] = V(36); step[rc, 2] = V(49); step[rc, 3] = -V(42); step[rc, 4] = V(47); step[rc, 5] = V(40); }
}
const string bases = "ACGT";
Console.WriteLine("Step parameters, the mean over each dinucleotide's 64 pentamer contexts (DNAshape Monte Carlo):");
Console.WriteLine("  step   tilt   roll  twist  shift  slide   rise");
var mean = new double[6];
for (int s = 0; s < 16; s++)
{
    var m = new double[6];
    for (int k = 0; k < 1024; k++) if (k / 4 % 16 == s) for (int j = 0; j < 6; j++) m[j] += step[k, j] / 64;
    Console.WriteLine($"  {bases[s / 4]}-{bases[s % 4]}  {m[0],6:F2} {m[1],6:F2} {m[2],6:F2} {m[3],6:F2} {m[4],6:F2} {m[5],6:F2}");
}
static double[] Mul(double[] a, double[] b)
{
    var r = new double[9];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) r[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
    return r;
}
static void MulInto(double[] a, double[] b, double[] r)
{
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) r[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
}
static double[] Rz(double deg) { double a = deg * Math.PI / 180, c = Math.Cos(a), s = Math.Sin(a); return [c, -s, 0, s, c, 0, 0, 0, 1]; }
static double[] Ry(double deg) { double a = deg * Math.PI / 180, c = Math.Cos(a), s = Math.Sin(a); return [c, 0, s, 0, 1, 0, -s, 0, c]; }
// The step's rotation and its displacement expressed in the first base's frame.
static (double[] T, double[] D) StepTransform(double tilt, double roll, double twist, double shift, double slide, double rise)
{
    double gamma = Math.Sqrt(tilt * tilt + roll * roll), phi = gamma < 1e-9 ? 0 : Math.Atan2(tilt, roll) * 180 / Math.PI;
    var t = Mul(Mul(Rz(twist / 2 - phi), Ry(gamma)), Rz(twist / 2 + phi));
    var mid = Mul(Mul(Rz(twist / 2 - phi), Ry(gamma / 2)), Rz(phi));
    double[] d = [mid[0] * shift + mid[1] * slide + mid[2] * rise, mid[3] * shift + mid[4] * slide + mid[5] * rise, mid[6] * shift + mid[7] * slide + mid[8] * rise];
    return (t, d);
}
var transforms = new (double[] T, double[] D)[1024];
for (int k = 0; k < 1024; k++) transforms[k] = StepTransform(step[k, 0], step[k, 1], step[k, 2], step[k, 3], step[k, 4], step[k, 5]);

static TextReader Open(string path) =>
    new StreamReader(new GZipStream(new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20), CompressionMode.Decompress), bufferSize: 1 << 20);
static byte[] ReverseComplement(byte[] s) { var r = new byte[s.Length]; for (int i = 0; i < s.Length; i++) r[s.Length - 1 - i] = s[i] == 4 ? (byte)4 : (byte)(3 - s[i]); return r; }
static bool IsGt(byte[] s, int i) => i >= 0 && i + 1 < s.Length && s[i] == 2 && s[i + 1] == 3;
static bool IsAg(byte[] s, int i) => i >= 0 && i + 1 < s.Length && s[i] == 0 && s[i + 1] == 2;
static int CodonAt(byte[] t, int i) => i < 0 || i + 2 >= t.Length || t[i] > 3 || t[i + 1] > 3 || t[i + 2] > 3 ? -1 : t[i] * 16 + t[i + 1] * 4 + t[i + 2];
static bool IsStopCodon(int c) => c is 48 or 50 or 56;
static int Pentamer(byte[] t, int i) { if (i < 2 || i + 2 >= t.Length) return -1; int k = 0; for (int j = i - 2; j <= i + 2; j++) { if (t[j] > 3) return -1; k = k * 4 + t[j]; } return k; }

// ---- genes, as check (f) keeps them -------------------------------------------------------------------------------------------------------
List<Gene> LoadGenes(string label, string prefix)
{
    var range = new Regex(@"(\d+)\.\.(\d+)"); var geneId = new Regex(@"GeneID:(\d+)");
    var longest = new Dictionary<string, (string SeqId, bool Complement, List<(long Start, long End)> Parts, long Length)>();
    using (var reader = Open(Path.Combine(folder, prefix + "_cds_from_genomic.fna.gz")))
    {
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] != '>') continue;
            int bar = line.IndexOf('|'), cds = line.IndexOf("_cds_", StringComparison.Ordinal), loc = line.IndexOf("[location=", StringComparison.Ordinal);
            if (bar < 0 || cds < 0 || loc < 0) continue;
            string location = line[loc..line.IndexOf(']', loc)];
            var parts = range.Matches(location).Select(m => (Start: long.Parse(m.Groups[1].Value) - 1, End: long.Parse(m.Groups[2].Value))).OrderBy(p => p.Start).ToList();
            if (parts.Count < 2 || location.IndexOfAny(['<', '>']) >= 0 || location.Contains("order(") || line.Contains("[partial=")) continue;
            string id = line[(bar + 1)..cds];
            var m2 = geneId.Match(line); string name = m2.Success ? m2.Groups[1].Value : id + ":" + parts[0].Start;
            long length = parts.Sum(p => p.End - p.Start);
            if (!longest.TryGetValue(name, out var current) || length > current.Length) longest[name] = (id, location.Contains("complement("), parts, length);
        }
    }
    var byStrand = longest.GroupBy(kv => (kv.Value.SeqId, kv.Value.Complement)).ToDictionary(x => x.Key, x => x.ToList());
    var genes = new List<Gene>(); int tried = 0;
    using (var reader = Open(Path.Combine(folder, prefix + "_genomic.fna.gz")))
        foreach (var (id, _, sequence) in FastaReader.ReadRecordStreams(reader))
        {
            if (!byStrand.ContainsKey((id, false)) && !byStrand.ContainsKey((id, true))) { sequence.ReadToEnd(); continue; }
            var forward = new byte[1 << 24]; int n = 0; var buffer = new char[1 << 16]; int read;
            while ((read = sequence.Read(buffer, 0, buffer.Length)) > 0)
                for (int i = 0; i < read; i++)
                {
                    byte b;
                    switch (buffer[i])
                    {
                        case 'A' or 'a': b = 0; break; case 'C' or 'c': b = 1; break; case 'G' or 'g': b = 2; break; case 'T' or 't': b = 3; break;
                        case '\n' or '\r' or ' ': continue;
                        default: b = 4; break;
                    }
                    if (n == forward.Length) Array.Resize(ref forward, forward.Length * 2);
                    forward[n++] = b;
                }
            Array.Resize(ref forward, n);
            foreach (bool reverse in new[] { false, true })
            {
                if (!byStrand.TryGetValue((id, reverse), out var list)) continue;
                var s = reverse ? ReverseComplement(forward) : forward; int length = s.Length;
                foreach (var (name, (_, _, gp, _)) in list)
                {
                    tried++;
                    var exons = gp.Select(p => reverse ? ((int)(length - p.End), (int)(length - p.Start)) : ((int)p.Start, (int)p.End)).OrderBy(e => e.Item1).ToArray();
                    bool ok = true;
                    for (int k = 0; k + 1 < exons.Length; k++)
                        if (exons[k + 1].Item1 - exons[k].Item2 < 4 || !IsGt(s, exons[k].Item2) || !IsAg(s, exons[k + 1].Item1 - 2)) ok = false;
                    if (!ok) continue;
                    var c = exons.SelectMany(e => s[e.Item1..e.Item2]).ToArray();
                    ok = c.Length % 3 == 0 && c.Length >= 6 && CodonAt(c, 0) == 14 && IsStopCodon(CodonAt(c, c.Length - 3)) && c.All(b => b < 4);
                    for (int i = 0; ok && i + 3 < c.Length; i += 3) if (IsStopCodon(CodonAt(c, i))) ok = false;
                    if (!ok) continue;
                    int sliceStart = Math.Max(0, exons[0].Item1 - pad), sliceEnd = Math.Min(length, exons[^1].Item2 + pad);
                    genes.Add(new Gene
                    {
                        Name = name, Slice = s[sliceStart..sliceEnd], G0 = exons[0].Item1 - sliceStart, G1 = exons[^1].Item2 - sliceStart,
                        Exons = exons.Select(e => (e.Item1 - sliceStart, e.Item2 - sliceStart)).ToArray(),
                    });
                }
            }
        }
    genes.Sort((a, b) => string.CompareOrdinal(a.Name, b.Name));
    Log($"{label}: {genes.Count:N0} of {tried:N0} spliced genes kept");
    return genes;
}

// ---- walking each gene in 3D and collecting the pair features ------------------------------------------------------------------------------
const int pose = 9 + 3 + 1;      // relative rotation, unit direction, log distance
const int width = 2 * pose + 1 + 2 + 1; // real pose, ideal-helix pose, log length, cos/sin summed twist, compaction
Samples Collect(List<Gene> genes)
{
    var introns = genes.SelectMany(g => Enumerable.Range(0, g.Exons.Length - 1).Select(i => g.Exons[i + 1].Start - g.Exons[i].End)).Order().ToArray();
    int minIntron = introns[(int)(0.01 * introns.Length)];
    // the ideal helix: the genome's mean step parameters over the gene spans
    var sum = new double[6]; long steps = 0;
    foreach (var g in genes) for (int i = g.G0; i < g.G1; i++) { int k = Pentamer(g.Slice, i); if (k < 0) continue; for (int j = 0; j < 6; j++) sum[j] += step[k, j]; steps++; }
    var ideal = sum.Select(x => x / steps).ToArray();
    var idealStep = StepTransform(ideal[0], ideal[1], ideal[2], ideal[3], ideal[4], ideal[5]);
    {
        double[] fs = [1, 0, 0, 0, 1, 0, 0, 0, 1], os = [0, 0, 0];
        for (int i = 0; i < 100; i++) { for (int r = 0; r < 3; r++) os[r] += fs[r * 3] * idealStep.D[0] + fs[r * 3 + 1] * idealStep.D[1] + fs[r * 3 + 2] * idealStep.D[2]; fs = Mul(fs, idealStep.T); }
        Log($"   sanity: the ideal helix spans {Math.Sqrt(os[0] * os[0] + os[1] * os[1] + os[2] * os[2]):F0} A over 100 bp, its axis turned {Math.Acos(Math.Clamp(fs[8], -1, 1)) * 180 / Math.PI:F1} deg");
    }
    var perGene = new Samples[genes.Count];
    Parallel.For(0, genes.Count, gi =>
    {
        var g = genes[gi]; var t = g.Slice; var local = new Samples(kinds.Length);
        var random = new Random(1013 + gi);
        var trueD = g.Exons[..^1].Select(e => e.End).ToList(); var trueA = g.Exons[1..].Select(e => e.Start).ToList();
        var trueDSet = trueD.ToHashSet(); var trueASet = trueA.ToHashSet();
        var gts = new List<int>(); var ags = new List<int>();
        for (int i = g.G0 + 1; i < g.G1; i++) { if (IsGt(t, i) && !trueDSet.Contains(i)) gts.Add(i); if (IsAg(t, i - 2) && !trueASet.Contains(i)) ags.Add(i); }
        List<int> Pick(List<int> from, Func<int, bool> keep, int count)
        {
            var x = from.Where(keep).ToList();
            for (int i = x.Count - 1; i > 0; i--) { int j = random.Next(i + 1); (x[i], x[j]) = (x[j], x[i]); }
            return x.Take(count).ToList();
        }
        var pairList = new List<(int Kind, bool Real, int From, int To)>();
        foreach (int d in trueD) pairList.Add((0, true, g.G0, d));
        foreach (int d in Pick(gts, _ => true, decoysPerSite * trueD.Count)) pairList.Add((0, false, g.G0, d));
        foreach (int a in trueA) pairList.Add((1, true, g.G0, a));
        foreach (int a in Pick(ags, _ => true, decoysPerSite * trueA.Count)) pairList.Add((1, false, g.G0, a));
        for (int i = 0; i + 1 < g.Exons.Length; i++)
        {
            int d = g.Exons[i].End, a = g.Exons[i + 1].Start, start = g.Exons[i].Start;
            pairList.Add((2, true, d, a));
            foreach (int x in Pick(ags, x => x >= d + minIntron, decoysPerSite)) pairList.Add((2, false, d, x));
            pairList.Add((3, true, start, d));
            foreach (int x in Pick(gts, x => x > start, decoysPerSite)) pairList.Add((3, false, start, x));
        }
        int stopAt = g.G1 - 3;
        var wanted = pairList.SelectMany(p => new[] { p.From, p.To }).Append(g.G0).Append(stopAt).ToHashSet();
        // walk from the ATG: frame F (3x3), origin o, summed twist; the same for the ideal helix
        var frames = new Dictionary<int, Pose>();
        double[] f = [1, 0, 0, 0, 1, 0, 0, 0, 1], o = [0, 0, 0]; double twist = 0;
        double[] fi = [1, 0, 0, 0, 1, 0, 0, 0, 1], oi = [0, 0, 0]; double twistIdeal = 0;
        var buffer = new double[9];
        int last = wanted.Max();
        for (int i = g.G0; i <= last; i++)
        {
            if (wanted.Contains(i)) frames[i] = new Pose((double[])f.Clone(), (double[])o.Clone(), twist, (double[])fi.Clone(), (double[])oi.Clone(), twistIdeal);
            int k = Pentamer(t, i);
            var (tr, d) = k >= 0 ? transforms[k] : idealStep;
            double tw = k >= 0 ? step[k, 2] : ideal[2];
            for (int r = 0; r < 3; r++) o[r] += f[r * 3] * d[0] + f[r * 3 + 1] * d[1] + f[r * 3 + 2] * d[2];
            MulInto(f, tr, buffer); (f, buffer) = (buffer, f); twist += tw;
            for (int r = 0; r < 3; r++) oi[r] += fi[r * 3] * idealStep.D[0] + fi[r * 3 + 1] * idealStep.D[1] + fi[r * 3 + 2] * idealStep.D[2];
            MulInto(fi, idealStep.T, buffer); (fi, buffer) = (buffer, fi); twistIdeal += ideal[2];
            if ((i - g.G0) % 64 == 63) { Orthonormalize(f); Orthonormalize(fi); }
        }
        foreach (var (kind, real, from, to) in pairList)
        {
            var p = frames[from]; var q = frames[to]; var x = new double[width];
            Relative(p.F, p.O, q.F, q.O, x, 0); Relative(p.Fi, p.Oi, q.Fi, q.Oi, x, pose);
            x[2 * pose] = Math.Log(to - from);
            double tw = (q.Twist - p.Twist) * Math.PI / 180; x[2 * pose + 1] = Math.Cos(tw); x[2 * pose + 2] = Math.Sin(tw);
            x[2 * pose + 3] = x[12] - x[pose + 12]; // log(real distance / ideal-helix distance)
            (real ? local.Real : local.Decoy)[kind].Add(x);
        }
        var s0 = frames[g.G0]; var s1 = frames[stopAt];
        local.StopTwist.Add(((s1.Twist - s0.Twist) % 360 + 360) % 360); local.StopTwistIdeal.Add(((s1.TwistIdeal - s0.TwistIdeal) % 360 + 360) % 360);
        perGene[gi] = local;
    });
    var samples = new Samples(kinds.Length) { MinIntron = minIntron, Ideal = ideal };
    foreach (var p in perGene)
    {
        for (int k = 0; k < kinds.Length; k++) { samples.Real[k].AddRange(p.Real[k]); samples.Decoy[k].AddRange(p.Decoy[k]); }
        samples.StopTwist.AddRange(p.StopTwist); samples.StopTwistIdeal.AddRange(p.StopTwistIdeal);
    }
    return samples;
}
static void Orthonormalize(double[] m)
{
    // columns are the frame's axes: Gram-Schmidt on x, y; z = x cross y
    double[] x = [m[0], m[3], m[6]], y = [m[1], m[4], m[7]];
    double nx = Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]); for (int i = 0; i < 3; i++) x[i] /= nx;
    double dot = x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; for (int i = 0; i < 3; i++) y[i] -= dot * x[i];
    double ny = Math.Sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]); for (int i = 0; i < 3; i++) y[i] /= ny;
    double[] z = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
    for (int i = 0; i < 3; i++) { m[i * 3] = x[i]; m[i * 3 + 1] = y[i]; m[i * 3 + 2] = z[i]; }
}
// q's frame and position as seen from p's frame: rotation Fp^T Fq (9), unit direction Fp^T (oq - op) (3), log distance (1).
static void Relative(double[] fp, double[] op, double[] fq, double[] oq, double[] x, int at)
{
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) x[at + i * 3 + j] = fp[i] * fq[j] + fp[3 + i] * fq[3 + j] + fp[6 + i] * fq[6 + j];
    double[] v = [oq[0] - op[0], oq[1] - op[1], oq[2] - op[2]];
    double[] local = [fp[0] * v[0] + fp[3] * v[1] + fp[6] * v[2], fp[1] * v[0] + fp[4] * v[1] + fp[7] * v[2], fp[2] * v[0] + fp[5] * v[1] + fp[8] * v[2]];
    double norm = Math.Sqrt(local[0] * local[0] + local[1] * local[1] + local[2] * local[2]);
    for (int i = 0; i < 3; i++) x[at + 9 + i] = norm < 1e-9 ? 0 : local[i] / norm;
    x[at + 12] = Math.Log(norm + 1);
}

static double[] Fit(List<(double[] X, bool Y)> rows)
{
    int d = rows[0].X.Length + 1;
    var w = new double[d];
    for (int iteration = 0; iteration < 25; iteration++)
    {
        var h = new double[d, d]; var g = new double[d];
        foreach (var (x, y) in rows)
        {
            double z = w[0]; for (int j = 1; j < d; j++) z += w[j] * x[j - 1];
            double p = 1 / (1 + Math.Exp(-z)), weight = Math.Max(p * (1 - p), 1e-9), err = (y ? 1 : 0) - p;
            for (int i = 0; i < d; i++) { double xi = i == 0 ? 1 : x[i - 1]; g[i] += err * xi; double wx = weight * xi; for (int j = i; j < d; j++) h[i, j] += wx * (j == 0 ? 1 : x[j - 1]); }
        }
        for (int i = 0; i < d; i++) for (int j = 0; j < i; j++) h[i, j] = h[j, i];
        for (int i = 0; i < d; i++) { h[i, i] += 1e-3 * rows.Count; g[i] -= 1e-3 * rows.Count * w[i]; }
        var stepW = Solve(h, g);
        double change = 0; for (int i = 0; i < d; i++) { w[i] += stepW[i]; change += Math.Abs(stepW[i]); }
        if (change < 1e-6) break;
    }
    return w;
}
static double[] Solve(double[,] m, double[] v)
{
    int n = v.Length; var a = new double[n, n + 1];
    for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) a[i, j] = m[i, j]; a[i, n] = v[i]; }
    for (int c = 0; c < n; c++)
    {
        int p = Enumerable.Range(c, n - c).MaxBy(r => Math.Abs(a[r, c]));
        for (int j = 0; j <= n; j++) (a[c, j], a[p, j]) = (a[p, j], a[c, j]);
        for (int r = 0; r < n; r++) if (r != c) { double k = a[r, c] / a[c, c]; for (int j = c; j <= n; j++) a[r, j] -= k * a[c, j]; }
    }
    return Enumerable.Range(0, n).Select(i => a[i, n] / a[i, i]).ToArray();
}
// Logistic model on the chosen columns: trained on one genome's pairs, AUC (Mann-Whitney) on the other's.
static double Auc(int[] columns, List<double[]> trainReal, List<double[]> trainDecoy, List<double[]> testReal, List<double[]> testDecoy)
{
    double[] Take(double[] x) => columns.Select(c => x[c]).ToArray();
    var rows = trainReal.Select(x => (Take(x), true)).Concat(trainDecoy.Select(x => (Take(x), false))).ToList();
    int d = columns.Length;
    var mu = Enumerable.Range(0, d).Select(j => rows.Average(r => r.Item1[j])).ToArray();
    var sd = Enumerable.Range(0, d).Select(j => { double v = Math.Sqrt(rows.Average(r => (r.Item1[j] - mu[j]) * (r.Item1[j] - mu[j]))); return v < 1e-12 ? 1 : v; }).ToArray();
    var w = Fit(rows.Select(r => (r.Item1.Select((v, j) => (v - mu[j]) / sd[j]).ToArray(), r.Item2)).ToList());
    double Score(double[] x) { double z = w[0]; for (int j = 0; j < d; j++) z += w[j + 1] * (x[columns[j]] - mu[j]) / sd[j]; return z; }
    var scored = testReal.Select(x => (Score(x), true)).Concat(testDecoy.Select(x => (Score(x), false))).OrderBy(s => s.Item1).ToList();
    double rankSum = 0; int i = 0;
    while (i < scored.Count)
    {
        int j = i; while (j + 1 < scored.Count && scored[j + 1].Item1 == scored[i].Item1) j++;
        double rank = (i + j) / 2.0 + 1; for (int k = i; k <= j; k++) if (scored[k].Item2) rankSum += rank;
        i = j + 1;
    }
    double nReal = testReal.Count, nDecoy = testDecoy.Count;
    return (rankSum - nReal * (nReal + 1) / 2) / (nReal * nDecoy);
}
static (double R, double Mean) Concentration(List<double> angles)
{
    double c = 0, s = 0; foreach (double a in angles) { c += Math.Cos(a * Math.PI / 180); s += Math.Sin(a * Math.PI / 180); }
    double mean = Math.Atan2(s, c) * 180 / Math.PI; if (mean < 0) mean += 360;
    return (Math.Sqrt(c * c + s * s) / angles.Count, mean);
}
static double Median(IEnumerable<double> values) { var v = values.Order().ToArray(); return v.Length == 0 ? double.NaN : v[v.Length / 2]; }

int[] control = [.. Enumerable.Range(pose, pose + 1)];                                 // ideal-helix pose + log length
int[] withTwist = [.. control, 2 * pose + 1, 2 * pose + 2];
int[] with3D = [.. control, .. Enumerable.Range(0, pose), 2 * pose + 3];
var rows = new List<string> { "pair,test,trained_on,kind,real,decoy,auc_control,auc_control_twist,auc_control_3d,median_bend_real,median_bend_decoy,median_compaction_real,median_compaction_decoy" };
var verdicts = new List<bool>();
foreach (var (pair, (labelA, prefixA), (labelB, prefixB)) in pairs)
{
    var a = Collect(LoadGenes(labelA, prefixA)); var b = Collect(LoadGenes(labelB, prefixB));
    Log($"{pair}: ideal helix (mean step) tilt {a.Ideal[0]:F2}, roll {a.Ideal[1]:F2}, twist {a.Ideal[2]:F2} deg, rise {a.Ideal[5]:F2} A " +
        $"(a turn every {360 / a.Ideal[2]:F2} bp); 1st-percentile intron {a.MinIntron} / {b.MinIntron} bp");
    foreach (var (name, s) in new[] { (labelA, a), (labelB, b) })
    {
        var (r, m) = Concentration(s.StopTwist); var (ri, _) = Concentration(s.StopTwistIdeal);
        Console.WriteLine($"   {name}: summed twist ATG -> stop over {s.StopTwist.Count:N0} genes: R {r:F3} (mean {m:F0} deg); ideal helix (length only) R {ri:F3}");
    }
    foreach (var (train, test, trainLabel, testLabel) in new[] { (a, b, labelA, labelB), (b, a, labelB, labelA) })
    {
        bool passes = false;
        Console.WriteLine($"== test {testLabel} (trained on {trainLabel})");
        Console.WriteLine($"   {"kind",-14} {"real",9} {"decoy",9} {"control",8} {"+twist",8} {"+3D",8} {"bend real/decoy",17} {"compaction real/decoy",22}");
        for (int k = 0; k < kinds.Length; k++)
        {
            double c = Auc(control, train.Real[k], train.Decoy[k], test.Real[k], test.Decoy[k]);
            double ct = Auc(withTwist, train.Real[k], train.Decoy[k], test.Real[k], test.Decoy[k]);
            double c3 = Auc(with3D, train.Real[k], train.Decoy[k], test.Real[k], test.Decoy[k]);
            if ((ct >= 0.55 && ct - c >= 0.02) || (c3 >= 0.55 && c3 - c >= 0.02)) passes = true;
            double Bend(double[] x) => Math.Acos(Math.Clamp(x[8], -1, 1)) * 180 / Math.PI; // angle between the two helix axes
            double br = Median(test.Real[k].Select(Bend)), bd = Median(test.Decoy[k].Select(Bend));
            double cr = Median(test.Real[k].Select(x => Math.Exp(x[2 * pose + 3]))), cd = Median(test.Decoy[k].Select(x => Math.Exp(x[2 * pose + 3])));
            Console.WriteLine($"   {kinds[k],-14} {test.Real[k].Count,9:N0} {test.Decoy[k].Count,9:N0} {c,8:F3} {ct,8:F3} {c3,8:F3} {br,8:F1}/{bd,-8:F1} {cr,11:F3}/{cd,-10:F3}");
            rows.Add(string.Create(CultureInfo.InvariantCulture,
                $"{pair},{testLabel},{trainLabel},{kinds[k]},{test.Real[k].Count},{test.Decoy[k].Count},{c:F4},{ct:F4},{c3:F4},{br:F2},{bd:F2},{cr:F4},{cd:F4}"));
        }
        if (pair == "Nelumbo") verdicts.Add(passes);
        Console.WriteLine($"   rule for this direction: {(passes ? "some test passes" : "no test passes")}");
    }
    File.WriteAllLines(output, rows);
}
Log(verdicts.Count == 2 && verdicts.All(v => v)
    ? "Rule: PASSES (Nelumbo, both directions: the angles reach AUC >= 0.55 and add >= 0.02 over the length-only control)"
    : "Rule: FAILS (Nelumbo: the angles don't reach AUC >= 0.55 while adding >= 0.02 over the length-only control in both directions)");
Log($"wrote {output}");

sealed record Pose(double[] F, double[] O, double Twist, double[] Fi, double[] Oi, double TwistIdeal);

sealed class Gene
{
    public string Name = "";
    public byte[] Slice = [];
    public int G0, G1;
    public (int Start, int End)[] Exons = [];
}

sealed class Samples
{
    public Samples(int kinds)
    {
        Real = new List<double[]>[kinds]; Decoy = new List<double[]>[kinds];
        for (int k = 0; k < kinds; k++) { Real[k] = []; Decoy[k] = []; }
    }
    public List<double[]>[] Real, Decoy;
    public List<double> StopTwist = [], StopTwistIdeal = [];
    public int MinIntron;
    public double[] Ideal = [];
}
