#:project ../../SequenceAnalise/SequenceAnalise.csproj
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.IO.Compression;
using System.Text.RegularExpressions;
using SequenceAnalyses;

// Run from the repo root (genomes in data/genomes/eukaryotes/, the DNAshape table in its shape/ folder); no database. Release:
//   dotnet run -c Release scripts/platform-checks/splice_join_check.cs [folder] [output.csv] [Podospora|Nelumbo]
// Check (f), 2026-11-12: once the candidate cut points are kept (check (e)'s model 3: site score + Z-curve bins + DNA shape bins, at the
// threshold keeping 95% of real sites), can the pieces be joined back into the right gene? Per annotated multi-exon gene (longest
// isoform, canonical GT..AG introns, ATG..stop, no internal stop), given its true start and end, a dynamic programme picks the best chain
// of pieces: every exon in an open reading frame, including the codons split across a join (the ORF rule), scored by the sites' log-odds
// plus w x a coding score: none (frame rule only), the Z-curve phase-specific composition (per codon position base frequencies, the 9
// phase-specific Z-curve parameters) or a 64-codon table (reference). w is tuned on the training genome, results are on its relative.
// "oracle" adds the true sites to the candidates: what's left is the joining's own error.
// Rule fixed before the run: continue only if in Nelumbo, both directions, kept sites, Z-curve score, intron sensitivity >= 50% and exact
// genes >= 20%; otherwise the puzzle-piece path is abandoned.
var folder = args.Length > 0 ? args[0] : "data/genomes/eukaryotes";
var output = args.Length > 1 ? args[1] : "results/SpliceJoin_Check.csv";
var only = args.Length > 2 ? args[2] : null;
var pairs = new[]
{
    ("Podospora", ("Podospora bellae-mahoneyi", "GCF_035222275.1_CBS112042p_1.0"), ("Podospora pseudocomata", "GCF_035222375.1_CBS415.72m_1.0")),
    ("Nelumbo", ("Nelumbo nucifera", "GCF_055504915.1_China_Antique-T2T"), ("Nelumbo lutea", "GCF_056320005.1_American_lotus-T2T")),
};
const int negativeSample = 200_000, chunk = 1_000_000, reach = 150, pad = reach + 2, tracks = 7, features = 1 + tracks * 6;
const int donorBefore = 3, donorAfter = 6, acceptorBefore = 14, acceptorAfter = 1, none = int.MaxValue;
const double recall = 0.95;
int[] binFrom = [-150, -75, -25, 5, 25, 75], binTo = [-75, -25, -5, 25, 75, 150];
double[] weights = [0.25, 0.5, 1, 2, 4];
var clock = Stopwatch.StartNew();
void Log(string message) => Console.WriteLine($"[{clock.Elapsed:hh\\:mm\\:ss}] {message}");

// ---- shape table (as check (e)) -----------------------------------------------------------------------------------------------------------
static int Code(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 };
var shape = new float[1024, 4]; var known = new bool[1024];
foreach (var line in File.ReadLines(Path.Combine(folder, "shape", "pentamers.txt")))
{
    var f = line.Split(' '); int k = 0; foreach (char c in f[0]) k = k * 4 + Code(c);
    float V(int i) => float.Parse(f[1 + i], CultureInfo.InvariantCulture);
    shape[k, 0] = V(0); shape[k, 1] = V(6); shape[k, 2] = V(18); shape[k, 3] = V(24); known[k] = true;
    int rc = 0; for (int i = 4; i >= 0; i--) rc = rc * 4 + (3 - Code(f[0][i]));
    if (!known[rc]) { shape[rc, 0] = V(0); shape[rc, 1] = V(6); shape[rc, 2] = V(15); shape[rc, 3] = V(21); }
}
var shapeMean = new float[4];
for (int t = 0; t < 4; t++) { double sum = 0; for (int k = 0; k < 1024; k++) sum += shape[k, t]; shapeMean[t] = (float)(sum / 1024); }

static TextReader Open(string path) =>
    new StreamReader(new GZipStream(new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20), CompressionMode.Decompress), bufferSize: 1 << 20);

Genome Load(string label, string prefix)
{
    var g = new Genome { Label = label };
    var range = new Regex(@"(\d+)\.\.(\d+)"); var geneId = new Regex(@"GeneID:(\d+)");
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
            if (parts.Count < 2) continue;
            string id = line[(bar + 1)..cds]; bool complement = location.Contains("complement(");
            if (!g.Cds.TryGetValue(id, out var list)) g.Cds[id] = list = [];
            list.Add((complement, parts));
            if (location.IndexOfAny(['<', '>']) >= 0 || location.Contains("order(") || line.Contains("[partial=")) continue;
            var m2 = geneId.Match(line); string name = m2.Success ? m2.Groups[1].Value : id + ":" + parts[0].Start;
            long length = parts.Sum(p => p.End - p.Start);
            if (!g.Longest.TryGetValue(name, out var current) || length > current.Length) g.Longest[name] = (id, complement, parts, length);
        }
    }
    using (var reader = Open(Path.Combine(folder, prefix + "_genomic.fna.gz")))
        foreach (var (id, _, sequence) in FastaReader.ReadRecordStreams(reader))
        {
            if (!g.Cds.ContainsKey(id)) { sequence.ReadToEnd(); continue; }
            var bytes = new byte[1 << 24]; int n = 0; var buffer = new char[1 << 16]; int read;
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
                    if (n == bytes.Length) Array.Resize(ref bytes, bytes.Length * 2);
                    bytes[n++] = b;
                }
            Array.Resize(ref bytes, n);
            g.Sequences[id] = bytes;
        }
    Log($"{label}: {g.Sequences.Count} sequences with spliced CDS, {g.Sequences.Values.Sum(s => (long)s.Length) / 1e6:F1} Mbp, {g.Longest.Count:N0} complete spliced genes");
    return g;
}

static byte[] ReverseComplement(byte[] s) { var r = new byte[s.Length]; for (int i = 0; i < s.Length; i++) r[s.Length - 1 - i] = s[i] == 4 ? (byte)4 : (byte)(3 - s[i]); return r; }
static bool IsGt(byte[] s, int i) => i + 1 < s.Length && s[i] == 2 && s[i + 1] == 3;
static bool IsAg(byte[] s, int i) => i + 1 < s.Length && s[i] == 0 && s[i + 1] == 2;
static int CodonAt(byte[] t, int i) => i < 0 || i + 2 >= t.Length || t[i] > 3 || t[i + 1] > 3 || t[i + 2] > 3 ? -1 : t[i] * 16 + t[i + 1] * 4 + t[i + 2];
static bool IsStopCodon(int c) => c is 48 or 50 or 56; // TAA, TAG, TGA
static int PendLength(int st) => st == 0 ? 0 : st <= 4 ? 1 : 2;
static int Extend(int st, int b) => st == 0 ? 1 + b : st <= 4 ? 5 + (st - 1) * 4 + b : throw new InvalidOperationException();
static int LowerBound(int[] a, int v) { int lo = 0, hi = a.Length; while (lo < hi) { int mid = (lo + hi) >> 1; if (a[mid] < v) lo = mid + 1; else hi = mid; } return lo; }
static double U(ulong k)
{
    k += 0x9E3779B97F4A7C15UL; k = (k ^ (k >> 30)) * 0xBF58476D1CE4E5B9UL; k = (k ^ (k >> 27)) * 0x94D049BB133111EBUL; k ^= k >> 31;
    return (k >> 11) * (1.0 / (1UL << 53));
}

(HashSet<int> Donors, HashSet<int> Acceptors) Truth(Genome g, string id, bool reverse, byte[] s)
{
    var donors = new HashSet<int>(); var acceptors = new HashSet<int>(); int length = s.Length;
    foreach (var (complement, parts) in g.Cds[id])
    {
        if (complement != reverse) continue;
        var exons = parts.Select(p => reverse ? ((int)(length - p.End), (int)(length - p.Start)) : ((int)p.Start, (int)p.End)).OrderBy(e => e.Item1).ToList();
        for (int k = 0; k + 1 < exons.Count; k++)
        {
            int b = exons[k].Item2, a = exons[k + 1].Item1;
            if (IsGt(s, b)) donors.Add(b);
            if (a >= 2 && IsAg(s, a - 2)) acceptors.Add(a);
        }
    }
    return (donors, acceptors);
}

static double[,] Pwm(IEnumerable<int[]> windows, int width, double[] background)
{
    var counts = new double[width, 4]; for (int i = 0; i < width; i++) for (int b = 0; b < 4; b++) counts[i, b] = 1;
    foreach (var w in windows) for (int i = 0; i < width; i++) if (w[i] < 4) counts[i, w[i]]++;
    var pwm = new double[width, 4];
    for (int i = 0; i < width; i++) { double total = Enumerable.Range(0, 4).Sum(b => counts[i, b]); for (int b = 0; b < 4; b++) pwm[i, b] = Math.Log(counts[i, b] / total / background[b]); }
    return pwm;
}
static double Score(double[,] pwm, byte[] s, int from)
{
    int width = pwm.GetLength(0); if (from < 0 || from + width > s.Length) return -20;
    double total = 0; for (int i = 0; i < width; i++) total += s[from + i] == 4 ? -3 : pwm[i, s[from + i]];
    return total;
}

// Every GT/AG cut in [rangeFrom, rangeTo) of one strand with check (e)'s model 3 features (site score, 18 Z-curve bins, 24 shape bins).
void Sites(byte[] s, int rangeFrom, int rangeTo, double[,] donorPwm, double[,] acceptorPwm, Action<bool, int, double[]> visit)
{
    int length = s.Length, size = Math.Min(chunk, Math.Max(1, rangeTo - rangeFrom)) + 2 * reach + 8;
    var prefix = new double[tracks][]; for (int t = 0; t < tracks; t++) prefix[t] = new double[size];
    var v = new double[tracks]; var x = new double[features];
    for (int start = rangeFrom; start < rangeTo; start += chunk)
    {
        int end = Math.Min(rangeTo, start + chunk), from = Math.Max(0, start - reach - 2), to = Math.Min(length, end + reach + 2), n = to - from;
        for (int i = 0; i < n; i++)
        {
            int p = from + i; byte b = s[p];
            Array.Clear(v);
            if (b < 4) { v[0] = b is 0 or 2 ? 1 : -1; v[1] = b is 0 or 1 ? 1 : -1; v[2] = b is 0 or 3 ? 1 : -1; }
            if (p >= 2 && p + 2 < length)
            {
                int k = 0; bool ok = true;
                for (int j = p - 2; j <= p + 2; j++) { if (s[j] == 4) { ok = false; break; } k = k * 4 + s[j]; }
                if (ok) for (int t = 0; t < 4; t++) v[3 + t] = shape[k, t] - shapeMean[t];
            }
            for (int t = 0; t < tracks; t++) prefix[t][i + 1] = prefix[t][i] + v[t];
        }
        for (int i = Math.Max(start - 2, 0); i < end; i++)
        {
            bool donor = IsGt(s, i); if (!donor && !IsAg(s, i)) continue;
            int cut = donor ? i : i + 2; if (cut < start || cut >= end) continue;
            x[0] = donor ? Score(donorPwm, s, cut - donorBefore) : Score(acceptorPwm, s, cut - acceptorBefore);
            for (int bin = 0; bin < 6; bin++)
            {
                int lo = Math.Clamp(cut + binFrom[bin], from, to), hi = Math.Clamp(cut + binTo[bin], from, to), w = Math.Max(1, hi - lo);
                for (int t = 0; t < tracks; t++) x[1 + t * 6 + bin] = (prefix[t][hi - from] - prefix[t][lo - from]) / w;
            }
            visit(donor, cut, x);
        }
    }
}

static double[] Fit(List<(double[] X, bool Y)> rows)
{
    int d = rows[0].X.Length + 1;
    var w = new double[d];
    for (int iteration = 0; iteration < 20; iteration++)
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
        var step = Solve(h, g);
        double change = 0; for (int i = 0; i < d; i++) { w[i] += step[i]; change += Math.Abs(step[i]); }
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

// ---- pass 1: genes, PWM windows, coding counts ----------------------------------------------------------------------------------------------
Trained Prepare(Genome g)
{
    var background = new double[4]; long total = 0;
    foreach (var s in g.Sequences.Values) foreach (var b in s) if (b < 4) { background[b]++; total++; }
    for (int b = 0; b < 4; b++) background[b] /= total;
    var ids = g.Sequences.Keys.OrderBy(id => id, StringComparer.Ordinal).ToList();
    var tasks = ids.SelectMany((id, index) => new[] { (Id: id, Index: index, Reverse: false), (Id: id, Index: index, Reverse: true) }).ToArray();
    var parts = new PrepareParts[tasks.Length];
    var byStrand = g.Longest.GroupBy(kv => (kv.Value.SeqId, kv.Value.Complement)).ToDictionary(x => x.Key, x => x.ToList());
    Parallel.For(0, tasks.Length, ti =>
    {
        var (id, _, reverse) = tasks[ti];
        var forward = g.Sequences[id]; var s = reverse ? ReverseComplement(forward) : forward; int length = s.Length;
        var part = new PrepareParts(); parts[ti] = part;
        var (donors, acceptors) = Truth(g, id, reverse, s);
        foreach (int b in donors) if (b - donorBefore >= 0 && b + donorAfter <= length) part.DonorWindows.Add(Enumerable.Range(b - donorBefore, donorBefore + donorAfter).Select(i => (int)s[i]).ToArray());
        foreach (int a in acceptors) if (a - acceptorBefore >= 0 && a + acceptorAfter <= length) part.AcceptorWindows.Add(Enumerable.Range(a - acceptorBefore, acceptorBefore + acceptorAfter).Select(i => (int)s[i]).ToArray());
        for (int i = 0; i + 1 < length; i++) { if (IsGt(s, i)) part.GtCount++; else if (IsAg(s, i)) part.AgCount++; }
        if (!byStrand.TryGetValue((id, reverse), out var genes)) return;
        foreach (var (name, (_, _, gp, _)) in genes)
        {
            part.Tried++;
            var exons = gp.Select(p => reverse ? ((int)(length - p.End), (int)(length - p.Start)) : ((int)p.Start, (int)p.End)).OrderBy(e => e.Item1).ToArray();
            bool canonical = true;
            for (int k = 0; k + 1 < exons.Length; k++)
                if (exons[k + 1].Item1 - exons[k].Item2 < 4 || !IsGt(s, exons[k].Item2) || !IsAg(s, exons[k + 1].Item1 - 2)) canonical = false;
            if (!canonical) { part.NonCanonical++; continue; }
            var cds = new List<byte>(); var ends = new List<int>();
            foreach (var (a, b) in exons) { for (int i = a; i < b; i++) cds.Add(s[i]); ends.Add(cds.Count); }
            var c = cds.ToArray();
            bool ok = c.Length % 3 == 0 && c.Length >= 6 && CodonAt(c, 0) == 0 * 16 + 3 * 4 + 2 && IsStopCodon(CodonAt(c, c.Length - 3)) && c.All(b => b < 4);
            for (int i = 0; ok && i + 3 < c.Length; i += 3) if (IsStopCodon(CodonAt(c, i))) ok = false;
            if (!ok) { part.BadFrame++; continue; }
            int sliceStart = Math.Max(0, exons[0].Item1 - pad), sliceEnd = Math.Min(length, exons[^1].Item2 + pad + 3);
            var pend = new int[exons.Length - 1];
            for (int k = 0; k + 1 < exons.Length; k++)
            {
                int lk = ends[k], r = lk % 3, st = 0;
                for (int i = lk - r; i < lk; i++) st = Extend(st, c[i]);
                pend[k] = st;
            }
            part.Genes.Add(new Gene
            {
                Name = name, Slice = s[sliceStart..sliceEnd], G0 = exons[0].Item1 - sliceStart, G1 = exons[^1].Item2 - sliceStart,
                Exons = exons.Select(e => (e.Item1 - sliceStart, e.Item2 - sliceStart)).ToArray(), PendAfter = pend,
            });
            for (int i = 0; i + 3 < c.Length; i += 3) { part.Codon[CodonAt(c, i)]++; for (int j = 0; j < 3; j++) part.Position[j * 4 + c[i + j]]++; }
            foreach (var (a, b) in exons) part.ExonLengths.Add(b - a);
            for (int k = 0; k + 1 < exons.Length; k++)
            {
                int from = exons[k].Item2, to = exons[k + 1].Item1; part.IntronLengths.Add(to - from);
                for (int i = from; i + 3 <= to; i++) { int t3 = CodonAt(s, i); if (t3 >= 0) part.Intron3[t3]++; }
                for (int i = from; i < to; i++) if (s[i] < 4) part.Intron1[s[i]]++;
            }
        }
    });
    var tr = new Trained { Label = g.Label, Background = background };
    tr.Genes = parts.SelectMany(p => p.Genes).OrderBy(x => x.Name, StringComparer.Ordinal).ToList();
    tr.DonorPwm = Pwm(parts.SelectMany(p => p.DonorWindows), donorBefore + donorAfter, background);
    tr.AcceptorPwm = Pwm(parts.SelectMany(p => p.AcceptorWindows), acceptorBefore + acceptorAfter, background);
    tr.GtCount = parts.Sum(p => p.GtCount); tr.AgCount = parts.Sum(p => p.AgCount);
    var codon = new double[64]; var position = new double[12]; var intron3 = new double[64]; var intron1 = new double[4];
    foreach (var p in parts) { for (int i = 0; i < 64; i++) { codon[i] += p.Codon[i]; intron3[i] += p.Intron3[i]; } for (int i = 0; i < 12; i++) position[i] += p.Position[i]; for (int i = 0; i < 4; i++) intron1[i] += p.Intron1[i]; }
    double codonTotal = codon.Sum() + 64, intron3Total = intron3.Sum() + 64, intron1Total = intron1.Sum() + 4;
    tr.CodonTable = new double[64]; tr.ZCurveTable = new double[64];
    for (int c = 0; c < 64; c++)
    {
        int[] b = [c / 16, c / 4 % 4, c % 4];
        tr.CodonTable[c] = Math.Log((codon[c] + 1) / codonTotal) - Math.Log((intron3[c] + 1) / intron3Total);
        for (int j = 0; j < 3; j++)
        {
            double positionTotal = position[j * 4] + position[j * 4 + 1] + position[j * 4 + 2] + position[j * 4 + 3] + 4;
            tr.ZCurveTable[c] += Math.Log((position[j * 4 + b[j]] + 1) / positionTotal) - Math.Log((intron1[b[j]] + 1) / intron1Total);
        }
    }
    var introns = parts.SelectMany(p => p.IntronLengths).Order().ToArray(); var exonLengths = parts.SelectMany(p => p.ExonLengths).Order().ToArray();
    tr.MinIntron = Math.Max(4, introns[(int)(0.01 * introns.Length)]);
    int Q(int[] a, double q) => a[(int)(q * (a.Length - 1))];
    Log($"{g.Label}: genes tried {parts.Sum(p => p.Tried):N0}, kept {tr.Genes.Count:N0} (non-canonical intron {parts.Sum(p => p.NonCanonical):N0}, " +
        $"frame/ATG/stop {parts.Sum(p => p.BadFrame):N0}); exons q1/median/q3 {Q(exonLengths, .25)}/{Q(exonLengths, .5)}/{Q(exonLengths, .75)} bp, " +
        $"introns q1/median/q3 {Q(introns, .25)}/{Q(introns, .5)}/{Q(introns, .75)} bp, 1st percentile {tr.MinIntron} bp, {introns.Count(i => i < 500) / (double)introns.Length:P0} under 500 bp");
    tr.Tasks = tasks.Select(t => (t.Id, t.Index, t.Reverse)).ToArray();
    return tr;
}

// ---- pass 2: model 3 over the whole genome (negatives: a uniform sample by hash, so the strands run in parallel) ----------------------------
void Train(Genome g, Trained tr)
{
    double[] keep = [Math.Min(1, 2.0 * negativeSample / tr.GtCount), Math.Min(1, 2.0 * negativeSample / tr.AgCount)];
    var results = new (List<double[]>[] Positives, List<(double U, double[] X)>[] Negatives, long[] Seen)[tr.Tasks.Length];
    Parallel.For(0, tr.Tasks.Length, ti =>
    {
        var (id, index, reverse) = tr.Tasks[ti];
        var forward = g.Sequences[id]; var s = reverse ? ReverseComplement(forward) : forward;
        var (donors, acceptors) = Truth(g, id, reverse, s);
        var positives = new[] { new List<double[]>(), new List<double[]>() }; var negatives = new[] { new List<(double, double[])>(), new List<(double, double[])>() }; var seen = new long[2];
        Sites(s, 0, s.Length, tr.DonorPwm, tr.AcceptorPwm, (donor, cut, x) =>
        {
            int t = donor ? 0 : 1;
            if ((donor ? donors : acceptors).Contains(cut)) { positives[t].Add((double[])x.Clone()); return; }
            seen[t]++;
            double u = U(((ulong)index << 40) ^ ((reverse ? 1UL : 0) << 39) ^ ((ulong)t << 38) ^ (ulong)cut);
            if (u < keep[t]) negatives[t].Add((u, (double[])x.Clone()));
        });
        results[ti] = (positives, negatives, seen);
    });
    for (int t = 0; t < 2; t++)
    {
        var positives = results.SelectMany(r => r.Positives[t]).ToList();
        var negatives = results.SelectMany(r => r.Negatives[t]).OrderBy(n => n.U).Take(negativeSample).Select(n => n.X).ToList();
        long seen = results.Sum(r => r.Seen[t]);
        var rows = positives.Select(x => (x, true)).Concat(negatives.Select(x => (x, false))).ToList();
        var mean = Enumerable.Range(0, features).Select(j => rows.Average(r => r.x[j])).ToArray();
        var sd = Enumerable.Range(0, features).Select(j => { double v = Math.Sqrt(rows.Average(r => (r.x[j] - mean[j]) * (r.x[j] - mean[j]))); return v < 1e-12 ? 1 : v; }).ToArray();
        var standardized = rows.Select(r => (r.x.Select((v, j) => (v - mean[j]) / sd[j]).ToArray(), r.Item2)).ToList();
        var model = new Model(mean, sd, Fit(standardized), Math.Log(negatives.Count / (double)seen));
        if (t == 0) tr.Donor = model; else tr.Acceptor = model;
        Log($"   model 3 {(t == 0 ? "donor" : "acceptor")} trained on {g.Label}: real {positives.Count:N0}, false {seen:N0} (sample {negatives.Count:N0})");
    }
}

static double Logit(Model m, double[] x) { double z = m.W[0]; for (int j = 0; j < x.Length; j++) z += m.W[j + 1] * (x[j] - m.Mean[j]) / m.Sd[j]; return z + m.LogRate; }

// ---- pass 3: site log-odds inside each gene span, under one genome's model -----------------------------------------------------------------
GeneSites[] ScoreGenes(List<Gene> genes, Trained model)
{
    var result = new GeneSites[genes.Count];
    Parallel.For(0, genes.Count, gi =>
    {
        var gene = genes[gi];
        var dp = new List<int>(); var ds = new List<double>(); var ap = new List<int>(); var @as = new List<double>();
        Sites(gene.Slice, gene.G0 + 1, gene.G1, model.DonorPwm, model.AcceptorPwm, (donor, cut, x) =>
        {
            if (donor) { dp.Add(cut); ds.Add(Logit(model.Donor, x)); } else { ap.Add(cut); @as.Add(Logit(model.Acceptor, x)); }
        });
        result[gi] = new GeneSites(dp.ToArray(), ds.ToArray(), ap.ToArray(), @as.ToArray());
    });
    return result;
}

// ---- the joining: best chain of pieces from the gene's start to its end ------------------------------------------------------------------
ParseResult Parse(Gene gene, int[] dPos, double[] dSc, int[] aPos, double[] aSc, double[] codon, int minIntron)
{
    var t = gene.Slice; int g0 = gene.G0, g1 = gene.G1, m = g1 - g0 + 3;
    var nextStop = new int[m + 3];
    for (int k = m + 2; k >= 0; k--) { int i = g0 + k; int c = CodonAt(t, i); nextStop[k] = c >= 0 && IsStopCodon(c) ? i : k + 3 <= m + 2 ? nextStop[k + 3] : none; }
    var cum = new double[m + 3];
    for (int k = 0; k < m; k++) { int c = CodonAt(t, g0 + k); cum[k + 3] = cum[k] + (c >= 0 && !IsStopCodon(c) ? codon[c] : 0); }
    int na = aPos.Length;
    var v = new double[na, 21]; var bestEnd = new int[na, 21];
    var sm = new double[na + 1, 21]; var smAt = new int[na + 1, 21];
    for (int st = 0; st < 21; st++) { sm[na, st] = double.NegativeInfinity; smAt[na, st] = -1; }
    (double, int) From(int e, int st) { int k = LowerBound(aPos, e + minIntron); return (sm[k, st], smAt[k, st]); }
    int PendAfter(int st, int b, int e)
    {
        int p = PendLength(st), need = (3 - p) % 3;
        if (e - b < need) { for (int i = b; i < e; i++) { if (t[i] > 3) return -1; st = Extend(st, t[i]); } return st; }
        int r = (e - b - need) % 3, ns = 0;
        for (int i = e - r; i < e; i++) { if (t[i] > 3) return -1; ns = Extend(ns, t[i]); }
        return ns;
    }
    (double, int) Best(int b, int st)
    {
        double best = double.NegativeInfinity; int bestE = -1;
        int p = PendLength(st), need = (3 - p) % 3;
        if (need > 0)
            for (int j = LowerBound(dPos, b + 1); j < dPos.Length && dPos[j] < b + need; j++)
            {
                int ns = PendAfter(st, b, dPos[j]); if (ns < 0) continue;
                double value = dSc[j] + From(dPos[j], ns).Item1;
                if (value > best) { best = value; bestE = dPos[j]; }
            }
        if (b + need > g1) return (best, bestE);
        double split = 0;
        if (p > 0)
        {
            int c = 0;
            if (st <= 4) c = st - 1; else c = (st - 5) / 4 * 4 + (st - 5) % 4;
            bool unknown = false;
            for (int i = b; i < b + need; i++) { if (t[i] > 3) unknown = true; else c = c * 4 + t[i]; }
            if (!unknown && IsStopCodon(c)) { if (b + need == g1 && 0 > best) { best = 0; bestE = g1; } return (best, bestE); }
            split = unknown ? 0 : codon[c];
        }
        int f = b + need, stop = nextStop[f - g0];
        if (f + 3 <= g1 && (g1 - f) % 3 == 0 && stop == g1 - 3)
        {
            double value = split + cum[g1 - 3 - g0] - cum[f - g0];
            if (value > best) { best = value; bestE = g1; }
        }
        int limit = stop == none ? g1 - 1 : Math.Min(g1 - 1, stop + 2);
        for (int j = LowerBound(dPos, Math.Max(f, b + 1)); j < dPos.Length && dPos[j] <= limit; j++)
        {
            int e = dPos[j], k3 = (e - f) / 3, ns = PendAfter(st, b, e); if (ns < 0) continue;
            double after = From(e, ns).Item1; if (double.IsNegativeInfinity(after)) continue;
            double value = split + cum[f + 3 * k3 - g0] - cum[f - g0] + dSc[j] + after;
            if (value > best) { best = value; bestE = e; }
        }
        return (best, bestE);
    }
    for (int k = na - 1; k >= 0; k--)
        for (int st = 0; st < 21; st++)
        {
            (v[k, st], bestEnd[k, st]) = Best(aPos[k], st);
            double value = aSc[k] + v[k, st];
            if (value > sm[k + 1, st]) { sm[k, st] = value; smAt[k, st] = k; } else { sm[k, st] = sm[k + 1, st]; smAt[k, st] = smAt[k + 1, st]; }
        }
    var (startValue, startEnd) = Best(g0, 0);
    var predicted = new List<(int, int)>();
    if (!double.IsNegativeInfinity(startValue))
    {
        int b = g0, st = 0, e = startEnd;
        while (true)
        {
            predicted.Add((b, e)); if (e == g1) break;
            int ns = PendAfter(st, b, e); var (_, k) = From(e, ns);
            b = aPos[k]; st = ns; e = bestEnd[k, st];
        }
    }
    // Given the true chain up to each true donor, is the true acceptor the best next piece?
    int acceptorFirst = 0;
    for (int i = 0; i + 1 < gene.Exons.Length; i++)
    {
        int d = gene.Exons[i].End, a = gene.Exons[i + 1].Start, st = gene.PendAfter[i], at = Array.BinarySearch(aPos, a);
        if (at < 0) continue;
        double truth = aSc[at] + v[at, st]; if (double.IsNegativeInfinity(truth)) continue;
        bool first = true;
        for (int k = LowerBound(aPos, d + minIntron); k < na && first; k++) if (k != at && aSc[k] + v[k, st] >= truth) first = false;
        if (first && a >= d + minIntron) acceptorFirst++;
    }
    return new ParseResult(predicted, acceptorFirst);
}

// ---- running one gene set under one model ---------------------------------------------------------------------------------------------------
List<Row> Evaluate(string pair, string set, Trained test, Trained model)
{
    var genes = test.Genes;
    var sites = ScoreGenes(genes, model);
    double Threshold(Func<Gene, GeneSites, IEnumerable<double>> pick)
    {
        var scores = genes.Zip(sites).SelectMany(x => pick(x.First, x.Second)).Order().ToArray();
        return scores[(int)((1 - recall) * scores.Length)];
    }
    double thresholdD = Threshold((g, s) => g.Exons[..^1].Select(e => s.DScore[Array.BinarySearch(s.DPos, e.End)]));
    double thresholdA = Threshold((g, s) => g.Exons[1..].Select(e => s.AScore[Array.BinarySearch(s.APos, e.Start)]));
    var configs = new List<(string Table, double W, double[] Codon)> { ("frame", 0, new double[64]) };
    foreach (double w in weights) { configs.Add(("zcurve", w, model.ZCurveTable.Select(c => c * w).ToArray())); configs.Add(("codon", w, model.CodonTable.Select(c => c * w).ToArray())); }
    var rows = new List<Row>();
    long spanBp = genes.Sum(g => (long)(g.G1 - g.G0)), allD = sites.Sum(s => (long)s.DPos.Length), allA = sites.Sum(s => (long)s.APos.Length);
    long keptD = sites.Sum(s => (long)s.DScore.Count(x => x >= thresholdD)), keptA = sites.Sum(s => (long)s.AScore.Count(x => x >= thresholdA));
    foreach (bool oracle in new[] { false, true })
    {
        var candidates = new (int[] D, double[] DS, int[] A, double[] AS, bool AllKept)[genes.Count];
        for (int gi = 0; gi < genes.Count; gi++)
        {
            var g = genes[gi]; var s = sites[gi];
            var trueD = g.Exons[..^1].Select(e => e.End).ToHashSet(); var trueA = g.Exons[1..].Select(e => e.Start).ToHashSet();
            var di = Enumerable.Range(0, s.DPos.Length).Where(i => s.DScore[i] >= thresholdD || (oracle && trueD.Contains(s.DPos[i]))).ToArray();
            var ai = Enumerable.Range(0, s.APos.Length).Where(i => s.AScore[i] >= thresholdA || (oracle && trueA.Contains(s.APos[i]))).ToArray();
            bool allKept = trueD.All(d => s.DScore[Array.BinarySearch(s.DPos, d)] >= thresholdD) && trueA.All(a => s.AScore[Array.BinarySearch(s.APos, a)] >= thresholdA);
            candidates[gi] = (di.Select(i => s.DPos[i]).ToArray(), di.Select(i => s.DScore[i]).ToArray(), ai.Select(i => s.APos[i]).ToArray(), ai.Select(i => s.AScore[i]).ToArray(), allKept);
        }
        foreach (var (table, w, codon) in configs)
        {
            var outcomes = new (bool Parsed, bool Exact, int Introns, int Predicted, int Hit, int AcceptorFirst)[genes.Count];
            Parallel.For(0, genes.Count, gi =>
            {
                var g = genes[gi]; var c = candidates[gi];
                var r = Parse(g, c.D, c.DS, c.A, c.AS, codon, model.MinIntron);
                var trueIntrons = Enumerable.Range(0, g.Exons.Length - 1).Select(i => (g.Exons[i].End, g.Exons[i + 1].Start)).ToHashSet();
                var predictedIntrons = Enumerable.Range(0, Math.Max(0, r.Exons.Count - 1)).Select(i => (r.Exons[i].Item2, r.Exons[i + 1].Item1)).ToList();
                bool exact = r.Exons.Count == g.Exons.Length && r.Exons.Zip(g.Exons).All(x => x.First.Item1 == x.Second.Start && x.First.Item2 == x.Second.End);
                outcomes[gi] = (r.Exons.Count > 0, exact, trueIntrons.Count, predictedIntrons.Count, predictedIntrons.Count(trueIntrons.Contains), r.AcceptorFirst);
            });
            int introns = outcomes.Sum(o => o.Introns), predicted = outcomes.Sum(o => o.Predicted), hit = outcomes.Sum(o => o.Hit);
            double ExactIn(int lo, int hi) { var x = outcomes.Where(o => o.Introns >= lo && o.Introns <= hi).ToList(); return x.Count == 0 ? double.NaN : x.Count(o => o.Exact) / (double)x.Count; }
            rows.Add(new Row(pair, test.Label, model.Label, set, oracle ? "oracle" : "kept", table, w, genes.Count, candidates.Count(c => c.AllKept), outcomes.Count(o => o.Parsed),
                outcomes.Count(o => o.Exact), outcomes.Count(o => o.Exact) / (double)genes.Count, introns, hit / (double)introns, predicted == 0 ? 0 : hit / (double)predicted,
                outcomes.Sum(o => o.AcceptorFirst) / (double)introns, ExactIn(1, 1), ExactIn(2, 4), ExactIn(5, 9), ExactIn(10, int.MaxValue),
                keptD / (double)allD, keptA / (double)allA, (keptD + keptA) / (spanBp / 1000.0)));
        }
        Log($"   {set} {test.Label} under {model.Label} model, {(oracle ? "oracle" : "kept")}: done");
    }
    return rows;
}

var allRows = new List<Row>();
var header = "pair,test,trained_on,set,mode,table,w,chosen,genes,genes_all_sites_kept,genes_parsed,exact_genes,exact_share,introns,intron_sn,intron_sp," +
             "acceptor_first_share,exact_1_intron,exact_2_4,exact_5_9,exact_10_plus,donor_kept_share,acceptor_kept_share,candidates_per_kbp";
void Write()
{
    Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(output))!);
    File.WriteAllLines(output, new[] { header }.Concat(allRows.Select(r => r.Csv())));
}
foreach (var (pair, (labelA, prefixA), (labelB, prefixB)) in pairs)
{
    if (only is not null && !pair.Equals(only, StringComparison.OrdinalIgnoreCase)) continue;
    var trained = new List<Trained>();
    foreach (var (label, prefix) in new[] { (labelA, prefixA), (labelB, prefixB) })
    {
        var genome = Load(label, prefix);
        var tr = Prepare(genome);
        Train(genome, tr);
        genome.Sequences.Clear(); trained.Add(tr);
        GC.Collect();
    }
    foreach (var (train, test) in new[] { (trained[0], trained[1]), (trained[1], trained[0]) })
    {
        var tune = Evaluate(pair, "tune", train, train);
        var result = Evaluate(pair, "test", test, train);
        foreach (var table in new[] { "zcurve", "codon" })
        {
            var best = tune.Where(r => r.Mode == "kept" && r.Table == table).OrderByDescending(r => r.ExactShare).ThenByDescending(r => r.IntronSn).First();
            foreach (var r in result.Where(r => r.Table == table && r.W == best.W)) r.Chosen = true;
        }
        foreach (var r in result.Where(r => r.Table == "frame")) r.Chosen = true;
        allRows.AddRange(tune); allRows.AddRange(result);
        Write();
        foreach (var r in result.Where(r => r.Chosen))
            Log($"== {r.Test} (trained on {r.TrainedOn}) {r.Mode,-6} {r.Table,-6} w={r.W}: exact {r.ExactShare:P1} of {r.Genes:N0} (all sites kept {r.GenesAllSitesKept / (double)r.Genes:P0}); " +
                $"introns Sn {r.IntronSn:P1} Sp {r.IntronSp:P1}; true acceptor first {r.AcceptorFirstShare:P1}; exact by introns 1 {r.Exact1:P0}, 2-4 {r.Exact2:P0}, 5-9 {r.Exact5:P0}, 10+ {r.Exact10:P0}; {r.CandidatesPerKbp:F1} cut points/kbp");
    }
}
Write();
var judged = allRows.Where(r => r.Pair == "Nelumbo" && r.Set == "test" && r.Mode == "kept" && r.Table == "zcurve" && r.Chosen).ToList();
if (judged.Count == 2)
    Log(judged.All(r => r.IntronSn >= 0.5 && r.ExactShare >= 0.2)
        ? "Rule: PASSES (Nelumbo, Z-curve, both directions: intron Sn >= 50% and exact genes >= 20%)"
        : "Rule: FAILS (Nelumbo, Z-curve: intron Sn < 50% or exact genes < 20% in a direction); the puzzle-piece path is abandoned");
Log($"wrote {output}");

sealed record Model(double[] Mean, double[] Sd, double[] W, double LogRate);
sealed record GeneSites(int[] DPos, double[] DScore, int[] APos, double[] AScore);
sealed record ParseResult(List<(int, int)> Exons, int AcceptorFirst);

sealed class Row(string pair, string test, string trainedOn, string set, string mode, string table, double w, int genes, int genesAllSitesKept, int genesParsed,
    int exactGenes, double exactShare, int introns, double intronSn, double intronSp, double acceptorFirstShare, double exact1, double exact2, double exact5, double exact10,
    double donorKept, double acceptorKept, double candidatesPerKbp)
{
    public string Pair = pair, Test = test, TrainedOn = trainedOn, Set = set, Mode = mode, Table = table;
    public double W = w, ExactShare = exactShare, IntronSn = intronSn, IntronSp = intronSp, AcceptorFirstShare = acceptorFirstShare, Exact1 = exact1, Exact2 = exact2, Exact5 = exact5, Exact10 = exact10, CandidatesPerKbp = candidatesPerKbp;
    public int Genes = genes, GenesAllSitesKept = genesAllSitesKept;
    public bool Chosen;
    public string Csv() => string.Create(CultureInfo.InvariantCulture,
        $"{Pair},{Test},{TrainedOn},{Set},{Mode},{Table},{W},{Chosen},{Genes},{GenesAllSitesKept},{genesParsed},{exactGenes},{ExactShare:F4},{introns},{IntronSn:F4},{IntronSp:F4}," +
        $"{AcceptorFirstShare:F4},{Exact1:F4},{Exact2:F4},{Exact5:F4},{Exact10:F4},{donorKept:F4},{acceptorKept:F4},{CandidatesPerKbp:F2}");
}

sealed class Gene
{
    public string Name = "";
    public byte[] Slice = [];
    public int G0, G1;
    public (int Start, int End)[] Exons = [];
    public int[] PendAfter = [];
}

sealed class PrepareParts
{
    public List<int[]> DonorWindows = [], AcceptorWindows = [];
    public List<Gene> Genes = [];
    public List<int> IntronLengths = [], ExonLengths = [];
    public long[] Codon = new long[64], Position = new long[12], Intron3 = new long[64], Intron1 = new long[4];
    public long GtCount, AgCount;
    public int Tried, NonCanonical, BadFrame;
}

sealed class Trained
{
    public string Label = "";
    public double[] Background = [];
    public List<Gene> Genes = [];
    public double[,] DonorPwm = new double[0, 0], AcceptorPwm = new double[0, 0];
    public Model Donor = null!, Acceptor = null!;
    public double[] CodonTable = [], ZCurveTable = [];
    public int MinIntron;
    public long GtCount, AgCount;
    public (string Id, int Index, bool Reverse)[] Tasks = [];
}

sealed class Genome
{
    public string Label = "";
    public Dictionary<string, List<(bool Complement, List<(long Start, long End)> Parts)>> Cds = [];
    public Dictionary<string, (string SeqId, bool Complement, List<(long Start, long End)> Parts, long Length)> Longest = [];
    public Dictionary<string, byte[]> Sequences = [];
}
