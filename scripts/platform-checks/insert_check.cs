#:project ../../SequenceAnalise/SequenceAnalise.csproj
#:property PublishAot=false
using System.Globalization;
using SequenceAnalyses;

// Cheap check (2026-11-04, after E3): can two genomes' stored 500 bp profiles show what one has that the other lacks? E3 asked whether a
// window looks unusual inside its own genome, a composition question; this asks whether a stretch has a counterpart in a relative,
// which is about order and identity, where the windows did well (localization). Synthetic and controlled, on bundled DNA:
// - Reference: E. coli K-12's profile. Query: K-12 with 16 inserts per copy (2, 5, 10, 20 kbp plus 0-499 bp, so every insert shifts
//   the 500 bp phase downstream; one of each length from each donor), placed at random >= 50 kbp apart; 10 copies (seeded).
// - Donors span the GC range and share no recent ancestry with E. coli: Vibrio cholerae (47.5 %), B. subtilis (43.5 %),
//   P. aeruginosa (66.6 %), S. aureus (32.9 %).
// - Comparative detector: every block of 8 query windows (4 kbp) is placed anywhere in the reference (BestAlignment); its expected
//   diagonal is the median diagonal of confidently placed blocks 8-40 windows upstream, and separately downstream; its local score is
//   the better RefineAlignment (radius 2) at those two diagonals. A window's score is the best local score of the blocks covering it.
// - Composition baseline: E3's GC detector (z residual from the genome mean, averaged over 8 windows, squared over its variance).
// - Thresholds from a control copy with no inserts, rotated by 250 bp so every window is out of phase: the 0.05th percentile of the
//   comparative score, the 99.95th of the GC score. Calls are runs of windows past the threshold.
// - Fixed before the run: pass if inserts of >= 10 kbp are found (a call overlaps them) in >= 80 % of cases with <= 1 false call per
//   Mbp. Writes <prefix>.csv.
const int window = 500, block = 8, near = 8, far = 40, radius = 2, copies = 10, seed = 20261104;
const int minGap = 50_000;
const double controlQuantile = 0.0005;
int[] lengths = [2_000, 5_000, 10_000, 20_000];
var prefix = args.FirstOrDefault() ?? "InsertCheck";
var genomes = Path.Combine("data", "genomes");
var reference = FastaReader.ReadPrimaryRecord(Path.Combine(genomes, "ecoli_k12_mg1655_GCF_000005845.2.fna.gz")).Sequence;
var donors = new (string Name, string File)[]
{
    ("Vibrio cholerae", "vibrio_cholerae_n16961_GCF_900205735.1.fna.gz"),
    ("Bacillus subtilis", "bsubtilis_168_GCF_000009045.1.fna.gz"),
    ("Pseudomonas aeruginosa", "pseudomonas_aeruginosa_pao1_GCF_000006765.1.fna.gz"),
    ("Staphylococcus aureus", "staphylococcus_aureus_nctc8325_GCF_000013425.1.fna.gz"),
}.Select(d => (d.Name, Sequence: FastaReader.ReadPrimaryRecord(Path.Combine(genomes, d.File)).Sequence)).ToArray();
var referenceProfile = new ZCurveWindowedProfile(reference, window);
Console.WriteLine($"K-12 {reference.Length:N0} bp, {referenceProfile.WindowCount} windows; donors " +
    string.Join(", ", donors.Select(d => $"{d.Name} GC {Gc(d.Sequence):P1}")));

var control = new ZCurveWindowedProfile(reference[250..] + reference[..250], window);
var controlComparative = Comparative(control, referenceProfile);
var controlGc = GcScores(control);
double comparativeThreshold = Quantile(controlComparative, controlQuantile);
double gcThreshold = Quantile(controlGc, 1 - controlQuantile);
Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
    $"Control: comparative score median {Quantile(controlComparative, 0.5):F4}, threshold {comparativeThreshold:F4}; GC threshold {gcThreshold:F2}."));

var found = new Dictionary<(string Detector, int Length, string Donor), (int Inserts, int Found)>();
var falseCalls = new Dictionary<string, int>();
double megabases = 0;
for (int copy = 0; copy < copies; copy++)
{
    var rng = new Random(seed + copy);
    var plan = Plan(rng, reference.Length, lengths, donors.Length);
    var (query, inserts) = Build(reference, plan, donors.Select(d => d.Sequence).ToArray(), rng);
    var profile = new ZCurveWindowedProfile(query, window);
    megabases += query.Length / 1e6;
    foreach (var (detector, called) in new[]
    {
        ("comparative", Comparative(profile, referenceProfile).Select(s => s < comparativeThreshold).ToArray()),
        ("gc", GcScores(profile).Select(s => s > gcThreshold).ToArray()),
    })
    {
        var calls = Runs(called);
        foreach (var insert in inserts)
        {
            int from = insert.Start / window, to = (insert.Start + insert.Length - 1) / window;
            bool hit = calls.Any(c => c.From <= to && c.To >= from);
            var key = (detector, insert.Nominal, donors[insert.Donor].Name);
            var (n, f) = found.GetValueOrDefault(key);
            found[key] = (n + 1, f + (hit ? 1 : 0));
        }
        falseCalls[detector] = falseCalls.GetValueOrDefault(detector)
            + calls.Count(c => !inserts.Any(i => c.From <= (i.Start + i.Length - 1) / window && c.To >= i.Start / window));
    }
    Console.WriteLine($"Copy {copy + 1}/{copies}: {query.Length:N0} bp.");
}

using var csv = File.CreateText($"{prefix}.csv");
csv.WriteLine("Detector, InsertKbp, Donor, Inserts, Found, Recall, FalseCalls, FalseCallsPerMbp, Threshold, Pass");
bool pass = false;
foreach (var detector in new[] { "comparative", "gc" })
{
    double perMbp = falseCalls[detector] / megabases;
    double threshold = detector == "comparative" ? comparativeThreshold : gcThreshold;
    foreach (var length in lengths)
    {
        foreach (var donor in donors)
        {
            var (n, f) = found[(detector, length, donor.Name)];
            csv.WriteLine(string.Create(CultureInfo.InvariantCulture, $"{detector}, {length / 1000}, {donor.Name}, {n}, {f}, {(double)f / n:F3}, , , , "));
        }
        var all = donors.Select(d => found[(detector, length, d.Name)]).Aggregate((a, b) => (a.Inserts + b.Inserts, a.Found + b.Found));
        var line = string.Create(CultureInfo.InvariantCulture,
            $"{detector}, {length / 1000}, all, {all.Inserts}, {all.Found}, {(double)all.Found / all.Inserts:F3}, {falseCalls[detector]}, {perMbp:F2}, {threshold:F4}, ");
        csv.WriteLine(line);
        Console.WriteLine(line);
    }
    var large = lengths.Where(l => l >= 10_000).SelectMany(l => donors.Select(d => found[(detector, l, d.Name)]))
        .Aggregate((a, b) => (a.Inserts + b.Inserts, a.Found + b.Found));
    bool detectorPass = (double)large.Found / large.Inserts >= 0.8 && perMbp <= 1;
    if (detector == "comparative") pass = detectorPass;
    var summary = string.Create(CultureInfo.InvariantCulture,
        $"{detector}, >=10, all, {large.Inserts}, {large.Found}, {(double)large.Found / large.Inserts:F3}, {falseCalls[detector]}, {perMbp:F2}, {threshold:F4}, {detectorPass}");
    csv.WriteLine(summary);
    Console.WriteLine(summary);
}
Console.WriteLine($"Verdict: {(pass ? "pass" : "fail")}.");

// Per query window, the best local score of the blocks covering it (NaN-free: a block that can't be placed scores 0).
static double[] Comparative(ZCurveWindowedProfile query, ZCurveWindowedProfile reference)
{
    var windows = query.Windows;
    int blocks = query.WindowCount - block + 1;
    var parts = new ZCurveWindowedProfile[blocks];
    var anywhere = new (double Score, int Offset)[blocks];
    Parallel.For(0, blocks, i =>
    {
        parts[i] = ZCurveWindowedProfile.FromPrecomputedWindows(window, windows.Skip(i).Take(block).ToArray());
        anywhere[i] = ZCurveWindowedProfile.BestAlignment(parts[i], reference);
    });
    double confident = Quantile(anywhere.Select(a => a.Score).ToArray(), 0.5);
    var local = new double[blocks];
    Parallel.For(0, blocks, i =>
    {
        double best = 0;
        foreach (var (from, to) in new[] { (i - far, i - near), (i + near, i + far) })
        {
            var diagonals = Enumerable.Range(Math.Max(0, from), Math.Max(0, Math.Min(blocks - 1, to) - Math.Max(0, from) + 1))
                .Where(j => anywhere[j].Score >= confident).Select(j => anywhere[j].Offset - j).Order().ToArray();
            if (diagonals.Length == 0) continue;
            var (score, offset) = ZCurveWindowedProfile.RefineAlignment(parts[i], reference, i + diagonals[diagonals.Length / 2], radius);
            if (offset >= 0 && score > best) best = score;
        }
        local[i] = best;
    });
    var result = new double[query.WindowCount];
    for (int w = 0; w < result.Length; w++)
        for (int i = Math.Max(0, w - block + 1); i <= Math.Min(blocks - 1, w); i++)
            result[w] = Math.Max(result[w], local[i]);
    return result;
}

// E3's GC detector against the genome mean: the z residual averaged over 8 windows, squared over its variance.
static double[] GcScores(ZCurveWindowedProfile profile)
{
    var z = profile.Windows.Select(w => w.Z).ToArray();
    double mean = z.Average();
    int n = z.Length;
    var smoothed = new double[n];
    for (int i = 0; i < n; i++)
    {
        int from = Math.Max(0, i - block / 2), to = Math.Min(n, i - block / 2 + block);
        double sum = 0;
        for (int j = from; j < to; j++) sum += z[j] - mean;
        smoothed[i] = sum / (to - from);
    }
    double variance = smoothed.Average(s => s * s);
    return smoothed.Select(s => s * s / variance).ToArray();
}

static List<(int From, int To)> Runs(bool[] called)
{
    var runs = new List<(int, int)>();
    for (int i = 0; i < called.Length; i++)
    {
        if (!called[i]) continue;
        int j = i;
        while (j + 1 < called.Length && called[j + 1]) j++;
        runs.Add((i, j));
        i = j;
    }
    return runs;
}

// Insert positions on the reference, sorted and >= minGap apart, each with a nominal length and a donor.
static List<(int Position, int Nominal, int Donor)> Plan(Random rng, int referenceLength, int[] lengths, int donorCount)
{
    var kinds = lengths.SelectMany(l => Enumerable.Range(0, donorCount).Select(d => (l, d))).ToArray();
    rng.Shuffle(kinds);
    var positions = new List<int>();
    while (positions.Count < kinds.Length)
    {
        int p = rng.Next(minGap, referenceLength - minGap);
        if (positions.All(q => Math.Abs(q - p) >= minGap)) positions.Add(p);
    }
    positions.Sort();
    return positions.Select((p, k) => (p, kinds[k].l, kinds[k].d)).ToList();
}

// The query DNA and where each insert landed in it.
static (string Query, List<(int Start, int Length, int Nominal, int Donor)> Inserts) Build(string reference,
    List<(int Position, int Nominal, int Donor)> plan, string[] donors, Random rng)
{
    var query = new System.Text.StringBuilder(reference.Length + plan.Sum(p => p.Nominal + window));
    var inserts = new List<(int, int, int, int)>();
    int at = 0;
    foreach (var (position, nominal, donor) in plan)
    {
        query.Append(reference, at, position - at);
        int length = nominal + rng.Next(window);
        int start = rng.Next(donors[donor].Length - length);
        inserts.Add((query.Length, length, nominal, donor));
        query.Append(donors[donor], start, length);
        at = position;
    }
    query.Append(reference, at, reference.Length - at);
    return (query.ToString(), inserts);
}

static double Quantile(double[] values, double q)
{
    var sorted = values.Order().ToArray();
    return sorted[Math.Clamp((int)(q * sorted.Length), 0, sorted.Length - 1)];
}

static double Gc(string sequence) => (double)sequence.Count(c => c is 'G' or 'C') / sequence.Count(c => c is 'A' or 'C' or 'G' or 'T');
