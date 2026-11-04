#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using SequencePlatform.Core;

// Check after E2 (2026-11-04): would k-mers longer than the stored dinucleotides predict a phage's host well enough to justify
// re-downloading the catalogue? E2's set and answer key (BV-BRC host_name's first word as a stored genus, provider lineage), a sample
// small enough to download: for each host genus with >= 5 phages, 3 phages and 3 host genomes (each the largest chromosome of its
// genome), seeded. Each sequence is streamed from NCBI (efetch, FASTA) and counted as it arrives: no file, no string, nothing kept
// but the counts, so the platform's no-DNA rule holds.
// Signatures, both strands (a k-mer and its reverse complement counted together), for k = 2..6:
// - markov: log of the observed frequency over its (k-2)-order Markov expectation f(w1..k-1) f(w2..k) / f(w2..k-1); at k = 2 this
//   is log rho*, E2's signature.
// - frequency: the k-mer frequencies themselves.
// Prediction as in E2: the genus of the nearest host genome (Euclidean); pool = the sampled genera. GC alone as context.
// Fixed before the run: a re-download is worth it only if the best k's top-1 genus accuracy beats k = 2 markov by >= 0.05 with the
// bootstrap 95 % lower bound of the difference above 0. Writes <prefix>.csv.
const int perGenus = 3, minPhages = 5, bootstraps = 1000, seed = 20261104, maxK = 6;
var prefix = args.FirstOrDefault() ?? "P9_KmerCheck";

await using var dataSource = ExperimentCorpus.OpenSnapshot();
var corpus = (await ExperimentCorpus.LoadAsync(dataSource)).Where(g => !g.Excluded).ToList();
var cellular = corpus.Where(g => g.Group == SignatureIndex.Cellular && g.Provider.Genus is not null).ToList();
var genera = cellular.Select(g => g.Provider.Genus!).ToHashSet(StringComparer.Ordinal);
string? HostGenus(CorpusGenome g) => g.Field("host_name")?.Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault();
var phages = corpus.Where(g => g.Group == SignatureIndex.Viral && HostGenus(g) is { } h && genera.Contains(h)).ToList();
var rng = new Random(seed);
var pool = phages.GroupBy(p => HostGenus(p)!).Where(g => g.Count() >= minPhages).OrderBy(g => g.Key, StringComparer.Ordinal).ToList();
var phageSample = pool.SelectMany(g => g.OrderBy(_ => rng.Next()).Take(perGenus).Select(p => (Genome: p, Genus: g.Key))).ToList();
var hostSample = pool.SelectMany(g => cellular.Where(c => c.Provider.Genus == g.Key).OrderBy(_ => rng.Next()).Take(perGenus)
    .Select(c => (Genome: c, Genus: g.Key))).ToList();
var accessions = await AccessionsAsync(phageSample.Select(p => p.Genome.Id).Concat(hostSample.Select(h => h.Genome.Id)).ToArray());
Console.WriteLine($"{pool.Count} host genera; {phageSample.Count} phages and {hostSample.Count} host genomes to stream.");

using var http = new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
var watch = Stopwatch.StartNew();
long bases = 0;
async Task<long[][]> CountAsync(long id)
{
    for (int attempt = 1; ; attempt++)
    {
        try
        {
            var url = $"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&tool=cat-platform-p9&id={accessions[id]}";
            await using var stream = await http.GetStreamAsync(url);
            var counts = await KmerCounts.FromFastaAsync(new StreamReader(stream), maxK);
            await Task.Delay(400);   // NCBI asks for at most 3 requests a second without a key
            return counts;
        }
        catch (Exception e) when (attempt < 4)
        {
            Console.WriteLine($"  retry {accessions[id]}: {e.Message}");
            await Task.Delay(2000 * attempt);
        }
    }
}
var phageCounts = new List<long[][]>();
foreach (var (p, i) in phageSample.Select((p, i) => (p, i)))
{
    phageCounts.Add(await CountAsync(p.Genome.Id));
    bases += phageCounts[^1][1].Sum();
}
var hostCounts = new List<long[][]>();
foreach (var (h, i) in hostSample.Select((h, i) => (h, i)))
{
    hostCounts.Add(await CountAsync(h.Genome.Id));
    bases += hostCounts[^1][1].Sum();
    if ((i + 1) % 25 == 0) Console.WriteLine($"  hosts {i + 1}/{hostSample.Count}, {watch.Elapsed.TotalMinutes:F1} min, {bases / 1e6:F0} Mbp");
}
Console.WriteLine($"Streamed {bases / 1e6:F0} Mbp in {watch.Elapsed.TotalMinutes:F1} min.");

// Sanity check: k = 2 from the stream should reproduce the stored dinucleotide counts' GC.
double maxGcGap = phageSample.Select(p => p.Genome).Concat(hostSample.Select(h => h.Genome))
    .Zip(phageCounts.Concat(hostCounts), (g, c) => Math.Abs(g.Gc - KmerCounts.Gc(c[1]))).Max();
Console.WriteLine($"Largest GC gap between streamed and stored counts: {maxGcGap:F5}.");

var truth = phageSample.Select(p => p.Genus).ToArray();
var hostGenus = hostSample.Select(h => h.Genus).ToArray();
var correct = new Dictionary<string, bool[]>();
var top5 = new Dictionary<string, double>();
foreach (var variant in new[] { "markov", "frequency" })
    for (int k = 2; k <= maxK; k++)
    {
        var hosts = hostCounts.Select(c => KmerCounts.Signature(c, k, variant)).ToArray();
        var ranks = phageCounts.Select(c => KmerCounts.Signature(c, k, variant))
            .Select((s, i) => Rank(hosts.Select(h => Squared(s, h)).ToArray(), hostGenus, truth[i])).ToArray();
        correct[$"{variant}{k}"] = ranks.Select(r => r == 1).ToArray();
        top5[$"{variant}{k}"] = ranks.Average(r => r <= 5 ? 1.0 : 0);
    }
{
    var ranks = phageCounts.Select((c, i) => Rank(hostCounts.Select(h => Math.Abs(KmerCounts.Gc(c[1]) - KmerCounts.Gc(h[1]))).ToArray(), hostGenus, truth[i])).ToArray();
    correct["gc"] = ranks.Select(r => r == 1).ToArray();
    top5["gc"] = ranks.Average(r => r <= 5 ? 1.0 : 0);
}

var baseline = correct["markov2"];
var draws = Enumerable.Range(0, bootstraps).Select(_ => Enumerable.Range(0, truth.Length).Select(_ => rng.Next(truth.Length)).ToArray()).ToArray();
var best = correct.Where(c => c.Key != "gc").MaxBy(c => c.Value.Average(v => v ? 1.0 : 0)).Key;
using var csv = File.CreateText($"{prefix}.csv");
csv.WriteLine("Signature, Phages, HostGenomes, Genera, Top1, Top5, DiffVsMarkov2, DiffLow, DiffHigh, Best");
foreach (var (name, hits) in correct)
{
    var differences = draws.Select(d => d.Average(i => (hits[i] ? 1.0 : 0) - (baseline[i] ? 1.0 : 0))).Order().ToArray();
    double difference = hits.Average(v => v ? 1.0 : 0) - baseline.Average(v => v ? 1.0 : 0);
    var line = string.Create(CultureInfo.InvariantCulture,
        $"{name}, {truth.Length}, {hostSample.Count}, {pool.Count}, {hits.Average(v => v ? 1.0 : 0):F4}, {top5[name]:F4}, {difference:F4}, " +
        $"{differences[(int)(0.025 * bootstraps)]:F4}, {differences[(int)(0.975 * bootstraps) - 1]:F4}, {name == best}");
    csv.WriteLine(line);
    Console.WriteLine(line);
}
var bestDiff = draws.Select(d => d.Average(i => (correct[best][i] ? 1.0 : 0) - (baseline[i] ? 1.0 : 0))).Order().ToArray();
double gain = correct[best].Average(v => v ? 1.0 : 0) - baseline.Average(v => v ? 1.0 : 0);
Console.WriteLine($"Best {best}: +{gain:F3} over markov2 (95 % {bestDiff[(int)(0.025 * bootstraps)]:F3} to {bestDiff[(int)(0.975 * bootstraps) - 1]:F3}). " +
    $"Re-download worth it: {gain >= 0.05 && bestDiff[(int)(0.025 * bootstraps)] > 0}.");

// 1-based rank of the true genus when genera are ordered by their nearest genome.
static int Rank(double[] distance, string[] genus, string truth)
{
    var nearest = new Dictionary<string, double>(StringComparer.Ordinal);
    for (int i = 0; i < distance.Length; i++)
        if (!nearest.TryGetValue(genus[i], out var d) || distance[i] < d) nearest[genus[i]] = distance[i];
    return 1 + nearest.Count(n => n.Value < nearest[truth]);
}

static double Squared(double[] a, double[] b)
{
    double sum = 0;
    for (int f = 0; f < a.Length; f++) { double d = a[f] - b[f]; sum += d * d; }
    return sum;
}

async Task<Dictionary<long, string>> AccessionsAsync(long[] ids)
{
    var map = new Dictionary<long, string>();
    await using var command = dataSource.CreateCommand("select id, record_id from sequence_record where id = any($1)");
    command.Parameters.AddWithValue(ids);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) map[reader.GetInt64(0)] = reader.GetString(1);
    return map;
}

static class KmerCounts
{
    // Counts of every k-mer for k = 1..maxK on one strand (index: base-4 number, A C G T = 0..3), read as the text arrives.
    public static async Task<long[][]> FromFastaAsync(TextReader reader, int maxK)
    {
        var counts = Enumerable.Range(0, maxK + 1).Select(k => new long[1 << (2 * k)]).ToArray();
        int code = 0, valid = 0, mask = (1 << (2 * maxK)) - 1;
        var buffer = new char[1 << 16];
        bool header = false, lineStart = true;
        int read;
        while ((read = await reader.ReadAsync(buffer)) > 0)
            for (int i = 0; i < read; i++)
            {
                char c = buffer[i];
                if (c == '\n') { lineStart = true; header = false; continue; }
                if (lineStart && c == '>') { header = true; valid = 0; }
                lineStart = false;
                if (header || c == '\r') continue;
                int b = c switch { 'A' or 'a' => 0, 'C' or 'c' => 1, 'G' or 'g' => 2, 'T' or 't' => 3, _ => -1 };
                if (b < 0) { valid = 0; continue; }
                code = ((code << 2) | b) & mask;
                valid++;
                for (int k = 1; k <= Math.Min(valid, maxK); k++) counts[k][code & ((1 << (2 * k)) - 1)]++;
            }
        return counts;
    }

    public static double Gc(long[] mono) => (double)(mono[1] + mono[2]) / mono.Sum();

    // Both strands: a k-mer's count plus its reverse complement's.
    private static double[] BothStrands(long[] counts, int k)
    {
        var both = new double[counts.Length];
        for (int w = 0; w < counts.Length; w++) both[w] = counts[w] + counts[ReverseComplement(w, k)];
        double total = both.Sum();
        return both.Select(c => (c + 0.5) / (total + 0.5 * both.Length)).ToArray();   // a half-count floor keeps logs finite
    }

    private static int ReverseComplement(int w, int k)
    {
        int r = 0;
        for (int i = 0; i < k; i++) { r = (r << 2) | (3 - (w & 3)); w >>= 2; }
        return r;
    }

    public static double[] Signature(long[][] counts, int k, string variant)
    {
        var f = BothStrands(counts[k], k);
        if (variant == "frequency") return f;
        var prefix = BothStrands(counts[k - 1], k - 1);
        var middle = k >= 3 ? BothStrands(counts[k - 2], k - 2) : null;
        int low = (1 << (2 * (k - 1))) - 1, inner = (1 << (2 * (k - 2))) - 1;
        return f.Select((value, w) =>
        {
            double expected = k == 2
                ? prefix[w >> 2] * prefix[w & 3]
                : prefix[w >> 2] * prefix[w & low] / middle![(w >> 2) & inner];
            return Math.Log(value / expected);
        }).ToArray();
    }
}
