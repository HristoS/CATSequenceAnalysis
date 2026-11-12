#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Text.RegularExpressions;
using SequencePlatform.Core;

// Checks before P9 S9.10 (E1b, misplaced genomes on amino-acid CVTree neighbours). Reads the round-2 corpus, writes nothing to the
// database. Fixed before the run (2026-11-12):
// 1. Keys. E1's keys as E1 defines them, over the genomes E1 scores (provider family with >= 5 genomes in the group). A positive in a
//    provider family where most scored genomes are positive is a wholesale rename, which no neighbour score can flag. A key with
//    >= 70 % of its positives in such families can't be tested as it stands; its minority positives (families where positives are
//    < 50 %) are what a fixed key would keep. Stop E1b if no key keeps >= 100 minority positives in >= 20 families.
// 2. Timing. 2,000 genomes per group with amino-acid 4-mers (seeded draw): read the payloads, build the hashed CVTree vectors (8,192
//    dimensions, 3-mers and 2-mers derived from the 4-mers, as CvTree_Derived3mers), kNN k = 10 by cosine on 4 threads. Extrapolated
//    to the group: reading and building linear, kNN quadratic. Stop (or decide a subsample first) if the total is over 1 h.
// Writes <prefix>_Keys.csv and <prefix>_Timing.csv.
const int minFamily = 5, sample = 2000, k = 10, dimensions = 8192, seed = 20261112;
var prefix = args.FirstOrDefault() ?? "P9_E1b_Check";
var suffix = new Regex(@"_[A-Z]+$", RegexOptions.Compiled);
string? Strip(string? name) => name is null ? null : suffix.Replace(name, "");

await using var dataSource = ExperimentCorpus.OpenLive();
var corpus = (await ExperimentCorpus.LoadAsync(dataSource, ExperimentCorpus.RoundTwoMaxId)).Where(g => !g.Excluded && g.Group != SignatureIndex.Eukaryotic).ToList();
var withKmers = new HashSet<long>();
await using (var command = dataSource.CreateCommand("select sequence_id from protein_kmer where k = 4 and status = 'done' and sequence_id <= $1"))
{
    command.Parameters.AddWithValue(ExperimentCorpus.RoundTwoMaxId);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) withKmers.Add(reader.GetInt64(0));
}
Console.WriteLine($"{corpus.Count} genomes (eukaryotes and partial records left out), {corpus.Count(g => withKmers.Contains(g.Id))} with amino-acid 4-mers.");

// 1. Keys.
using (var csv = File.CreateText($"{prefix}_Keys.csv"))
{
    csv.WriteLine("Group, Key, Scored, WithKey, Positives, WholesalePositives, WholesaleShare, MinorityPositives, MinorityFamilies, " +
        "MinorityPositivesWithKmers, Testable, IsolatedPositives, IsolatedFamilies");
    foreach (var group in new[] { SignatureIndex.Cellular, SignatureIndex.Viral })
    {
        var members = corpus.Where(g => g.Group == group).ToList();
        var size = members.Where(g => g.Provider.Family is not null).GroupBy(g => g.Provider.Family!).ToDictionary(f => f.Key, f => f.Count());
        var scored = members.Where(g => g.Provider.Family is { } f && size[f] >= minFamily).ToList();
        var keys = group == SignatureIndex.Cellular
            ? new (string Name, Func<CorpusGenome, bool?> Positive)[]
            {
                ("gtdbFamilyDiffers", g => g.Gtdb?.Family is { } gf ? Strip(gf) != g.Provider.Family : null),
                ("checkmPoor", g => Number(g.Field("checkm_contamination")) is { } c && Number(g.Field("checkm_completeness")) is { } m ? c >= 5 || m < 90 : null),
            }
            : [("ncbiFamilyDiffers", g => g.Ncbi?.Family is { } nf ? nf != g.Provider.Family : null)];
        foreach (var (name, positive) in keys)
        {
            var withKey = scored.Where(g => positive(g) is not null).ToList();
            var rate = withKey.GroupBy(g => g.Provider.Family!).ToDictionary(f => f.Key, f => f.Count(g => positive(g) == true) / (double)f.Count());
            var positives = withKey.Where(g => positive(g) == true).ToList();
            int wholesale = positives.Count(g => rate[g.Provider.Family!] > 0.5);
            var minority = positives.Where(g => rate[g.Provider.Family!] < 0.5).ToList();
            int minorityFamilies = minority.Select(g => g.Provider.Family).Distinct().Count();
            double wholesaleShare = positives.Count == 0 ? 0 : wholesale / (double)positives.Count;
            bool testable = minority.Count >= 100 && minorityFamilies >= 20;
            // Added after the first run, for designing E1b only (not part of the rule above): a minority positive whose provider genus
            // is also mostly negative, so a renamed genus doesn't count as a misplaced genome.
            var genusRate = withKey.Where(g => g.Provider.Genus is not null).GroupBy(g => g.Provider.Genus!)
                .ToDictionary(f => f.Key, f => f.Count(g => positive(g) == true) / (double)f.Count());
            var isolated = minority.Where(g => g.Provider.Genus is not { } genus || genusRate[genus] < 0.5).ToList();
            var line = string.Create(CultureInfo.InvariantCulture,
                $"{group}, {name}, {scored.Count}, {withKey.Count}, {positives.Count}, {wholesale}, {wholesaleShare:F3}, {minority.Count}, {minorityFamilies}, " +
                $"{minority.Count(g => withKmers.Contains(g.Id))}, {testable}, {isolated.Count}, {isolated.Select(g => g.Provider.Family).Distinct().Count()}");
            csv.WriteLine(line);
            Console.WriteLine(line);
            foreach (var family in positives.GroupBy(g => g.Provider.Family!).OrderByDescending(f => f.Count()).Take(5))
                Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
                    $"    {family.Key}: {family.Count()} positives, rate {rate[family.Key]:F2}, e.g. {Strip(family.First().Gtdb?.Family) ?? family.First().Ncbi?.Family}"));
        }
    }
}

// 2. Timing.
using (var csv = File.CreateText($"{prefix}_Timing.csv"))
{
    csv.WriteLine("Group, GenomesWithKmers, Sample, ReadSeconds, BuildSeconds, KnnSeconds, MegabytesRead, ExtrapolatedMinutes");
    double totalMinutes = 0;
    foreach (var group in new[] { SignatureIndex.Cellular, SignatureIndex.Viral })
    {
        var eligible = corpus.Where(g => g.Group == group && withKmers.Contains(g.Id)).Select(g => g.Id).Order().ToArray();
        var ids = eligible.OrderBy(_ => 0).ToArray();
        new Random(seed).Shuffle(ids);
        ids = ids.Take(sample).ToArray();

        var watch = Stopwatch.StartNew();
        var payloads = new List<byte[]>();
        await using (var command = dataSource.CreateCommand("select payload from protein_kmer where k = 4 and status = 'done' and sequence_id = any($1)"))
        {
            command.Parameters.AddWithValue(ids);
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) payloads.Add(reader.GetFieldValue<byte[]>(0));
        }
        double read = watch.Elapsed.TotalSeconds;
        long bytes = payloads.Sum(p => (long)p.Length);

        watch.Restart();
        var vectors = payloads.AsParallel().AsOrdered().WithDegreeOfParallelism(4).Select(p => Normalize(Hash(Composition(ProteinKmers.Decode(p, 4))))).ToArray();
        double build = watch.Elapsed.TotalSeconds;

        watch.Restart();
        var neighbours = Neighbours(vectors);
        double knn = watch.Elapsed.TotalSeconds;

        double scale = eligible.Length / (double)ids.Length;
        double minutes = ((read + build) * scale + knn * scale * scale) / 60;
        totalMinutes += minutes;
        var line = string.Create(CultureInfo.InvariantCulture,
            $"{group}, {eligible.Length}, {ids.Length}, {read:F1}, {build:F1}, {knn:F1}, {bytes / 1e6:F0}, {minutes:F1}");
        csv.WriteLine(line);
        Console.WriteLine(line + $"  (first neighbours of genome 0: {string.Join(" ", neighbours[0].Take(3))})");
    }
    Console.WriteLine($"Extrapolated total {totalMinutes:F1} min (limit 60). Within: {totalMinutes <= 60}.");
}

static double? Number(string? text) => double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out var value) ? value : null;

// Hao's composition vector of the 4-mers, with the 3-mers and 2-mers derived from the 4-mer counts (ProteinKmers.Shorter).
static float[] Composition(int[] words)
{
    var shorterCounts = ProteinKmers.Shorter(words);
    var shortestCounts = ProteinKmers.Shorter(shorterCounts);
    double total = words.Sum(w => (double)w), shorterTotal = shorterCounts.Sum(w => (double)w), shortestTotal = shortestCounts.Sum(w => (double)w);
    var vector = new float[words.Length];
    int suffixSpace = shorterCounts.Length, middleSpace = shortestCounts.Length;
    for (int w = 0; w < words.Length; w++)
    {
        int prefix = w / 20, suffix = w % suffixSpace, middle = prefix % middleSpace;
        double shortest = shortestCounts[middle] / shortestTotal;
        double expected = shortest > 0 ? (shorterCounts[prefix] / shorterTotal) * (shorterCounts[suffix] / shorterTotal) / shortest : 0;
        vector[w] = expected > 0 ? (float)((words[w] / total - expected) / expected) : 0;
    }
    return vector;
}

// Feature hashing exactly as cvtree_compression_check.cs (commit f5cb655), so the 8,192 result applies.
static float[] Hash(float[] vector)
{
    var hashed = new float[dimensions];
    for (int w = 0; w < vector.Length; w++)
    {
        ulong h = (ulong)w * 0x9E3779B97F4A7C15UL;
        h ^= h >> 29;
        h *= 0xBF58476D1CE4E5B9UL;
        h ^= h >> 32;
        hashed[(int)(h % dimensions)] += (h & (1UL << 63)) == 0 ? vector[w] : -vector[w];
    }
    return hashed;
}

static float[] Normalize(float[] vector)
{
    double norm = Math.Sqrt(vector.Sum(x => (double)x * x));
    if (norm > 0) for (int i = 0; i < vector.Length; i++) vector[i] = (float)(vector[i] / norm);
    return vector;
}

// k nearest by cosine (vectors are unit length, so the largest dot product), itself excluded.
static int[][] Neighbours(float[][] vectors)
{
    int n = vectors.Length;
    var result = new int[n][];
    Parallel.For(0, n, new ParallelOptions { MaxDegreeOfParallelism = 4 }, i =>
    {
        var best = new int[k];
        var bestDot = new float[k];
        Array.Fill(bestDot, float.NegativeInfinity);
        var a = vectors[i];
        for (int j = 0; j < n; j++)
        {
            if (j == i) continue;
            float dot = Dot(a, vectors[j]);
            if (dot <= bestDot[k - 1]) continue;
            int slot = k - 1;
            while (slot > 0 && bestDot[slot - 1] < dot) { bestDot[slot] = bestDot[slot - 1]; best[slot] = best[slot - 1]; slot--; }
            bestDot[slot] = dot;
            best[slot] = j;
        }
        result[i] = best;
    });
    return result;
}

static float Dot(float[] a, float[] b)
{
    var sum = Vector<float>.Zero;
    int i = 0;
    for (; i <= a.Length - Vector<float>.Count; i += Vector<float>.Count)
        sum += new Vector<float>(a, i) * new Vector<float>(b, i);
    float total = Vector.Sum(sum);
    for (; i < a.Length; i++) total += a[i] * b[i];
    return total;
}
