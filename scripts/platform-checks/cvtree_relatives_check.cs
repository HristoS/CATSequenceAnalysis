#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// CVTree relatives card, check (2026-11-22): would the relatives service (SignatureIndex, 10 nearest) gain from amino-acid CVTree over
// log rho*? Fixed before the run:
// - Genomes: the public genomes SignatureIndex holds (largest chromosome of a bacterium or archaeon, every virus record) that have stored
//   amino-acid 4-mers; both arms rank the same candidates. Eukaryotes left out (131 with 4-mers).
// - Arms: hashed CVTree exactly as p9_e1b_cvtree.cs (4-mers, 3- and 2-mers derived, 8,192 dimensions, cosine), the same rounded to 16-bit
//   halves, and log rho* (Euclidean, SignatureIndex.Signature). Added the same day, after the go-ahead: the platform's vectors
//   (CvTreeSignature, 16-bit fixed point, what the index holds), checked against float32.
// - Queries: 500 per group (seed 20261122) with a known genus; the 10 nearest, itself excluded. Keys: GTDB genus and family for bacteria
//   and archaea (primary), NCBI's current genus and family for both groups (primary for viruses). A candidate without the key counts
//   as a different taxon in every arm. Context: the same with the query's own species removed from the candidates.
// - Rule (card): per group, CVTree's mean same-genus share among the 10 nearest beats log rho*'s by >= 5 points; also build <= 5 min on
//   the server (4 local threads x 2, extrapolated) and 16-bit memory <= 1 GB. Bootstrap 95% CI of the paired difference as context.
// Arg: output path (default CvTreeRelatives_Check.csv).
const int k = 10, queries = 500, dimensions = 8192, batch = 1000, threads = 4, bootstraps = 1000, seed = 20261122;
var output = args.Length > 0 ? args[0] : "CvTreeRelatives_Check.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var genomes = new List<Genome>();
await using (var command = dataSource.CreateCommand(
    "select r.id, r.kind, r.dinucleotide_counts, g.tax_family, g.tax_genus, g.tax_species, t.tax_family, t.tax_genus, t.tax_species " +
    "from (select r.*, row_number() over (partition by r.source_id, r.kind in ('bacterium', 'archaeon') order by r.length_bp desc, r.id) as n" +
    "      from sequence_record r join sequence_source s on s.id = r.source_id" +
    "      where s.owner_id is null and r.dinucleotide_counts is not null and r.kind in ('bacterium', 'archaeon', 'virus')) r " +
    "join protein_kmer p on p.sequence_id = r.id and p.k = 4 and p.status = 'done' " +
    "left join taxon t on t.taxid = r.ncbi_taxid " +
    "left join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "  and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind = 'virus' or r.n = 1 order by r.id"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
    {
        string? Text(int i) => reader.IsDBNull(i) || reader.GetString(i) == "" ? null : reader.GetString(i);
        var kind = reader.GetString(1);
        genomes.Add(new Genome(reader.GetInt64(0), kind == "virus" ? SignatureIndex.Viral : SignatureIndex.Cellular,
            SignatureIndex.Signature(reader.GetFieldValue<long[]>(2)),
            new Taxa(Text(3), Text(4), Text(5)), new Taxa(Text(6), Text(7), Text(8))));
    }

var rows = new List<string> { "group,key,arm,ownSpecies,queries,genusShare,familyShare,genusMinusRho,genusMinusRhoLow,genusMinusRhoHigh" };
var notes = new List<string>();
foreach (var group in new[] { SignatureIndex.Cellular, SignatureIndex.Viral })
{
    var members = genomes.Where(g => g.Group == group).ToArray();
    var watch = Stopwatch.StartNew();
    var vectors = await CvTreeVectorsAsync(dataSource, members.Select(g => g.Id).ToArray());
    double build = watch.Elapsed.TotalSeconds;
    double memory16 = (double)members.Length * dimensions * 2 / 1e9;
    var keys = group == SignatureIndex.Cellular
        ? new (string Name, Func<Genome, Taxa> Of)[] { ("gtdb", g => g.Gtdb), ("ncbi", g => g.Ncbi) }
        : new (string Name, Func<Genome, Taxa> Of)[] { ("ncbi", g => g.Ncbi) };
    var primary = keys[0];
    var random = new Random(seed);
    var picked = Enumerable.Range(0, members.Length).Where(i => primary.Of(members[i]).Genus is not null)
        .OrderBy(_ => random.Next()).Take(queries).ToArray();
    Console.Error.WriteLine($"{group}: {members.Length:N0} genomes, {members.Count(g => primary.Of(g).Genus is not null):N0} with a {primary.Name} genus; " +
        $"vectors built in {build:F1} s on {threads} threads; {picked.Length} queries");

    // One query, single thread, over the whole group: what a relatives request would cost.
    var times = new List<double>();
    foreach (var q in picked.Take(20))
    {
        watch.Restart();
        Nearest(q, members.Length, (a, b) => -Dot(vectors[a], vectors[b]), _ => false);
        times.Add(watch.Elapsed.TotalMilliseconds);
    }
    double queryMs = times.Order().ElementAt(times.Count / 2);

    foreach (bool ownSpecies in new[] { true, false })
    {
        var results = new Dictionary<string, int[][]>();
        foreach (var key in keys)
        {
            Func<int, int, bool> skip = (q, j) => !ownSpecies && key.Of(members[q]).Species is { } s && key.Of(members[j]).Species == s;
            results[$"cvtree|{key.Name}"] = NearestAll(picked, members.Length, (a, b) => -Dot(vectors[a], vectors[b]), skip);
            results[$"rho|{key.Name}"] = NearestAll(picked, members.Length, (a, b) => Euclidean(members[a].Rho, members[b].Rho), skip);
        }
        foreach (var key in keys)
            foreach (var arm in new[] { "cvtree", "rho" })
                rows.Add(Row(group, key, arm, ownSpecies, results[$"{arm}|{key.Name}"], results[$"rho|{key.Name}"]));
    }

    // 16-bit: round every vector through a half, rank again (own species kept), and compare with float32. Then the platform's own
    // vectors (CvTreeSignature: 16-bit fixed point), built and scored by the code the index runs.
    var full = NearestAll(picked, members.Length, (a, b) => -Dot(vectors[a], vectors[b]), (_, _) => false);
    foreach (var v in vectors)
        for (int i = 0; i < v.Length; i++) v[i] = (float)(Half)v[i];
    var half = NearestAll(picked, members.Length, (a, b) => -Dot(vectors[a], vectors[b]), (_, _) => false);
    vectors = null!;
    watch.Restart();
    var stored = await CvTreeSignature.ReadAsync(dataSource, members.Select(g => g.Id).ToArray(), threads);
    double platformBuild = watch.Elapsed.TotalSeconds;
    // A stored row with no word gives no vector (the index leaves such a genome out); here it scores 0.5 to all, as its zero float vector.
    int withoutVector = members.Count(g => !stored.ContainsKey(g.Id));
    var fixedPoint = members.Select(g => stored.GetValueOrDefault(g.Id) ?? new short[CvTreeSignature.Dimensions]).ToArray();
    stored = null!;
    var platform = NearestAll(picked, members.Length, (a, b) => CvTreeSignature.Distance(fixedPoint[a], fixedPoint[b]), (_, _) => false);
    times.Clear();
    foreach (var q in picked.Take(20))
    {
        watch.Restart();
        double nearest = double.MaxValue;
        for (int j = 0; j < fixedPoint.Length; j++) nearest = Math.Min(nearest, CvTreeSignature.Distance(fixedPoint[q], fixedPoint[j]));
        times.Add(watch.Elapsed.TotalMilliseconds);
    }
    double platformMs = times.Order().ElementAt(times.Count / 2);
    var rhoAgain = NearestAll(picked, members.Length, (a, b) => Euclidean(members[a].Rho, members[b].Rho), (_, _) => false);
    foreach (var key in keys)
    {
        rows.Add(Row(group, key, "cvtree16", true, half, rhoAgain));
        rows.Add(Row(group, key, "platform", true, platform, rhoAgain));
    }
    string Same(int[][] other) => $"identical {Enumerable.Range(0, picked.Length).Count(i => full[i].SequenceEqual(other[i]))}/{picked.Length}, " +
        $"mean overlap {Enumerable.Range(0, picked.Length).Average(i => full[i].Intersect(other[i]).Count() / (double)k):P1}";
    notes.Add($"{group}: {members.Length} genomes; build {build:F1} s on {threads} threads (server, 2 vCPU: about {build * 2:F0} s, extrapolated); " +
        $"16-bit index {memory16:F2} GB; one query {queryMs:F1} ms single-thread; vs float32 10-nearest: half {Same(half)}; " +
        $"platform fixed point {Same(platform)}, its build {platformBuild:F1} s, its scan {platformMs:F1} ms single-thread, " +
        $"{withoutVector} genomes whose stored 4-mers hold no word");

    string Row(string g, (string Name, Func<Genome, Taxa> Of) key, string arm, bool own, int[][] lists, int[][] rho)
    {
        double Share(int[][] l, int i, Func<Taxa, string?> rank) =>
            rank(key.Of(members[picked[i]])) is { } name ? l[i].Count(j => rank(key.Of(members[j])) == name) / (double)k : double.NaN;
        var genus = Enumerable.Range(0, picked.Length).Select(i => Share(lists, i, t => t.Genus)).ToArray();
        var family = Enumerable.Range(0, picked.Length).Select(i => Share(lists, i, t => t.Family)).Where(x => !double.IsNaN(x)).ToArray();
        var rhoGenus = Enumerable.Range(0, picked.Length).Select(i => Share(rho, i, t => t.Genus)).ToArray();
        var paired = Enumerable.Range(0, picked.Length).Where(i => !double.IsNaN(genus[i])).Select(i => genus[i] - rhoGenus[i]).ToArray();
        var boot = new Random(seed);
        var means = Enumerable.Range(0, bootstraps).Select(_ => Enumerable.Range(0, paired.Length).Average(_ => paired[boot.Next(paired.Length)])).Order().ToArray();
        return string.Join(",", g, key.Name, arm, own ? "kept" : "removed", paired.Length, F(genus.Where(x => !double.IsNaN(x)).Average()),
            F(family.Average()), F(paired.Average()), F(means[(int)(0.025 * bootstraps)]), F(means[(int)(0.975 * bootstraps)]));
    }
}

await File.WriteAllLinesAsync(output, rows);
foreach (var row in rows) Console.WriteLine(row);
foreach (var note in notes) Console.WriteLine(note);

static string F(double x) => x.ToString("F4", CultureInfo.InvariantCulture);

static int[][] NearestAll(int[] picked, int n, Func<int, int, double> distance, Func<int, int, bool> skip)
{
    var result = new int[picked.Length][];
    Parallel.For(0, picked.Length, new ParallelOptions { MaxDegreeOfParallelism = threads },
        i => result[i] = Nearest(picked[i], n, distance, j => skip(picked[i], j)));
    return result;
}

// The k nearest of q (smallest distance; ties by index), q itself and skipped candidates excluded.
static int[] Nearest(int q, int n, Func<int, int, double> distance, Func<int, bool> skip)
{
    var best = new int[k];
    var bestDistance = new double[k];
    Array.Fill(bestDistance, double.PositiveInfinity);
    for (int j = 0; j < n; j++)
    {
        if (j == q || skip(j)) continue;
        double d = distance(q, j);
        if (d >= bestDistance[k - 1]) continue;
        int slot = k - 1;
        while (slot > 0 && bestDistance[slot - 1] > d) { bestDistance[slot] = bestDistance[slot - 1]; best[slot] = best[slot - 1]; slot--; }
        bestDistance[slot] = d;
        best[slot] = j;
    }
    return best;
}

static double Euclidean(double[] a, double[] b)
{
    double sum = 0;
    for (int i = 0; i < a.Length; i++) sum += (a[i] - b[i]) * (a[i] - b[i]);
    return Math.Sqrt(sum);
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

// Hashed CVTree vectors (unit length) in the order of ids, read in batches, built on `threads` threads; as p9_e1b_cvtree.cs.
static async Task<float[][]> CvTreeVectorsAsync(NpgsqlDataSource dataSource, long[] ids)
{
    var vectors = new Dictionary<long, float[]>();
    foreach (var chunk in ids.Chunk(batch))
    {
        var payloads = new List<(long Id, byte[] Payload)>();
        await using (var command = dataSource.CreateCommand("select sequence_id, payload from protein_kmer where k = 4 and status = 'done' and sequence_id = any($1)"))
        {
            command.Parameters.AddWithValue(chunk);
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) payloads.Add((reader.GetInt64(0), reader.GetFieldValue<byte[]>(1)));
        }
        foreach (var (id, vector) in payloads.AsParallel().WithDegreeOfParallelism(threads)
            .Select(p => (p.Id, Normalize(Hash(Composition(ProteinKmers.Decode(p.Payload, 4)))))).ToList())
            vectors[id] = vector;
    }
    return ids.Select(id => vectors[id]).ToArray();
}

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

record Taxa(string? Family, string? Genus, string? Species);

record Genome(long Id, string Group, double[] Rho, Taxa Gtdb, Taxa Ncbi);
