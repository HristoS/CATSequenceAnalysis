#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Buffers;
using System.Collections.Concurrent;
using System.Data;
using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Search breakdown (2026-11-09, plan "faster protein search and localization"): where a protein search and a streamed localization
// spend their time on the local catalogue, stage by stage, before any change. Run with `dotnet run -c Release`.
// Protein: seeded random proteins at 60, 120, 250, 400 and 700 aa (cost depends only on length, which picks the buckets; no network).
// Stages as ProteinSearch runs them, timed one at a time over the same rows: fixed overhead (visible ids, protein null), the database
// read as now (buffered reader, a byte[] per row) and as a pure transfer (SequentialAccess into one pooled buffer), SHA-256 of every
// row, FromBinary of every row, Score of the candidates; then ProteinSearch.SearchAsync end to end.
// Localization: fragments of 4, 8, 20, 50 and 128 kbp from the bundled O157:H7 Sakai (stored in the catalogue), both strands, as
// Localization.LocalizeAsync runs without a cache: ids, the coarse level streamed (read only, then with hash and decode), the scan
// over every phase and strand, the refine's read and compute; then LocalizeAsync streamed and cached end to end.
// Each case: one cold run dropped, median of 3 warm. Baseline hits go to the second argument for the after-change equality check.
// Args: output CSV (default SearchBreakdown_Before.csv), hits file (default search_breakdown_hits.txt).
const int seed = 20261109, warmRuns = 3;
int[] proteinLengths = [60, 120, 250, 400, 700];
int[] fragmentLengths = [4_000, 8_000, 20_000, 50_000, 128_000];
var output = args.Length > 0 ? args[0] : "SearchBreakdown_Before.csv";
var hitsPath = args.Length > 1 ? args[1] : "search_breakdown_hits.txt";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var rows = new List<string> { "Search,Query,Stage,MedianMs,Bytes,Rows,Items,MBPerS,Gen2Collections" };
var hits = new List<string>();

double Median(List<double> values) => values.Order().ElementAt(values.Count / 2);

// One cold run, then the median of warmRuns; the stage returns (bytes, rows, items) of its last run.
async Task Stage(string search, string query, string stage, Func<Task<(long Bytes, long Rows, long Items)>> run)
{
    await run();
    var times = new List<double>();
    (long Bytes, long Rows, long Items) last = default;
    int gen2 = GC.CollectionCount(2);
    for (int i = 0; i < warmRuns; i++)
    {
        var watch = Stopwatch.StartNew();
        last = await run();
        times.Add(watch.Elapsed.TotalMilliseconds);
    }
    double ms = Median(times);
    string mbps = last.Bytes > 0 ? (last.Bytes / 1e6 / (ms / 1000)).ToString("F0", CultureInfo.InvariantCulture) : "";
    var line = string.Join(",", search, query, stage.Replace(", ", "; "), ms.ToString("F0", CultureInfo.InvariantCulture), last.Bytes, last.Rows, last.Items, mbps,
        (GC.CollectionCount(2) - gen2) / warmRuns);
    rows.Add(line);
    Console.WriteLine(line);
}

// ---------------- Protein search ----------------
var random = new Random(seed);
const string aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
foreach (var length in proteinLengths)
{
    var protein = new string(Enumerable.Range(0, length).Select(_ => aminoAcids[random.Next(aminoAcids.Length)]).ToArray());
    var query = ProteinQuery.FromProtein(protein);
    var name = $"{length} aa";
    var buckets = GeneProfiles.BucketsFor(query.LengthAa, ProteinSearch.LengthBand).Select(b => (short)b).ToArray();
    var candidate = ProteinSearch.CandidateRule(query);
    const string sql = "select sequence_id, payload, sha256 from gene_profile where bucket = any($1) and window_count > 0";

    await Stage("protein", name, "fixed: visible ids + protein null", async () =>
    {
        var visible = await Catalogue.IdsAsync(dataSource, new CatalogueFilter());
        await ProteinScoreNull.LoadAsync(dataSource);
        return (0, 0, visible.Count);
    });

    List<(long Id, byte[] Payload, byte[] Sha256)> fetched = [];
    await Stage("protein", name, "read: buffered, byte[] per row (as now)", async () =>
    {
        var list = new List<(long, byte[], byte[])>();
        long bytes = 0;
        await using var command = dataSource.CreateCommand(sql);
        command.Parameters.AddWithValue(buckets);
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync())
        {
            var payload = reader.GetFieldValue<byte[]>(1);
            bytes += payload.Length;
            list.Add((reader.GetInt64(0), payload, reader.GetFieldValue<byte[]>(2)));
        }
        fetched = list;
        return (bytes, list.Count, 0);
    });

    await Stage("protein", name, "read: sequential, pooled buffer (transfer only)", async () =>
    {
        long bytes = 0, count = 0;
        var buffer = ArrayPool<byte>.Shared.Rent(4 << 20);
        await using var command = dataSource.CreateCommand(sql);
        command.Parameters.AddWithValue(buckets);
        await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess);
        while (await reader.ReadAsync())
        {
            reader.GetInt64(0);
            int length = (int)reader.GetBytes(1, 0, null, 0, 0);
            reader.GetBytes(1, 0, buffer, 0, length);
            bytes += length;
            count++;
        }
        ArrayPool<byte>.Shared.Return(buffer);
        return (bytes, count, 0);
    });

    long totalBytes = fetched.Sum(r => (long)r.Payload.Length);
    await Stage("protein", name, "sha256: every row, all cores", () =>
    {
        Parallel.ForEach(fetched, r => { if (!SHA256.HashData(r.Payload).AsSpan().SequenceEqual(r.Sha256)) throw new InvalidDataException(); });
        return Task.FromResult((totalBytes, (long)fetched.Count, 0L));
    });

    ZCurveGeneProfile[] decoded = [];
    await Stage("protein", name, "decode: FromBinary every row, all cores", () =>
    {
        var profiles = new ZCurveGeneProfile[fetched.Count];
        Parallel.For(0, fetched.Count, i => profiles[i] = ZCurveGeneProfile.FromBinary(fetched[i].Payload));
        decoded = profiles;
        return Task.FromResult((totalBytes, (long)fetched.Count, (long)profiles.Sum(p => p.Genes.Count)));
    });

    await Stage("protein", name, "score: candidates, all cores", () =>
    {
        long scored = 0;
        Parallel.For(0, decoded.Length, () => 0L, (i, _, local) =>
        {
            var profile = decoded[i];
            for (int g = 0; g < profile.Genes.Count; g++)
                if (candidate(profile.Genes[g])) { profile.Score(query.Axes, g); local++; }
            return local;
        }, local => Interlocked.Add(ref scored, local));
        return Task.FromResult((0L, (long)decoded.Length, scored));
    });
    await Stage("protein", name, "score in place: ScoreRow every row, all cores (since 2026-11-09)", () =>
    {
        long scored = 0;
        Parallel.For(0, fetched.Count, () => (Genes: Array.Empty<ZCurveGeneProfile.Gene>(), Count: 0L), (i, _, local) =>
        {
            var payload = fetched[i].Payload;
            int geneCount = ZCurveGeneProfile.GeneCountOf(payload);
            var genes = local.Genes.Length < geneCount ? new ZCurveGeneProfile.Gene[geneCount] : local.Genes;
            long count = local.Count;
            ZCurveGeneProfile.ScoreRow(payload, query.Axes, genes, candidate, (_, _) => count++);
            return (genes, count);
        }, local => Interlocked.Add(ref scored, local.Count));
        return Task.FromResult((totalBytes, (long)fetched.Count, scored));
    });
    fetched = [];
    decoded = [];

    ProteinSearchResult? result = null;
    await Stage("protein", name, "end to end: ProteinSearch.SearchAsync", async () =>
    {
        result = await ProteinSearch.SearchAsync(dataSource, query);
        return (totalBytes, 0, result.Candidates);
    });
    foreach (var (hit, rank) in result!.Hits.Select((h, i) => (h, i + 1)))
        hits.Add(string.Join("\t", "protein", name, rank, hit.SequenceId, hit.GeneSequence, hit.GeneStart, hit.GeneEnd, hit.Reverse,
            hit.Score.ToString("R", CultureInfo.InvariantCulture)));
}

// ---------------- Localization (streamed, as the worker runs it today) ----------------
var o157 = FastaReader.ReadPrimaryRecord(Path.Combine("data/genomes", "ecoli_o157h7_sakai_GCF_000008865.2.fna.gz")).Sequence;
var streamed = new Localization(dataSource, useCache: false);
var cached = new Localization(dataSource, useCache: true);
foreach (var length in fragmentLengths)
{
    int start = random.Next(o157.Length - length);
    var fragment = QueryProfile.FromSequence(o157.Substring(start, length));
    var name = $"{length / 1000} kbp at {start}";
    int factor = ZCurveProfileLadder.SelectQueryFactor(fragment.Forward, ZCurveProfileLadder.DefaultMinFragmentWindows)!.Value;
    int windowSize = RecordStore.FineWindowSize * factor;
    var lengths = new LengthRange(fragment.LengthBp, null);
    HashSet<long> candidates = [];

    await Stage("localize", name, $"fixed: candidate ids (level {windowSize} bp)", async () =>
    {
        candidates = await Catalogue.IdsAsync(dataSource, new CatalogueFilter { MinLengthBp = fragment.LengthBp });
        return (0, 0, candidates.Count);
    });

    await Stage("localize", name, "read: level rows, sequential (transfer only)", async () =>
    {
        long bytes = 0, count = 0;
        await using var command = dataSource.CreateCommand(
            "select p.payload from profile_level p join sequence_record r on r.id = p.sequence_id " +
            "where p.method = $1 and p.method_version = $2 and p.window_size = $3 and p.id > 0 and r.length_bp >= $4 and r.length_bp < $5 order by p.id");
        command.Parameters.AddWithValue(RecordStore.Method);
        command.Parameters.AddWithValue(RecordStore.MethodVersion);
        command.Parameters.AddWithValue(windowSize);
        command.Parameters.AddWithValue(lengths.Min);
        command.Parameters.AddWithValue(lengths.MaxOrUnbounded);
        await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess);
        while (await reader.ReadAsync())
        {
            bytes += reader.GetBytes(0, 0, null, 0, 0);
            count++;
        }
        return (bytes, count, 0);
    });

    var profiles = new Dictionary<long, ZCurveWindowedProfile>();
    await Stage("localize", name, "read: StreamAfterAsync (hash + decode, as now)", async () =>
    {
        profiles = new Dictionary<long, ZCurveWindowedProfile>();
        long bytes = 0;
        await foreach (var (_, id, level) in ProfileLevels.StreamAfterAsync(dataSource, RecordStore.Method, RecordStore.MethodVersion, windowSize, 0, lengths))
        {
            profiles[id] = level;
            bytes += 16 + 6L * level.WindowCount;
        }
        return (bytes, profiles.Count, 0);
    });

    var comparable = candidates.Where(profiles.ContainsKey).ToArray();
    ZCurveWindowedProfile[] strands = fragment.Reverse is null ? [fragment.Forward.Fine] : [fragment.Forward.Fine, fragment.Reverse.Fine];
    var phased = strands.Select(s => Enumerable.Range(0, factor).Select(p => s.Downsample(factor, p).ToHalfPrecision()).ToArray()).ToArray();
    var placements = new (double Score, int Strand, int Phase, int Offset)[comparable.Length];
    await Stage("localize", name, $"scan: {phased.Length} strands x {factor} phases, all cores", () =>
    {
        Parallel.For(0, comparable.Length, i =>
        {
            var candidate = profiles[comparable[i]];
            var best = (Score: double.NegativeInfinity, Strand: 0, Phase: 0, Offset: -1);
            for (int strand = 0; strand < phased.Length; strand++)
                for (int phase = 0; phase < phased[strand].Length; phase++)
                {
                    var query = phased[strand][phase];
                    if (query.WindowCount == 0 || query.WindowCount > candidate.WindowCount) continue;
                    var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, candidate);
                    if (score > best.Score) best = (score, strand, phase, offset);
                }
            placements[i] = best;
        });
        return Task.FromResult((0L, (long)comparable.Length, (long)comparable.Length * phased.Length * factor));
    });

    var best = Enumerable.Range(0, comparable.Length).Where(i => placements[i].Offset >= 0)
        .OrderByDescending(i => placements[i].Score).ThenBy(i => comparable[i]).Take(Localization.RefineDepth).ToArray();
    Dictionary<long, ZCurveWindowedProfile> fine = [];
    if (factor > 1)
    {
        await Stage("localize", name, "refine read: 100 fine levels", async () =>
        {
            fine = await ProfileLevels.ReadManyAsync(dataSource, best.Select(i => comparable[i]).ToArray(), RecordStore.Method,
                RecordStore.MethodVersion, RecordStore.FineWindowSize);
            return (fine.Values.Sum(f => 16 + 6L * f.WindowCount), fine.Count, 0);
        });
        await Stage("localize", name, "refine compute", () =>
        {
            foreach (var i in best)
                if (fine.TryGetValue(comparable[i], out var reference))
                    ZCurveWindowedProfile.RefineAlignment(strands[placements[i].Strand], reference,
                        placements[i].Offset * factor - placements[i].Phase, radius: factor);
            return Task.FromResult((0L, (long)best.Length, 0L));
        });
    }
    profiles.Clear();

    LocalizationResult? result = null;
    await Stage("localize", name, "end to end: LocalizeAsync streamed", async () =>
    {
        result = await streamed.LocalizeAsync(fragment);
        return (0, 0, result.Candidates);
    });
    await Stage("localize", name, "end to end: LocalizeAsync cached", async () =>
    {
        var again = await cached.LocalizeAsync(fragment);
        return (0, 0, again.Candidates);
    });
    foreach (var (hit, rank) in result!.Hits.Select((h, i) => (h, i + 1)))
        hits.Add(string.Join("\t", "localize", name, rank, hit.SequenceId, hit.OffsetBp, hit.Strand,
            hit.Score.ToString("R", CultureInfo.InvariantCulture), hit.CoarseScore.ToString("R", CultureInfo.InvariantCulture)));
}

File.WriteAllLines(output, rows);
File.WriteAllLines(hitsPath, hits);
Console.WriteLine($"Wrote {output} and {hits.Count} baseline hits to {hitsPath}.");
