#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Rotation check B: does rotating every circular bacterial chromosome to its replication origin change find-similar's ranking?
// The 150 bacterial queries of neighbour_eval.cs (same pool, seed and skip rule; precision@10 0.131 and MRR 0.260 by score in
// results/S34b_NeighbourRanking.csv), scored by SimilaritySearch's rule re-implemented here: records of length [L/2, 2L] except the
// query's own source, at SelectBandFactor's level, the better query strand, by score. Once as stored, once with every record
// that is kind 'bacterium', topology circular, >= 500 kbp and has a clear skew (origin_catalogue.cs's rule) rotated at profile
// level to its skew-minimum origin, the query included; and "best of", where the query and every such candidate are tried both as
// stored and rotated and the best score is kept (so an unrotated query still meets its candidates as stored). Every level is
// downsampled from the stored 500 bp profile in every run.
// Pass an output path (default S32_RotationRanking.csv); a summary goes to the console.
const int perKind = 150, top = 100;
const double ratio = SimilaritySearch.DefaultLengthRatio;
int[] coarseFactors = { 4, 16, 64 };
var output = args.Length > 0 ? args[0] : "S32_RotationRanking.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var meta = new List<(long Id, long SourceId, long Length, string Kind, string? Genus, string? Topology)>();
await using (var cmd = dataSource.CreateCommand(
    "select r.id, r.source_id, r.length_bp, r.kind, t.tax_genus, r.annotations->'contig'->>'topology' from sequence_record r " +
    "left join taxon t on t.taxid = r.ncbi_taxid order by r.id"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
        meta.Add((r.GetInt64(0), r.GetInt64(1), r.GetInt64(2), r.GetString(3), r.IsDBNull(4) ? null : r.GetString(4), r.IsDBNull(5) ? null : r.GetString(5)));
var byId = meta.ToDictionary(m => m.Id);

var cache = new ProfileCache(dataSource, RecordStore.FineWindowSize);
await cache.RefreshAsync();
var fine = cache.Profiles;

// Levels as stored (downsampled from the stored fine profile) and, for rotated records, from the rotated fine profile.
var stored = new Dictionary<long, Dictionary<int, ZCurveWindowedProfile>>();
var rotated = new Dictionary<long, Dictionary<int, ZCurveWindowedProfile>>();
var rotatedFine = new Dictionary<long, ZCurveWindowedProfile>();
foreach (var m in meta)
{
    if (!fine.TryGetValue(m.Id, out var profile)) continue;
    stored[m.Id] = coarseFactors.ToDictionary(f => f, f => profile.Downsample(f));
    if (m.Kind != "bacterium" || m.Topology != "circular" || m.Length < 500_000 || !ClearSkew(profile, m.Length)) continue;
    var turned = ToOrigin(profile);
    rotatedFine[m.Id] = turned;
    rotated[m.Id] = coarseFactors.ToDictionary(f => f, f => turned.Downsample(f));
}
Console.Error.WriteLine($"{fine.Count:N0} fine profiles, {rotated.Count:N0} chromosomes rotated");

ZCurveWindowedProfile Level(long id, int factor, bool rotate) =>
    rotate && rotated.TryGetValue(id, out var r) ? (factor == 1 ? rotatedFine[id] : r[factor]) : (factor == 1 ? fine[id] : stored[id][factor]);

// The versions of a record a mode compares: as stored, rotated (when it qualifies), or both for best-of.
IEnumerable<bool> Versions(long id, Mode mode) => mode switch
{
    Mode.Stored => [false],
    Mode.Rotated => [rotated.ContainsKey(id)],
    _ => rotated.ContainsKey(id) ? [false, true] : [false],
};

// One search, as SimilaritySearch.SearchRecordAsync without neighbours: ids by score, best first.
List<long> Search(long query, Mode mode)
{
    var q = byId[query];
    if (ZCurveProfileLadder.SelectBandFactor(new ZCurveProfileLadder(fine[query]), ratio) is not int factor) return [];
    var queries = Versions(query, mode).SelectMany(rotate =>
    {
        var queryFine = rotate ? rotatedFine[query] : fine[query];
        return new[] { new ZCurveProfileLadder(queryFine).Level(factor), new ZCurveProfileLadder(Flip(queryFine)).Level(factor) };
    }).ToArray();
    var band = SimilaritySearch.Band(q.Length, ratio);
    var candidates = meta.Where(m => band.Contains(m.Length) && m.SourceId != q.SourceId && fine.ContainsKey(m.Id)
        && Level(m.Id, factor, false).WindowCount > 0).Select(m => m.Id).ToArray();
    var scores = new double[candidates.Length];
    Parallel.For(0, candidates.Length, i =>
    {
        double best = double.NegativeInfinity;
        foreach (var rotate in Versions(candidates[i], mode))
        {
            var c = Level(candidates[i], factor, rotate);
            foreach (var queryLevel in queries) best = Math.Max(best, ZCurveWindowedProfile.Compare(queryLevel, c));
        }
        scores[i] = best;
    });
    return Enumerable.Range(0, candidates.Length).OrderByDescending(i => scores[i]).ThenBy(i => candidates[i]).Select(i => candidates[i]).ToList();
}

var rng = new Random(20261022);
var pool = meta.Where(m => m.Kind == "bacterium" && m.Genus is not null && m.Length >= 8_000).Select(m => m.Id).OrderBy(_ => rng.Next()).ToList();
var rows = new List<(long Query, bool QueryRotated, double P10, double P10Rotated, double Rr, double RrRotated, double P10Best, double RrBest)>();
using var report = File.CreateText(output);
report.WriteLine("QueryId, LengthBp, QueryRotated, Candidates, RotatedCandidates, PrecisionAt10Stored, PrecisionAt10Rotated, " +
    "ReciprocalRankStored, ReciprocalRankRotated, PrecisionAt10BestOf, ReciprocalRankBestOf");
foreach (var query in pool)
{
    if (rows.Count == perKind) break;
    if (!fine.ContainsKey(query)) continue;   // no stored 500 bp level
    string genus = byId[query].Genus!;
    var asStored = Search(query, Mode.Stored);
    if (!asStored.Any(id => byId[id].Genus == genus)) continue;   // no relative in the band
    var turned = Search(query, Mode.Rotated);
    var bestOf = Search(query, Mode.BestOf);

    (double P10, double Rr) Metrics(List<long> ids)
    {
        int first = ids.Take(top).ToList().FindIndex(id => byId[id].Genus == genus);
        return (ids.Take(10).Count(id => byId[id].Genus == genus) / 10.0, 1.0 / (first < 0 ? top + 1 : first + 1));
    }
    var (p10, rr) = Metrics(asStored);
    var (p10Rotated, rrRotated) = Metrics(turned);
    var (p10Best, rrBest) = Metrics(bestOf);
    rows.Add((query, rotated.ContainsKey(query), p10, p10Rotated, rr, rrRotated, p10Best, rrBest));
    report.WriteLine($"{query}, {byId[query].Length}, {rotated.ContainsKey(query)}, {asStored.Count}, {asStored.Count(rotated.ContainsKey)}, " +
        $"{p10:F1}, {p10Rotated:F1}, {rr:F4}, {rrRotated:F4}, {p10Best:F1}, {rrBest:F4}");
}

Console.WriteLine("Queries, Subset, Variant, MeanPrecisionAt10Stored, MeanPrecisionAt10Variant, MrrStored, MrrVariant, Better, Worse, SignTestP");
foreach (var (subset, list) in new[] { ("all", rows), ("query rotated", rows.Where(r => r.QueryRotated).ToList()),
    ("query not rotated", rows.Where(r => !r.QueryRotated).ToList()) })
    foreach (var (variant, p10Of, rrOf) in new (string, Func<(long, bool, double, double, double, double, double, double), double>,
        Func<(long, bool, double, double, double, double, double, double), double>)[]
        { ("rotated", r => r.Item4, r => r.Item6), ("best of", r => r.Item7, r => r.Item8) })
    {
        int better = list.Count(r => p10Of(r) > r.P10), worse = list.Count(r => p10Of(r) < r.P10);
        Console.WriteLine($"{list.Count}, {subset}, {variant}, {list.Average(r => r.P10):F3}, {list.Average(p10Of):F3}, {list.Average(r => r.Rr):F3}, " +
            $"{list.Average(rrOf):F3}, {better}, {worse}, {SignTest(better, worse):F4}");
    }

static bool ClearSkew(ZCurveWindowedProfile profile, long lengthBp)
{
    var windows = profile.Windows;
    double sum = 0, min = 0, max = 0;
    int origin = 0, terminus = 0;
    for (int i = 0; i < windows.Count; i++)
    {
        sum += windows[i].X - windows[i].Y;
        if (sum < min) (min, origin) = (sum, i + 1);
        if (sum > max) (max, terminus) = (sum, i + 1);
    }
    long n = windows.Count, d = Math.Abs(origin - terminus) % n;
    return (double)Math.Min(d, n - d) / n >= 0.4 && (max - min) * 250 / (lengthBp / 1000.0) >= 3;
}

// As StrandAndRotationBenchmark.ToOrigin.
static ZCurveWindowedProfile ToOrigin(ZCurveWindowedProfile profile)
{
    var windows = profile.Windows;
    double sum = 0, min = 0;
    int origin = 0;
    for (int i = 0; i < windows.Count; i++)
    {
        sum += windows[i].X - windows[i].Y;
        if (sum < min) (min, origin) = (sum, (i + 1) % windows.Count);
    }
    return ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize, windows.Skip(origin).Concat(windows.Take(origin)).ToList());
}

static ZCurveWindowedProfile Flip(ZCurveWindowedProfile p) =>
    ZCurveWindowedProfile.FromPrecomputedWindows(p.WindowSize, p.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

// Two-sided exact sign test over the queries that changed.
static double SignTest(int better, int worse)
{
    int n = better + worse, k = Math.Min(better, worse);
    if (n == 0) return 1;
    double tail = 0;
    for (int i = 0; i <= k; i++) tail += Math.Exp(LogChoose(n, i) - n * Math.Log(2));
    return Math.Min(1, 2 * tail);
}

static double LogChoose(int n, int k) => LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k);

static double LogFactorial(int n)
{
    double sum = 0;
    for (int i = 2; i <= n; i++) sum += Math.Log(i);
    return sum;
}

enum Mode { Stored, Rotated, BestOf }
