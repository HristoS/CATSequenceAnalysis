#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Would find-similar gain from wrap-around scoring? Phage genomes are often circularly permuted: the same genome deposited from a
// different start (Pantoea phage ROV85 scores 0.870 against Escherichia phage Lidtsur as deposited, 0.924 from a start 24.5 kbp on).
// find-similar slides the shorter profile along the longer and never wraps. Here each query is scored as find-similar scores it
// (level chosen by length at ratio 0.5, candidates of length [0.5 L, 2 L], both query strands), once plain and once against each
// candidate's windows followed by themselves (every start). Viruses with an NCBI genus, 10–400 kbp; 400 queries with a same-genus
// record in the band. Rule fixed before the run: stop (search stays as it is) if wrap-around gains < 0.005 mean per-query AUC
// (same genus vs other family) and changes < 2% of top-10 hits. Run from the repo root with SEQUENCE_PLATFORM_DB set.
const int queries = 400;
var output = args.Length > 0 ? args[0] : "results/CircularPermutation_Check.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var meta = new Dictionary<long, (long Length, string Genus, string? Family)>();
await using (var command = dataSource.CreateCommand(
    "select r.id, r.length_bp, t.tax_genus, t.tax_family from sequence_record r join taxon t on t.taxid = r.ncbi_taxid " +
    "where r.kind = 'virus' and t.tax_genus is not null and r.length_bp between 10000 and 400000"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
        meta[reader.GetInt64(0)] = (reader.GetInt64(1), reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetString(3));

var fine = new Dictionary<long, ZCurveWindowedProfile>();
await foreach (var (id, level) in ProfileLevels.StreamAsync(dataSource, RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize))
    if (meta.ContainsKey(id)) fine[id] = level;
var pool = fine.Keys.Order().ToArray();
Console.WriteLine($"{pool.Length:N0} viruses with a genus and a 500 bp level");

var factors = ZCurveProfileLadder.DefaultFactors;
var levels = factors.Select(_ => new System.Collections.Concurrent.ConcurrentDictionary<long, (ZCurveWindowedProfile Plain, ZCurveWindowedProfile Wrapped)>()).ToArray();
(ZCurveWindowedProfile Plain, ZCurveWindowedProfile Wrapped) At(long id, int level) => levels[level].GetOrAdd(id, _ =>
{
    var plain = factors[level] == 1 ? fine[id] : fine[id].Downsample(factors[level]);
    var windows = plain.Windows.Select(w => (w.X, w.Y, w.Z)).ToList();
    return (plain, ZCurveWindowedProfile.FromPrecomputedWindows(plain.WindowSize, windows.Concat(windows.Take(Math.Max(0, windows.Count - 1))).ToList()));
});
static ZCurveWindowedProfile Flip(ZCurveWindowedProfile p) =>
    ZCurveWindowedProfile.FromPrecomputedWindows(p.WindowSize, p.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

var rng = new Random(20261112);
var order = pool.OrderBy(_ => rng.Next()).ToList();
var perQuery = new List<(long Query, int Candidates, double AucPlain, double AucWrapped, int SameGenusTopPlain, int SameGenusTopWrapped,
    int TopChanged, int GenusGains, int GenusPairs, int OtherGains, int OtherPairs, int ReciprocalRankPlain1000, int ReciprocalRankWrapped1000)>();
foreach (var q in order)
{
    if (perQuery.Count == queries) break;
    var (length, genus, family) = meta[q];
    var band = SimilaritySearch.Band(length, 0.5);
    var candidates = pool.Where(c => c != q && meta[c].Length >= band.Min && meta[c].Length < band.Max).ToArray();
    if (!candidates.Any(c => meta[c].Genus == genus)) continue;
    int level = SimilaritySearch.BandLevel(length, 0.5);
    var query = At(q, level).Plain;
    var forms = new[] { query, Flip(query) };
    var scores = new (long Id, double Plain, double Wrapped)[candidates.Length];
    Parallel.For(0, candidates.Length, i =>
    {
        var (plain, wrapped) = At(candidates[i], level);
        scores[i] = (candidates[i], forms.Max(f => ZCurveWindowedProfile.Compare(f, plain)), forms.Max(f => ZCurveWindowedProfile.Compare(f, wrapped)));
    });

    bool Same(long id) => meta[id].Genus == genus;
    bool Other(long id) => family is not null && meta[id].Family is not null && meta[id].Family != family;
    double Auc(Func<(long Id, double Plain, double Wrapped), double> score)
    {
        var positives = scores.Where(s => Same(s.Id)).Select(score).ToArray();
        var negatives = scores.Where(s => Other(s.Id)).Select(score).Order().ToArray();
        if (positives.Length == 0 || negatives.Length == 0) return double.NaN;
        double wins = 0;
        foreach (var p in positives)
        {
            int below = LowerBound(negatives, p), upTo = UpperBound(negatives, p);
            wins += below + 0.5 * (upTo - below);
        }
        return wins / ((double)positives.Length * negatives.Length);
    }
    var topPlain = scores.OrderByDescending(s => s.Plain).ThenBy(s => s.Id).Take(10).Select(s => s.Id).ToList();
    var topWrapped = scores.OrderByDescending(s => s.Wrapped).ThenBy(s => s.Id).Take(10).Select(s => s.Id).ToList();
    int Rank(Func<(long Id, double Plain, double Wrapped), double> score) =>
        scores.OrderByDescending(score).ThenBy(s => s.Id).Select((s, i) => (s, i)).First(x => Same(x.s.Id)).i + 1;
    var genusPairs = scores.Where(s => Same(s.Id)).ToList();
    var otherPairs = scores.Where(s => Other(s.Id)).ToList();
    perQuery.Add((q, candidates.Length, Auc(s => s.Plain), Auc(s => s.Wrapped), topPlain.Count(Same), topWrapped.Count(Same),
        topWrapped.Except(topPlain).Count(), genusPairs.Count(s => s.Wrapped - s.Plain > 0.03), genusPairs.Count,
        otherPairs.Count(s => s.Wrapped - s.Plain > 0.03), otherPairs.Count, (int)(1000.0 / Rank(s => s.Plain)), (int)(1000.0 / Rank(s => s.Wrapped))));
    if (perQuery.Count % 50 == 0) Console.WriteLine($"  {perQuery.Count} queries");
}

var withAuc = perQuery.Where(p => !double.IsNaN(p.AucPlain)).ToList();
double aucPlain = withAuc.Average(p => p.AucPlain), aucWrapped = withAuc.Average(p => p.AucWrapped);
double topChanged = perQuery.Sum(p => p.TopChanged) / (10.0 * perQuery.Count);
var lines = new List<string>
{
    "measure,plain,wrapped,note",
    $"queries,{perQuery.Count},{perQuery.Count},viruses with a same-genus record in the length band; {pool.Length} in the pool",
    $"median_candidates,{Median(perQuery.Select(p => (double)p.Candidates))},,",
    $"mean_auc_genus_vs_other_family,{aucPlain:F4},{aucWrapped:F4},{withAuc.Count} queries with both classes",
    $"queries_auc_better,,{withAuc.Count(p => p.AucWrapped > p.AucPlain + 1e-9)},",
    $"queries_auc_worse,,{withAuc.Count(p => p.AucWrapped < p.AucPlain - 1e-9)},",
    $"mean_same_genus_in_top10,{perQuery.Average(p => p.SameGenusTopPlain):F3},{perQuery.Average(p => p.SameGenusTopWrapped):F3},",
    $"mean_reciprocal_rank_first_same_genus,{perQuery.Average(p => p.ReciprocalRankPlain1000) / 1000:F4},{perQuery.Average(p => p.ReciprocalRankWrapped1000) / 1000:F4},",
    $"top10_hits_changed,,{topChanged:F4},share of top-10 hits not in the plain top 10",
    $"same_genus_pairs_gaining_over_0.03,,{(double)perQuery.Sum(p => p.GenusGains) / perQuery.Sum(p => p.GenusPairs):F4},of {perQuery.Sum(p => p.GenusPairs)} pairs",
    $"other_family_pairs_gaining_over_0.03,,{(double)perQuery.Sum(p => p.OtherGains) / Math.Max(1, perQuery.Sum(p => p.OtherPairs)):F4},of {perQuery.Sum(p => p.OtherPairs)} pairs",
    $"verdict,,{(aucWrapped - aucPlain < 0.005 && topChanged < 0.02 ? "stop" : "go")},rule: stop if AUC gain < 0.005 and top-10 change < 2%",
};
File.WriteAllLines(output, lines);
foreach (var line in lines) Console.WriteLine(line);

static int LowerBound(double[] sorted, double value)
{
    int lo = 0, hi = sorted.Length;
    while (lo < hi) { int mid = (lo + hi) / 2; if (sorted[mid] < value) lo = mid + 1; else hi = mid; }
    return lo;
}
static int UpperBound(double[] sorted, double value)
{
    int lo = 0, hi = sorted.Length;
    while (lo < hi) { int mid = (lo + hi) / 2; if (sorted[mid] <= value) lo = mid + 1; else hi = mid; }
    return lo;
}
static double Median(IEnumerable<double> values) { var v = values.Order().ToArray(); return v.Length == 0 ? double.NaN : v[v.Length / 2]; }
