#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Check (2026-10-22): does re-ranking by shared neighbours beat raw score over many queries, not three?
// 150 random bacteria and 150 random viruses with an NCBI genus are queries (stored profile, profile-level flip for
// the other strand). Each query's top 100 (other assemblies only) is re-ranked by |top-20(q) ∩ top-20(c)|, then by
// score; also by mutual rank. Metrics per query: precision@10 (same NCBI genus) and reciprocal rank of the first
// same-genus hit; averaged over queries that have at least one same-genus record in their band.
const int topN = 100, shareK = 20, perKind = 150;
const double ratio = 0.5;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var search = new SimilaritySearch(dataSource);

var meta = new Dictionary<long, (long Source, long Length, string Kind, string? Genus)>();
await using (var cmd = dataSource.CreateCommand(
    "select r.id, r.source_id, r.length_bp, r.kind, t.tax_genus from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync()) meta[r.GetInt64(0)] = (r.GetInt64(1), r.GetInt64(2), r.GetString(3), r.IsDBNull(4) ? null : r.GetString(4));

ZCurveWindowedProfile Flip(ZCurveWindowedProfile p) =>
    ZCurveWindowedProfile.FromPrecomputedWindows(p.WindowSize, p.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

var profiles = new Dictionary<long, QueryProfile>();
async Task<QueryProfile?> Own(long id)
{
    if (profiles.TryGetValue(id, out var cached)) return cached;
    var fine = await ProfileLevels.ReadAsync(dataSource, id, RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
    if (fine is null) return null;
    return profiles[id] = QueryProfile.FromProfiles(fine, Flip(fine)) with { LengthBp = meta[id].Length };
}

// A record's list without its own assembly (its other chromosomes and plasmids aren't "relatives").
async Task<List<SimilarityHit>> ListOf(long id, QueryProfile query) =>
    (await search.SearchAsync(query, lengthRatio: ratio, top: 100_000)).Ranked.Where(h => meta[h.SequenceId].Source != meta[id].Source).ToList();

var rng = new Random(20261022);
var results = new List<(string Kind, string Ranking, double Precision, double ReciprocalRank)>();
var watch = Stopwatch.StartNew();
foreach (var kind in new[] { "bacterium", "virus" })
{
    var pool = meta.Where(m => m.Value.Kind == kind && m.Value.Genus is not null && m.Value.Length >= 8_000).Select(m => m.Key).OrderBy(_ => rng.Next()).ToList();
    int used = 0;
    foreach (var q in pool)
    {
        if (used == perKind) break;
        var query = await Own(q);
        if (query is null) continue;
        var list = await ListOf(q, query);
        string genus = meta[q].Genus!;
        if (!list.Any(h => meta[h.SequenceId].Genus == genus)) continue;   // no relative in the band: nothing to rank
        used++;

        var top = list.Take(topN).ToList();
        var queryTopK = top.Take(shareK).Select(h => h.SequenceId).ToHashSet();
        var keyed = new List<(long Id, double Score, int Rank, int Mutual, int Shared)>();
        for (int i = 0; i < top.Count; i++)
        {
            long c = top[i].SequenceId;
            var own = await Own(c);
            var ownList = own is null ? new List<SimilarityHit>() : await ListOf(c, own);
            int back = ownList.FindIndex(h => h.SequenceId == q);
            keyed.Add((c, top[i].Score, i + 1, Math.Max(i + 1, back < 0 ? int.MaxValue : back + 1),
                ownList.Take(shareK).Count(h => queryTopK.Contains(h.SequenceId))));
        }

        void Score(string ranking, IEnumerable<long> order)
        {
            var ids = order.ToList();
            double precision = ids.Take(10).Count(id => meta[id].Genus == genus) / 10.0;
            int first = ids.FindIndex(id => meta[id].Genus == genus);
            // A relative past the top 100 counts as rank 101 for every ranking alike.
            results.Add((kind, ranking, precision, 1.0 / (first < 0 ? topN + 1 : first + 1)));
        }
        Score("score", keyed.Select(k => k.Id));
        Score("mutual rank", keyed.OrderBy(k => k.Mutual).ThenByDescending(k => k.Score).Select(k => k.Id));
        Score($"shared top-{shareK}", keyed.OrderByDescending(k => k.Shared).ThenByDescending(k => k.Score).Select(k => k.Id));
    }
    Console.Error.WriteLine($"{kind}: {used} queries, {watch.Elapsed.TotalSeconds:F0} s");
}

Console.WriteLine("Kind, Ranking, Queries, MeanPrecisionAt10, MeanReciprocalRank");
foreach (var g in results.GroupBy(r => (r.Kind, r.Ranking)))
    Console.WriteLine($"{g.Key.Kind}, {g.Key.Ranking}, {g.Count()}, {g.Average(r => r.Precision):F3}, {g.Average(r => r.ReciprocalRank):F3}");
