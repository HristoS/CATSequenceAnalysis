#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Check (2026-10-22): re-rank a query's top N by mutual-neighbour agreement instead of raw score.
// mutual rank = max(rank of c in q's list, rank of q in c's list); shared = |top-K(q) ∩ top-K(c)|.
// Stored records have no DNA, so a candidate's own search uses its stored 500 bp profile and a profile-level
// flip for the other strand (approximate, results/C3_ProfileFlip.csv), with its stored length.
const int topN = 300, shareK = 20;
const double ratio = 0.5;
var root = "data/genomes/";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var search = new SimilaritySearch(dataSource);

var meta = new Dictionary<long, (long Length, string Genus, string Name)>();
await using (var cmd = dataSource.CreateCommand("select id, length_bp, coalesce(tax_genus, ''), coalesce(name, '') from sequence_record"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync()) meta[r.GetInt64(0)] = (r.GetInt64(1), r.GetString(2), r.GetString(3));

ZCurveWindowedProfile Flip(ZCurveWindowedProfile p) =>
    ZCurveWindowedProfile.FromPrecomputedWindows(p.WindowSize, p.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

Console.WriteLine("Query, Target, Ranking, FirstRanksWithinTop300, TargetsInTop300");
var watch = Stopwatch.StartNew();
foreach (var (name, file, target) in new (string, string, Func<long, bool>)[]
{
    ("E. coli K-12", "ecoli_k12_mg1655_GCF_000005845.2.fna.gz", id => meta[id].Genus is "Escherichia" or "Shigella"),
    ("SARS-CoV-2", "viral/sars_cov_2_NC_045512.2.fna.gz", id => meta[id].Genus == "Betacoronavirus"),
    ("HPV16", "viral/hpv16_NC_001526.4.fna.gz", id => meta[id].Name.Contains("papillomavirus", StringComparison.OrdinalIgnoreCase)),
})
{
    var query = QueryProfile.FromSequence(FastaReader.ReadPrimaryRecord(root + file).Sequence);
    var list = (await search.SearchAsync(query, lengthRatio: ratio, top: 100_000)).Ranked;
    var top = list.Take(topN).ToList();
    var queryTopK = top.Take(shareK).Select(h => h.SequenceId).ToHashSet();

    var keys = new List<(long Id, double Score, int OwnRank, int MutualRank, int Shared)>();
    for (int i = 0; i < top.Count; i++)
    {
        long id = top[i].SequenceId;
        var fine = await ProfileLevels.ReadAsync(dataSource, id, RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
        var own = QueryProfile.FromProfiles(fine!, Flip(fine!)) with { LengthBp = meta[id].Length };
        var ownList = (await search.SearchAsync(own, lengthRatio: ratio, top: 100_000)).Ranked;

        // Where the query would sit in the candidate's list: the candidate's score for the query, at the candidate's level.
        var level = ZCurveProfileLadder.SelectBandFactor(own.Forward, ratio) ?? 1;
        double back = Math.Max(ZCurveWindowedProfile.Compare(own.Forward.Level(level), query.Forward.Level(level)),
                               ZCurveWindowedProfile.Compare(own.Reverse!.Level(level), query.Forward.Level(level)));
        int reverseRank = 1 + ownList.Count(h => h.SequenceId != id && h.Score > back);
        int shared = ownList.Take(shareK).Count(h => queryTopK.Contains(h.SequenceId));
        keys.Add((id, top[i].Score, i + 1, Math.Max(i + 1, reverseRank), shared));
    }

    void Report(string ranking, IEnumerable<(long Id, double Score, int OwnRank, int MutualRank, int Shared)> order)
    {
        var ranks = order.Select((k, i) => (k.Id, Rank: i + 1)).Where(x => target(x.Id)).Select(x => x.Rank).ToList();
        Console.WriteLine($"{name}, {(name == "E. coli K-12" ? "Escherichia or Shigella" : name == "HPV16" ? "any papillomavirus" : "Betacoronavirus")}, " +
            $"{ranking}, \"{string.Join("; ", ranks.Take(8))}\", {ranks.Count}");
    }
    Report("score", keys);
    Report("mutual rank", keys.OrderBy(k => k.MutualRank).ThenByDescending(k => k.Score));
    Report($"shared top-{shareK}", keys.OrderByDescending(k => k.Shared).ThenByDescending(k => k.Score));
}
Console.Error.WriteLine($"{watch.Elapsed.TotalSeconds:F0} s");
