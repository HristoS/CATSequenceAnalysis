#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S3.4 pre-check: what one find-similar costs on the local catalogue, from the in-memory caches. It measured the
// engine before the length band (window-count short list, 500 bp re-rank of the top 50), so it scores unbanded.
const int minWindows = 8;
var root = "data/genomes/";
int[] factors = { 1, 4, 16, 64 };
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var caches = new Dictionary<int, IReadOnlyDictionary<long, ZCurveWindowedProfile>>();
foreach (var f in factors)
{
    var cache = new ProfileCache(dataSource, 500 * f);
    await cache.RefreshAsync();
    caches[f] = cache.Profiles;
}
var kinds = new Dictionary<long, string>();
await using (var cmd = dataSource.CreateCommand("select id, kind from sequence_record"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync()) kinds[r.GetInt64(0)] = r.GetString(1);

string Rc(string s)
{
    var c = new char[s.Length];
    for (int i = 0; i < s.Length; i++)
        c[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c', var o => o };
    return new string(c);
}

double Both(ZCurveWindowedProfile f, ZCurveWindowedProfile r, ZCurveWindowedProfile c) =>
    Math.Max(ZCurveWindowedProfile.Compare(f, c), ZCurveWindowedProfile.Compare(r, c));

Console.WriteLine("Query, LengthBp, QueryWindowSize, Filter, Ranked, ShortList, RankedMs, ShortListMs, Rerank50Ms, Threads");
foreach (var (name, file) in new[] { ("E. coli K-12", "ecoli_k12_mg1655_GCF_000005845.2.fna.gz"), ("SARS-CoV-2", "viral/sars_cov_2_NC_045512.2.fna.gz"), ("HPV16", "viral/hpv16_NC_001526.4.fna.gz") })
{
    var dna = FastaReader.ReadPrimaryRecord(root + file).Sequence;
    var fwd = new ZCurveProfileLadder(new ZCurveWindowedProfile(dna));
    var rev = new ZCurveProfileLadder(new ZCurveWindowedProfile(Rc(dna)));
    int qf = ZCurveProfileLadder.SelectQueryFactor(fwd)!.Value;
    string queryKind = file.StartsWith("viral/") ? "virus" : "bacterium";

    foreach (var filter in new[] { "all", "same kind" })
    {
        var candidates = caches[1].Keys.Where(id => filter == "all" || kinds[id] == queryKind).ToArray();
        int threshold = Math.Min(minWindows, fwd.Level(qf).WindowCount);
        var ranked = candidates.Where(id => caches[qf].TryGetValue(id, out var p) && p.WindowCount >= threshold).ToArray();
        var shortList = candidates.Except(ranked).ToArray();

        // Warm-up, then timed.
        var rankedScores = new double[ranked.Length];
        for (int pass = 0; pass < 2; pass++)
        {
            var watch = Stopwatch.StartNew();
            Parallel.For(0, ranked.Length, i => rankedScores[i] = Both(fwd.Level(qf), rev.Level(qf), caches[qf][ranked[i]]));
            var rankedMs = watch.Elapsed.TotalMilliseconds;

            watch.Restart();
            Parallel.For(0, shortList.Length, i =>
            {
                var fine = caches[1][shortList[i]];
                int f = factors.Last(f => f == 1 || (caches[f].TryGetValue(shortList[i], out var p) && p.WindowCount >= minWindows));
                _ = Both(fwd.Level(f), rev.Level(f), caches[f][shortList[i]]);
            });
            var shortMs = watch.Elapsed.TotalMilliseconds;

            var top = Enumerable.Range(0, ranked.Length).OrderByDescending(i => rankedScores[i]).Take(50).Select(i => ranked[i]).ToArray();
            watch.Restart();
            Parallel.For(0, top.Length, i => _ = Both(fwd.Fine, rev.Fine, caches[1][top[i]]));
            var rerankMs = watch.Elapsed.TotalMilliseconds;

            if (pass == 1)
                Console.WriteLine($"{name}, {dna.Length}, {500 * qf}, {filter}, {ranked.Length}, {shortList.Length}, {rankedMs:F0}, {shortMs:F0}, {rerankMs:F0}, {Environment.ProcessorCount}");
        }
    }
}
