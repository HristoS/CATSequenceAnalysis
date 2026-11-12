#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Diagnostics;
using Npgsql;
using SequencePlatform.Core;

// Stop criterion for wrap-around scoring in find-similar (results/CircularPermutation_Check.csv measured it outside the search):
// through SimilaritySearch.SearchRecordAsync itself, 300 virus queries with an NCBI genus and a same-genus record of another source
// in the length band, WrapViruses off and on. Rule fixed before the run: stop if the production path keeps less than half the gain
// (mean reciprocal rank or precision@10 gain under half of the check's) or virus searches get more than 1.5x slower.
//   score:      both settings ranked by score alone (no neighbour graph), and the time per search, cache warm.
//   production: the setting given (plain | wrapped) as users get it, re-ranked by the stored neighbour graph; run plain before
//               rebuilding build-null and build-neighbours, and wrapped after.
// Rows are appended to results/WrapAround_Search.csv. Run from the repo root with SEQUENCE_PLATFORM_DB set.
const int queries = 300;
var mode = args.Length > 0 ? args[0] : "score";
var output = "results/WrapAround_Search.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var meta = new Dictionary<long, (long Length, long Source, string Genus)>();
await using (var command = dataSource.CreateCommand(
    "select r.id, r.length_bp, r.source_id, t.tax_genus from sequence_record r join taxon t on t.taxid = r.ncbi_taxid " +
    "join sequence_source s on s.id = r.source_id where r.kind = 'virus' and t.tax_genus is not null and s.owner_id is null and r.length_bp >= 8000"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync()) meta[reader.GetInt64(0)] = (reader.GetInt64(1), reader.GetInt64(2), reader.GetString(3));
var byGenus = meta.GroupBy(m => m.Value.Genus).ToDictionary(g => g.Key, g => g.Select(m => m.Key).ToArray());

var rng = new Random(20261112);
var picked = new List<long>();
foreach (var id in meta.Keys.Order().OrderBy(_ => rng.Next()))
{
    if (picked.Count == queries) break;
    var (length, source, genus) = meta[id];
    var band = SimilaritySearch.Band(length, SimilaritySearch.DefaultLengthRatio);
    if (byGenus[genus].Any(o => o != id && meta[o].Source != source && meta[o].Length >= band.Min && meta[o].Length < band.MaxOrUnbounded))
        picked.Add(id);
}
Console.WriteLine($"{meta.Count:N0} public viruses with a genus; {picked.Count} queries");

var settings = mode switch
{
    "score" => new[] { (Name: "plain", Wrap: false, Neighbours: false), ("wrapped", true, false) },
    "production" => new[] { (Name: args[1], Wrap: args[1] == "wrapped", Neighbours: true) },
    _ => throw new ArgumentException("mode: score | production plain|wrapped"),
};
var rows = new List<string>();
// Every setting warmed with a full pass first (caches, wrapped forms, graph: a long-running worker has them), then timed alternating
// query by query, so neither pays the other's cold start or a background load.
var searches = new List<(string Name, bool Neighbours, SimilaritySearch Search, Stopwatch Watch, double[] Totals)>();
foreach (var (name, wrap, neighbours) in settings)
{
    var search = new SimilaritySearch(dataSource) { WrapViruses = wrap };
    await search.LoadAllAsync();
    foreach (var q in picked) await search.SearchRecordAsync(q, top: 100, useNeighbours: neighbours);
    searches.Add((name, neighbours, search, new Stopwatch(), new double[2]));
}
foreach (var q in picked)
    foreach (var (_, neighbours, search, watch, totals) in searches)
    {
        watch.Start();
        var result = await search.SearchRecordAsync(q, top: 100, useNeighbours: neighbours);
        watch.Stop();
        var ids = result.Ranked.Select(h => h.SequenceId).ToList();
        string genus = meta[q].Genus;
        bool Same(long id) => meta.TryGetValue(id, out var m) && m.Genus == genus;
        int first = ids.FindIndex(Same);
        totals[0] += ids.Take(10).Count(Same) / 10.0;
        totals[1] += 1.0 / (first < 0 ? 101 : first + 1);
    }
foreach (var (name, neighbours, _, watch, totals) in searches)
{
    var row = $"{DateTime.UtcNow:yyyy-MM-dd},{mode},{name},{(neighbours ? "neighbours" : "score")},{picked.Count},{totals[0] / picked.Count:F4}," +
        $"{totals[1] / picked.Count:F4},{watch.Elapsed.TotalMilliseconds / picked.Count:F1}";
    rows.Add(row);
    Console.WriteLine(row);
}
if (!File.Exists(output)) File.WriteAllLines(output, ["date,mode,setting,ranking,queries,precision_at_10_same_genus,mean_reciprocal_rank,ms_per_search"]);
File.AppendAllLines(output, rows);
