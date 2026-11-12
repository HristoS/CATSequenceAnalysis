#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using Npgsql;
using SequencePlatform.Core;

// Should find-similar skip the neighbour re-ranking for viruses now that it scores them wrapped (results/WrapAround_Search.csv: on
// 300 virus queries re-ranking moved precision@10 0.668 -> 0.670 but mean reciprocal rank 0.908 -> 0.869)? Fresh queries (another
// seed): 600 public viruses and 300 public bacteria with an NCBI genus and a same-genus record of another source in the band, each
// searched by score alone and re-ranked by the stored graph (built with wrapping). Three rankings: score, neighbours (today), and the
// rule "score alone when most of the query's top 20 by score are viruses", which also works for an upload of unknown kind.
// Rule fixed before the run: skip for viruses if, on the virus queries, the rule gains >= 0.02 mean reciprocal rank with a 95%
// bootstrap interval above 0 and loses <= 0.01 precision@10; and it changes the ranking of < 1% of the bacterial queries.
// Run from the repo root with SEQUENCE_PLATFORM_DB set, after build-null and build-neighbours.
var output = args.Length > 0 ? args[0] : "results/VirusRerank_Check.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var kinds = new Dictionary<long, string>();
var meta = new Dictionary<long, (long Length, long Source, string Kind, string Genus)>();
await using (var command = dataSource.CreateCommand(
    "select r.id, r.length_bp, r.source_id, r.kind, t.tax_genus, s.owner_id is null from sequence_record r " +
    "join sequence_source s on s.id = r.source_id left join taxon t on t.taxid = r.ncbi_taxid"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
    {
        kinds[reader.GetInt64(0)] = reader.GetString(3);
        if (!reader.IsDBNull(4) && reader.GetBoolean(5) && reader.GetInt64(1) >= 8000)
            meta[reader.GetInt64(0)] = (reader.GetInt64(1), reader.GetInt64(2), reader.GetString(3), reader.GetString(4));
    }
var byGenus = meta.GroupBy(m => m.Value.Genus).ToDictionary(g => g.Key, g => g.Select(m => m.Key).ToArray());

var search = new SimilaritySearch(dataSource);
await search.LoadAllAsync();
var rows = new List<string> { "kind,queries,ranking,precision_at_10,mean_reciprocal_rank,top1,reciprocal_rank_gain_vs_neighbours,ci95_low,ci95_high,queries_rule_skips" };
foreach (var (kind, count) in new[] { ("virus", 600), ("bacterium", 300) })
{
    var rng = new Random(20261113);
    var picked = new List<long>();
    foreach (var id in meta.Where(m => m.Value.Kind == kind).Select(m => m.Key).Order().OrderBy(_ => rng.Next()))
    {
        if (picked.Count == count) break;
        var (length, source, _, genus) = meta[id];
        var band = SimilaritySearch.Band(length, SimilaritySearch.DefaultLengthRatio);
        if (byGenus[genus].Any(o => o != id && meta[o].Source != source && meta[o].Length >= band.Min && meta[o].Length < band.MaxOrUnbounded))
            picked.Add(id);
    }

    var measures = new Dictionary<string, List<(double Precision, double Reciprocal, double Top1)>> { ["score"] = [], ["neighbours"] = [], ["rule"] = [] };
    int skips = 0;
    foreach (var q in picked)
    {
        var byScore = (await search.SearchRecordAsync(q, top: 100, useNeighbours: false)).Ranked.Select(h => h.SequenceId).ToList();
        var byNeighbours = (await search.SearchRecordAsync(q, top: 100, useNeighbours: true)).Ranked.Select(h => h.SequenceId).ToList();
        bool viral = byScore.Take(SimilaritySearch.NeighbourCount).Count(id => kinds.GetValueOrDefault(id) == "virus") * 2 > Math.Min(SimilaritySearch.NeighbourCount, byScore.Count);
        if (viral) skips++;
        string genus = meta[q].Genus;
        bool Same(long id) => meta.TryGetValue(id, out var m) && m.Genus == genus;
        (double, double, double) Measure(List<long> ids)
        {
            int first = ids.FindIndex(Same);
            return (ids.Take(10).Count(Same) / 10.0, 1.0 / (first < 0 ? 101 : first + 1), first == 0 ? 1 : 0);
        }
        measures["score"].Add(Measure(byScore));
        measures["neighbours"].Add(Measure(byNeighbours));
        measures["rule"].Add(Measure(viral ? byScore : byNeighbours));
    }

    foreach (var (ranking, list) in measures)
    {
        var gains = list.Zip(measures["neighbours"], (a, b) => a.Reciprocal - b.Reciprocal).ToArray();
        var (low, high) = Bootstrap(gains, 2000, new Random(7));
        var row = $"{kind},{picked.Count},{ranking},{list.Average(m => m.Precision):F4},{list.Average(m => m.Reciprocal):F4},{list.Average(m => m.Top1):F4}," +
            $"{gains.Average():F4},{low:F4},{high:F4},{skips}";
        rows.Add(row);
        Console.WriteLine(row);
    }
}
File.WriteAllLines(output, rows);

static (double Low, double High) Bootstrap(double[] values, int rounds, Random rng)
{
    var means = new double[rounds];
    for (int r = 0; r < rounds; r++)
    {
        double sum = 0;
        for (int i = 0; i < values.Length; i++) sum += values[rng.Next(values.Length)];
        means[r] = sum / values.Length;
    }
    Array.Sort(means);
    return (means[(int)(0.025 * rounds)], means[(int)(0.975 * rounds) - 1]);
}
