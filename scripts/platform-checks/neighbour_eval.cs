#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using Npgsql;
using SequencePlatform.Core;

// S3.4b stop criterion: the 300-query evaluation of results/Check_MutualNeighboursSample.csv, through the production
// path (SimilaritySearch.SearchRecordAsync with and without the stored neighbour graph). Same sampling rule and seed;
// the pool is ordered by id this time, so the queries differ from the check's.
const int perKind = 150;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var search = new SimilaritySearch(dataSource);

var meta = new Dictionary<long, (long Length, string Kind, string? Genus)>();
await using (var cmd = dataSource.CreateCommand(
    "select r.id, r.length_bp, r.kind, t.tax_genus from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid order by r.id"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync()) meta[r.GetInt64(0)] = (r.GetInt64(1), r.GetString(2), r.IsDBNull(3) ? null : r.GetString(3));

var rng = new Random(20261022);
var results = new List<(string Kind, string Ranking, double Precision, double ReciprocalRank)>();
foreach (var kind in new[] { "bacterium", "virus" })
{
    var pool = meta.Where(m => m.Value.Kind == kind && m.Value.Genus is not null && m.Value.Length >= 8_000)
        .Select(m => m.Key).OrderBy(_ => rng.Next()).ToList();
    int used = 0;
    foreach (var q in pool)
    {
        if (used == perKind) break;
        SimilarityResult plain;
        try { plain = await search.SearchRecordAsync(q, top: 100, useNeighbours: false); }
        catch (ArgumentException) { continue; }   // no stored 500 bp level
        string genus = meta[q].Genus!;
        var all = await search.SearchRecordAsync(q, top: 100_000, useNeighbours: false);
        if (!all.Ranked.Any(h => meta[h.SequenceId].Genus == genus)) continue;   // no relative in the band
        used++;

        var ranked = await search.SearchRecordAsync(q, top: 100, useNeighbours: true);
        foreach (var (name, list) in new[] { ("score", plain.Ranked), ("shared top-20 (production)", ranked.Ranked) })
        {
            var ids = list.Select(h => h.SequenceId).ToList();
            int first = ids.FindIndex(id => meta[id].Genus == genus);
            results.Add((kind, name, ids.Take(10).Count(id => meta[id].Genus == genus) / 10.0, 1.0 / (first < 0 ? 101 : first + 1)));
        }
    }
}

Console.WriteLine("Kind, Ranking, Queries, MeanPrecisionAt10, MeanReciprocalRank, QueriesBetter, QueriesWorse");
foreach (var kind in new[] { "bacterium", "virus" })
{
    var score = results.Where(r => r.Kind == kind && r.Ranking == "score").ToList();
    var shared = results.Where(r => r.Kind == kind && r.Ranking != "score").ToList();
    int better = score.Zip(shared).Count(p => p.Second.Precision > p.First.Precision);
    int worse = score.Zip(shared).Count(p => p.Second.Precision < p.First.Precision);
    Console.WriteLine($"{kind}, score, {score.Count}, {score.Average(r => r.Precision):F3}, {score.Average(r => r.ReciprocalRank):F3}, , ");
    Console.WriteLine($"{kind}, shared top-20 (production), {shared.Count}, {shared.Average(r => r.Precision):F3}, {shared.Average(r => r.ReciprocalRank):F3}, {better}, {worse}");
}
