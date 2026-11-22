using System.Text;
using System.Text.Json.Nodes;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// neighbour-ranking (P3 S3.4b), ported from scripts/platform-checks/neighbour_eval.cs unchanged in what it computes: does re-ranking
    /// find-similar's top 100 by shared top-20 neighbours (the stored neighbour graph) find more same-genus records? 150 bacteria and 150
    /// viruses of at least 8 kbp with an NCBI genus and a same-genus record in their length band, searched through the production path
    /// (<see cref="SimilaritySearch.SearchRecordAsync"/>) with and without the neighbours; precision at 10 and the reciprocal rank of the first
    /// same-genus hit. Pass, as the article run read it: re-ranking raises bacteria's precision at 10 in more queries than it lowers it
    /// (sign test, one-sided p &lt; 0.05). The search is the production one, so a rerun measures today's search: rotation and wrap-around
    /// scoring have been added since the article run.
    /// </summary>
    public sealed class NeighbourRanking : IExperiment
    {
        public const string ExperimentName = "neighbour-ranking";
        public const string FileName = "S34b_NeighbourRanking.csv";

        private const int perKind = 150, minLength = 8_000, seed = 20261022;

        public string Name => ExperimentName;

        public async Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken)
        {
            var snapshot = await Snapshots.LiveAsync(dataSource, cancellationToken);
            var search = new SimilaritySearch(dataSource);
            var meta = new Dictionary<long, (long Length, string Kind, string? Genus)>();
            await using (var command = dataSource.CreateCommand(
                "select r.id, r.length_bp, r.kind, t.tax_genus from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid order by r.id"))
            await using (var reader = await command.ExecuteReaderAsync(cancellationToken))
                while (await reader.ReadAsync(cancellationToken))
                    meta[reader.GetInt64(0)] = (reader.GetInt64(1), reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetString(3));

            var rng = new Random(seed);
            var results = new List<(string Kind, string Ranking, double Precision, double ReciprocalRank)>();
            foreach (var kind in new[] { "bacterium", "virus" })
            {
                var pool = meta.Where(m => m.Value.Kind == kind && m.Value.Genus is not null && m.Value.Length >= minLength)
                    .Select(m => m.Key).OrderBy(_ => rng.Next()).ToList();
                int used = 0;
                foreach (var q in pool)
                {
                    if (used == perKind) break;
                    SimilarityResult plain;
                    try { plain = await search.SearchRecordAsync(q, top: 100, useNeighbours: false, cancellationToken: cancellationToken); }
                    catch (ArgumentException) { continue; }   // no stored 500 bp level
                    string genus = meta[q].Genus!;
                    var all = await search.SearchRecordAsync(q, top: 100_000, useNeighbours: false, cancellationToken: cancellationToken);
                    if (!all.Ranked.Any(h => meta[h.SequenceId].Genus == genus)) continue;   // no relative in the band
                    used++;

                    var ranked = await search.SearchRecordAsync(q, top: 100, useNeighbours: true, cancellationToken: cancellationToken);
                    foreach (var (name, list) in new[] { ("score", plain.Ranked), ("shared top-20 (production)", ranked.Ranked) })
                    {
                        var ids = list.Select(h => h.SequenceId).ToList();
                        int first = ids.FindIndex(id => meta[id].Genus == genus);
                        results.Add((kind, name, ids.Take(10).Count(id => meta[id].Genus == genus) / 10.0, 1.0 / (first < 0 ? 101 : first + 1)));
                    }
                    if (used % 50 == 0) context.Progress($"{kind}: {used} queries");
                }
            }

            // The script prints with the current culture, which on the author's machine writes these as invariant ones do; the server's may not.
            var csv = new StringBuilder("Kind, Ranking, Queries, MeanPrecisionAt10, MeanReciprocalRank, QueriesBetter, QueriesWorse\n");
            var keyNumbers = new JsonArray();
            bool pass = false;
            foreach (var kind in new[] { "bacterium", "virus" })
            {
                var score = results.Where(r => r.Kind == kind && r.Ranking == "score").ToList();
                var shared = results.Where(r => r.Kind == kind && r.Ranking != "score").ToList();
                int better = score.Zip(shared).Count(p => p.Second.Precision > p.First.Precision);
                int worse = score.Zip(shared).Count(p => p.Second.Precision < p.First.Precision);
                csv.Append(FormattableString.Invariant($"{kind}, score, {score.Count}, {score.Average(r => r.Precision):F3}, {score.Average(r => r.ReciprocalRank):F3}, , \n"));
                csv.Append(FormattableString.Invariant(
                    $"{kind}, shared top-20 (production), {shared.Count}, {shared.Average(r => r.Precision):F3}, {shared.Average(r => r.ReciprocalRank):F3}, {better}, {worse}\n"));
                if (kind == "bacterium")
                {
                    pass = better > worse && Habitat.SignTest(better, worse) < 0.05;
                    keyNumbers.Add(KeyNumber("bacteria precision@10, by score", score.Average(r => r.Precision)));
                    keyNumbers.Add(KeyNumber("bacteria precision@10, re-ranked", shared.Average(r => r.Precision)));
                }
                else
                {
                    keyNumbers.Add(KeyNumber("viruses precision@10, by score", score.Average(r => r.Precision)));
                    keyNumbers.Add(KeyNumber("viruses precision@10, re-ranked", shared.Average(r => r.Precision)));
                }
            }
            var summary = new JsonObject
            {
                ["keyNumbers"] = keyNumbers,
                ["analysed"] = results.Count / 2,
                ["analysedNote"] = $"queries, {perKind} bacteria and {perKind} viruses with a same-genus record in their band, through today's search",
            };
            return new ExperimentResult(pass ? ExperimentStore.Pass : ExperimentStore.Fail, snapshot, "SequencePlatform.Experiments/NeighbourRanking.cs",
                summary.ToJsonString(), [new ExperimentFile(FileName, Encoding.UTF8.GetBytes(csv.ToString()))], []);
        }

        private static JsonObject KeyNumber(string label, double value) => new() { ["label"] = label, ["value"] = Math.Round(value, 3) };
    }
}
