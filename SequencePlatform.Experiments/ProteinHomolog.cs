using System.Collections.Concurrent;
using System.Globalization;
using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// protein-homolog (2026-11-12), moved here from scripts/platform-checks/protein_homolog_check.cs: does the right protein family come
    /// first without the query's own gene, genome, species or genus? The script drew its proteins from BV-BRC; here they come from the
    /// stored query set (<see cref="ExperimentQueries"/>, filled by the Ingest's seed-experiment-queries with the script's rule and seed),
    /// as their expected Z-curve axes, never as sequences. Compared as ProteinSearch compares them (<see cref="ZCurveGeneProfile.Score"/>,
    /// weighted since 2026-11-20, so a rerun is weighted where the article run was not). Per query and exclusion tier, over the genes left:
    /// the best hit's PGfam, and whether the query's PGfam is in the top 10. Tiers: gene (only the exact gene removed), genome (every record
    /// of the query's source), species, genus (records with the same NCBI tax_species / tax_genus; a query without that name is left out of
    /// the tier). A query counts in a tier only if some gene of its family is left there ("reachable").
    /// Rule fixed before the article run: family top-1 with the own species removed under 50% of reachable queries in any band fails.
    /// </summary>
    public sealed class ProteinHomolog : IExperiment
    {
        public const string ExperimentName = "protein-homolog";
        public const string FileName = "ProteinHomolog_Check.csv";

        private const int topK = 10;
        private static readonly string[] bandNames = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
        private static readonly string[] tierNames = ["gene", "genome", "species", "genus"];

        public string Name => ExperimentName;

        public static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;

        /// <summary>A stored query's fields: the gene it came from, located in its record's gene profile, and the protein's length.</summary>
        public sealed record QueryFields(string Kind, int Residues, int Start, int End, bool Reverse);

        public static ExperimentQuery ToQuery(long sequenceId, QueryFields fields, double[] axes) =>
            new(sequenceId, JsonSerializer.SerializeToElement(fields), axes);

        private sealed record Taxon(long Source, string? Species, string? Genus);

        private sealed record Query(int Band, int Residues, long SequenceId, Taxon Taxon, ZCurveGeneProfile.Gene Truth, double[] Axes, HashSet<int> Buckets);

        private sealed class Tally
        {
            public int SameFamily;
            public readonly PriorityQueue<int, double> Top = new();
        }

        public async Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken)
        {
            var snapshot = await Snapshots.LiveAsync(dataSource, cancellationToken);
            var stored = await ExperimentQueries.ReadAsync(dataSource, ExperimentName, cancellationToken);
            if (stored.Count == 0)
                throw new ExperimentNotRunnableException($"No query set is stored for {ExperimentName}: run the Ingest's seed-experiment-queries {ExperimentName} first.");

            var taxa = new Dictionary<long, Taxon>();
            await using (var command = dataSource.CreateCommand("select id, source_id, tax_species, tax_genus from sequence_record"))
            await using (var reader = await command.ExecuteReaderAsync(cancellationToken))
                while (await reader.ReadAsync(cancellationToken))
                    taxa[reader.GetInt64(0)] = new Taxon(reader.GetInt64(1), reader.IsDBNull(2) ? null : reader.GetString(2),
                        reader.IsDBNull(3) ? null : reader.GetString(3));

            // Each query's own gene as stored today, which gives its family; a query whose gene is no longer stored drops out.
            var queries = new List<Query>();
            var profiles = new Dictionary<long, ZCurveGeneProfile?>();
            foreach (var query in stored)
            {
                var fields = query.Fields.Deserialize<QueryFields>()!;
                if (!profiles.TryGetValue(query.SequenceId, out var profile))
                    profiles[query.SequenceId] = profile = await GeneProfiles.ReadAsync(dataSource, query.SequenceId);
                int gene = profile is null ? -1 : Enumerable.Range(0, profile.Genes.Count).FirstOrDefault(g => profile.Genes[g].Start == fields.Start
                    && profile.Genes[g].End == fields.End && profile.Genes[g].Reverse == fields.Reverse, -1);
                if (gene < 0 || !taxa.TryGetValue(query.SequenceId, out var taxon)) continue;
                queries.Add(new Query(BandOf(fields.Residues), fields.Residues, query.SequenceId, taxon, profile!.Genes[gene], query.Profile,
                    GeneProfiles.BucketsFor(fields.Residues, ProteinSearch.LengthBand).ToHashSet()));
            }
            var scoredQueries = queries.Where(q => q.Truth.Family != ZCurveGeneProfile.NoFamily).ToList();
            context.Progress($"{stored.Count} stored proteins, {queries.Count} with their gene still stored, {scoredQueries.Count} with a PGfam");

            var perWorker = await ScoreAsync(dataSource, scoredQueries, taxa, context, cancellationToken);

            // Per query and tier: merged top 10 over the workers, and whether any same-family gene was left.
            var outcomes = new List<(Query Query, int Tier, bool Reachable, bool Top1, bool Top10)>();
            for (int q = 0; q < scoredQueries.Count; q++)
                for (int t = 0; t < tierNames.Length; t++)
                {
                    if (!TierApplies(scoredQueries[q], t)) continue;
                    var merged = perWorker.SelectMany(w => w[q][t].Top.UnorderedItems).OrderByDescending(e => e.Priority).Take(topK).ToList();
                    int family = scoredQueries[q].Truth.Family;
                    outcomes.Add((scoredQueries[q], t, perWorker.Sum(w => w[q][t].SameFamily) > 0,
                        merged.Count > 0 && merged[0].Element == family, merged.Any(e => e.Element == family)));
                }

            static string Percent(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
            var lines = new List<string> { "Band,Tier,Proteins,Reachable,Family top-1,Family top-10" };
            var speciesShares = new Dictionary<int, double>();
            double Share(int tier, Func<(Query Query, int Tier, bool Reachable, bool Top1, bool Top10), bool> hit)
            {
                var reachable = outcomes.Where(o => o.Tier == tier && o.Reachable).ToList();
                return reachable.Count == 0 ? double.NaN : reachable.Count(hit) / (double)reachable.Count;
            }
            foreach (int band in Enumerable.Range(0, bandNames.Length).Append(-1))
                for (int t = 0; t < tierNames.Length; t++)
                {
                    var cell = outcomes.Where(o => o.Tier == t && (band < 0 || o.Query.Band == band)).ToList();
                    if (cell.Count == 0) continue;
                    var reachable = cell.Where(o => o.Reachable).ToList();
                    if (t == 2 && band >= 0 && reachable.Count > 0) speciesShares[band] = reachable.Count(o => o.Top1) / (double)reachable.Count;
                    lines.Add(string.Join(",", band < 0 ? "all" : bandNames[band], tierNames[t], cell.Count, reachable.Count,
                        Percent(reachable.Count(o => o.Top1), reachable.Count), Percent(reachable.Count(o => o.Top10), reachable.Count)));
                }
            lines.Add("");
            var failing = speciesShares.Where(s => s.Value < 0.5).Select(s => bandNames[s.Key]).ToList();
            lines.Add("Rule: species-removed family top-1 under 50% in any band reframes the claim," +
                (failing.Count == 0 ? "no band under 50%: identity curve goes ahead" : "under 50% in " + string.Join(" ", failing) + ": reframe"));

            var summary = new JsonObject
            {
                ["keyNumbers"] = new JsonArray(
                    KeyNumber("own gene removed, top-1", Share(0, o => o.Top1)),
                    KeyNumber("own genus removed, top-1", Share(3, o => o.Top1)),
                    KeyNumber("own genus removed, top-10", Share(3, o => o.Top10))),
                ["rule"] = failing.Count == 0
                    ? "Under 50% family top-1 with the own species removed in any band reframes the claim. No band was under 50%."
                    : $"Under 50% family top-1 with the own species removed in any band reframes the claim. Under 50% in {string.Join(", ", failing)}.",
                ["analysed"] = scoredQueries.Count,
                ["analysedNote"] = "proteins with a PGfam from the stored query set, weighted score",
            };
            var csv = string.Join("\n", lines) + "\n";
            return new ExperimentResult(failing.Count == 0 ? ExperimentStore.Pass : ExperimentStore.Fail, snapshot,
                "SequencePlatform.Experiments/ProteinHomolog.cs", summary.ToJsonString(), [new ExperimentFile(FileName, Encoding.UTF8.GetBytes(csv))], []);
        }

        private static JsonObject KeyNumber(string label, double value) =>
            new() { ["label"] = label, ["value"] = double.IsNaN(value) ? null : Math.Round(value, 4) };

        // Which tiers a candidate record stays in for a query: tier t keeps it if it isn't excluded at t.
        private static bool Kept(Query query, Taxon candidate, int tier) => tier switch
        {
            0 => true,
            1 => candidate.Source != query.Taxon.Source,
            2 => candidate.Source != query.Taxon.Source && (candidate.Species is null || candidate.Species != query.Taxon.Species),
            _ => candidate.Source != query.Taxon.Source && (candidate.Genus is null || candidate.Genus != query.Taxon.Genus),
        };

        private static bool TierApplies(Query query, int tier) =>
            tier switch { 2 => query.Taxon.Species is not null, 3 => query.Taxon.Genus is not null, _ => true };

        private static bool Comparable(Query query, ZCurveGeneProfile.Gene gene) => query.Residues >= ZCurveGeneProfile.LongResidues
            ? gene.Residues >= ZCurveGeneProfile.LongResidues
            : gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - query.Residues) <= ProteinSearch.LengthBand * query.Residues;

        // One pass over the gene rows the queries need, scored on every worker: per worker, query and tier, the top 10 families and how
        // many same-family genes were left.
        private static async Task<Tally[][][]> ScoreAsync(NpgsqlDataSource dataSource, List<Query> scoredQueries, Dictionary<long, Taxon> taxa,
            ExperimentContext context, CancellationToken cancellationToken)
        {
            var needed = scoredQueries.SelectMany(q => q.Buckets).Distinct().Select(b => (short)b).ToArray();
            var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
            var watch = System.Diagnostics.Stopwatch.StartNew();
            var workers = Enumerable.Range(0, context.Parallelism).Select(_ => Task.Run(() =>
            {
                var tallies = scoredQueries.Select(_ => tierNames.Select(_ => new Tally()).ToArray()).ToArray();
                foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable(cancellationToken))
                {
                    var candidate = taxa[id];
                    var profile = ZCurveGeneProfile.FromBinary(payload);
                    for (int q = 0; q < scoredQueries.Count; q++)
                    {
                        var query = scoredQueries[q];
                        if (!query.Buckets.Contains(bucket)) continue;
                        for (int g = 0; g < profile.Genes.Count; g++)
                        {
                            var gene = profile.Genes[g];
                            if (!Comparable(query, gene) || (id == query.SequenceId && gene.Start == query.Truth.Start && gene.End == query.Truth.End
                                && gene.Reverse == query.Truth.Reverse)) continue;
                            double score = profile.Score(query.Axes, g);
                            for (int t = 0; t < tierNames.Length; t++)
                            {
                                if (!Kept(query, candidate, t)) continue;
                                var tally = tallies[q][t];
                                if (gene.Family == query.Truth.Family) tally.SameFamily++;
                                if (tally.Top.Count < topK) tally.Top.Enqueue(gene.Family, score);
                                else if (tally.Top.TryPeek(out _, out double lowest) && score > lowest) tally.Top.EnqueueDequeue(gene.Family, score);
                            }
                        }
                    }
                }
                return tallies;
            }, cancellationToken)).ToArray();
            try
            {
                await using var command = dataSource.CreateCommand("select sequence_id, bucket, payload from gene_profile where bucket = any($1) and window_count > 0");
                command.Parameters.AddWithValue(needed);
                await using var reader = await command.ExecuteReaderAsync(cancellationToken);
                while (await reader.ReadAsync(cancellationToken)) rows.Add((reader.GetInt64(0), reader.GetInt16(1), reader.GetFieldValue<byte[]>(2)), cancellationToken);
            }
            finally
            {
                rows.CompleteAdding();
            }
            var perWorker = await Task.WhenAll(workers);
            context.Progress($"Scoring pass: {watch.Elapsed.TotalSeconds:F0} s");
            return perWorker;
        }
    }
}
