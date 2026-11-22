using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// gc-stratum (S5.1, P5 R2), moved here from scripts/platform-checks/gc_stratum_check.cs unchanged in what it computes: does 1 - Compare
    /// separate taxa beyond GC content on pairs whose GC is nearly equal? R1 predicts the added value has to show in the low-|ΔGC| stratum,
    /// where GC-only is near chance. Genomes per D9: the largest bacterial or archaeal record of each source, labelled by GTDB
    /// (gtdb_taxon), at least 256 kbp so even the 32 kbp level keeps 8 windows; the floor doesn't depend on the level, so every level
    /// scores the same sample. One genome per GTDB species; up to 150 genera with at least 2 species, up to 4 species each, so same-genus
    /// pairs are different species. All pairs compared at one fixed level (default 32 kbp, the level find-similar uses for chromosomes),
    /// as stored (no strand flip or origin rotation: those are R2 candidates, not the baseline). GC from the stored dinucleotide counts
    /// (first base of each pair). AUC by rank sum per rank (genus, family, phylum: positive = same taxon at that rank) and stratum of |ΔGC|;
    /// 95% CIs from 200 bootstrap resamples of genomes, and the paired CI of AUC(Z-curve) - AUC(GC). No pass rule: a context experiment.
    /// </summary>
    public sealed class GcStratum(int windowSize = GcStratum.DefaultWindowSize) : IExperiment
    {
        public const string ExperimentName = "gc-stratum";
        public const string FileName = "S51_GcStratumCheck.csv";
        public const int DefaultWindowSize = 32_000;

        private const int maxGenera = 150, maxSpeciesPerGenus = 4, resamples = 200, minLength = 256_000, seed = 20261029;
        private const string header = "WindowSize, Rank, Stratum, Positives, Negatives, AucZCurve, ZCurveCiLow, ZCurveCiHigh, AucGc, GcCiLow, GcCiHigh, " +
            "AucDifference, DifferenceCiLow, DifferenceCiHigh";

        private static readonly (string Name, Func<Genome, string> Of)[] ranks =
            [("genus", g => g.Genus), ("family", g => g.Family), ("phylum", g => g.Phylum)];
        private static readonly (string Name, double Min, double Max)[] strata =
            [("all", 0, 1), ("dGC<0.02", 0, 0.02), ("dGC 0.02-0.10", 0.02, 0.10), ("dGC>=0.10", 0.10, 1)];

        public string Name => ExperimentName;

        public sealed record Genome(long Id, long Length, double Gc, string Phylum, string Family, string Genus, string Species);

        /// <param name="Auc">Rank-sum AUC; NaN when a class is empty.</param>
        public sealed record Row(int WindowSize, string Rank, string Stratum, int Positives, int Negatives, double AucZCurve, double ZCurveCiLow,
            double ZCurveCiHigh, double AucGc, double GcCiLow, double GcCiHigh, double AucDifference, double DifferenceCiLow, double DifferenceCiHigh)
        {
            public string ToCsv() => string.Format(CultureInfo.InvariantCulture,
                "{0}, {1}, {2}, {3}, {4}, {5:F4}, {6:F4}, {7:F4}, {8:F4}, {9:F4}, {10:F4}, {11:F4}, {12:F4}, {13:F4}", WindowSize, Rank, Stratum,
                Positives, Negatives, AucZCurve, ZCurveCiLow, ZCurveCiHigh, AucGc, GcCiLow, GcCiHigh, AucDifference, DifferenceCiLow, DifferenceCiHigh);
        }

        public async Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken)
        {
            var snapshot = await Snapshots.LiveAsync(dataSource, cancellationToken);
            var (sample, levels) = await SampleAsync(dataSource, context, cancellationToken);
            var rows = Tabulate(sample, levels, windowSize, context, cancellationToken);

            var csv = new StringBuilder(header).Append('\n');
            foreach (var row in rows) csv.Append(row.ToCsv()).Append('\n');
            Row Cell(string rank, string stratum) => rows.Single(r => r.Rank == rank && r.Stratum == stratum);
            var genus = Cell("genus", "dGC<0.02");
            var family = Cell("family", "dGC<0.02");
            var summary = new JsonObject
            {
                ["keyNumbers"] = new JsonArray(
                    KeyNumber("same-GC genus AUC, Z-curve", genus.AucZCurve),
                    KeyNumber("same-GC genus AUC, GC", genus.AucGc),
                    KeyNumber("same-GC family, Z-curve minus GC", family.AucDifference)),
                ["analysed"] = sample.Count,
                ["analysedNote"] = $"one per GTDB species, {sample.Select(g => g.Genus).Distinct().Count()} genera, compared at {windowSize / 1000} kbp",
            };
            return new ExperimentResult(ExperimentStore.Context, snapshot, "SequencePlatform.Experiments/GcStratum.cs", summary.ToJsonString(),
                [new ExperimentFile(FileName, Encoding.UTF8.GetBytes(csv.ToString()))], []);
        }

        private static JsonObject KeyNumber(string label, double value) =>
            new() { ["label"] = label, ["value"] = double.IsNaN(value) ? null : Math.Round(value, 4) };

        private async Task<(List<Genome> Sample, IReadOnlyDictionary<long, ZCurveWindowedProfile> Levels)> SampleAsync(NpgsqlDataSource dataSource,
            ExperimentContext context, CancellationToken cancellationToken)
        {
            var all = new List<Genome>();
            await using (var command = dataSource.CreateCommand(
                "select distinct on (r.source_id) r.id, r.length_bp, r.dinucleotide_counts, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species " +
                "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
                "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
                "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null " +
                "order by r.source_id, r.length_bp desc, r.id"))
            await using (var reader = await command.ExecuteReaderAsync(cancellationToken))
                while (await reader.ReadAsync(cancellationToken))
                {
                    if (reader.GetInt64(1) < minLength || reader.IsDBNull(6) || reader.IsDBNull(5) || reader.IsDBNull(4) || reader.IsDBNull(3)) continue;
                    var counts = reader.GetFieldValue<long[]>(2);
                    double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
                    all.Add(new Genome(reader.GetInt64(0), reader.GetInt64(1), gc, reader.GetString(3), reader.GetString(4), reader.GetString(5),
                        reader.GetString(6)));
                }

            var random = new Random(seed);
            var perSpecies = all.GroupBy(g => g.Species).Select(s => s.OrderBy(_ => random.Next()).First()).ToList();
            var sample = perSpecies.GroupBy(g => g.Genus).Where(s => s.Count() >= 2).OrderBy(_ => random.Next()).Take(maxGenera)
                .SelectMany(s => s.OrderBy(_ => random.Next()).Take(maxSpeciesPerGenus)).ToList();

            var levels = await ProfileLevels.ReadManyAsync(dataSource, sample.Select(g => g.Id).ToArray(), RecordStore.Method, RecordStore.MethodVersion,
                windowSize, cancellationToken);
            sample = sample.Where(g => levels.ContainsKey(g.Id)).ToList();
            context.Progress($"{all.Count:N0} labelled genomes, {perSpecies.Count:N0} species; sample {sample.Count} genomes in " +
                $"{sample.Select(g => g.Genus).Distinct().Count()} genera, {sample.Select(g => g.Family).Distinct().Count()} families, " +
                $"{sample.Select(g => g.Phylum).Distinct().Count()} phyla");
            return (sample, levels);
        }

        /// <summary>Every rank and stratum's AUCs with bootstrap CIs, in the CSV's row order.</summary>
        public static IReadOnlyList<Row> Tabulate(IReadOnlyList<Genome> sample, IReadOnlyDictionary<long, ZCurveWindowedProfile> levels, int windowSize,
            ExperimentContext context, CancellationToken cancellationToken)
        {
            var parallel = new ParallelOptions { MaxDegreeOfParallelism = context.Parallelism, CancellationToken = cancellationToken };
            int n = sample.Count;
            var zCurve = new double[n, n];
            Parallel.For(0, n, parallel, a =>
            {
                for (int b = a + 1; b < n; b++)
                    zCurve[a, b] = zCurve[b, a] = ZCurveWindowedProfile.Compare(levels[sample[a].Id], levels[sample[b].Id]);
            });
            context.Progress($"{n * (n - 1) / 2:N0} pairs compared");

            // AUCs of one draw of genomes (indexes into sample; a genome drawn twice is not paired with itself): [rank, stratum, method].
            (double Auc, int Positives, int Negatives)[,,] Aucs(int[] drawn)
            {
                var result = new (double, int, int)[ranks.Length, strata.Length, 2];
                for (int k = 0; k < ranks.Length; k++)
                    for (int s = 0; s < strata.Length; s++)
                    {
                        var zPairs = new List<(double Score, bool Positive)>();
                        var gcPairs = new List<(double Score, bool Positive)>();
                        for (int i = 0; i < drawn.Length; i++)
                            for (int j = i + 1; j < drawn.Length; j++)
                            {
                                int a = drawn[i], b = drawn[j];
                                if (a == b) continue;
                                double dGc = Math.Abs(sample[a].Gc - sample[b].Gc);
                                if (dGc < strata[s].Min || dGc >= strata[s].Max) continue;
                                bool positive = ranks[k].Of(sample[a]) == ranks[k].Of(sample[b]);
                                zPairs.Add((zCurve[a, b], positive));
                                gcPairs.Add((-dGc, positive));
                            }
                        result[k, s, 0] = Auc(zPairs);
                        result[k, s, 1] = Auc(gcPairs);
                    }
                return result;
            }

            var observed = Aucs(Enumerable.Range(0, n).ToArray());
            var draws = new (double Auc, int Positives, int Negatives)[resamples][,,];
            Parallel.For(0, resamples, parallel, r =>
            {
                var local = new Random(seed + r);
                draws[r] = Aucs(Enumerable.Range(0, n).Select(_ => local.Next(n)).ToArray());
            });

            var rows = new List<Row>();
            for (int k = 0; k < ranks.Length; k++)
                for (int s = 0; s < strata.Length; s++)
                {
                    var (zAuc, positives, negatives) = observed[k, s, 0];
                    double gcAuc = observed[k, s, 1].Auc;
                    var zDraws = draws.Select(d => d[k, s, 0].Auc).ToList();
                    var gcDraws = draws.Select(d => d[k, s, 1].Auc).ToList();
                    var differences = draws.Select(d => d[k, s, 0].Auc - d[k, s, 1].Auc).ToList();
                    rows.Add(new Row(windowSize, ranks[k].Name, strata[s].Name, positives, negatives,
                        zAuc, Percentile(zDraws, 0.025), Percentile(zDraws, 0.975),
                        gcAuc, Percentile(gcDraws, 0.025), Percentile(gcDraws, 0.975),
                        zAuc - gcAuc, Percentile(differences, 0.025), Percentile(differences, 0.975)));
                }
            return rows;
        }

        private static double Percentile(IEnumerable<double> values, double p)
        {
            var sorted = values.Where(v => !double.IsNaN(v)).Order().ToArray();
            return sorted.Length == 0 ? double.NaN : sorted[(int)Math.Clamp(Math.Round(p * (sorted.Length - 1)), 0, sorted.Length - 1)];
        }

        // Rank-sum AUC with average ranks for ties, as RetrievalClassificationBenchmark.ComputeAuc; NaN when a class is empty.
        public static (double Auc, int Positives, int Negatives) Auc(List<(double Score, bool Positive)> pairs)
        {
            int positives = pairs.Count(p => p.Positive), negatives = pairs.Count - positives;
            if (positives == 0 || negatives == 0) return (double.NaN, positives, negatives);
            var sorted = pairs.OrderBy(p => p.Score).ToArray();
            double rankSum = 0;
            for (int i = 0; i < sorted.Length;)
            {
                int j = i;
                while (j + 1 < sorted.Length && sorted[j + 1].Score == sorted[i].Score) j++;
                double averageRank = (i + j + 2) / 2.0;
                for (int t = i; t <= j; t++) if (sorted[t].Positive) rankSum += averageRank;
                i = j + 1;
            }
            return ((rankSum - positives * (positives + 1) / 2.0) / ((double)positives * negatives), positives, negatives);
        }
    }
}
