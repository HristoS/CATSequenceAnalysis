using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// e2-phage-host (P9 S9.2, E2), predictor A only, ported from scripts/platform-checks/p9_e2_phage_host.cs unchanged in what it computes;
    /// predictor B (localization, --localize) stays in the script. The set: viruses whose BV-BRC host_name starts with a stored bacterial or
    /// archaeal genus (provider lineage). Cellular taxa are ranked by the distance from the phage's signature (log rho*, Euclidean) to the
    /// nearest genome of each taxon, at genus, family and order, over every stored taxon (and, at genus, among the host genera only); GC
    /// alone (|ΔGC|), the most common host and a shuffle of host labels across phages are the baselines. Pass, as fixed before the article
    /// run: top-1 genus accuracy minus GC's has a bootstrap 95% lower bound above 0, overall and macro (mean over host genera with at least 5
    /// phages), on the whole set; the set without records under 5 kbp is a check.
    /// </summary>
    public sealed class PhageHost : IExperiment
    {
        public const string ExperimentName = "e2-phage-host";
        public const string FileName = "P9_E2_PhageHost.csv";

        private const int minHostPhages = 5, shuffles = 200, bootstraps = 1000, seed = 20261104;
        private const long minLengthCheck = 5_000;
        private const string header = "Set, Pool, Rank, Predictor, Phages, Top1, Top5, MacroTop1, MacroTop5, MacroGenera, DiffLow, DiffHigh, MacroDiffLow, MacroDiffHigh, Pass";

        private static readonly string[] ranks = ["genus", "family", "order"];

        public string Name => ExperimentName;

        private sealed record Phage(Corpus.Genome Genome, string?[] Host);

        private sealed record Ranking(int[][][] RankOf, int[] PoolRank);

        public async Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken)
        {
            var snapshot = await Snapshots.LiveAsync(dataSource, cancellationToken);
            var corpus = (await Corpus.LoadAsync(dataSource, cancellationToken)).Where(g => !g.Excluded).ToList();
            var cellular = corpus.Where(g => g.Group == SignatureIndex.Cellular && g.Provider.Genus is not null).ToList();
            var taxa = ranks.Select((_, r) => new TaxonIndex(cellular.Select(g => Taxon(g.Provider, r)).ToArray())).ToArray();

            // A host genus's family and order: the majority over its stored genomes.
            var hostLineage = cellular.GroupBy(g => g.Provider.Genus!).ToDictionary(g => g.Key, g => new[]
            {
                g.Key,
                g.Select(x => x.Provider.Family).Where(f => f is not null).GroupBy(f => f).OrderByDescending(f => f.Count()).ThenBy(f => f.Key, StringComparer.Ordinal).FirstOrDefault()?.Key,
                g.Select(x => x.Provider.Order).Where(o => o is not null).GroupBy(o => o).OrderByDescending(o => o.Count()).ThenBy(o => o.Key, StringComparer.Ordinal).FirstOrDefault()?.Key,
            });
            var phages = corpus.Where(g => g.Group == SignatureIndex.Viral && HostGenus(g) is { } h && hostLineage.ContainsKey(h))
                .Select(g => new Phage(g, hostLineage[HostGenus(g)!])).ToList();
            var hostGenera = phages.Select(p => taxa[0].Of(p.Host[0])).ToHashSet();
            context.Progress($"{phages.Count} phages, {hostGenera.Count} host genera; {cellular.Count} cellular genomes in {taxa[0].Count} genera");

            var parallel = new ParallelOptions { MaxDegreeOfParallelism = context.Parallelism, CancellationToken = cancellationToken };
            var signature = Rankings(phages, cellular, taxa, hostGenera, (p, g) => Squared(p.Genome.Signature, g.Signature), parallel);
            var gc = Rankings(phages, cellular, taxa, hostGenera, (p, g) => Math.Abs(p.Genome.Gc - g.Gc), parallel);
            context.Progress("rankings done");

            var truth = ranks.Select((_, r) => phages.Select(p => taxa[r].Of(p.Host[r])).ToArray()).ToArray();
            var hostOf = phages.Select(p => p.Host[0]!).ToArray();
            var macroGenera = hostOf.GroupBy(h => h).Where(h => h.Count() >= minHostPhages).Select(h => h.Key).ToHashSet();
            var rows = new List<string>();
            var top1 = new Dictionary<string, double>();
            bool pass = false;

            foreach (var (set, members) in new[]
            {
                ("all", Enumerable.Range(0, phages.Count).ToArray()),
                ($"atLeast{minLengthCheck / 1000}kbp", Enumerable.Range(0, phages.Count).Where(i => phages[i].Genome.LengthBp >= minLengthCheck).ToArray()),
            })
            {
                void Add(string pool, string rank, string predictor, int[] scored, Func<int, int> rankOfTruth)
                {
                    var hit1 = scored.Select(i => rankOfTruth(i) <= 1).ToArray();
                    var hit5 = scored.Select(i => rankOfTruth(i) <= 5).ToArray();
                    var macroScored = scored.Where(i => macroGenera.Contains(hostOf[i])).ToArray();
                    double Macro(Func<int, bool> hit) => macroScored.GroupBy(i => hostOf[i]).Average(g => g.Average(i => hit(i) ? 1.0 : 0));
                    var position = scored.Select((i, j) => (i, j)).ToDictionary(x => x.i, x => x.j);
                    double share = hit1.Average(h => h ? 1.0 : 0);
                    if (set == "all" && pool == "all" && rank == "genus") top1[predictor] = share;
                    rows.Add(Row(set, pool, rank, predictor, scored.Length, share, hit5.Average(h => h ? 1.0 : 0),
                        Macro(i => hit1[position[i]]), Macro(i => hit5[position[i]]), macroScored.GroupBy(i => hostOf[i]).Count()));
                }

                for (int r = 0; r < ranks.Length; r++)
                {
                    var scored = members.Where(i => truth[r][i] >= 0).ToArray();
                    Add("all", ranks[r], "signature", scored, i => signature.RankOf[r][i][truth[r][i]]);
                    Add("all", ranks[r], "gc", scored, i => gc.RankOf[r][i][truth[r][i]]);
                    var common = scored.Select(i => truth[r][i]).GroupBy(t => t).OrderByDescending(t => t.Count()).ThenBy(t => t.Key).Select(t => t.Key).ToArray();
                    Add("all", ranks[r], "mostCommon", scored, i => Array.IndexOf(common, truth[r][i]) is var at && at >= 0 ? at + 1 : int.MaxValue);

                    // Host labels shuffled across phages: what the signature's ranking scores against a label it has nothing to do with.
                    var rng = new Random(seed + r);
                    var shuffledTop1 = new double[shuffles];
                    var shuffledTop5 = new double[shuffles];
                    var labels = scored.Select(i => truth[r][i]).ToArray();
                    for (int s = 0; s < shuffles; s++)
                    {
                        cancellationToken.ThrowIfCancellationRequested();
                        rng.Shuffle(labels);
                        shuffledTop1[s] = scored.Select((i, j) => signature.RankOf[r][i][labels[j]] <= 1 ? 1.0 : 0).Average();
                        shuffledTop5[s] = scored.Select((i, j) => signature.RankOf[r][i][labels[j]] <= 5 ? 1.0 : 0).Average();
                    }
                    Array.Sort(shuffledTop1);
                    Array.Sort(shuffledTop5);
                    rows.Add(Row(set, "all", ranks[r], "shuffledMean", scored.Length, shuffledTop1.Average(), shuffledTop5.Average(), null, null, null));
                    rows.Add(Row(set, "all", ranks[r], "shuffled97.5", scored.Length, shuffledTop1[(int)(0.975 * shuffles) - 1],
                        shuffledTop5[(int)(0.975 * shuffles) - 1], null, null, null));
                }
                Add("hostGenera", "genus", "signature", members, i => signature.PoolRank[i]);
                Add("hostGenera", "genus", "gc", members, i => gc.PoolRank[i]);

                var (overall, overallLow, overallHigh, macro, macroLow, macroHigh) = BootstrapDifference(members,
                    members.Select(i => signature.RankOf[0][i][truth[0][i]] == 1).ToArray(), members.Select(i => gc.RankOf[0][i][truth[0][i]] == 1).ToArray(),
                    members.Select(i => hostOf[i]).ToArray(), macroGenera);
                bool setPass = overallLow > 0 && macroLow > 0;
                if (set == "all") pass = setPass;
                rows.Add(string.Create(CultureInfo.InvariantCulture,
                    $"{set}, all, genus, signatureMinusGc, {members.Length}, {overall:F4}, , {macro:F4}, , {macroGenera.Count}, {overallLow:F4}, {overallHigh:F4}, {macroLow:F4}, {macroHigh:F4}, {setPass}"));
            }

            var csv = new StringBuilder(header).Append('\n');
            foreach (var row in rows) csv.Append(row).Append('\n');
            var summary = new JsonObject
            {
                ["keyNumbers"] = new JsonArray(
                    KeyNumber("top-1 host genus, signature", top1["signature"]),
                    KeyNumber("GC", top1["gc"]),
                    KeyNumber("always the commonest host", top1["mostCommon"])),
                ["analysed"] = phages.Count,
                ["analysedNote"] = $"phages with a stored host genus, {hostGenera.Count} host genera, {taxa[0].Count} stored genera; predictor A only",
            };
            return new ExperimentResult(pass ? ExperimentStore.Pass : ExperimentStore.Fail, snapshot, "SequencePlatform.Experiments/PhageHost.cs",
                summary.ToJsonString(), [new ExperimentFile(FileName, Encoding.UTF8.GetBytes(csv.ToString()))], []);
        }

        private static JsonObject KeyNumber(string label, double value) => new() { ["label"] = label, ["value"] = Math.Round(value, 4) };

        public static string? HostGenus(Corpus.Genome g) => g.Field("host_name")?.Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault();

        private static string? Taxon(Lineage lineage, int rank) => rank switch { 0 => lineage.Genus, 1 => lineage.Family, _ => lineage.Order };

        private static double Squared(double[] a, double[] b)
        {
            double sum = 0;
            for (int f = 0; f < a.Length; f++) { double diff = a[f] - b[f]; sum += diff * diff; }
            return sum;
        }

        private static string Row(string set, string pool, string rank, string predictor, int phages, double top1, double top5, double? macro1, double? macro5,
            int? macroGenera) =>
            string.Create(CultureInfo.InvariantCulture,
                $"{set}, {pool}, {rank}, {predictor}, {phages}, {top1:F4}, {top5:F4}, {macro1:F4}, {macro5:F4}, {macroGenera}, , , , , ");

        // Per phage and rank: every taxon's 1-based position when taxa are ordered by their nearest genome (ties in order of taxon index),
        // and the true genus's position among the host genera only.
        private static Ranking Rankings(List<Phage> phages, List<Corpus.Genome> cellular, TaxonIndex[] taxa, HashSet<int> hostGenera,
            Func<Phage, Corpus.Genome, double> distance, ParallelOptions parallel)
        {
            var rankOf = taxa.Select(_ => new int[phages.Count][]).ToArray();
            var poolRank = new int[phages.Count];
            Parallel.For(0, phages.Count, parallel, i =>
            {
                var nearest = taxa.Select(t => Enumerable.Repeat(double.PositiveInfinity, t.Count).ToArray()).ToArray();
                for (int g = 0; g < cellular.Count; g++)
                {
                    double d = distance(phages[i], cellular[g]);
                    for (int r = 0; r < taxa.Length; r++)
                        if (taxa[r].OfGenome[g] is var t and >= 0 && d < nearest[r][t]) nearest[r][t] = d;
                }
                for (int r = 0; r < taxa.Length; r++)
                {
                    var order = Enumerable.Range(0, taxa[r].Count).OrderBy(t => nearest[r][t]).ThenBy(t => t).ToArray();
                    var positions = new int[order.Length];
                    for (int j = 0; j < order.Length; j++) positions[order[j]] = j + 1;
                    rankOf[r][i] = positions;
                }
                int host = taxa[0].Of(phages[i].Host[0]);
                poolRank[i] = 1 + hostGenera.Count(u => rankOf[0][i][u] < rankOf[0][i][host]);
            });
            return new Ranking(rankOf, poolRank);
        }

        // Signature minus GC in top-1 genus accuracy, overall and macro, with a percentile bootstrap over phages.
        private static (double Overall, double OverallLow, double OverallHigh, double Macro, double MacroLow, double MacroHigh) BootstrapDifference(int[] members,
            bool[] a, bool[] b, string[] host, HashSet<string> macroGenera)
        {
            var draw = new int[members.Length];
            (double Overall, double Macro) Difference(int[] picks)
            {
                double overall = picks.Average(j => (a[j] ? 1.0 : 0) - (b[j] ? 1.0 : 0));
                double macro = picks.Where(j => macroGenera.Contains(host[j])).GroupBy(j => host[j])
                    .Average(g => g.Average(j => (a[j] ? 1.0 : 0) - (b[j] ? 1.0 : 0)));
                return (overall, macro);
            }
            var (overallAll, macroAll) = Difference(Enumerable.Range(0, members.Length).ToArray());
            var rng = new Random(seed);
            var overalls = new double[bootstraps];
            var macros = new double[bootstraps];
            for (int s = 0; s < bootstraps; s++)
            {
                for (int j = 0; j < draw.Length; j++) draw[j] = rng.Next(members.Length);
                (overalls[s], macros[s]) = Difference(draw);
            }
            Array.Sort(overalls);
            Array.Sort(macros);
            int low = (int)(0.025 * bootstraps), high = (int)(0.975 * bootstraps) - 1;
            return (overallAll, overalls[low], overalls[high], macroAll, macros[low], macros[high]);
        }

        private sealed class TaxonIndex
        {
            private readonly Dictionary<string, int> index = new(StringComparer.Ordinal);
            private readonly List<string> names = [];

            public int[] OfGenome { get; }
            public int Count => names.Count;

            public TaxonIndex(string?[] genomeTaxa)
            {
                OfGenome = genomeTaxa.Select(t =>
                {
                    if (t is null) return -1;
                    if (!index.TryGetValue(t, out var i)) { i = names.Count; index[t] = i; names.Add(t); }
                    return i;
                }).ToArray();
            }

            public int Of(string? name) => name is not null && index.TryGetValue(name, out var i) ? i : -1;
        }
    }
}
