using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// e5-habitat (P9 S9.5, E5), ported from scripts/platform-checks/p9_e5_habitat.cs unchanged in what it computes: within a genus, do a
    /// genome's nearest neighbours by signature share its habitat more than a shuffle of habitat labels inside that genus predicts?
    /// Shuffling only inside the genus keeps the genus's ancestry in the null. Genomes: <see cref="Corpus"/>'s bacteria and archaea, partial
    /// ones left out. Habitat from BV-BRC's
    /// isolation_source by <see cref="Classify"/>. Per genus (provider lineage) with at least 10 labelled genomes and 2 habitats: each
    /// genome's 5 nearest labelled genomes of the genus (log rho*, Euclidean), the mean share with the same habitat, against 1,000 shuffles
    /// of the genus's labels; the same with |ΔGC|, and the same per species, as context. Pass: more genera above the shuffle than below,
    /// one-sided sign test p &lt; 0.05. Stops (fail) when fewer than 20 genera qualify.
    /// </summary>
    public sealed class Habitat : IExperiment
    {
        public const string ExperimentName = "e5-habitat";
        public const string HabitatFileName = "P9_E5_Habitat.csv";
        public const string RulesFileName = "P9_E5_Rules.csv";

        private const int k = 5, minGenomes = 10, minGenera = 20, shuffles = 1000, seed = 20261104;

        private static readonly (string Habitat, Regex Pattern)[] rules =
        [
            ("none", new(@"^(missing|not available.*|not applicable|not collected|not provided|unknown|n/?a|-|none|environmental|environment|water|r2a|envo:\d+|other|sample)$")),
            ("clinical", new(@"blood|sputum|urine|wound|patient|hospital|clinical|abscess|cerebrospinal|\bcsf\b|bronch|trache|catheter|\bpus\b|swab|peritoneal|pleural|aspirate|biopsy")),
            ("food", new(@"\bfood|milk|cheese|dairy|meat|sausage|salami|kimchi|fermented|yog(h)?urt|kefir|\bwine|\bbeer|vinegar|sourdough|pickle|sauerkraut")),
            ("host-associated", new(@"^[a-z ]+:|fec(es|al)|faec|stool|\bgut\b|intestin|ca?ecal|cecum|caecum|rumen|oral|skin|nasal|saliva|nodule|\broot|\bleaf|phyllosphere|endophyt|insect|beetle|\btick|mosquito|larva|\bhost\b")),
            ("marine", new(@"marine|\bsea\b|seawater|ocean|coastal|hydrothermal|cold seep|estuar|mangrove|coral|deep-sea")),
            ("freshwater", new(@"freshwater|fresh water|\blake|\briver|\bpond|stream|groundwater|drinking water")),
            ("soil", new(@"soil|rhizosphere|compost|sediment|\bmud\b|permafrost|desert")),
        ];
        public string Name => ExperimentName;

        public sealed record Genome(long Id, double Gc, double[] Signature, string? Genus, string? Species, string? IsolationSource);

        /// <summary>A group's scores: the mean same-habitat share of the k nearest, its mean over the shuffles, and the one-sided p.</summary>
        public sealed record Score(double Observed, double Shuffled, double P)
        {
            public double Excess => Observed - Shuffled;
        }

        public async Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken)
        {
            var snapshot = await Snapshots.LiveAsync(dataSource, cancellationToken);
            var genomes = (await Corpus.LoadAsync(dataSource, cancellationToken)).Where(g => !g.Excluded && g.Group == SignatureIndex.Cellular)
                .Select(g => new Genome(g.Id, g.Gc, g.Signature, g.Provider.Genus, g.Provider.Species, g.Field("isolation_source"))).ToList();

            var labelled = new List<(Genome Genome, string Habitat)>();
            var coverage = new Dictionary<string, int>();
            int withSource = 0;
            foreach (var g in genomes)
            {
                if (g.IsolationSource is not { } source) continue;
                withSource++;
                var habitat = Classify(source);
                coverage[habitat] = coverage.GetValueOrDefault(habitat) + 1;
                if (habitat is not ("none" or "other")) labelled.Add((g, habitat));
            }
            var rulesCsv = new StringBuilder("Habitat, Genomes, ShareOfWithSource\n");
            foreach (var (habitat, count) in coverage.OrderByDescending(c => c.Value))
                rulesCsv.Append(string.Create(CultureInfo.InvariantCulture, $"{habitat}, {count}, {(double)count / withSource:F4}")).Append('\n');
            context.Progress($"{genomes.Count} cellular genomes, {withSource} with an isolation source, {labelled.Count} labelled");

            var rows = new List<string>();
            var summary = new JsonObject();
            string verdict = ExperimentStore.Fail;
            foreach (var rank in new[] { "genus", "species" })
            {
                var groups = labelled.Where(l => (rank == "genus" ? l.Genome.Genus : l.Genome.Species) is not null)
                    .GroupBy(l => rank == "genus" ? l.Genome.Genus! : l.Genome.Species!)
                    .Where(g => g.Count() >= minGenomes && g.Select(l => l.Habitat).Distinct().Count() >= 2)
                    .OrderBy(g => g.Key, StringComparer.Ordinal).ToList();
                context.Progress($"{rank}: {groups.Count} qualify");
                if (rank == "genus" && groups.Count < minGenera)
                {
                    summary["stopped"] = $"{groups.Count} genera qualify, under {minGenera}";
                    break;
                }
                if (groups.Count == 0) continue;   // species: the catalogue holds about one reference genome per bacterial species

                var results = groups.Select((group, index) => (group, index)).AsParallel().WithDegreeOfParallelism(Math.Max(1, context.Parallelism))
                    .WithCancellation(cancellationToken).Select(item =>
                    {
                        var (group, index) = item;
                        var members = group.Select(m => m.Genome).ToArray();
                        var labels = group.Select(m => m.Habitat).ToArray();
                        var signature = ScoreGroup(Neighbours(members, (a, b) => Squared(a.Signature, b.Signature)), labels, seed + 2 * index);
                        var gc = ScoreGroup(Neighbours(members, (a, b) => Math.Abs(a.Gc - b.Gc)), labels, seed + 2 * index + 1);
                        return (Name: group.Key, Genomes: members.Length, Habitats: labels.Distinct().Count(),
                            Majority: labels.GroupBy(l => l).Max(l => l.Count()) / (double)members.Length, Signature: signature, Gc: gc);
                    }).OrderBy(r => r.Name, StringComparer.Ordinal).ToList();

                foreach (var r in results)
                    rows.Add(string.Create(CultureInfo.InvariantCulture,
                        $"{rank}, \"{r.Name}\", {r.Genomes}, {r.Habitats}, {r.Majority:F3}, {r.Signature.Observed:F4}, {r.Signature.Shuffled:F4}, " +
                        $"{r.Signature.Excess:F4}, {r.Signature.P:F4}, {r.Gc.Observed:F4}, {r.Gc.Shuffled:F4}, {r.Gc.Excess:F4}, , "));
                int positive = results.Count(r => r.Signature.Excess > 1e-9), negative = results.Count(r => r.Signature.Excess < -1e-9);
                int better = results.Count(r => r.Signature.Excess > r.Gc.Excess + 1e-9), worse = results.Count(r => r.Signature.Excess < r.Gc.Excess - 1e-9);
                double p = SignTest(positive, negative), versusGc = SignTest(better, worse);
                bool rankPass = positive > negative && p < 0.05;
                if (rank == "genus")
                {
                    verdict = rankPass ? ExperimentStore.Pass : ExperimentStore.Fail;
                    summary["keyNumbers"] = new JsonArray(
                        KeyNumber("genera above the shuffle", positive),
                        KeyNumber("genera below the shuffle", negative),
                        KeyNumber("genera where the signature beats GC", better));
                    summary["analysed"] = results.Sum(r => r.Genomes);
                    summary["analysedNote"] = string.Create(CultureInfo.InvariantCulture,
                        $"labelled genomes in {results.Count} genera; sign test p {p:G3}, against GC {better} / {worse} (p {versusGc:G3})");
                }
                rows.Add(string.Create(CultureInfo.InvariantCulture,
                    $"{rank}, all, {results.Sum(r => r.Genomes)}, {results.Count}, , {results.Average(r => r.Signature.Observed):F4}, " +
                    $"{results.Average(r => r.Signature.Shuffled):F4}, {results.Average(r => r.Signature.Excess):F4}, {p:G3}, " +
                    $"{results.Average(r => r.Gc.Observed):F4}, {results.Average(r => r.Gc.Shuffled):F4}, {results.Average(r => r.Gc.Excess):F4}, " +
                    $"{versusGc:G3}, {(rank == "genus" ? rankPass.ToString() : "")}"));
            }

            var habitatCsv = new StringBuilder("Rank, Taxon, Genomes, Habitats, MajorityShare, SignatureShare, SignatureShuffled, SignatureExcess, " +
                "SignatureP, GcShare, GcShuffled, GcExcess, SignVsGcP, Pass\n");
            foreach (var row in rows) habitatCsv.Append(row).Append('\n');
            return new ExperimentResult(verdict, snapshot, "SequencePlatform.Experiments/Habitat.cs", summary.ToJsonString(),
                [new ExperimentFile(HabitatFileName, Encoding.UTF8.GetBytes(habitatCsv.ToString())),
                 new ExperimentFile(RulesFileName, Encoding.UTF8.GetBytes(rulesCsv.ToString()))], []);
        }

        /// <summary>The habitat of an isolation_source by the keyword rules, first match wins; "none" for vague values, "other" if no rule matches.</summary>
        public static string Classify(string isolationSource)
        {
            var text = isolationSource.Trim().ToLowerInvariant();
            return rules.FirstOrDefault(r => r.Pattern.IsMatch(text)).Habitat ?? "other";
        }

        private static JsonObject KeyNumber(string label, double value) => new() { ["label"] = label, ["value"] = value };

        private static double Squared(double[] a, double[] b)
        {
            double sum = 0;
            for (int f = 0; f < a.Length; f++) { double d = a[f] - b[f]; sum += d * d; }
            return sum;
        }

        /// <summary>The k nearest of each member within its group, nearest first, itself excluded; ties by index.</summary>
        public static int[][] Neighbours(Genome[] members, Func<Genome, Genome, double> distance) =>
            members.Select((m, i) => Enumerable.Range(0, members.Length).Where(j => j != i)
                .OrderBy(j => distance(m, members[j])).ThenBy(j => j).Take(k).ToArray()).ToArray();

        public static Score ScoreGroup(int[][] neighbours, string[] labels, int seed)
        {
            double Share(string[] l) => neighbours.Select((n, i) => n.Count(j => l[j] == l[i]) / (double)n.Length).Average();
            double observed = Share(labels);
            var rng = new Random(seed);
            var shuffled = labels.ToArray();
            double sum = 0;
            int atLeast = 0;
            for (int s = 0; s < shuffles; s++)
            {
                rng.Shuffle(shuffled);
                double value = Share(shuffled);
                sum += value;
                if (value >= observed - 1e-12) atLeast++;
            }
            return new Score(observed, sum / shuffles, (atLeast + 1.0) / (shuffles + 1));
        }

        /// <summary>One-sided binomial sign test: P(X &gt;= positive) with X ~ Binomial(positive + negative, 1/2).</summary>
        public static double SignTest(int positive, int negative)
        {
            int n = positive + negative;
            if (n == 0) return 1;
            double p = 0;
            for (int x = positive; x <= n; x++) p += Math.Exp(LogChoose(n, x) - n * Math.Log(2));
            return Math.Min(1, p);
        }

        private static double LogChoose(int n, int x) => LogFactorial(n) - LogFactorial(x) - LogFactorial(n - x);

        private static double LogFactorial(int n)
        {
            double sum = 0;
            for (int i = 2; i <= n; i++) sum += Math.Log(i);
            return sum;
        }
    }
}
