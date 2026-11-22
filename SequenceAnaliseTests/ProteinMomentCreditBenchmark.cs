using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    using G9 = ProteinHomologRankingBenchmark;
    using P = ProteinBackTranslationFeasibility;

    /// <summary>
    /// Data-set credits for the gene moments (2026-11-22, after the per-moment check, protein_family_diagnostic.cs --moments at commit f7c796e): instead
    /// of the weighted Euclidean distance, each of the 9 window features (P0, P1, P2 of x, y, z) scores a difference by its log-odds,
    /// log P(d | homolog) / P(d | unrelated), learned from pairs, as BLOSUM scores residue pairs. If both distributions are normal the
    /// log-odds is a weighted squared distance, today's form; a learned table can only beat it where they aren't.
    ///
    /// Pool, homolog key and length rule as ProteinHomologRankingBenchmark (G9): genes called on the 30 bundled bacteria from 50 aa, gene
    /// windows from the DNA rounded to float16 (the stored profile), queries from the protein; a homolog is an alignment at Rost's HSSP
    /// curve + 5 points over 60% of the shorter protein. 20 queries per genome; genomes split by a seeded shuffle into a training and a
    /// test half. Training pairs (training queries only, keyed by alignment, so independent of any Z-curve score): every aligned homolog
    /// as a positive, 200 random unaligned candidates that aren't among its aligned homologs as negatives (about 1% of random genes are
    /// homologs: a small contamination). Differences are in noise units, d * sqrt(W), W the query's window size in codon positions 1 and
    /// 2, since on random DNA every moment's noise has standard deviation 1 / sqrt(W).
    /// Arms on the test half: A today (1/variance weighted L2, ZCurveGeneProfile.Score); T1 1-D credit tables, mean over windows; T1s the
    /// same summed over windows (the plain log-likelihood ratio); T2 2-D tables over (query, gene) values in quantile bins, which learn on
    /// their own whether the same sign or the composition level earns credit beyond the difference; G Gaussian log-odds per feature.
    /// Rule fixed before the run: the best of T1, T1s, T2 beats A by >= 3 points of homolog top-1 with the paired bootstrap CI above 0 ->
    /// a card for a database check on G4's sets; otherwise Euclidean stays. Writes ProteinMomentCredits.csv and
    /// ProteinMomentCredits_Tables.csv.
    /// </summary>
    public class ProteinMomentCreditBenchmark
    {
        private const int minResidues = 50, queriesPerGenome = 20, prefilterTop = 50, topK = 10, negativesPerQuery = 200,
            bootstrapRounds = 2000, seed = 20261122, zBins = 64, quantileBins = 12;
        private const double zLimit = 8, minCoverage = 0.6, hsspMargin = 5, histogramPseudocount = 1, cellPseudocount = 0.5;
        private const int F = ZCurveGeneProfile.FeatureCount;
        private static readonly string[] arms = { "A today", "T1 credits", "T1s credits summed", "T2 credits 2-D", "G Gaussian" };
        private static readonly string[] featureNames = { "P0 x", "P0 y", "P0 z", "P1 x", "P1 y", "P1 z", "P2 x", "P2 y", "P2 z" };
        private static readonly string[] bandNames = { "50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa" };
        private static readonly double[] weights = ZCurveGeneProfile.FeatureWeights.ToArray();

        private sealed record Gene(int Genome, string Protein, double[] Windows, double[] QueryWindows, int[] Trimers);

        private sealed record Outcome(int Band, bool Reachable, double BestIdentity, bool[] Top1, bool[] Top10);

        // Per class (0 homolog, 1 unrelated) and feature: the query's and the gene's values in noise units, window by window.
        private sealed class Observations
        {
            public readonly List<double>[,] Query = New(), Gene = New();
            public readonly long[] Pairs = new long[2];

            private static List<double>[,] New()
            {
                var lists = new List<double>[2, F];
                for (int c = 0; c < 2; c++) for (int f = 0; f < F; f++) lists[c, f] = new List<double>();
                return lists;
            }
        }

        private sealed class Tables
        {
            public readonly double[][] Credit1 = new double[F][], Credit2 = new double[F][], Edges = new double[F][];
            public readonly (double Mean, double Variance)[,] Normal = new (double, double)[2, F];
        }

        // The pool, the genome split and the query sets both tests use, with the candidate rule and the alignment key.
        private sealed class Corpus(List<Gene> pool, string[] genusOf, List<int> trainQueries, List<int> testQueries)
        {
            public List<Gene> Pool => pool;
            public List<int> TrainQueries => trainQueries;
            public List<int> TestQueries => testQueries;

            public List<int> Candidates(Gene query)
            {
                int length = query.Protein.Length;
                bool longRule = length >= ZCurveGeneProfile.LongResidues;
                var candidates = new List<int>();
                for (int g = 0; g < pool.Count; g++)
                {
                    var gene = pool[g];
                    if (genusOf[gene.Genome] == genusOf[query.Genome]) continue;
                    int residues = gene.Protein.Length;
                    if (longRule ? residues >= ZCurveGeneProfile.LongResidues
                        : residues < ZCurveGeneProfile.LongResidues && Math.Abs(residues - length) <= 0.1 * length) candidates.Add(g);
                }
                return candidates;
            }

            public Dictionary<int, (bool Is, double Identity)> Key(Gene query, IEnumerable<int> toAlign)
            {
                var key = new Dictionary<int, (bool Is, double Identity)>();
                foreach (int c in toAlign)
                {
                    if (key.ContainsKey(c)) continue;
                    var (identity, aligned) = ProteinAlignment.Align(query.Protein, pool[c].Protein);
                    key[c] = (aligned >= minCoverage * Math.Min(query.Protein.Length, pool[c].Protein.Length)
                        && identity * 100 >= G9.Hssp(aligned) + hsspMargin, identity);
                }
                return key;
            }

            public IEnumerable<int> Prefiltered(Gene query, List<int> candidates) =>
                candidates.Select(c => (c, G9.Shared(query.Trimers, pool[c].Trimers))).OrderByDescending(x => x.Item2).Take(prefilterTop).Select(x => x.c);
        }

        private static Corpus LoadCorpus()
        {
            var files = P.Corpus().Where(c => c.Group == "Bacteria").Select(c => c.FileName).ToList();
            var genusOf = files.Select(f => GenomeTaxonomy.Genomes.Single(g => g.FileName == Path.GetFileName(f)).Genus).ToArray();
            var pool = new List<Gene>();
            var byGenome = new List<List<int>>();
            for (int genome = 0; genome < files.Count; genome++)
            {
                var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(files[genome])).Sequence;
                var code = Path.GetFileName(files[genome]).StartsWith("mycoplasma_", StringComparison.Ordinal) ? P.mycoplasmaCode : P.standardCode;
                var indices = new List<int>();
                foreach (var dna in GeneLevelSignatureBenchmark.CallGenes(sequence, code, 3 * minResidues))
                {
                    var protein = P.Translate(dna, code).TrimEnd('*');
                    if (protein.Length < minResidues || protein.Contains('*')) continue;
                    var windows = ZCurveGeneProfile.QueryFromGene(dna, code).Select(v => (double)(Half)v).ToArray();
                    if (genome == 0 && indices.Count < 20) CheckAgainstLibrary(dna, code, windows, pool.Count > 0 ? pool[0].QueryWindows : windows);
                    indices.Add(pool.Count);
                    pool.Add(new Gene(genome, protein, windows, ZCurveGeneProfile.QueryFromProtein(protein, null, code), G9.Trimers(protein)));
                }
                byGenome.Add(indices);
            }
            var rng = new Random(seed);
            var training = Enumerable.Range(0, files.Count).OrderBy(_ => rng.Next()).Take(files.Count / 2).ToHashSet();
            var queries = byGenome.SelectMany(indices => indices.OrderBy(_ => rng.Next()).Take(queriesPerGenome)).ToList();
            var trainQueries = queries.Where(q => training.Contains(pool[q].Genome)).ToList();
            var testQueries = queries.Where(q => !training.Contains(pool[q].Genome)).ToList();
            TestContext.Out.WriteLine($"{pool.Count:N0} genes in {files.Count} genomes; {trainQueries.Count} training and {testQueries.Count} test queries.");
            return new Corpus(pool, genusOf, trainQueries, testQueries);
        }

        // Components of one pair for the combination arms: 0 today's score; per cutoff c (1, 2, 4 noise units) the axis-weighted mean
        // credit max(0, 1 - |d| sqrt(W) / c) of P0 (1-3), of P1 (4-6), of P1 with full credit for the same sign (7-9), of P2 (10-12);
        // 13 the share of P1 and P2 values whose signs agree. Axes are weighted by FeatureWeights within each moment.
        private static readonly double[] cutoffs = { 1, 2, 4 };
        private const int componentCount = 14;
        private static readonly double[] withinMoment = Enumerable.Range(0, F).Select(f => weights[f] / (weights[3 * (f / 3)] + weights[3 * (f / 3) + 1] + weights[3 * (f / 3) + 2])).ToArray();

        private static bool Components(double[] query, double[] gene, double unit, double[] components)
        {
            int count = Math.Min(query.Length, gene.Length) / F;
            if (count == 0) return false;
            Array.Clear(components);
            double today = 0;
            int agree = 0;
            for (int w = 0; w < count; w++)
            {
                double squared = 0;
                for (int f = 0; f < F; f++)
                {
                    double q = query[w * F + f], g = gene[w * F + f], d = q - g;
                    squared += weights[f] * (d * d);
                    int moment = f / 3;
                    bool same = Math.Sign(q) == Math.Sign(g);
                    if (moment > 0 && same) agree++;
                    double z = Math.Abs(d) * unit, share = withinMoment[f];
                    for (int k = 0; k < cutoffs.Length; k++)
                    {
                        double credit = Math.Max(0, 1 - z / cutoffs[k]);
                        if (moment == 0) components[1 + k] += share * credit;
                        else if (moment == 1) { components[4 + k] += share * credit; components[7 + k] += share * (same ? 1 : credit); }
                        else components[10 + k] += share * credit;
                    }
                }
                today += Math.Sqrt(squared);
            }
            for (int i = 1; i < 13; i++) components[i] /= count;
            components[0] = 1 - today / count;
            components[13] = agree / (6.0 * count);
            return true;
        }

        // Every combination arm as coefficients over the components, grouped in families; arm 0 is today's score.
        private static List<(string Family, string Name, double[] Coefficients)> CombinationArms()
        {
            var list = new List<(string Family, string Name, double[] Coefficients)>();
            var today = new double[componentCount];
            today[0] = 1;
            list.Add(("A", "A today", today));
            double[] betas = { 0.02, 0.05, 0.1, 0.2, 0.5, 1 };
            foreach (bool sign in new[] { false, true })
                for (int k = 0; k < cutoffs.Length; k++)
                {
                    int trend = (sign ? 7 : 4) + k;
                    string rule = sign ? "sign" : "distance";
                    for (int a = 0; a <= 4; a++)
                        for (int b = 0; a + b <= 4; b++)
                        {
                            var coefficients = new double[componentCount];
                            (coefficients[1 + k], coefficients[trend], coefficients[10 + k]) = (a / 4.0, b / 4.0, (4 - a - b) / 4.0);
                            list.Add(($"credits ({rule})", $"credits ({rule}) c={cutoffs[k]} lambda={a / 4.0}/{b / 4.0}/{(4 - a - b) / 4.0}", coefficients));
                        }
                    foreach (double beta in betas)
                    {
                        var coefficients = new double[componentCount];
                        coefficients[0] = 1;
                        coefficients[1 + k] = coefficients[trend] = coefficients[10 + k] = beta / 3;
                        list.Add(("today + credits", $"today + {beta} x credits ({rule}) c={cutoffs[k]}", coefficients));
                    }
                }
            foreach (double beta in betas)
            {
                var coefficients = new double[componentCount];
                (coefficients[0], coefficients[13]) = (1, beta);
                list.Add(("today + sign", $"today + {beta} x sign agreement", coefficients));
            }
            return list;
        }

        /// <summary>
        /// The author's hand-made per-moment credits, combined (2026-11-22, after LogOddsCreditsAgainstEuclidean). Per window and axis,
        /// mean and curvature earn max(0, 1 - |d| sqrt(W) / c); the trend the same, or full credit when the signs agree. The three moment
        /// credits are combined with weights lambda (simplex in steps of 0.25), alone or as a bonus on today's score (today + beta x the
        /// equal-weight credit), and sign agreement alone as a bonus (today + beta x sign share). c, lambda, beta and the trend rule are
        /// chosen per family on the training half by homolog top-1; each family's choice is judged once on the test half. Same corpus,
        /// split and key as LogOddsCreditsAgainstEuclidean. Rule fixed before the run: a family's choice beats today by >= 3 points of
        /// homolog top-1 on the test half with the paired bootstrap CI above 0 -> a card for a database check on G4's sets; otherwise
        /// Euclidean stays. Writes ProteinMomentCombinations.csv.
        /// </summary>
        [Test]
        public void HandCreditCombinations()
        {
            var corpus = LoadCorpus();
            var pool = corpus.Pool;
            var arms = CombinationArms();
            int R = arms.Count;

            // Every arm's top hit per query (and the top 10 of a few arms), judged by the alignment key.
            (bool Reachable, int Band, double BestIdentity, bool[] Top1, bool[] Top10)[] Evaluate(List<int> queries, int[] topTenArms)
            {
                var results = new (bool, int, double, bool[], bool[])[queries.Count];
                Parallel.For(0, queries.Count, n =>
                {
                    var query = pool[queries[n]];
                    var candidates = corpus.Candidates(query);
                    double unit = Math.Sqrt(ZCurveGeneProfile.WindowOf(query.Protein.Length));
                    var components = new double[componentCount];
                    var best = Enumerable.Repeat(double.NegativeInfinity, R).ToArray();
                    var bestGene = Enumerable.Repeat(-1, R).ToArray();
                    var tops = topTenArms.Select(_ => new PriorityQueue<int, double>()).ToArray();
                    foreach (int c in candidates)
                    {
                        if (!Components(query.QueryWindows, pool[c].Windows, unit, components)) continue;
                        for (int a = 0; a < R; a++)
                        {
                            var coefficients = arms[a].Coefficients;
                            double score = 0;
                            for (int i = 0; i < componentCount; i++) score += coefficients[i] * components[i];
                            if (score > best[a]) (best[a], bestGene[a]) = (score, c);
                        }
                        for (int t = 0; t < topTenArms.Length; t++)
                        {
                            double score = 0;
                            var coefficients = arms[topTenArms[t]].Coefficients;
                            for (int i = 0; i < componentCount; i++) score += coefficients[i] * components[i];
                            if (tops[t].Count < topK) tops[t].Enqueue(c, score);
                            else if (tops[t].TryPeek(out _, out double lowest) && score > lowest) tops[t].EnqueueDequeue(c, score);
                        }
                    }
                    var ranked = tops.Select(t => t.UnorderedItems.OrderByDescending(e => e.Priority).Select(e => e.Element).ToArray()).ToArray();
                    var key = corpus.Key(query, corpus.Prefiltered(query, candidates).Concat(bestGene.Where(g => g >= 0)).Concat(ranked.SelectMany(r => r)));
                    results[n] = (key.Values.Any(h => h.Is), BandOf(query.Protein.Length),
                        key.Values.Where(h => h.Is).Select(h => h.Identity).DefaultIfEmpty(double.NaN).Max(),
                        bestGene.Select(g => g >= 0 && key[g].Is).ToArray(), ranked.Select(r => r.Any(c => key[c].Is)).ToArray());
                });
                return results;
            }

            var train = Evaluate(corpus.TrainQueries, []).Where(r => r.Reachable).ToList();
            double TrainTop1(int a) => train.Count(r => r.Top1[a]) / (double)train.Count;
            var families = arms.Select(a => a.Family).Distinct().Where(f => f != "A").ToList();
            var chosen = families.Select(family => Enumerable.Range(0, R).Where(a => arms[a].Family == family).MaxBy(TrainTop1)).ToList();
            var judged = new[] { 0 }.Concat(chosen).ToArray();
            var test = Evaluate(corpus.TestQueries, judged).Where(r => r.Reachable).ToList();

            (double Delta, double Low, double High) Paired(int after)
            {
                double Delta(IEnumerable<(bool, int, double, bool[] Top1, bool[])> set) => set.Average(r => (r.Top1[after] ? 1 : 0) - (r.Top1[0] ? 1 : 0));
                var random = new Random(seed);
                var deltas = Enumerable.Range(0, bootstrapRounds).Select(_ => Delta(Enumerable.Range(0, test.Count).Select(_ => test[random.Next(test.Count)]))).Order().ToArray();
                return (Delta(test), deltas[(int)(0.025 * bootstrapRounds)], deltas[(int)(0.975 * bootstrapRounds) - 1]);
            }
            string Points((double Delta, double Low, double High) d) =>
                $"{100 * d.Delta:+0.0;-0.0} points (95% CI {100 * d.Low:+0.0;-0.0} to {100 * d.High:+0.0;-0.0})";

            var lines = new List<string> { $"Judged arms (chosen on the training half: {train.Count} queries with a homolog; test half: {test.Count}),Training top-1,Test top-1,Test top-10,Test vs today" };
            for (int j = 0; j < judged.Length; j++)
            {
                int a = judged[j];
                lines.Add(string.Join(",", arms[a].Name, Format(TrainTop1(a)), Format(test.Count(r => r.Top1[a]) / (double)test.Count),
                    Format(test.Count(r => r.Top10[j]) / (double)test.Count), a == 0 ? "" : Points(Paired(a))));
            }
            int winner = chosen.MaxBy(a => test.Count(r => r.Top1[a]));
            var rule = Paired(winner);
            lines.Add($"Rule (a family's choice >= 3 points over today with CI above 0),{arms[winner].Name}: " +
                (rule.Delta >= 0.03 && rule.Low > 0 ? "go: card for a database check on G4's sets" : "stop: Euclidean stays"));
            lines.Add("");
            lines.Add("Every arm (test top-1 shown for information; the rule uses only the training choice),Family,Training top-1,Test top-1");
            for (int a = 0; a < R; a++)
                lines.Add(string.Join(",", arms[a].Name, arms[a].Family, Format(TrainTop1(a)), Format(test.Count(r => r.Top1[a]) / (double)test.Count)));

            File.WriteAllLines("ProteinMomentCombinations.csv", lines);
            foreach (var line in lines.Take(judged.Length + 2)) TestContext.Out.WriteLine(line);
            Assert.That(test.Count, Is.GreaterThan(50), "too few test queries with a homolog to judge");
        }

        [Test]
        public void LogOddsCreditsAgainstEuclidean()
        {
            var corpus = LoadCorpus();
            var pool = corpus.Pool;
            var (trainQueries, testQueries) = (corpus.TrainQueries, corpus.TestQueries);
            List<int> Candidates(Gene query) => corpus.Candidates(query);
            Dictionary<int, (bool Is, double Identity)> Key(Gene query, IEnumerable<int> toAlign) => corpus.Key(query, toAlign);
            IEnumerable<int> Prefiltered(Gene query, List<int> candidates) => corpus.Prefiltered(query, candidates);

            // Training pairs.
            var perQuery = new Observations[trainQueries.Count];
            Parallel.For(0, trainQueries.Count, n =>
            {
                var query = pool[trainQueries[n]];
                var candidates = Candidates(query);
                var key = Key(query, Prefiltered(query, candidates));
                var observations = new Observations();
                double unit = Math.Sqrt(ZCurveGeneProfile.WindowOf(query.Protein.Length));
                void Add(int pairClass, Gene gene)
                {
                    int count = Math.Min(query.QueryWindows.Length, gene.Windows.Length) / F;
                    if (count == 0) return;
                    observations.Pairs[pairClass]++;
                    for (int i = 0; i < count * F; i++)
                    {
                        observations.Query[pairClass, i % F].Add(query.QueryWindows[i] * unit);
                        observations.Gene[pairClass, i % F].Add(gene.Windows[i] * unit);
                    }
                }
                foreach (var (c, entry) in key) if (entry.Is) Add(0, pool[c]);
                var random = new Random(seed + n);
                int taken = 0;
                for (int attempt = 0; attempt < 4 * negativesPerQuery && taken < negativesPerQuery && candidates.Count > 0; attempt++)
                {
                    int c = candidates[random.Next(candidates.Count)];
                    if (key.TryGetValue(c, out var entry) && entry.Is) continue;
                    Add(1, pool[c]);
                    taken++;
                }
                perQuery[n] = observations;
            });
            var all = new Observations();
            foreach (var o in perQuery)
                for (int c = 0; c < 2; c++)
                {
                    all.Pairs[c] += o.Pairs[c];
                    for (int f = 0; f < F; f++) { all.Query[c, f].AddRange(o.Query[c, f]); all.Gene[c, f].AddRange(o.Gene[c, f]); }
                }
            var tables = Build(all, out var tableLines);
            TestContext.Out.WriteLine($"Training: {all.Pairs[0]:N0} homolog pairs ({all.Query[0, 0].Count:N0} windows), {all.Pairs[1]:N0} unrelated pairs ({all.Query[1, 0].Count:N0} windows).");

            // Test queries: every candidate under every arm; the 3-mer prefilter's top 50 and every arm's top 10 aligned.
            int R = arms.Length;
            var results = new Outcome[testQueries.Count];
            Parallel.For(0, testQueries.Count, n =>
            {
                var query = pool[testQueries[n]];
                var candidates = Candidates(query);
                double unit = Math.Sqrt(ZCurveGeneProfile.WindowOf(query.Protein.Length));
                var tops = Enumerable.Range(0, R).Select(_ => new PriorityQueue<int, double>()).ToArray();
                var scores = new double[R];
                foreach (int c in candidates)
                {
                    if (!Scores(query.QueryWindows, pool[c].Windows, unit, tables, scores)) continue;
                    for (int a = 0; a < R; a++)
                    {
                        var top = tops[a];
                        if (top.Count < topK) top.Enqueue(c, scores[a]);
                        else if (top.TryPeek(out _, out double lowest) && scores[a] > lowest) top.EnqueueDequeue(c, scores[a]);
                    }
                }
                var ranked = tops.Select(t => t.UnorderedItems.OrderByDescending(e => e.Priority).Select(e => e.Element).ToArray()).ToArray();
                var key = Key(query, Prefiltered(query, candidates).Concat(ranked.SelectMany(r => r)));
                results[n] = new Outcome(BandOf(query.Protein.Length), key.Values.Any(h => h.Is),
                    key.Values.Where(h => h.Is).Select(h => h.Identity).DefaultIfEmpty(double.NaN).Max(),
                    ranked.Select(r => r.Length > 0 && key[r[0]].Is).ToArray(), ranked.Select(r => r.Any(c => key[c].Is)).ToArray());
            });

            var reach = results.Where(r => r.Reachable).ToList();
            var lines = new List<string> { "Arm,Band,Reachable,Top1,Top10" };
            foreach (int band in Enumerable.Range(0, bandNames.Length).Prepend(-1))
            {
                var cell = reach.Where(r => band < 0 || r.Band == band).ToList();
                if (cell.Count == 0) continue;
                for (int a = 0; a < R; a++)
                    lines.Add(string.Join(",", arms[a], band < 0 ? "all" : bandNames[band], cell.Count,
                        Format(cell.Count(r => r.Top1[a]) / (double)cell.Count), Format(cell.Count(r => r.Top10[a]) / (double)cell.Count)));
            }
            int bestTable = new[] { 1, 2, 3 }.MaxBy(a => reach.Count(r => r.Top1[a]));
            lines.Add("");
            lines.Add($"Best homolog identity,Queries,A top-1,{arms[bestTable]} top-1");
            foreach (var (low, high) in new[] { (0.0, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01) })
            {
                var bin = reach.Where(r => r.BestIdentity >= low && r.BestIdentity < high).ToList();
                if (bin.Count > 0)
                    lines.Add($"{low:F1}-{Math.Min(high, 1):F1},{bin.Count},{Format(bin.Count(r => r.Top1[0]) / (double)bin.Count)},{Format(bin.Count(r => r.Top1[bestTable]) / (double)bin.Count)}");
            }

            (double Delta, double Low, double High) Paired(int after, int before)
            {
                double Delta(IEnumerable<Outcome> set) => set.Average(r => (r.Top1[after] ? 1 : 0) - (r.Top1[before] ? 1 : 0));
                var random = new Random(seed);
                var deltas = Enumerable.Range(0, bootstrapRounds).Select(_ => Delta(Enumerable.Range(0, reach.Count).Select(_ => reach[random.Next(reach.Count)]))).Order().ToArray();
                return (Delta(reach), deltas[(int)(0.025 * bootstrapRounds)], deltas[(int)(0.975 * bootstrapRounds) - 1]);
            }
            string Points((double Delta, double Low, double High) d) =>
                $"{100 * d.Delta:+0.0;-0.0} points (95% CI {100 * d.Low:+0.0;-0.0} to {100 * d.High:+0.0;-0.0})";
            lines.Add("");
            lines.Add($"Test queries with a homolog,{reach.Count} of {results.Length}");
            lines.Add($"Training pairs,{all.Pairs[0]:N0} homolog / {all.Pairs[1]:N0} unrelated ({all.Query[0, 0].Count:N0} / {all.Query[1, 0].Count:N0} windows)");
            for (int a = 1; a < R; a++) lines.Add($"{arms[a]} vs {arms[0]},{Points(Paired(a, 0))}");
            lines.Add($"{arms[1]} vs {arms[4]} (do the tails matter),{Points(Paired(1, 4))}");
            lines.Add($"{arms[3]} vs {arms[1]} (sign or level beyond the difference),{Points(Paired(3, 1))}");
            var rule = Paired(bestTable, 0);
            lines.Add($"Rule (best of T1 / T1s / T2 vs A: >= 3 points with CI above 0),{arms[bestTable]}: " +
                (rule.Delta >= 0.03 && rule.Low > 0 ? "go: card for a database check on G4's sets" : "stop: Euclidean stays"));
            lines.Add("");
            lines.AddRange(tableLines.Where(l => l.StartsWith("Credit at", StringComparison.Ordinal)));

            File.WriteAllLines("ProteinMomentCredits.csv", lines);
            File.WriteAllLines("ProteinMomentCredits_Tables.csv", tableLines.Where(l => !l.StartsWith("Credit at", StringComparison.Ordinal)));
            foreach (var line in lines) TestContext.Out.WriteLine(line);
            Assert.That(reach.Count, Is.GreaterThan(50), "too few test queries with a homolog to judge");
        }

        // The stored profile and today's score are the library's: windows equal FromCodingSequences' stored halves, and arm A equals Score.
        private static void CheckAgainstLibrary(string dna, string code, double[] windows, double[] query)
        {
            var profile = ZCurveGeneProfile.FromCodingSequences(new[] { (0, 0, dna.Length, false, ZCurveGeneProfile.NoFamily, dna) }, code);
            Assert.That(profile.WindowsOf(0), Is.EqualTo(windows), "stored windows");
            var scores = new double[arms.Length];
            if (Scores(query, windows, 1, null, scores))
                Assert.That(Math.Abs(scores[0] - profile.Score(query, 0)), Is.LessThan(1e-12), "arm A against ZCurveGeneProfile.Score");
        }

        // Every arm's score of one pair over the common windows; false when there's none. Without tables only arm A is filled.
        private static bool Scores(double[] query, double[] gene, double unit, Tables? tables, double[] scores)
        {
            int count = Math.Min(query.Length, gene.Length) / F;
            if (count == 0) return false;
            double today = 0, credit1 = 0, credit2 = 0, normal = 0;
            for (int w = 0; w < count; w++)
            {
                double squared = 0;
                for (int f = 0; f < F; f++)
                {
                    double q = query[w * F + f], g = gene[w * F + f], d = q - g;
                    squared += weights[f] * (d * d);
                    if (tables is null) continue;
                    double z = d * unit;
                    credit1 += tables.Credit1[f][ZBin(z)];
                    credit2 += tables.Credit2[f][QuantileBin(tables.Edges[f], q * unit) * quantileBins + QuantileBin(tables.Edges[f], g * unit)];
                    normal += LogNormal(z, tables.Normal[0, f]) - LogNormal(z, tables.Normal[1, f]);
                }
                today += Math.Sqrt(squared);
            }
            scores[0] = 1 - today / count;
            scores[1] = credit1 / count;
            scores[2] = credit1;
            scores[3] = credit2 / count;
            scores[4] = normal / count;
            return true;
        }

        private static Tables Build(Observations all, out List<string> lines)
        {
            var tables = new Tables();
            lines = new List<string> { "Table,Feature,Bin,From,To,Homolog windows,Unrelated windows,Credit" };
            for (int f = 0; f < F; f++)
            {
                // 1-D: the difference in noise units.
                var counts = new double[2, zBins + 2];
                var moments = new (double Sum, double Squares, long N)[2];
                for (int c = 0; c < 2; c++)
                    for (int i = 0; i < all.Query[c, f].Count; i++)
                    {
                        double z = all.Query[c, f][i] - all.Gene[c, f][i];
                        counts[c, ZBin(z)]++;
                        moments[c] = (moments[c].Sum + z, moments[c].Squares + z * z, moments[c].N + 1);
                    }
                for (int c = 0; c < 2; c++)
                {
                    double mean = moments[c].Sum / moments[c].N;
                    tables.Normal[c, f] = (mean, moments[c].Squares / moments[c].N - mean * mean);
                    Assert.That(Enumerable.Range(0, zBins + 2).Sum(b => counts[c, b]), Is.EqualTo(all.Query[c, f].Count), "1-D table total");
                }
                tables.Credit1[f] = LogOdds(counts, zBins + 2, histogramPseudocount);
                for (int b = 0; b < zBins + 2; b++)
                {
                    double from = b == 0 ? double.NegativeInfinity : -zLimit + (b - 1) * 2 * zLimit / zBins;
                    double to = b == zBins + 1 ? double.PositiveInfinity : -zLimit + b * 2 * zLimit / zBins;
                    lines.Add(string.Join(",", "1-D", featureNames[f], b, Format(from), Format(to), counts[0, b], counts[1, b], Format(tables.Credit1[f][b])));
                }
                int zero = ZBin(0), best = Enumerable.Range(0, zBins + 2).MaxBy(b => tables.Credit1[f][b]);
                lines.Add($"Credit at d = 0 / its maximum,{featureNames[f]}: {Format(tables.Credit1[f][zero])} / {Format(tables.Credit1[f][best])} " +
                    $"at bin {best}; Gaussian homolog sd {Format(Math.Sqrt(tables.Normal[0, f].Variance))} vs unrelated {Format(Math.Sqrt(tables.Normal[1, f].Variance))}");

                // 2-D: the query's and the gene's values, in quantile bins of the pooled values.
                var pooled = new List<double>();
                for (int c = 0; c < 2; c++) { pooled.AddRange(all.Query[c, f]); pooled.AddRange(all.Gene[c, f]); }
                pooled.Sort();
                tables.Edges[f] = Enumerable.Range(1, quantileBins - 1).Select(k => pooled[(int)((long)k * pooled.Count / quantileBins)]).ToArray();
                var cells = new double[2, quantileBins * quantileBins];
                for (int c = 0; c < 2; c++)
                    for (int i = 0; i < all.Query[c, f].Count; i++)
                        cells[c, QuantileBin(tables.Edges[f], all.Query[c, f][i]) * quantileBins + QuantileBin(tables.Edges[f], all.Gene[c, f][i])]++;
                tables.Credit2[f] = LogOdds(cells, quantileBins * quantileBins, cellPseudocount);
                for (int cell = 0; cell < quantileBins * quantileBins; cell++)
                    lines.Add(string.Join(",", "2-D", featureNames[f], $"{cell / quantileBins} {cell % quantileBins}",
                        Format(cell / quantileBins == 0 ? double.NegativeInfinity : tables.Edges[f][cell / quantileBins - 1]),
                        Format(cell % quantileBins == 0 ? double.NegativeInfinity : tables.Edges[f][cell % quantileBins - 1]),
                        cells[0, cell], cells[1, cell], Format(tables.Credit2[f][cell])));
            }
            return tables;
        }

        // log of each bin's homolog share over its unrelated share, with a pseudocount per bin.
        private static double[] LogOdds(double[,] counts, int bins, double pseudocount)
        {
            var totals = new double[2];
            for (int c = 0; c < 2; c++) for (int b = 0; b < bins; b++) totals[c] += counts[c, b] + pseudocount;
            return Enumerable.Range(0, bins).Select(b => Math.Log((counts[0, b] + pseudocount) / totals[0]) - Math.Log((counts[1, b] + pseudocount) / totals[1])).ToArray();
        }

        // Bins of 2 * zLimit / zBins over [-zLimit, zLimit), with a tail bin at each end.
        private static int ZBin(double z) => z < -zLimit ? 0 : z >= zLimit ? zBins + 1 : 1 + Math.Min(zBins - 1, (int)((z + zLimit) * zBins / (2 * zLimit)));

        private static int QuantileBin(double[] edges, double value)
        {
            int bin = 0;
            while (bin < edges.Length && value >= edges[bin]) bin++;
            return bin;
        }

        private static double LogNormal(double z, (double Mean, double Variance) normal) =>
            -0.5 * Math.Log(normal.Variance) - (z - normal.Mean) * (z - normal.Mean) / (2 * normal.Variance);

        private static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;

        private static string Format(double value) => double.IsInfinity(value) ? (value > 0 ? "inf" : "-inf") : value.ToString("F4", CultureInfo.InvariantCulture);
    }
}
