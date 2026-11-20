using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    using P = ProteinBackTranslationFeasibility;

    /// <summary>
    /// Homolog ranking on the bundled bacteria (2026-11-20, after protein_offset_check.cs, G8). G8 found that missed proteins' right-family
    /// genes are distant homologs (median identity 48%), that a shift mostly repairs the moments (P1, P2), and that 45.7% of misses beat
    /// the winning hit once both are re-measured from the protein, but flipping against one hit is not a ranking. This ranks every gene.
    ///
    /// Pool: genes called on the 30 bundled bacteria from 50 aa (GeneLevelSignatureBenchmark.CallGenes, longest ORF from its first ATG, so
    /// starts can be early: offsets found here are an upper bound for annotated genes). Queries: 10 proteins per genome (seeded), searched
    /// against the genes of other genera with ProteinSearch's length rule (500 aa and up: genes 500 aa and up; shorter: within 10%).
    /// Arms, each unweighted and with G4's 1/variance weights:
    ///   A  today: the query's expected windows against the gene's windows from its DNA (the stored profile);
    ///   B  the gene's windows from its translation instead (expected axes, as the query), each with its own window rule;
    ///   B0 as B but over the common length from the start, with the query's window rule (G8's comparison at offset 0);
    ///   C  as B0, best over shifts of -30..30 aa in steps of 3 (d > 0 drops the gene's first d residues, d < 0 the query's).
    /// Homolog key from alignment (global, free end gaps, BLOSUM62, gap 11 + k): a gene is a homolog when its identity is at least Rost's
    /// HSSP curve + 5 points for the aligned length and the alignment covers 60% of the shorter protein. Aligned per query: the 50 genes
    /// sharing the most amino-acid 3-mers with it, and every arm's top 10, so every arm's top hit is judged directly. A query counts when
    /// any aligned gene is a homolog. Measures: top-1 and top-10 homolog, per arm, with paired bootstrap CIs over queries.
    /// Rule fixed before the run: B beats A by >= 3 points of top-1 with the CI above 0 -> a card for profiling genes from their
    /// translation (a gene-layer rebuild); C beats the better of B and B0 by >= 3 points with the CI above 0 -> a card for an offset search
    /// in the catalogue; otherwise neither. Writes ProteinHomologRanking.csv.
    /// </summary>
    public class ProteinHomologRankingBenchmark
    {
        private const int minResidues = 50, queriesPerGenome = 10, prefilterTop = 50, topK = 10, maxShift = 30, shiftStep = 3,
            bootstrapRounds = 2000, seed = 20261120;
        private const double minCoverage = 0.6, hsspMargin = 5;
        private const int F = ZCurveGeneProfile.FeatureCount;
        private static readonly double[] varianceWeights = { 0.741, 0.838, 0.333, 1.17, 1.16, 1.15, 1.19, 1.21, 1.2 };
        private static readonly string[] bandNames = { "50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa" };
        private static readonly string[] arms = { "A", "B", "B0", "C" };

        private sealed record Gene(int Genome, string Protein, string Code, double[] DnaWindows, double[] ProteinWindows, int[] Trimers);

        [Test]
        public void HomologRankingWithOffsetsAndTranslatedTargets()
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
                    indices.Add(pool.Count);
                    pool.Add(new Gene(genome, protein, code, ZCurveGeneProfile.QueryFromGene(dna, code), ZCurveGeneProfile.QueryFromProtein(protein, null, code), Trimers(protein)));
                }
                byGenome.Add(indices);
            }
            var rng = new Random(seed);
            var queries = byGenome.SelectMany(indices => indices.OrderBy(_ => rng.Next()).Take(queriesPerGenome)).ToList();
            TestContext.Out.WriteLine($"{pool.Count:N0} genes in {files.Count} genomes, {queries.Count} queries.");

            // Self-checks: the windows rebuilt from prefix sums equal the library's, and arm A's unweighted score is ZCurveGeneProfile.Score's.
            foreach (var q in queries.Take(20))
            {
                var gene = pool[q];
                var rebuilt = Windows(Prefix(gene.Protein, gene.Code), 0, gene.Protein.Length, gene.Protein.Length >= ZCurveGeneProfile.LongResidues);
                Assert.That(rebuilt.Length, Is.EqualTo(gene.ProteinWindows.Length), "prefix-sum window count");
                Assert.That(rebuilt.Zip(gene.ProteinWindows).Max(p => Math.Abs(p.First - p.Second)), Is.LessThan(1e-7), "prefix-sum windows");
            }

            int W = 2, R = arms.Length;
            var results = new Outcome[queries.Count];
            Parallel.For(0, queries.Count, n =>
            {
                var query = pool[queries[n]];
                int length = query.Protein.Length;
                bool longRule = length >= ZCurveGeneProfile.LongResidues;
                var queryPrefix = Prefix(query.Protein, query.Code);
                var geneBuffer = Enumerable.Range(0, 9).Select(_ => new double[1]).ToArray();
                var candidates = new List<int>();
                for (int g = 0; g < pool.Count; g++)
                {
                    var gene = pool[g];
                    if (genusOf[gene.Genome] == genusOf[query.Genome]) continue;
                    int residues = gene.Protein.Length;
                    if (longRule ? residues >= ZCurveGeneProfile.LongResidues
                        : residues < ZCurveGeneProfile.LongResidues && Math.Abs(residues - length) <= 0.1 * length) candidates.Add(g);
                }

                // Scores per arm and weighting; the top 10 of each kept.
                var tops = new PriorityQueue<int, double>[R, W];
                for (int a = 0; a < R; a++) for (int w = 0; w < W; w++) tops[a, w] = new PriorityQueue<int, double>();
                void Keep(int arm, int weighting, int candidate, double score)
                {
                    var top = tops[arm, weighting];
                    if (top.Count < topK) top.Enqueue(candidate, score);
                    else if (top.TryPeek(out _, out double lowest) && score > lowest) top.EnqueueDequeue(candidate, score);
                }
                foreach (int c in candidates)
                {
                    var gene = pool[c];
                    var genePrefix = Prefix(gene.Protein, gene.Code, geneBuffer);
                    for (int w = 0; w < W; w++)
                    {
                        var weights = w == 0 ? null : varianceWeights;
                        Keep(0, w, c, Score(query.ProteinWindows, gene.DnaWindows, weights));
                        Keep(1, w, c, Score(query.ProteinWindows, gene.ProteinWindows, weights));
                        double best = double.NegativeInfinity, atZero = double.NaN;
                        for (int d = -maxShift; d <= maxShift; d += shiftStep)
                        {
                            int qStart = Math.Max(0, -d), gStart = Math.Max(0, d);
                            int common = Math.Min(length - qStart, gene.Protein.Length - gStart);
                            if (common < minResidues) continue;
                            double score = Score(Windows(queryPrefix, qStart, common, longRule), Windows(genePrefix, gStart, common, longRule), weights);
                            if (d == 0) atZero = score;
                            best = Math.Max(best, score);
                        }
                        if (!double.IsNaN(atZero)) Keep(2, w, c, atZero);
                        if (best > double.NegativeInfinity) Keep(3, w, c, best);
                    }
                }
                var ranked = new int[R, W][];
                for (int a = 0; a < R; a++)
                    for (int w = 0; w < W; w++)
                        ranked[a, w] = tops[a, w].UnorderedItems.OrderByDescending(e => e.Priority).Select(e => e.Element).ToArray();

                // Homolog key: align the 3-mer prefilter's top 50 and every arm's top 10.
                var toAlign = candidates.Select(c => (c, Shared(query.Trimers, pool[c].Trimers))).OrderByDescending(x => x.Item2).Take(prefilterTop).Select(x => x.c)
                    .Concat(Enumerable.Range(0, R).SelectMany(a => Enumerable.Range(0, W).SelectMany(w => ranked[a, w]))).ToHashSet();
                var homolog = new Dictionary<int, (bool Is, double Identity)>();
                foreach (int c in toAlign)
                {
                    var (identity, aligned) = ProteinAlignment.Align(query.Protein, pool[c].Protein);
                    bool isHomolog = aligned >= minCoverage * Math.Min(length, pool[c].Protein.Length)
                        && identity * 100 >= Hssp(aligned) + hsspMargin;
                    homolog[c] = (isHomolog, identity);
                }
                bool reachable = homolog.Values.Any(h => h.Is);
                double bestIdentity = homolog.Values.Where(h => h.Is).Select(h => h.Identity).DefaultIfEmpty(double.NaN).Max();
                var top1 = new bool[R, W];
                var top10 = new bool[R, W];
                for (int a = 0; a < R; a++)
                    for (int w = 0; w < W; w++)
                    {
                        top1[a, w] = ranked[a, w].Length > 0 && homolog[ranked[a, w][0]].Is;
                        top10[a, w] = ranked[a, w].Any(c => homolog[c].Is);
                    }
                results[n] = new Outcome(BandOf(length), reachable, bestIdentity, top1, top10);
            });

            var reach = results.Where(r => r.Reachable).ToList();
            var lines = new List<string> { "Arm,Weighting,Band,Reachable,Top1,Top10" };
            foreach (int band in Enumerable.Range(0, bandNames.Length).Prepend(-1))
            {
                var cell = reach.Where(r => band < 0 || r.Band == band).ToList();
                if (cell.Count == 0) continue;
                for (int a = 0; a < R; a++)
                    for (int w = 0; w < W; w++)
                        lines.Add(string.Join(",", arms[a], w == 0 ? "unweighted" : "1/variance", band < 0 ? "all" : bandNames[band], cell.Count,
                            Format(cell.Count(r => r.Top1[a, w]) / (double)cell.Count), Format(cell.Count(r => r.Top10[a, w]) / (double)cell.Count)));
            }
            lines.Add("");
            lines.Add("Best homolog identity,Queries,A top-1 unweighted,C top-1 unweighted");
            foreach (var (low, high) in new[] { (0.0, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01) })
            {
                var bin = reach.Where(r => r.BestIdentity >= low && r.BestIdentity < high).ToList();
                if (bin.Count > 0)
                    lines.Add($"{low:F1}-{Math.Min(high, 1):F1},{bin.Count},{Format(bin.Count(r => r.Top1[0, 0]) / (double)bin.Count)},{Format(bin.Count(r => r.Top1[3, 0]) / (double)bin.Count)}");
            }

            // Paired differences of top-1 with bootstrap CIs over queries.
            (double Delta, double Low, double High) Paired(Func<Outcome, bool> after, Func<Outcome, bool> before)
            {
                double delta = reach.Average(r => (after(r) ? 1 : 0) - (before(r) ? 1 : 0));
                var random = new Random(seed);
                var deltas = Enumerable.Range(0, bootstrapRounds).Select(_ =>
                    Enumerable.Range(0, reach.Count).Select(_ => reach[random.Next(reach.Count)]).Average(r => (after(r) ? 1 : 0) - (before(r) ? 1 : 0))).Order().ToArray();
                return (delta, deltas[(int)(0.025 * bootstrapRounds)], deltas[(int)(0.975 * bootstrapRounds) - 1]);
            }
            string Points((double Delta, double Low, double High) d) =>
                $"{100 * d.Delta:+0.0;-0.0} points (95% CI {100 * d.Low:+0.0;-0.0} to {100 * d.High:+0.0;-0.0})";
            lines.Add("");
            lines.Add($"Queries with a homolog,{reach.Count} of {results.Length}");
            for (int w = 0; w < W; w++)
            {
                string name = w == 0 ? "unweighted" : "1/variance";
                int ww = w;
                var translation = Paired(r => r.Top1[1, ww], r => r.Top1[0, ww]);
                var overlap = Paired(r => r.Top1[2, ww], r => r.Top1[1, ww]);
                lines.Add($"{name}: B (translation) vs A (DNA),{Points(translation)}");
                lines.Add($"{name}: B0 (common length) vs B,{Points(overlap)}");
                lines.Add($"{name}: C (shifts) vs the better of B and B0,{Points(Paired(r => r.Top1[3, ww], r => BetterOf(r, ww)))}");
            }
            lines.Add($"Weights: A 1/variance vs A unweighted,{Points(Paired(r => r.Top1[0, 1], r => r.Top1[0, 0]))}");
            var ruleTranslation = Paired(r => r.Top1[1, 0], r => r.Top1[0, 0]);
            var ruleShift = Paired(r => r.Top1[3, 0], r => BetterOf(r, 0));
            lines.Add("Rule (unweighted; >= 3 points with CI above 0)," +
                $"translation: {(ruleTranslation.Delta >= 0.03 && ruleTranslation.Low > 0 ? "card for profiling genes from their translation" : "no")}; " +
                $"shifts: {(ruleShift.Delta >= 0.03 && ruleShift.Low > 0 ? "card for an offset search" : "no")}");

            File.WriteAllLines("ProteinHomologRanking.csv", lines);
            foreach (var line in lines) TestContext.Out.WriteLine(line);
            Assert.That(reach.Count, Is.GreaterThan(50), "too few queries with a homolog to judge");
        }

        // The better of B and B0 per query: B0's comparison or B's, whichever ranks the homolog first (the arm C is measured against).
        private static bool BetterOf(Outcome r, int w) => r.Top1[1, w] || r.Top1[2, w];

        private sealed record Outcome(int Band, bool Reachable, double BestIdentity, bool[,] Top1, bool[,] Top10);

        private static int BandOf(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;

        // Rost (1999) HSSP curve: the identity (%) above which an alignment of L residues is in the safe zone of homology.
        private static double Hssp(int length) => length <= 11 ? 100 : length > 450 ? 19.5 : 480 * Math.Pow(length, -0.32 * (1 + Math.Exp(-length / 1000.0)));

        private static double Score(double[] query, double[] gene, double[]? weights)
        {
            int count = Math.Min(query.Length, gene.Length) / F;
            if (count == 0) return double.NegativeInfinity;
            double sum = 0;
            for (int w = 0; w < count; w++)
            {
                double s = 0;
                for (int f = 0; f < F; f++)
                {
                    double d = query[w * F + f] - gene[w * F + f];
                    s += (weights?[f] ?? 1) * d * d;
                }
                sum += Math.Sqrt(s);
            }
            return 1 - sum / count;
        }

        private static readonly System.Collections.Concurrent.ConcurrentDictionary<string, (double X, double Y, double Z)[]> expectedTables = new();

        // Prefix sums over codon positions 1 and 2 of a protein's expected axes: for each axis, sums of a_i, a_i * i and a_i * i^2.
        // With a buffer, its arrays are reused (grown when too short) and returned.
        private static double[][] Prefix(string protein, string code, double[][]? buffer = null)
        {
            var table = expectedTables.GetOrAdd(code, c => ZCurveGeneProfile.ExpectedAxes(null, c));
            int bases = 2 * protein.Length;
            var sums = buffer ?? new double[9][];
            for (int k = 0; k < 9; k++)
                if (sums[k] is null || sums[k].Length < bases + 1) sums[k] = new double[Math.Max(bases + 1, 2 * (sums[k]?.Length ?? 0))];
            for (int i = 0; i < bases; i++)
            {
                int letter = char.ToUpperInvariant(protein[i / 2]) - 'A';
                var (x, y, z) = letter is >= 0 and < 26 ? table[2 * letter + i % 2] : (0, 0, 0);
                double[] a = { x, y, z };
                for (int axis = 0; axis < 3; axis++)
                {
                    sums[3 * axis][i + 1] = sums[3 * axis][i] + a[axis];
                    sums[3 * axis + 1][i + 1] = sums[3 * axis + 1][i] + a[axis] * i;
                    sums[3 * axis + 2][i + 1] = sums[3 * axis + 2][i] + a[axis] * (double)i * i;
                }
            }
            return sums;
        }

        // The windows of `residues` residues from residue `start`, by the window rule (long: 333-base windows; else 5 equal windows),
        // as ZCurveGeneProfile lays them out (window by window, 3 * order + axis), from the prefix sums in O(1) per window.
        private static double[] Windows(double[][] prefix, int start, int residues, bool longRule)
        {
            int bases = 2 * residues, window = longRule ? ZCurveGeneProfile.WindowBases : 2 * residues / ZCurveGeneProfile.SplitWindows;
            int count = window == 0 ? 0 : bases / window;
            var result = new double[count * F];
            for (int w = 0; w < count; w++)
            {
                int from = 2 * start + w * window, to = from + window;
                double shift = from - 0.5;
                for (int axis = 0; axis < 3; axis++)
                {
                    double s0 = prefix[3 * axis][to] - prefix[3 * axis][from];
                    double s1 = prefix[3 * axis + 1][to] - prefix[3 * axis + 1][from];
                    double s2 = prefix[3 * axis + 2][to] - prefix[3 * axis + 2][from];
                    // u = i - from + 0.5 and t = 2u / W - 1, so sums of a u and a u^2 follow from those of a i and a i^2.
                    double u1 = s1 - shift * s0, u2 = s2 - 2 * shift * s1 + shift * shift * s0;
                    double p0 = s0, p1 = 2 * u1 / window - s0, p2 = 6 * u2 / ((double)window * window) - 6 * u1 / window + s0;
                    result[w * F + axis] = p0 / window;
                    result[w * F + 3 + axis] = p1 * Math.Sqrt(3) / window;
                    result[w * F + 6 + axis] = p2 * Math.Sqrt(5) / window;
                }
            }
            return result;
        }

        // Distinct amino-acid 3-mers (20 letters), sorted, for the alignment prefilter.
        private static int[] Trimers(string protein)
        {
            const string letters = "ACDEFGHIKLMNPQRSTVWY";
            var set = new HashSet<int>();
            for (int i = 0; i + 3 <= protein.Length; i++)
            {
                int a = letters.IndexOf(protein[i]), b = letters.IndexOf(protein[i + 1]), c = letters.IndexOf(protein[i + 2]);
                if (a >= 0 && b >= 0 && c >= 0) set.Add(400 * a + 20 * b + c);
            }
            var sorted = set.ToArray();
            Array.Sort(sorted);
            return sorted;
        }

        private static int Shared(int[] a, int[] b)
        {
            int i = 0, j = 0, shared = 0;
            while (i < a.Length && j < b.Length)
            {
                if (a[i] == b[j]) { shared++; i++; j++; }
                else if (a[i] < b[j]) i++;
                else j++;
            }
            return shared;
        }

        private static string Format(double value) => value.ToString("F4", CultureInfo.InvariantCulture);
    }

    /// <summary>
    /// Global protein alignment with free end gaps (Gotoh, BLOSUM62, a gap of k residues costs 11 + k), for homolog keys in experiments.
    /// Returns the identity over aligned pairs and how many pairs there are.
    /// </summary>
    internal static class ProteinAlignment
    {
        private const string letters = "ARNDCQEGHILKMFPSTWYVBZX*";
        private static readonly int[] blosum = (
            " 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 " +
            "-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 " +
            "-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 " +
            "-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 " +
            " 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 " +
            "-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 " +
            "-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 " +
            " 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 " +
            "-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 " +
            "-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 " +
            "-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 " +
            "-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 " +
            "-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 " +
            "-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 " +
            "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 " +
            " 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 " +
            " 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 " +
            "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 " +
            "-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 " +
            " 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 " +
            "-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 " +
            "-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 " +
            " 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 " +
            "-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1")
            .Split(' ', StringSplitOptions.RemoveEmptyEntries).Select(s => int.Parse(s, CultureInfo.InvariantCulture)).ToArray();

        private static int Index(char c)
        {
            int i = letters.IndexOf(char.ToUpperInvariant(c));
            return i >= 0 ? i : letters.IndexOf('X');
        }

        public static (double Identity, int Aligned) Align(string query, string gene)
        {
            const int open = 12, extend = 1, negative = int.MinValue / 4;
            int n = query.Length, m = gene.Length, width = letters.Length;
            var trace = new byte[(n + 1) * (m + 1)];
            var mPrev = new int[m + 1]; var xPrev = new int[m + 1]; var yPrev = new int[m + 1];
            var mCur = new int[m + 1]; var xCur = new int[m + 1]; var yCur = new int[m + 1];
            for (int j = 0; j <= m; j++) { mPrev[j] = j == 0 ? 0 : negative; xPrev[j] = negative; yPrev[j] = j == 0 ? negative : 0; }
            var geneIndex = gene.Select(Index).ToArray();
            int bestScore = negative, bestI = 0, bestJ = 0, bestState = 0;
            for (int i = 1; i <= n; i++)
            {
                mCur[0] = negative; xCur[0] = 0; yCur[0] = negative;
                int row = Index(query[i - 1]) * width;
                for (int j = 1; j <= m; j++)
                {
                    byte t = 0;
                    int diagonal = mPrev[j - 1];
                    if (xPrev[j - 1] > diagonal) { diagonal = xPrev[j - 1]; t = 1; }
                    if (yPrev[j - 1] > diagonal) { diagonal = yPrev[j - 1]; t = 2; }
                    mCur[j] = diagonal + blosum[row + geneIndex[j - 1]];
                    int xOpen = mPrev[j] - open, xExtend = xPrev[j] - extend;
                    if (xExtend > xOpen) { xCur[j] = xExtend; t |= 4; } else xCur[j] = xOpen;
                    int yOpen = mCur[j - 1] - open, yExtend = yCur[j - 1] - extend;
                    if (yExtend > yOpen) { yCur[j] = yExtend; t |= 8; } else yCur[j] = yOpen;
                    trace[i * (m + 1) + j] = t;
                    if (i == n || j == m)
                    {
                        if (mCur[j] > bestScore) (bestScore, bestI, bestJ, bestState) = (mCur[j], i, j, 0);
                        if (xCur[j] > bestScore) (bestScore, bestI, bestJ, bestState) = (xCur[j], i, j, 1);
                        if (yCur[j] > bestScore) (bestScore, bestI, bestJ, bestState) = (yCur[j], i, j, 2);
                    }
                }
                (mPrev, mCur) = (mCur, mPrev); (xPrev, xCur) = (xCur, xPrev); (yPrev, yCur) = (yCur, yPrev);
            }
            int ci = bestI, cj = bestJ, state = bestState, pairs = 0, identical = 0;
            while (ci > 0 && cj > 0)
            {
                byte t = trace[ci * (m + 1) + cj];
                if (state == 0)
                {
                    pairs++;
                    if (char.ToUpperInvariant(query[ci - 1]) == char.ToUpperInvariant(gene[cj - 1])) identical++;
                    state = t & 3; ci--; cj--;
                }
                else if (state == 1) { state = (t & 4) != 0 ? 1 : 0; ci--; }
                else { state = (t & 8) != 0 ? 2 : 0; cj--; }
            }
            return (pairs == 0 ? 0 : identical / (double)pairs, pairs);
        }
    }
}
