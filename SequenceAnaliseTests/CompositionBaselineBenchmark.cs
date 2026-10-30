using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// P5, R1's baselines on the bundled genomes (2026-10-30). E4: are the signatures the catalogue can only derive from its stored
    /// profile (Natural Vector, Zheng-style curve, Z-curve mean) faithful to their exact definitions from the DNA? Per signature, the
    /// Spearman correlation of exact and derived pairwise distances and the largest feature error in units of that feature's SD
    /// across genomes. E1: every signature and Z-curve distance through one evaluation, with Mash for the bacteria: AUC per rank and
    /// clades recovered by NJ and UPGMA trees with a 1,000-permutation label null. Distances: compare (1 - Compare; bacteria at
    /// 32 kbp, viruses at the ladder level), summary (window statistics at 32 kbp, viruses 500 bp), zMean (whole-genome mean x, y,
    /// z), zheng, naturalVector (standardized, Euclidean), karlin (delta*), gc, gcLength, mash. Exact signatures from the DNA here.
    /// Writes S52b_ProfileDerivedCheck.csv and S52b_Baselines_Bundled.csv. About a minute.
    /// </summary>
    public class CompositionBaselineBenchmark
    {
        private const int permutations = 1000;

        private sealed record Genome(string File, string Name, string Genus, string Higher);

        private static List<Genome> Corpus(string corpus) => corpus == "bacteria"
            ? GenomeTaxonomy.Genomes.Select(g => new Genome(g.FileName, g.OrganismName, g.Genus, g.Phylum)).ToList()
            : ViralTaxonomy.Genomes.Select(g => new Genome(g.FileName, g.Name, g.Genus, g.Family)).ToList();

        [Test]
        public void ProfileDerivedSignaturesOnBundledGenomes()
        {
            using var report = File.CreateText("S52b_ProfileDerivedCheck.csv");
            report.WriteLine("Corpus, Signature, Genomes, Pairs, SpearmanExactVsDerived, MaxFeatureErrorInSds");
            foreach (var corpus in new[] { "bacteria", "viruses" })
            {
                var exact = new Dictionary<string, List<double[]>> { ["naturalVector"] = new(), ["zheng"] = new(), ["zMean"] = new() };
                var derived = new Dictionary<string, List<double[]>> { ["naturalVector"] = new(), ["zheng"] = new(), ["zMean"] = new() };
                foreach (var g in Corpus(corpus))
                {
                    var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.File)).Sequence;
                    var fine = new ZCurveWindowedProfile(sequence);
                    exact["naturalVector"].Add(GenomeSignatures.NaturalVector(sequence));
                    derived["naturalVector"].Add(GenomeSignatures.NaturalVector(fine, sequence.Length));
                    exact["zheng"].Add(GenomeSignatures.ZhengCurve(sequence));
                    derived["zheng"].Add(GenomeSignatures.ZhengCurve(fine));
                    exact["zMean"].Add(GenomeSignatures.ZCurveMean(sequence));
                    derived["zMean"].Add(GenomeSignatures.ZCurveMean(fine));
                }
                foreach (var signature in exact.Keys)
                {
                    var e = exact[signature].ToArray();
                    var d = derived[signature].ToArray();
                    // Both standardized with the exact signature's means and SDs, so the distances are on one scale.
                    var (eStd, dStd) = (StandardizeBy(e, e), StandardizeBy(d, e));
                    double rho = Spearman(Pairs(ZCurveTreeBenchmark.Euclidean(eStd)), Pairs(ZCurveTreeBenchmark.Euclidean(dStd)));
                    double maxError = eStd.Zip(dStd, (x, y) => x.Zip(y, (p, q) => Math.Abs(p - q)).Max()).Max();
                    string line = string.Create(CultureInfo.InvariantCulture,
                        $"{corpus}, {signature}, {e.Length}, {e.Length * (e.Length - 1) / 2}, {rho:F4}, {maxError:F4}");
                    report.WriteLine(line);
                    TestContext.Out.WriteLine(line);
                }
            }
            Assert.Pass();
        }

        [Test]
        public void BaselinesOnBundledGenomes()
        {
            using var report = File.CreateText("S52b_Baselines_Bundled.csv");
            report.WriteLine("Corpus, Distance, Rank, Positives, Negatives, Auc, Taxa, NjClades, UpgmaClades, NjNullMean, UpgmaNullMean, NjP, UpgmaP");
            foreach (var corpus in new[] { "bacteria", "viruses" })
            {
                var genomes = Corpus(corpus);
                int n = genomes.Count;
                var sequences = genomes.Select(g => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.File)).Sequence).ToArray();
                var fine = sequences.Select(s => new ZCurveWindowedProfile(s)).ToArray();
                var gc = sequences.Select(ZCurveTreeBenchmark.Gc).ToArray();
                var length = sequences.Select(s => (double)s.Length).ToArray();
                var rho = sequences.Select(s => GenomeSignatures.KarlinRho(GenomeSignatures.DinucleotideCounts(s))).ToArray();

                var compare = new double[n, n];
                var ladders = fine.Select(f => new ZCurveProfileLadder(f)).ToArray();
                var coarse = fine.Select(f => f.Downsample(64)).ToArray();
                var karlin = new double[n, n];
                for (int a = 0; a < n; a++)
                    for (int b = a + 1; b < n; b++)
                    {
                        double score = corpus == "bacteria" ? ZCurveWindowedProfile.Compare(coarse[a], coarse[b])
                            : ZCurveProfileLadder.Compare(ladders[a], ladders[b]) ?? 0;
                        compare[a, b] = compare[b, a] = 1 - score;
                        karlin[a, b] = karlin[b, a] = GenomeSignatures.KarlinDelta(rho[a], rho[b]);
                    }

                var distances = new List<(string Name, double[,] Matrix)>
                {
                    ("compare", compare),
                    ("summary", Std((corpus == "bacteria" ? coarse : fine).Select(GenomeSignatures.WindowSummary))),
                    ("zMean", ZCurveTreeBenchmark.Euclidean(sequences.Select(GenomeSignatures.ZCurveMean).ToArray())),
                    ("zheng", Std(sequences.Select(GenomeSignatures.ZhengCurve))),
                    ("naturalVector", Std(sequences.Select(GenomeSignatures.NaturalVector))),
                    ("karlin", karlin),
                    ("gc", ZCurveTreeBenchmark.Euclidean(gc.Select(v => new[] { v }).ToArray())),
                    ("gcLength", Std(gc.Select((v, i) => new[] { v, Math.Log10(length[i]) }))),
                };
                if (corpus == "bacteria")
                {
                    var mash = ExternalBaselineData.ParseMashDistances(ExternalBaselineData.RawResultsPath("mash_dist_raw.tsv"));
                    var matrix = new double[n, n];
                    for (int a = 0; a < n; a++)
                        for (int b = 0; b < n; b++)
                            if (a != b) matrix[a, b] = mash[ExternalBaselineData.PairKey(genomes[a].File, genomes[b].File)];
                    distances.Add(("mash", matrix));
                }

                var labels = genomes.Select(g => $"{g.Name}|{g.Higher}").ToArray();
                var ranks = new[] { ("genus", genomes.Select(g => g.Genus).ToArray()),
                    (corpus == "bacteria" ? "phylum" : "family", genomes.Select(g => g.Higher).ToArray()) };
                foreach (var (name, matrix) in distances)
                {
                    var nj = Phylogeny.NeighbourJoining(matrix, labels);
                    var upgma = Phylogeny.Upgma(matrix, labels);
                    foreach (var (rank, taxa) in ranks)
                    {
                        var positives = new List<double>();
                        var negatives = new List<double>();
                        for (int a = 0; a < n; a++)
                            for (int b = a + 1; b < n; b++)
                                (taxa[a] == taxa[b] ? positives : negatives).Add(-matrix[a, b]);
                        double auc = RetrievalClassificationBenchmark.ComputeAuc(positives, negatives);
                        var (njClades, njNull, njP, total) = Clades(nj, taxa);
                        var (upgmaClades, upgmaNull, upgmaP, _) = Clades(upgma, taxa);
                        string line = string.Create(CultureInfo.InvariantCulture,
                            $"{corpus}, {name}, {rank}, {positives.Count}, {negatives.Count}, {auc:F4}, {total}, {njClades}, {upgmaClades}, " +
                            $"{njNull:F2}, {upgmaNull:F2}, {njP:F4}, {upgmaP:F4}");
                        report.WriteLine(line);
                        TestContext.Out.WriteLine(line);
                    }
                }
            }
            Assert.Pass();
        }

        private static (int Clades, double NullMean, double P, int Taxa) Clades(PhyloTree tree, string[] taxa)
        {
            var recovery = tree.Recovery(taxa);
            int observed = recovery.Count(r => r.Monophyletic);
            var rng = new Random(20261030);
            var nullCounts = Enumerable.Range(0, permutations)
                .Select(_ => tree.Recovery(taxa.OrderBy(_ => rng.Next()).ToArray()).Count(r => r.Monophyletic)).ToList();
            return (observed, nullCounts.Average(), (1 + nullCounts.Count(c => c >= observed)) / (permutations + 1.0), recovery.Count);
        }

        private static double[,] Std(IEnumerable<double[]> rows) => ZCurveTreeBenchmark.Euclidean(ZCurveTreeBenchmark.Standardize(rows.ToArray()));

        private static double[][] StandardizeBy(double[][] rows, double[][] reference)
        {
            int width = reference[0].Length;
            var means = Enumerable.Range(0, width).Select(f => reference.Average(r => r[f])).ToArray();
            var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(reference.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
            return rows.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray();
        }

        private static List<double> Pairs(double[,] matrix)
        {
            var values = new List<double>();
            for (int a = 0; a < matrix.GetLength(0); a++)
                for (int b = a + 1; b < matrix.GetLength(0); b++)
                    values.Add(matrix[a, b]);
            return values;
        }

        private static double Spearman(List<double> xs, List<double> ys)
        {
            var (rx, ry) = (Ranks(xs), Ranks(ys));
            double mx = rx.Average(), my = ry.Average();
            double sxy = rx.Zip(ry, (x, y) => (x - mx) * (y - my)).Sum();
            return sxy / Math.Sqrt(rx.Sum(x => (x - mx) * (x - mx)) * ry.Sum(y => (y - my) * (y - my)));
        }

        private static double[] Ranks(List<double> values)
        {
            var index = Enumerable.Range(0, values.Count).OrderBy(i => values[i]).ToArray();
            var ranks = new double[values.Count];
            for (int i = 0; i < index.Length;)
            {
                int j = i;
                while (j + 1 < index.Length && values[index[j + 1]] == values[index[i]]) j++;
                for (int t = i; t <= j; t++) ranks[index[t]] = (i + j) / 2.0;
                i = j + 1;
            }
            return ranks;
        }
    }
}
