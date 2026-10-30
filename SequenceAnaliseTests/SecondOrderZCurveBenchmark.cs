using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// P5, second-order Z-curve, step 2 (2026-10-30): does windowing the second-order Z-curve (lag-1 correlations between
    /// neighbouring bases' axes, <see cref="GenomeSignatures.SecondOrder"/>) add anything over the global one? On the bundled
    /// genomes, which have DNA. Per window the 9 correlations over both strands; two comparisons: windowSummary (mean and SD of each
    /// correlation across windows, standardized, Euclidean) and sliding (<see cref="GenomeSignatures.SlidingDistance"/>, either
    /// window order). Bacteria at 32 kbp and 8 kbp windows, laid out from the replication origin where OriginRotation's rule places
    /// one; viruses at 2 kbp, as read. Context: karlin (delta*), secondOrder (global), the first-order window summary, Mash for the
    /// bacteria. AUC per rank and NJ / UPGMA clades with a 1,000-permutation label null. The bundled bacteria have 8 same-genus pairs,
    /// so this only points the way; the decision is the catalogue's. Writes S52c_SecondOrder_Bundled.csv. A few minutes.
    /// </summary>
    public class SecondOrderZCurveBenchmark
    {
        [Test]
        public void SecondOrderOnBundledGenomes()
        {
            using var report = File.CreateText("S52c_SecondOrder_Bundled.csv");
            report.WriteLine("Corpus, Distance, Rank, Positives, Negatives, Auc, Taxa, NjClades, UpgmaClades, NjNullMean, UpgmaNullMean, NjP, UpgmaP");

            var bacteria = GenomeTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.OrganismName}|{g.Phylum}", g.Genus, Higher: g.Phylum)).ToList();
            var viruses = ViralTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.Name}|{g.Family}", g.Genus, Higher: g.Family)).ToList();
            foreach (var (corpus, genomes, windowSizes, higher) in new[] { ("bacteria", bacteria, new[] { 32_000, 8_000 }, "phylum"),
                ("viruses", viruses, new[] { 2_000 }, "family") })
            {
                int n = genomes.Count;
                var sequences = genomes.Select(g => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence).ToArray();
                var fine = sequences.Select(s => new ZCurveWindowedProfile(s)).ToArray();
                if (corpus == "bacteria")
                    for (int i = 0; i < n; i++)
                        if (ZCurveTreeBenchmark.OriginWindow(fine[i], sequences[i].Length) is int origin)
                        {
                            int start = origin * fine[i].WindowSize;
                            sequences[i] = sequences[i][start..] + sequences[i][..start];
                        }
                var counts = sequences.Select(GenomeSignatures.DinucleotideCounts).ToArray();

                var distances = new List<(string Name, double[,] Matrix)>
                {
                    ("karlin", Matrix(n, (a, b) => GenomeSignatures.KarlinDelta(GenomeSignatures.KarlinRho(counts[a]), GenomeSignatures.KarlinRho(counts[b])))),
                    ("secondOrder", ZCurveTreeBenchmark.Euclidean(counts.Select(c => GenomeSignatures.SecondOrder(c)).ToArray())),
                    ("summary", CompositionBaselineBenchmark.Std((corpus == "bacteria" ? fine.Select(f => f.Downsample(64)) : fine).Select(GenomeSignatures.WindowSummary))),
                };
                foreach (int size in windowSizes)
                {
                    var windows = sequences.Select(s => GenomeSignatures.SecondOrderWindows(s, size)).ToArray();
                    distances.Add(($"windowSummary{size / 1000}k", CompositionBaselineBenchmark.Std(windows.Select(MeanAndSd))));
                    var sliding = new double[n, n];
                    Parallel.For(0, n * n, cell =>
                    {
                        int a = cell / n, b = cell % n;
                        if (b > a) sliding[a, b] = sliding[b, a] = GenomeSignatures.SlidingDistance(windows[a], windows[b]);
                    });
                    double worst = Pairs(sliding).Where(double.IsFinite).DefaultIfEmpty(1).Max();
                    for (int a = 0; a < n; a++)
                        for (int b = 0; b < n; b++)
                            if (!double.IsFinite(sliding[a, b])) sliding[a, b] = worst;   // no windows: as unrelated as anything here
                    distances.Add(($"sliding{size / 1000}k", sliding));
                }
                if (corpus == "bacteria")
                {
                    var mash = ExternalBaselineData.ParseMashDistances(ExternalBaselineData.RawResultsPath("mash_dist_raw.tsv"));
                    distances.Add(("mash", Matrix(n, (a, b) => mash[ExternalBaselineData.PairKey(genomes[a].FileName, genomes[b].FileName)])));
                }

                var labels = genomes.Select(g => g.Label).ToArray();
                var ranks = new[] { ("genus", genomes.Select(g => g.Genus).ToArray()), (higher, genomes.Select(g => g.Higher).ToArray()) };
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
                        var (njClades, njNull, njP, total) = CompositionBaselineBenchmark.Clades(nj, taxa);
                        var (upgmaClades, upgmaNull, upgmaP, _) = CompositionBaselineBenchmark.Clades(upgma, taxa);
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

        private static double[] MeanAndSd(double[][] windows)
        {
            if (windows.Length == 0) return new double[18];
            var result = new double[18];
            for (int k = 0; k < 9; k++)
            {
                double mean = windows.Average(w => w[k]);
                (result[k], result[9 + k]) = (mean, Math.Sqrt(windows.Average(w => (w[k] - mean) * (w[k] - mean))));
            }
            return result;
        }

        private static double[,] Matrix(int n, Func<int, int, double> distance)
        {
            var result = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    result[a, b] = result[b, a] = distance(a, b);
            return result;
        }

        private static IEnumerable<double> Pairs(double[,] matrix)
        {
            for (int a = 0; a < matrix.GetLength(0); a++)
                for (int b = a + 1; b < matrix.GetLength(0); b++)
                    yield return matrix[a, b];
        }
    }
}
