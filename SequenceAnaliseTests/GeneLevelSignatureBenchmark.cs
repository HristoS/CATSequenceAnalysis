using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// R6e applied to neighbours and trees (2026-11-04): do signatures built from the genes follow taxonomy better than whole-genome
    /// dinucleotides (Karlin)? Codon positions 1 and 2 carry the protein, while position 3 carries most of the GC drift that makes
    /// whole-genome composition trees GC trees (S5.2). Bundled genomes, which have DNA.
    /// Genes are stop-free stretches on both strands from their first in-frame ATG, at least minGeneLength long, taken longest first,
    /// and a candidate overlapping accepted genes by more than maxOverlap bases is dropped (antisense ORFs, long in GC-rich genomes
    /// because stops are AT-rich, would otherwise add a GC-driven second proteome).
    /// Signatures, fixed before the run: pos12Mean (the genome mean of R6e's profile from the DNA, x, y, z at codon positions 1 and 2),
    /// expectedPos12Mean (the same from the proteins with R6e's uniform expected table), aminoAcids (20 frequencies) and codons (64
    /// frequencies; stops are zero except TGA in Mycoplasma), all standardized Euclidean; pos12Rho (Karlin's delta on the within-codon
    /// position-1-2 dinucleotides, coding strand, position-specific marginals); cvTree3 and cvTree4 (amino-acid k-mers with Hao's
    /// Markov correction, distance (1 - cosine) / 2); geneBestHit (R6e's search used directly: each gene's first two positions-1-and-2
    /// windows against the other genome's genes, 1 minus the mean best score, averaged both ways). Context: karlin, gc, for the
    /// bacteria Mash, and dnaCv6 (added after the first run, outside the rule). AUC per rank and NJ / UPGMA clades with the S52c label null.
    /// Rule fixed before the run: a gene-level signature passes if on the bacteria it beats karlin's phylum AUC and its genus clades
    /// (NJ + UPGMA); if none does, stop. The bundled bacteria have 8 same-genus pairs, so a pass only earns a catalogue sample.
    /// Writes R6e_Trees_Bundled.csv. About a minute.
    /// </summary>
    public class GeneLevelSignatureBenchmark
    {
        private const int minGeneLength = 300;
        private const int maxOverlap = 60;
        private const int bestHitGeneLength = 999;
        private const string aminoAcidLetters = "ACDEFGHIKLMNPQRSTVWY";
        private static readonly string[] geneLevel =
            { "pos12Mean", "expectedPos12Mean", "aminoAcids", "codons", "pos12Rho", "cvTree3", "cvTree4", "geneBestHit" };

        [Test]
        public void GeneLevelSignaturesOnBundledGenomes()
        {
            using var report = File.CreateText("R6e_Trees_Bundled.csv");
            report.WriteLine("Corpus, Distance, Rank, Positives, Negatives, Auc, Taxa, NjClades, UpgmaClades, NjNullMean, UpgmaNullMean, NjP, UpgmaP");

            var bacteria = GenomeTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.OrganismName}|{g.Phylum}", g.Genus, Higher: g.Phylum)).ToList();
            var viruses = ViralTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.Name}|{g.Family}", g.Genus, Higher: g.Family)).ToList();
            var results = new Dictionary<(string Corpus, string Distance, string Rank), (double Auc, int Clades)>();
            foreach (var (corpus, genomes, higher) in new[] { ("bacteria", bacteria, "phylum"), ("viruses", viruses, "family") })
            {
                int n = genomes.Count;
                var sequences = genomes.Select(g => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence).ToArray();
                var codes = genomes.Select(g => Path.GetFileName(g.FileName).StartsWith("mycoplasma_", StringComparison.Ordinal)
                    ? ProteinBackTranslationFeasibility.mycoplasmaCode : ProteinBackTranslationFeasibility.standardCode).ToArray();
                var genes = new List<string>[n];
                Parallel.For(0, n, i => genes[i] = CallGenes(sequences[i], codes[i]));
                var proteins = genes.Select((g, i) => g.Select(x => ProteinBackTranslationFeasibility.Translate(x, codes[i])).ToList()).ToArray();
                for (int i = 0; i < n; i++)
                    TestContext.Out.WriteLine(string.Create(CultureInfo.InvariantCulture,
                        $"{corpus}, {Path.GetFileName(genomes[i].FileName)}: {genes[i].Count} genes, {genes[i].Sum(g => (long)g.Length) / (double)sequences[i].Length:F3} of the length"));
                Assert.That(genes.All(g => g.Count >= (corpus == "bacteria" ? 100 : 1)), "every genome needs genes called (bacteria at least 100)");

                var counts = sequences.Select(GenomeSignatures.DinucleotideCounts).ToArray();
                var gc = sequences.Select(Gc).ToArray();
                var cvTree3 = proteins.Select(p => CvTreeVector(p, 3)).ToArray();
                var cvTree4 = proteins.Select(p => CvTreeVector(p, 4)).ToArray();
                var pos12Rho = genes.Select(Pos12Rho).ToArray();
                var dnaCv6 = sequences.Select(s => CvTreeVector(new[] { s, ProteinBackTranslationFeasibility.ReverseComplement(s) }, 6, "ACGT")).ToArray();
                var distances = new List<(string Name, double[,] Matrix)>
                {
                    ("karlin", Matrix(n, (a, b) => GenomeSignatures.KarlinDelta(GenomeSignatures.KarlinRho(counts[a]), GenomeSignatures.KarlinRho(counts[b])))),
                    ("gc", Matrix(n, (a, b) => Math.Abs(gc[a] - gc[b]))),
                    ("pos12Mean", CompositionBaselineBenchmark.Std(genes.Select(Pos12Mean))),
                    ("expectedPos12Mean", CompositionBaselineBenchmark.Std(proteins.Select((p, i) => ExpectedPos12Mean(p, codes[i])))),
                    ("aminoAcids", CompositionBaselineBenchmark.Std(proteins.Select(AminoAcidFrequencies))),
                    ("codons", CompositionBaselineBenchmark.Std(genes.Select(CodonFrequencies))),
                    ("pos12Rho", Matrix(n, (a, b) => GenomeSignatures.KarlinDelta(pos12Rho[a], pos12Rho[b]))),
                    ("cvTree3", Matrix(n, (a, b) => CosineDistance(cvTree3[a], cvTree3[b]))),
                    ("cvTree4", Matrix(n, (a, b) => CosineDistance(cvTree4[a], cvTree4[b]))),
                    ("geneBestHit", GeneBestHitMatrix(genes)),
                    // Context added after the first run: the same Markov correction on both-strand DNA 6-mers, which the platform
                    // stores per record since migration 0013, so no gene calls would be needed. Not part of the rule.
                    ("dnaCv6", Matrix(n, (a, b) => CosineDistance(dnaCv6[a], dnaCv6[b]))),
                };
                // Added 2026-11-04 at your suggestion, after the first runs: GC combined with the amino-acid signals. gcAminoAcids
                // standardizes GC with the 20 frequencies; gcCvTree4 adds the two distances, each over its median, so they weigh equally.
                var aminoAcidsWithGc = proteins.Select((p, i) => AminoAcidFrequencies(p).Append(gc[i]).ToArray());
                distances.Add(("gcAminoAcids", CompositionBaselineBenchmark.Std(aminoAcidsWithGc)));
                var cvTree4Matrix = distances.Single(d => d.Name == "cvTree4").Matrix;
                var gcMatrix = distances.Single(d => d.Name == "gc").Matrix;
                var (cvMedian, gcMedian) = (MedianOffDiagonal(cvTree4Matrix), MedianOffDiagonal(gcMatrix));
                distances.Add(("gcCvTree4", Matrix(n, (a, b) => cvTree4Matrix[a, b] / cvMedian + gcMatrix[a, b] / gcMedian)));
                if (corpus == "bacteria")
                {
                    var mash = ExternalBaselineData.ParseMashDistances(ExternalBaselineData.RawResultsPath("mash_dist_raw.tsv"));
                    distances.Add(("mash", Matrix(n, (a, b) => mash[ExternalBaselineData.PairKey(genomes[a].FileName, genomes[b].FileName)])));
                }

                ScoreTrees(corpus, genomes, higher, distances, "R6e", report, results);
            }
            ReportRule(geneLevel, results);
        }

        /// <summary>
        /// Gene classes (2026-11-04, your idea): isolate the signal per gene, group genes into classes, then describe a genome by how
        /// many of its genes fall in each class. Genes as in GeneLevelSignaturesOnBundledGenomes. Per-gene features, fixed before the
        /// run: zcurve (the gene as one dynamic window: mean x, y, z at codon positions 1, 2 and 3, 9 values) and codons (the gene's 64
        /// codon frequencies), standardized over all genes of the corpus. Classes by k-means (k-means++ start, seeded, kMeansRounds
        /// rounds) over all genes of the corpus, K = 64 and 256. Genome distance: Hellinger between class histograms. Same scoring and
        /// rule as the test above (beat karlin on bacterial phylum AUC and genus clades); karlin and cvTree4 as context.
        /// Writes R6e_GeneClasses_Bundled.csv. A few minutes.
        /// </summary>
        [Test]
        public void GeneClassHistogramsOnBundledGenomes()
        {
            const int kMeansRounds = 25;
            using var report = File.CreateText("R6e_GeneClasses_Bundled.csv");
            report.WriteLine("Corpus, Distance, Rank, Positives, Negatives, Auc, Taxa, NjClades, UpgmaClades, NjNullMean, UpgmaNullMean, NjP, UpgmaP");

            var bacteria = GenomeTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.OrganismName}|{g.Phylum}", g.Genus, Higher: g.Phylum)).ToList();
            var viruses = ViralTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.Name}|{g.Family}", g.Genus, Higher: g.Family)).ToList();
            var results = new Dictionary<(string Corpus, string Distance, string Rank), (double Auc, int Clades)>();
            var classSignatures = new List<string>();
            foreach (var (corpus, genomes, higher) in new[] { ("bacteria", bacteria, "phylum"), ("viruses", viruses, "family") })
            {
                int n = genomes.Count;
                var sequences = genomes.Select(g => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence).ToArray();
                var codes = genomes.Select(g => Path.GetFileName(g.FileName).StartsWith("mycoplasma_", StringComparison.Ordinal)
                    ? ProteinBackTranslationFeasibility.mycoplasmaCode : ProteinBackTranslationFeasibility.standardCode).ToArray();
                var genes = new List<string>[n];
                Parallel.For(0, n, i => genes[i] = CallGenes(sequences[i], codes[i]));
                var proteins = genes.Select((g, i) => g.Select(x => ProteinBackTranslationFeasibility.Translate(x, codes[i])).ToList()).ToArray();
                var genomeOfGene = genes.SelectMany((g, i) => g.Select(_ => i)).ToArray();
                var allGenes = genes.SelectMany(g => g).ToArray();

                var counts = sequences.Select(GenomeSignatures.DinucleotideCounts).ToArray();
                var cvTree4 = proteins.Select(p => CvTreeVector(p, 4)).ToArray();
                var distances = new List<(string Name, double[,] Matrix)>
                {
                    ("karlin", Matrix(n, (a, b) => GenomeSignatures.KarlinDelta(GenomeSignatures.KarlinRho(counts[a]), GenomeSignatures.KarlinRho(counts[b])))),
                    ("cvTree4", Matrix(n, (a, b) => CosineDistance(cvTree4[a], cvTree4[b]))),
                };
                foreach (var (feature, extract) in new (string, Func<string, double[]>)[] { ("zcurve", GeneZCurve), ("codons", GeneCodonFrequencies) })
                {
                    var points = ZCurveTreeBenchmark.Standardize(allGenes.Select(extract).ToArray());
                    foreach (var k in new[] { 64, 256 })
                    {
                        var assignment = KMeans(points, k, 20261104, kMeansRounds);
                        var histograms = new double[n][];
                        for (int a = 0; a < n; a++)
                            histograms[a] = new double[k];
                        for (int g = 0; g < assignment.Length; g++)
                            histograms[genomeOfGene[g]][assignment[g]] += 1.0 / genes[genomeOfGene[g]].Count;
                        var name = $"geneClass_{feature}_{k}";
                        if (!classSignatures.Contains(name)) classSignatures.Add(name);
                        distances.Add((name, Matrix(n, (a, b) => Hellinger(histograms[a], histograms[b]))));
                    }
                }
                ScoreTrees(corpus, genomes, higher, distances, "R6eClasses", report, results);
            }
            ReportRule(classSignatures, results);
        }

        internal static void ScoreTrees(
            string corpus, List<(string FileName, string Label, string Genus, string Higher)> genomes, string higher,
            List<(string Name, double[,] Matrix)> distances, string treePrefix, StreamWriter report,
            Dictionary<(string Corpus, string Distance, string Rank), (double Auc, int Clades)> results)
        {
            int n = genomes.Count;
            var labels = genomes.Select(g => g.Label).ToArray();
            var ranks = new[] { ("genus", genomes.Select(g => g.Genus).ToArray()), (higher, genomes.Select(g => g.Higher).ToArray()) };
            Directory.CreateDirectory("trees");
            foreach (var (name, matrix) in distances)
            {
                var nj = Phylogeny.NeighbourJoining(matrix, labels);
                var upgma = Phylogeny.Upgma(matrix, labels);
                foreach (var (method, tree) in new[] { ("NJ", nj), ("UPGMA", upgma) })
                    File.WriteAllText(Path.Combine("trees", $"{treePrefix}_{corpus}_{name}_{method}.nwk"), tree.ToNewick() + Environment.NewLine);
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
                    results[(corpus, name, rank)] = (auc, njClades + upgmaClades);
                    string line = string.Create(CultureInfo.InvariantCulture,
                        $"{corpus}, {name}, {rank}, {positives.Count}, {negatives.Count}, {auc:F4}, {total}, {njClades}, {upgmaClades}, " +
                        $"{njNull:F2}, {upgmaNull:F2}, {njP:F4}, {upgmaP:F4}");
                    report.WriteLine(line);
                    TestContext.Out.WriteLine(line);
                }
            }
        }

        internal static void ReportRule(IEnumerable<string> candidates, Dictionary<(string Corpus, string Distance, string Rank), (double Auc, int Clades)> results)
        {
            var karlinPhylum = results[("bacteria", "karlin", "phylum")].Auc;
            var karlinGenusClades = results[("bacteria", "karlin", "genus")].Clades;
            var passing = candidates.Where(s => results[("bacteria", s, "phylum")].Auc > karlinPhylum
                && results[("bacteria", s, "genus")].Clades > karlinGenusClades).ToList();
            var verdict = passing.Count > 0 ? "pass for " + string.Join(", ", passing) : "stop, no candidate beats karlin on both";
            TestContext.Out.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"Rule (karlin phylum AUC {karlinPhylum:F4}, genus clades NJ + UPGMA {karlinGenusClades}): {verdict}"));
        }

        // Mean x, y, z at codon positions 1, 2 and 3: the gene as one window.
        private static double[] GeneZCurve(string gene)
        {
            var sums = new double[9];
            int codons = gene.Length / 3;
            for (int i = 0; i < codons * 3; i++)
            {
                var (x, y, z) = ProteinBackTranslationFeasibility.BaseAxes(gene[i]);
                sums[3 * (i % 3)] += x;
                sums[3 * (i % 3) + 1] += y;
                sums[3 * (i % 3) + 2] += z;
            }
            return sums.Select(s => s / codons).ToArray();
        }

        private static double[] GeneCodonFrequencies(string gene) => CodonFrequencies(new List<string> { gene });

        // Lloyd's k-means from a k-means++ start; an emptied cluster keeps its previous centre.
        private static int[] KMeans(double[][] points, int k, int seed, int rounds)
        {
            int n = points.Length, width = points[0].Length;
            var rng = new Random(seed);
            var centres = new List<double[]> { points[rng.Next(n)] };
            var nearest = points.Select(p => SquaredDistance(p, centres[0])).ToArray();
            while (centres.Count < k)
            {
                double draw = rng.NextDouble() * nearest.Sum();
                int chosen = 0;
                for (; chosen < n - 1 && draw > nearest[chosen]; chosen++)
                    draw -= nearest[chosen];
                centres.Add(points[chosen]);
                var latest = centres[^1];
                Parallel.For(0, n, i => nearest[i] = Math.Min(nearest[i], SquaredDistance(points[i], latest)));
            }

            var assignment = new int[n];
            for (int round = 0; round < rounds; round++)
            {
                Parallel.For(0, n, i =>
                {
                    int best = 0;
                    double bestDistance = double.PositiveInfinity;
                    for (int c = 0; c < k; c++)
                    {
                        var distance = SquaredDistance(points[i], centres[c]);
                        if (distance < bestDistance) (best, bestDistance) = (c, distance);
                    }
                    assignment[i] = best;
                });
                var sums = new double[k, width];
                var sizes = new int[k];
                for (int i = 0; i < n; i++)
                {
                    sizes[assignment[i]]++;
                    for (int f = 0; f < width; f++)
                        sums[assignment[i], f] += points[i][f];
                }
                for (int c = 0; c < k; c++)
                    if (sizes[c] > 0)
                        centres[c] = Enumerable.Range(0, width).Select(f => sums[c, f] / sizes[c]).ToArray();
            }
            return assignment;
        }

        private static double SquaredDistance(double[] a, double[] b)
        {
            double sum = 0;
            for (int f = 0; f < a.Length; f++)
                sum += (a[f] - b[f]) * (a[f] - b[f]);
            return sum;
        }

        private static double Hellinger(double[] p, double[] q)
        {
            double sum = 0;
            for (int i = 0; i < p.Length; i++)
                sum += (Math.Sqrt(p[i]) - Math.Sqrt(q[i])) * (Math.Sqrt(p[i]) - Math.Sqrt(q[i]));
            return Math.Sqrt(sum / 2);
        }

        // Longest first, so at a locus the real gene usually wins over shorter antisense or nested ORFs. Returned in coding orientation.
        internal static List<string> CallGenes(string forward, string code, int minLength = minGeneLength)
        {
            var reverse = ProteinBackTranslationFeasibility.ReverseComplement(forward);
            int length = forward.Length;
            var candidates = ProteinBackTranslationFeasibility.AtgGenes(forward, code, minLength).Select(g => (g.Start, g.End, Reverse: false))
                .Concat(ProteinBackTranslationFeasibility.AtgGenes(reverse, code, minLength).Select(g => (Start: length - g.End, End: length - g.Start, Reverse: true)))
                .OrderByDescending(g => g.End - g.Start)
                .ToList();
            var covered = new bool[length];
            var genes = new List<string>();
            foreach (var (start, end, isReverse) in candidates)
            {
                int overlap = 0;
                for (int i = start; i < end; i++)
                    if (covered[i]) overlap++;
                if (overlap > maxOverlap) continue;
                for (int i = start; i < end; i++)
                    covered[i] = true;
                genes.Add(isReverse ? reverse.Substring(length - end, end - start) : forward.Substring(start, end - start));
            }
            return genes;
        }

        private static double Gc(string sequence)
        {
            long strong = 0, all = 0;
            foreach (var c in sequence)
            {
                if (c is 'G' or 'C') strong++;
                if (c is 'A' or 'C' or 'G' or 'T') all++;
            }
            return strong / (double)all;
        }

        private static double[] Pos12Mean(List<string> genes)
        {
            var sums = new double[6];
            long codons = 0;
            foreach (var gene in genes)
                for (int i = 0; i + 3 <= gene.Length; i += 3, codons++)
                    for (int position = 0; position < 2; position++)
                    {
                        var (x, y, z) = ProteinBackTranslationFeasibility.BaseAxes(gene[i + position]);
                        sums[3 * position] += x;
                        sums[3 * position + 1] += y;
                        sums[3 * position + 2] += z;
                    }
            return sums.Select(s => s / codons).ToArray();
        }

        private static double[] ExpectedPos12Mean(List<string> proteins, string code)
        {
            var table = ProteinBackTranslationFeasibility.ExpectedAxesTable(
                ProteinBackTranslationFeasibility.Synonyms(code), Enumerable.Repeat(1L, 64).ToArray());
            var sums = new double[6];
            long residues = 0;
            foreach (var protein in proteins)
                foreach (var aminoAcid in protein)
                {
                    residues++;
                    for (int position = 0; position < 2; position++)
                    {
                        var axes = table[2 * (aminoAcid - 'A') + position];
                        sums[3 * position] += axes.X;
                        sums[3 * position + 1] += axes.Y;
                        sums[3 * position + 2] += axes.Z;
                    }
                }
            return sums.Select(s => s / residues).ToArray();
        }

        private static double[] AminoAcidFrequencies(List<string> proteins)
        {
            var counts = new double[aminoAcidLetters.Length];
            foreach (var protein in proteins)
                foreach (var aminoAcid in protein)
                    counts[aminoAcidLetters.IndexOf(aminoAcid)]++;
            var total = counts.Sum();
            return counts.Select(c => c / total).ToArray();
        }

        private static double[] CodonFrequencies(List<string> genes)
        {
            var counts = new double[64];
            foreach (var gene in genes)
                for (int i = 0; i + 3 <= gene.Length; i += 3)
                    counts[ProteinBackTranslationFeasibility.CodonIndex(gene, i)]++;
            var total = counts.Sum();
            return counts.Select(c => c / total).ToArray();
        }

        // rho(xy) = f(xy) / (f1(x) f2(y)) over each codon's first two bases, with each position's own base frequencies as marginals.
        private static double[] Pos12Rho(List<string> genes)
        {
            const string bases = "ACGT";
            var pairs = new double[16];
            foreach (var gene in genes)
                for (int i = 0; i + 3 <= gene.Length; i += 3)
                    pairs[4 * bases.IndexOf(gene[i]) + bases.IndexOf(gene[i + 1])]++;
            var total = pairs.Sum();
            var first = Enumerable.Range(0, 4).Select(x => Enumerable.Range(0, 4).Sum(y => pairs[4 * x + y]) / total).ToArray();
            var second = Enumerable.Range(0, 4).Select(y => Enumerable.Range(0, 4).Sum(x => pairs[4 * x + y]) / total).ToArray();
            return Enumerable.Range(0, 16).Select(p => pairs[p] / total / (first[p / 4] * second[p % 4])).ToArray();
        }

        // Hao's composition vector (CVTree): a(w) = (f(w) - f0(w)) / f0(w), where f0 predicts w from its two (k-1)-mers and their
        // shared (k-2)-mer, so what remains is what the shorter words don't explain. Words don't cross protein boundaries.
        internal static double[] CvTreeVector(IReadOnlyList<string> texts, int k, string alphabet = aminoAcidLetters)
        {
            var words = KmerFrequencies(texts, k, alphabet);
            var shorter = KmerFrequencies(texts, k - 1, alphabet);
            var shortest = KmerFrequencies(texts, k - 2, alphabet);
            int suffixSpace = shorter.Length;
            int middleSpace = shortest.Length;
            var vector = new double[words.Length];
            for (int w = 0; w < words.Length; w++)
            {
                int prefix = w / alphabet.Length;
                int suffix = w % suffixSpace;
                int middle = prefix % middleSpace;
                double expected = shortest[middle] > 0 ? shorter[prefix] * shorter[suffix] / shortest[middle] : 0;
                vector[w] = expected > 0 ? (words[w] - expected) / expected : 0;
            }
            return vector;
        }

        // A word holding a letter outside the alphabet (N and other IUPAC codes in DNA) isn't counted.
        private static double[] KmerFrequencies(IReadOnlyList<string> texts, int k, string alphabet)
        {
            int space = 1;
            for (int i = 0; i < k; i++)
                space *= alphabet.Length;
            var counts = new double[space];
            long total = 0;
            foreach (var text in texts)
                for (int start = 0; start + k <= text.Length; start++)
                {
                    int word = 0;
                    for (int i = 0; i < k && word >= 0; i++)
                    {
                        int letter = alphabet.IndexOf(text[start + i]);
                        word = letter < 0 ? -1 : word * alphabet.Length + letter;
                    }
                    if (word < 0) continue;
                    counts[word]++;
                    total++;
                }
            for (int w = 0; w < space; w++)
                counts[w] /= total;
            return counts;
        }

        internal static double CosineDistance(double[] a, double[] b)
        {
            double dot = 0, normA = 0, normB = 0;
            for (int i = 0; i < a.Length; i++)
            {
                dot += a[i] * b[i];
                normA += a[i] * a[i];
                normB += b[i] * b[i];
            }
            return (1 - dot / Math.Sqrt(normA * normB)) / 2;
        }

        // A genome with no gene long enough gets the largest distance seen, as unrelated as anything here (S52c's convention).
        private static double[,] GeneBestHitMatrix(List<string>[] genes)
        {
            int n = genes.Length;
            int windowCount = bestHitGeneLength / 3 * 2 / ProteinBackTranslationFeasibility.pos12WindowSize;
            var windows = genes.Select(g => g.Where(x => x.Length >= bestHitGeneLength)
                .Select(x => ProteinBackTranslationFeasibility.Pos12Windows(x[..bestHitGeneLength])).ToArray()).ToArray();
            var matrix = new double[n, n];
            Parallel.For(0, n * n, cell =>
            {
                int a = cell / n, b = cell % n;
                if (b <= a) return;
                if (windows[a].Length == 0 || windows[b].Length == 0)
                {
                    matrix[a, b] = matrix[b, a] = double.NaN;
                    return;
                }
                var bestA = Enumerable.Repeat(double.NegativeInfinity, windows[a].Length).ToArray();
                var bestB = Enumerable.Repeat(double.NegativeInfinity, windows[b].Length).ToArray();
                for (int i = 0; i < windows[a].Length; i++)
                    for (int j = 0; j < windows[b].Length; j++)
                    {
                        var score = ProteinBackTranslationFeasibility.LeadingWindowsScore(windows[a][i], windows[b][j], windowCount);
                        if (score > bestA[i]) bestA[i] = score;
                        if (score > bestB[j]) bestB[j] = score;
                    }
                matrix[a, b] = matrix[b, a] = 1 - (bestA.Average() + bestB.Average()) / 2;
            });
            double worst = 0;
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (double.IsFinite(matrix[a, b])) worst = Math.Max(worst, matrix[a, b]);
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                    if (a != b && !double.IsFinite(matrix[a, b])) matrix[a, b] = worst;
            return matrix;
        }

        private static double MedianOffDiagonal(double[,] matrix)
        {
            var values = new List<double>();
            for (int a = 0; a < matrix.GetLength(0); a++)
                for (int b = a + 1; b < matrix.GetLength(0); b++)
                    values.Add(matrix[a, b]);
            values.Sort();
            return values[values.Count / 2];
        }

        internal static double[,] Matrix(int n, Func<int, int, double> distance)
        {
            var result = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    result[a, b] = result[b, a] = distance(a, b);
            return result;
        }
    }
}
