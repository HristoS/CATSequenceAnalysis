using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// P5 S5.2 (R3) on the bundled genomes: neighbour-joining and UPGMA trees built from Z-curve distances only, so what the
    /// trees show is the profile's own signal; GC, GC + length and Mash trees are built beside them as context. Distances:
    ///  - compare: 1 - Compare of the profiles as read.
    ///  - strandOrigin: 1 - the best Compare over both strands (exact, from the reverse complement) and, for bacterial chromosomes
    ///    with a clear GC skew (OriginRotation's rule in SequencePlatform.Core), rotation to the replication origin. This is how
    ///    the platform compares.
    ///  - summary: Euclidean distance of standardized per-genome window statistics (|mean x|, |mean y|, mean z, SD and lag-1
    ///    autocorrelation per axis), a fixed-length alternative that is a metric.
    /// Bacteria at 32 kbp and 500 bp; viruses (no rotation) at the per-pair ladder level, the coarsest keeping 8 windows.
    /// Per tree: clades recovered per rank (taxa with at least 2 leaves), mean consistency, and a 1,000-permutation null of the
    /// labels over the same tree. Writes S52_Trees_Bundled.csv, S52_LeafOrder_Bundled.csv and trees/S52_*.nwk. About 2 minutes.
    /// </summary>
    public class ZCurveTreeBenchmark
    {
        private const int permutations = 1000;
        // OriginRotation (SequencePlatform.Core) rotates only chromosomes it can place an origin on; the same rule here.
        private const long minRotationLength = 500_000;
        private const double minOriginToTerminusFraction = 0.4, minSkewBpPerKbp = 3;

        [Test]
        public void TreesFromBundledGenomes()
        {
            Directory.CreateDirectory("trees");
            using var summary = File.CreateText("S52_Trees_Bundled.csv");
            summary.WriteLine("Corpus, Level, Distance, Method, Rank, Taxa, Monophyletic, MeanConsistency, NullMeanMonophyletic, NullP, NjVsUpgmaNrf");
            using var order = File.CreateText("S52_LeafOrder_Bundled.csv");
            order.WriteLine("Corpus, Level, Distance, Method, Position, Organism, Genus, Higher");

            Bacteria(summary, order);
            Viruses(summary, order);
            Assert.Pass();
        }

        private static void Bacteria(StreamWriter summary, StreamWriter order)
        {
            var genomes = GenomeTaxonomy.Genomes;
            int n = genomes.Length;
            var versions = new ZCurveWindowedProfile[n][];
            var gc = new double[n];
            var length = new double[n];
            for (int i = 0; i < n; i++)
            {
                var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genomes[i].FileName)).Sequence;
                (gc[i], length[i]) = (Gc(sequence), sequence.Length);
                var forward = new ZCurveWindowedProfile(sequence);
                var reverse = new ZCurveWindowedProfile(ReverseComplement(sequence));
                var list = new List<ZCurveWindowedProfile> { forward, reverse };
                if (ToOrigin(forward, sequence.Length) is { } turned && ToOrigin(reverse, sequence.Length) is { } turnedReverse)
                {
                    list.Add(turned);
                    list.Add(turnedReverse);
                }
                versions[i] = list.ToArray();
            }
            TestContext.Out.WriteLine($"bacteria: {versions.Count(v => v.Length > 2)} of {n} rotated to their origin");

            var labels = genomes.Select(g => $"{g.OrganismName}|{g.Phylum}").ToArray();
            var ranks = new (string Name, string[] Taxa)[] { ("genus", genomes.Select(g => g.Genus).ToArray()), ("phylum", genomes.Select(g => g.Phylum).ToArray()) };
            void Order(string level, string distance, string method, PhyloTree tree)
            {
                int position = 0;
                foreach (int leaf in tree.LeafOrder())
                    order.WriteLine($"bacteria, {level}, {distance}, {method}, {position++}, \"{genomes[leaf].OrganismName}\", {genomes[leaf].Genus}, {genomes[leaf].Phylum}");
            }

            foreach (var (level, factor) in new[] { ("32kbp", 64), ("500bp", 1) })
            {
                var atLevel = versions.Select(v => v.Select(p => factor == 1 ? p : p.Downsample(factor)).ToArray()).ToArray();
                var plain = new double[n, n];
                var best = new double[n, n];
                Parallel.For(0, n * n, cell =>
                {
                    int a = cell / n, b = cell % n;
                    if (b <= a) return;
                    plain[a, b] = plain[b, a] = 1 - ZCurveWindowedProfile.Compare(atLevel[a][0], atLevel[b][0]);
                    double top = double.NegativeInfinity;
                    foreach (var x in atLevel[a])
                        foreach (var y in atLevel[b])
                            top = Math.Max(top, ZCurveWindowedProfile.Compare(x, y));
                    best[a, b] = best[b, a] = 1 - top;
                });
                var features = Standardize(atLevel.Select(v => Features(v[0])).ToArray());
                Evaluate("bacteria", level, "compare", plain, labels, ranks, summary, Order);
                Evaluate("bacteria", level, "strandOrigin", best, labels, ranks, summary, Order);
                Evaluate("bacteria", level, "summary", Euclidean(features), labels, ranks, summary, Order);
            }

            var mash = ExternalBaselineData.ParseMashDistances(ExternalBaselineData.RawResultsPath("mash_dist_raw.tsv"));
            var mashDistances = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                    if (a != b) mashDistances[a, b] = mash[ExternalBaselineData.PairKey(genomes[a].FileName, genomes[b].FileName)];
            Evaluate("bacteria", "-", "mash", mashDistances, labels, ranks, summary, Order);
            Evaluate("bacteria", "-", "gc", Euclidean(gc.Select(v => new[] { v }).ToArray()), labels, ranks, summary, Order);
            Evaluate("bacteria", "-", "gcLength", Euclidean(Standardize(gc.Select((v, i) => new[] { v, Math.Log10(length[i]) }).ToArray())),
                labels, ranks, summary, Order);
        }

        private static void Viruses(StreamWriter summary, StreamWriter order)
        {
            var genomes = ViralTaxonomy.Genomes;
            int n = genomes.Length;
            var forward = new ZCurveWindowedProfile[n];
            var ladders = new (ZCurveProfileLadder Forward, ZCurveProfileLadder Reverse)[n];
            var gc = new double[n];
            var length = new double[n];
            for (int i = 0; i < n; i++)
            {
                var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genomes[i].FileName)).Sequence;
                (gc[i], length[i]) = (Gc(sequence), sequence.Length);
                forward[i] = new ZCurveWindowedProfile(sequence);
                ladders[i] = (new ZCurveProfileLadder(forward[i]), new ZCurveProfileLadder(new ZCurveWindowedProfile(ReverseComplement(sequence))));
            }

            // Not comparable (fewer than 8 windows at 500 bp) counts as unrelated, as Compare's 0 does.
            var plain = new double[n, n];
            var best = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                {
                    double same = ZCurveProfileLadder.Compare(ladders[a].Forward, ladders[b].Forward) ?? 0;
                    double flipped = ZCurveProfileLadder.Compare(ladders[a].Reverse, ladders[b].Forward) ?? 0;
                    plain[a, b] = plain[b, a] = 1 - same;
                    best[a, b] = best[b, a] = 1 - Math.Max(same, flipped);
                }

            var labels = genomes.Select(g => $"{g.Name}|{g.Family}").ToArray();
            var ranks = new (string Name, string[] Taxa)[] { ("genus", genomes.Select(g => g.Genus).ToArray()), ("family", genomes.Select(g => g.Family).ToArray()) };
            void Order(string level, string distance, string method, PhyloTree tree)
            {
                int position = 0;
                foreach (int leaf in tree.LeafOrder())
                    order.WriteLine($"viruses, {level}, {distance}, {method}, {position++}, \"{genomes[leaf].Name}\", {genomes[leaf].Genus}, {genomes[leaf].Family}");
            }

            Evaluate("viruses", "ladder", "compare", plain, labels, ranks, summary, Order);
            Evaluate("viruses", "ladder", "strandOrigin", best, labels, ranks, summary, Order);
            Evaluate("viruses", "500bp", "summary", Euclidean(Standardize(forward.Select(Features).ToArray())), labels, ranks, summary, Order);
            Evaluate("viruses", "-", "gc", Euclidean(gc.Select(v => new[] { v }).ToArray()), labels, ranks, summary, Order);
            Evaluate("viruses", "-", "gcLength", Euclidean(Standardize(gc.Select((v, i) => new[] { v, Math.Log10(length[i]) }).ToArray())),
                labels, ranks, summary, Order);
        }

        private static void Evaluate(string corpus, string level, string distance, double[,] distances, string[] labels,
            (string Name, string[] Taxa)[] ranks, StreamWriter summary, Action<string, string, string, PhyloTree> order)
        {
            var trees = new[] { ("NJ", Phylogeny.NeighbourJoining(distances, labels)), ("UPGMA", Phylogeny.Upgma(distances, labels)) };
            double disagreement = trees[0].Item2.NormalizedRobinsonFoulds(trees[1].Item2);
            foreach (var (method, tree) in trees)
            {
                File.WriteAllText(Path.Combine("trees", $"S52_{corpus}_{level}_{distance}_{method}.nwk"), tree.ToNewick() + Environment.NewLine);
                order(level, distance, method, tree);
                foreach (var (rank, taxa) in ranks)
                {
                    var recovery = tree.Recovery(taxa);
                    int observed = recovery.Count(r => r.Monophyletic);
                    var rng = new Random(20261029);
                    var nullCounts = Enumerable.Range(0, permutations)
                        .Select(_ => tree.Recovery(taxa.OrderBy(_ => rng.Next()).ToArray()).Count(r => r.Monophyletic)).ToList();
                    double p = (1 + nullCounts.Count(c => c >= observed)) / (permutations + 1.0);
                    string line = string.Create(CultureInfo.InvariantCulture,
                        $"{corpus}, {level}, {distance}, {method}, {rank}, {recovery.Count}, {observed}, {recovery.Average(r => r.Consistency):F3}, " +
                        $"{nullCounts.Average():F2}, {p:F4}, {disagreement:F3}");
                    summary.WriteLine(line);
                    TestContext.Out.WriteLine(line);
                }
            }
        }

        private static ZCurveWindowedProfile? ToOrigin(ZCurveWindowedProfile fine, long lengthBp)
        {
            if (lengthBp < minRotationLength || fine.WindowCount == 0) return null;
            var skew = fine.GcSkewExtremes();
            int count = fine.WindowCount, apart = Math.Abs(skew.OriginWindow - skew.TerminusWindow);
            if ((double)Math.Min(apart, count - apart) / count < minOriginToTerminusFraction) return null;
            if (skew.AmplitudeBp / (lengthBp / 1000.0) < minSkewBpPerKbp) return null;
            return fine.RotateTo(skew.OriginWindow);
        }

        private static double[] Features(ZCurveWindowedProfile level)
        {
            var result = new List<double>();
            foreach (var (values, signed) in new[] { (level.Windows.Select(w => w.X).ToArray(), false),
                (level.Windows.Select(w => w.Y).ToArray(), false), (level.Windows.Select(w => w.Z).ToArray(), true) })
            {
                double mean = values.Average();
                double variance = values.Sum(v => (v - mean) * (v - mean)) / values.Length;
                double lag1 = 0;
                for (int i = 1; i < values.Length; i++) lag1 += (values[i] - mean) * (values[i - 1] - mean);
                result.Add(signed ? mean : Math.Abs(mean));   // x and y change sign with the strand; z is GC and doesn't
                result.Add(Math.Sqrt(variance));
                result.Add(variance == 0 ? 0 : lag1 / ((values.Length - 1) * variance));
            }
            return result.ToArray();
        }

        internal static double[][] Standardize(double[][] rows)
        {
            int width = rows[0].Length;
            var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
            var sds = Enumerable.Range(0, width).Select(f => Math.Sqrt(rows.Average(r => (r[f] - means[f]) * (r[f] - means[f])))).ToArray();
            return rows.Select(r => r.Select((v, f) => sds[f] == 0 ? 0 : (v - means[f]) / sds[f]).ToArray()).ToArray();
        }

        internal static double[,] Euclidean(double[][] points)
        {
            int n = points.Length;
            var result = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    result[a, b] = result[b, a] = Math.Sqrt(points[a].Zip(points[b], (p, q) => (p - q) * (p - q)).Sum());
            return result;
        }

        internal static double Gc(string sequence)
        {
            long strong = 0, all = 0;
            foreach (char c in sequence)
                switch (char.ToUpperInvariant(c))
                {
                    case 'G' or 'C': strong++; all++; break;
                    case 'A' or 'T': all++; break;
                }
            return (double)strong / all;
        }

        // IUPAC codes other than the four bases are kept unchanged; both profile types score them as zero anyway.
        internal static string ReverseComplement(string sequence)
        {
            var chars = new char[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
                chars[sequence.Length - 1 - i] = sequence[i] switch
                {
                    'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
                    'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c',
                    var other => other,
                };
            return new string(chars);
        }
    }
}
