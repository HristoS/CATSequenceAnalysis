using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    using P = ProteinBackTranslationFeasibility;

    /// <summary>
    /// Genome-profile check (2026-11-04, your question): would the per-CDS winner, 500 bp windows with Legendre moments, also improve the
    /// platform's whole-genome 500 bp profile? Genome windows aren't anchored to a gene start, so a fragment's windows sit at an arbitrary
    /// phase against the stored ones, and a trend can change under a shift where a mean barely does. Cells: axes x y z or x y, times means
    /// (P0) or moments (P0-P2, each scaled by sqrt(2k + 1)); all bases, forward strand, 500 bp windows from the record start. Scoring is
    /// BestAlignment's: slide the shorter profile along the longer, 1 - mean per-window Euclidean distance, best offset.
    /// Retrieval: the 30 bundled genomes and RetrievalClassificationBenchmark's two thresholds (xyz means must reproduce its AUCs 0.835 /
    /// 0.789). Localization: per bacterium, 20 fragments per condition, aligned (start a multiple of 500) or arbitrary phase, 4 or 10 kbp,
    /// exact or with 2% random substitutions, slid along the genome's own profile; a hit is a best offset within one window of the start.
    /// Rule fixed before the run: a cell replaces xyz means only if both retrieval AUCs are higher and its mean hit rate over the arbitrary-
    /// phase conditions is not more than 2 points lower; among qualifying cells, the cheapest (numbers per window) wins.
    /// Writes GenomeMoments_Retrieval.csv and GenomeMoments_Localization.csv.
    /// </summary>
    public class GenomeProfileMomentsBenchmark
    {
        private const int window = 500;
        private const int fragmentsPerCondition = 20;
        private static readonly (string Name, int[] Axes)[] axisSets = { ("xyz", new[] { 0, 1, 2 }), ("xy", new[] { 0, 1 }) };

        [Test]
        public void MomentsOnTheGenomeProfile()
        {
            var genomes = GenomeTaxonomy.Genomes;
            var sequences = genomes.Select(g => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence).ToArray();
            int n = genomes.Length;
            using var retrieval = File.CreateText("GenomeMoments_Retrieval.csv");
            using var localization = File.CreateText("GenomeMoments_Localization.csv");
            retrieval.WriteLine("Axes,Features,NumbersPerWindow,Threshold,Positives,Negatives,Auc");
            localization.WriteLine("Axes,Features,NumbersPerWindow,Phase,FragmentBp,Substitutions,Fragments,HitRate");
            var summary = new List<(string Cell, int Numbers, double GenusAuc, double PhylumAuc, double ArbitraryHits)>();

            foreach (var (axisName, axes) in axisSets)
                foreach (var orders in new[] { 1, 3 })
                {
                    var cell = $"{axisName},{(orders == 1 ? "means" : "moments")},{axes.Length * orders}";
                    var profiles = sequences.Select(s => Profile(s, axes, orders)).ToArray();

                    var scores = new double[n, n];
                    Parallel.For(0, n * n, pair =>
                    {
                        int a = pair / n, b = pair % n;
                        if (b <= a) return;
                        var (shorter, longer) = profiles[a][0].Length <= profiles[b][0].Length ? (profiles[a], profiles[b]) : (profiles[b], profiles[a]);
                        scores[a, b] = BestAlignment(shorter, longer).Score;
                    });
                    var aucs = new List<double>();
                    foreach (var (threshold, positive) in new (string, Func<RelationshipLevel, bool>)[]
                    {
                        ("SameGenusOrCloser_vs_Rest", level => level is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus),
                        ("SamePhylumOrCloser_vs_DifferentPhylum", level => level != RelationshipLevel.DifferentPhylum),
                    })
                    {
                        var positives = new List<double>();
                        var negatives = new List<double>();
                        for (int a = 0; a < n; a++)
                            for (int b = a + 1; b < n; b++)
                                (positive(GenomeTaxonomy.Classify(genomes[a], genomes[b])) ? positives : negatives).Add(scores[a, b]);
                        var auc = RetrievalClassificationBenchmark.ComputeAuc(positives, negatives);
                        aucs.Add(auc);
                        retrieval.WriteLine(string.Create(CultureInfo.InvariantCulture, $"{cell},{threshold},{positives.Count},{negatives.Count},{auc:F4}"));
                    }

                    var arbitrary = new List<double>();
                    foreach (var phase in new[] { "aligned", "arbitrary" })
                        foreach (var fragmentLength in new[] { 4_000, 10_000 })
                            foreach (var substitutions in new[] { 0.0, 0.02 })
                            {
                                var hits = new double[n];
                                Parallel.For(0, n, g =>
                                {
                                    var rng = new Random(P.StableSeed(genomes[g].FileName) ^ fragmentLength ^ (phase == "aligned" ? 1 : 2) ^ (int)(substitutions * 1000));
                                    int found = 0;
                                    for (int f = 0; f < fragmentsPerCondition; f++)
                                    {
                                        int start = rng.Next(0, sequences[g].Length - fragmentLength);
                                        if (phase == "aligned") start -= start % window;
                                        var fragment = Mutate(sequences[g].Substring(start, fragmentLength), substitutions, rng);
                                        var (_, offset) = BestAlignment(Profile(fragment, axes, orders), profiles[g]);
                                        if (Math.Abs(offset * window - start) <= window) found++;
                                    }
                                    hits[g] = found / (double)fragmentsPerCondition;
                                });
                                var rate = hits.Average();
                                if (phase == "arbitrary") arbitrary.Add(rate);
                                localization.WriteLine(string.Create(CultureInfo.InvariantCulture,
                                    $"{cell},{phase},{fragmentLength},{substitutions},{n * fragmentsPerCondition},{rate:F4}"));
                            }
                    summary.Add((cell, axes.Length * orders, aucs[0], aucs[1], arbitrary.Average()));
                    TestContext.Out.WriteLine(string.Create(CultureInfo.InvariantCulture,
                        $"{cell}: genus AUC {aucs[0]:F4}, phylum AUC {aucs[1]:F4}, arbitrary-phase hits {arbitrary.Average():F4}"));
                }

            var baseline = summary.Single(s => s.Cell == "xyz,means,3");
            var qualifying = summary.Where(s => s.Cell != baseline.Cell && s.GenusAuc > baseline.GenusAuc && s.PhylumAuc > baseline.PhylumAuc
                && s.ArbitraryHits >= baseline.ArbitraryHits - 0.02).OrderBy(s => s.Numbers).ToList();
            TestContext.Out.WriteLine(qualifying.Count > 0 ? $"Rule: replace xyz means with {qualifying[0].Cell}" : "Rule: keep xyz means");
            Assert.That(baseline.GenusAuc, Is.EqualTo(0.8352).Within(0.002), "xyz means must reproduce RetrievalClassificationBenchmark");
            Assert.That(baseline.PhylumAuc, Is.EqualTo(0.7891).Within(0.002), "xyz means must reproduce RetrievalClassificationBenchmark");
        }

        // Structure of arrays: features[f][w], with f = axes.Length * k + a and each value stored as float, as the profile stores means.
        private static float[][] Profile(string sequence, int[] axes, int orders)
        {
            int count = sequence.Length / window;
            var features = Enumerable.Range(0, axes.Length * orders).Select(_ => new float[count]).ToArray();
            var sums = new double[axes.Length * orders];
            for (int w = 0; w < count; w++)
            {
                Array.Clear(sums);
                for (int i = 0; i < window; i++)
                {
                    var (x, y, z) = Axes(sequence[w * window + i]);
                    double t = -1 + 2 * (i + 0.5) / window;
                    for (int k = 0; k < orders; k++)
                    {
                        double weight = k switch { 0 => 1, 1 => t, _ => (3 * t * t - 1) / 2 };
                        for (int a = 0; a < axes.Length; a++)
                            sums[axes.Length * k + a] += (axes[a] switch { 0 => x, 1 => y, _ => z }) * weight;
                    }
                }
                for (int k = 0; k < orders; k++)
                    for (int a = 0; a < axes.Length; a++)
                        features[axes.Length * k + a][w] = (float)(sums[axes.Length * k + a] * Math.Sqrt(2 * k + 1) / window);
            }
            return features;
        }

        // ZCurveProfile.Axes (internal to the library), with ambiguity codes contributing nothing as there.
        private static (int X, int Y, int Z) Axes(char @base) => @base switch
        {
            'A' or 'a' => (1, 1, -1),
            'C' or 'c' => (-1, 1, 1),
            'G' or 'g' => (1, -1, 1),
            'T' or 't' => (-1, -1, -1),
            _ => (0, 0, 0),
        };

        private static (double Score, int Offset) BestAlignment(float[][] query, float[][] reference)
        {
            int m = query[0].Length, n = reference[0].Length, width = query.Length;
            if (m == 0 || n < m) return (0, -1);
            double best = double.NegativeInfinity;
            int bestOffset = -1;
            for (int offset = 0; offset + m <= n; offset++)
            {
                double sum = 0;
                for (int w = 0; w < m; w++)
                {
                    double squared = 0;
                    for (int f = 0; f < width; f++)
                    {
                        double difference = query[f][w] - reference[f][offset + w];
                        squared += difference * difference;
                    }
                    sum += Math.Sqrt(squared);
                }
                double score = 1 - sum / m;
                if (score > best) (best, bestOffset) = (score, offset);
            }
            return (best, bestOffset);
        }

        private static string Mutate(string fragment, double rate, Random rng)
        {
            if (rate == 0) return fragment;
            var bases = fragment.ToCharArray();
            for (int i = 0; i < bases.Length; i++)
                if (rng.NextDouble() < rate)
                {
                    var others = "ACGT".Where(c => c != bases[i]).ToArray();
                    bases[i] = others[rng.Next(others.Length)];
                }
            return new string(bases);
        }
    }
}
