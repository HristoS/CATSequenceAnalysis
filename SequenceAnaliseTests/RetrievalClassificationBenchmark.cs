using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// §8's mean-similarity-and-margin reporting is not how alignment-free methods are conventionally
    /// benchmarked in the literature - a threshold-based retrieval/classification framing (ROC-AUC,
    /// precision/recall) is. The §8.6 corpus already carries the exact ground-truth relationship
    /// labels needed for this; only the analysis was missing. This reuses the same 30-genome corpus,
    /// taxonomy labels, and ZCurveWindowedProfile.Compare scores as
    /// <see cref="RealGenomeDistanceMatrixBenchmark"/>, framed as two binary retrieval questions at
    /// the two thresholds that split the four ordinal relationship levels into a meaningful positive
    /// class of reasonable size (SameSpecies alone has only one pair, too few for a stable AUC):
    ///
    ///  - "same genus or closer" (SameSpecies + SameGenus, n=8) vs. everything else (n=427)
    ///  - "same phylum or closer" (SameSpecies + SameGenus + SamePhylum, n=105) vs. DifferentPhylum
    ///    (n=330) - the same split ExternalBaselineComparisonBenchmark already uses for its
    ///    "ExcludingDifferentPhylum" correlation, so the two analyses partition the corpus the same way.
    ///
    /// Mash is scored at the identical thresholds on the identical pairs, which is the head-to-head
    /// comparison §8.7's mean-similarity table cannot make: Mash covers all 435 pairs (FastANI covers
    /// 14, too few for an AUC, so it is not included here). Note that a saturated score is not the same
    /// as an uninformative one for ranking purposes - Mash's floor-valued pairs are overwhelmingly the
    /// distant ones, so its AUC can stay high even where §8.7 shows it has no graded signal left; the
    /// two measurements answer different questions and are reported side by side rather than merged.
    ///
    /// AUC is computed via the Mann-Whitney-U / rank-sum identity (with average ranks for ties), which
    /// is exactly equivalent to the area under the ROC curve without needing to sweep thresholds
    /// explicitly: AUC = (rankSum(positives) - nPos*(nPos+1)/2) / (nPos*nNeg).
    /// </summary>
    public class RetrievalClassificationBenchmark
    {
        private static readonly GenomeMeta[] genomes = GenomeTaxonomy.Genomes;

        private static double ComputeAuc(List<double> positiveScores, List<double> negativeScores)
        {
            var combined = positiveScores.Select(s => (Score: s, IsPositive: true))
                .Concat(negativeScores.Select(s => (Score: s, IsPositive: false)))
                .OrderBy(x => x.Score)
                .ToList();

            var ranks = new double[combined.Count];
            int i = 0;
            while (i < combined.Count)
            {
                int j = i;
                while (j + 1 < combined.Count && combined[j + 1].Score == combined[i].Score) j++;
                double averageRank = (i + 1 + j + 1) / 2.0; // 1-based rank, averaged across the tie block
                for (int k = i; k <= j; k++) ranks[k] = averageRank;
                i = j + 1;
            }

            double rankSumPositive = 0;
            for (int k = 0; k < combined.Count; k++)
                if (combined[k].IsPositive) rankSumPositive += ranks[k];

            int nPos = positiveScores.Count;
            int nNeg = negativeScores.Count;
            return (rankSumPositive - nPos * (nPos + 1) / 2.0) / (nPos * (double)nNeg);
        }

        /// <summary>
        /// Sweeps every distinct score as a candidate "predict positive if score >= threshold" cutoff
        /// and returns the one maximizing F1, with its precision and recall. Reporting the F1-optimal
        /// operating point is the standard way to summarize a single point on the curve without
        /// hand-picking a cutoff; it is not a claim that this threshold generalizes beyond this corpus.
        /// </summary>
        private static (double Threshold, double Precision, double Recall, double F1) BestF1OperatingPoint(
            List<double> positiveScores, List<double> negativeScores)
        {
            var candidates = positiveScores.Concat(negativeScores).Distinct().OrderByDescending(s => s).ToList();
            double bestThreshold = double.NaN, bestPrecision = 0, bestRecall = 0, bestF1 = -1;

            foreach (var threshold in candidates)
            {
                int truePositives = positiveScores.Count(s => s >= threshold);
                int falsePositives = negativeScores.Count(s => s >= threshold);
                int falseNegatives = positiveScores.Count - truePositives;

                if (truePositives == 0) continue;
                double precision = truePositives / (double)(truePositives + falsePositives);
                double recall = truePositives / (double)(truePositives + falseNegatives);
                double f1 = 2 * precision * recall / (precision + recall);

                if (f1 > bestF1)
                    (bestF1, bestPrecision, bestRecall, bestThreshold) = (f1, precision, recall, threshold);
            }

            return (bestThreshold, bestPrecision, bestRecall, bestF1);
        }

        [Test]
        public void RetrievalAucAtTwoRelationshipThresholds()
        {
            var mash = ExternalBaselineData.ParseMashDistances(ExternalBaselineData.RawResultsPath("mash_dist_raw.tsv"));

            var profiles = genomes.Select(g =>
            {
                var primary = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName));
                return (Meta: g, Profile: new ZCurveWindowedProfile(primary.Sequence));
            }).ToList();

            var rows = new List<(RelationshipLevel Level, double ZCurve, double MashSimilarity)>();
            for (int a = 0; a < profiles.Count; a++)
                for (int b = a + 1; b < profiles.Count; b++)
                {
                    var score = ZCurveWindowedProfile.Compare(profiles[a].Profile, profiles[b].Profile);
                    var mashSimilarity = 1 - mash[ExternalBaselineData.PairKey(profiles[a].Meta.FileName, profiles[b].Meta.FileName)];
                    rows.Add((GenomeTaxonomy.Classify(profiles[a].Meta, profiles[b].Meta), score, mashSimilarity));
                }

            Assert.That(rows.Count, Is.EqualTo(genomes.Length * (genomes.Length - 1) / 2));

            bool SameGenusOrCloser(RelationshipLevel l) => l is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus;
            bool SamePhylumOrCloser(RelationshipLevel l) => l is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus or RelationshipLevel.SamePhylum;

            var thresholds = new (string Name, Func<RelationshipLevel, bool> IsPositive)[]
            {
                ("SameGenusOrCloser_vs_Rest", SameGenusOrCloser),
                ("SamePhylumOrCloser_vs_DifferentPhylum", SamePhylumOrCloser),
            };

            using var aucReport = File.CreateText("ReportRetrievalClassification.csv");
            aucReport.WriteLine("Threshold, Method, PositiveCount, NegativeCount, AUC");

            using var operatingPointReport = File.CreateText("ReportRetrievalOperatingPoint.csv");
            operatingPointReport.WriteLine("Threshold, Method, ScoreCutoff, Precision, Recall, F1");

            var zCurveAucs = new List<double>();
            foreach (var (name, isPositive) in thresholds)
            {
                var zPositive = rows.Where(r => isPositive(r.Level)).Select(r => r.ZCurve).ToList();
                var zNegative = rows.Where(r => !isPositive(r.Level)).Select(r => r.ZCurve).ToList();
                var mPositive = rows.Where(r => isPositive(r.Level)).Select(r => r.MashSimilarity).ToList();
                var mNegative = rows.Where(r => !isPositive(r.Level)).Select(r => r.MashSimilarity).ToList();

                var zAuc = ComputeAuc(zPositive, zNegative);
                var mAuc = ComputeAuc(mPositive, mNegative);
                zCurveAucs.Add(zAuc);

                aucReport.WriteLine($"{name}, ZCurveWindowedProfile, {zPositive.Count}, {zNegative.Count}, {zAuc}");
                aucReport.WriteLine($"{name}, Mash, {mPositive.Count}, {mNegative.Count}, {mAuc}");

                var zPoint = BestF1OperatingPoint(zPositive, zNegative);
                var mPoint = BestF1OperatingPoint(mPositive, mNegative);
                operatingPointReport.WriteLine($"{name}, ZCurveWindowedProfile, {zPoint.Threshold}, {zPoint.Precision}, {zPoint.Recall}, {zPoint.F1}");
                operatingPointReport.WriteLine($"{name}, Mash, {mPoint.Threshold}, {mPoint.Precision}, {mPoint.Recall}, {mPoint.F1}");

                TestContext.Out.WriteLine(
                    $"{name} (n={zPositive.Count}/{zNegative.Count}):  " +
                    $"AUC ZCurve={zAuc:F4} Mash={mAuc:F4}  |  " +
                    $"best-F1 ZCurve: cutoff={zPoint.Threshold:F4} P={zPoint.Precision:F3} R={zPoint.Recall:F3} F1={zPoint.F1:F3}  |  " +
                    $"Mash: cutoff={mPoint.Threshold:F4} P={mPoint.Precision:F3} R={mPoint.Recall:F3} F1={mPoint.F1:F3}");
            }

            // Central claim: ZCurveWindowedProfile's score ranks positives above negatives at a rate
            // meaningfully better than chance (AUC = 0.5) at both thresholds, on the full unselected
            // 435-pair corpus - a stronger bar than the mean-ordering assertions already made
            // elsewhere, since AUC is sensitive to the group overlap §8.6 documents rather than only
            // to the group means. No assertion is made about beating Mash: Mash's own AUC is reported
            // alongside because the honest finding is what it is, not because this method must win.
            foreach (var auc in zCurveAucs)
                Assert.That(auc, Is.GreaterThan(0.5));
        }
    }
}
