using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// §2.3's motivating claim for the whole manuscript: two unrelated bacteriophage genomes of
    /// similar composition score CATProfile.Compare ~ 0.99, indistinguishable in practice from a
    /// true near-duplicate. Until this test, that number had no backing anywhere in the repository
    /// - no phage sequence under data/, no test computing it, only a code comment
    /// (ZCurveProfileTests.UnrelatedRandomSequencesOfSameLengthScoreNearOne) referencing it (see
    /// docs/zcurve-windowed-profile-factcheck.md, discrepancy 3). This computes the real number
    /// against the real genomes, following the same real-data convention already used for §8.6's
    /// 30-genome corpus (RealGenomeDistanceMatrixBenchmark.cs).
    /// </summary>
    public class PhageMotivatingExampleTests
    {
        [Test]
        public void UnrelatedPhageGenomesScoreNearOneUnderCATProfile()
        {
            var phiX174 = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo("phix174_GCF_000819615.1.fna.gz"));
            var lambda = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo("lambda_phage_GCF_000840245.1.fna.gz"));

            var score = CATProfile.Compare(new CATProfile(phiX174.Sequence), new CATProfile(lambda.Sequence));

            TestContext.Out.WriteLine(
                $"phiX174 ({phiX174.Length} bp) vs lambda phage ({lambda.Length} bp), CATProfile.Compare: {score}");

            using (var report = File.CreateText("ReportPhageMotivatingExample.csv"))
            {
                report.WriteLine("GenomeA, LengthA, GenomeB, LengthB, CATProfileCompare");
                report.WriteLine($"phiX174, {phiX174.Length}, lambda phage, {lambda.Length}, {score}");
            }

            Assert.That(score, Is.GreaterThan(0.9),
                "two unrelated bacteriophage genomes of similar composition should still score close " +
                "to a near-duplicate under CATProfile's single-point compression - the motivating " +
                "false positive for the whole manuscript");
        }
    }
}
