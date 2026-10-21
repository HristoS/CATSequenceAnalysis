using SequenceAnalise.Benchmarks;

namespace SequenceAnaliseTests
{
    public class BenchmarkSummaryTests
    {
        [Test]
        public void DivergenceSummaryAggregatesByMethodAndRate()
        {
            var directory = Path.Combine(Path.GetTempPath(), $"benchmark_summary_{Guid.NewGuid():N}");
            Directory.CreateDirectory(directory);
            try
            {
                var rawPath = Path.Combine(directory, "divergence.csv");
                File.WriteAllLines(rawPath, new[]
                {
                    "ReferenceFileName,Accession,LocusId,QueryId,Method,DivergenceRate,WindowSizeBp,QueryLengthBp,TrueStartBp,BestOffsetBp,AbsoluteErrorBp,ExactHit,TrueHitRank,TrueScore,BestFalseScore,ScoreMargin,BestScore",
                    "ref.fa,GCF_000000001.1,locus-01,query-1,ZCurve,0.01,500,2000,1000,1000,0,true,1,42,30,12,42",
                    "ref.fa,GCF_000000001.1,locus-01,query-2,ZCurve,0.01,500,2000,1000,1500,500,false,2,40,30,10,40",
                    "ref.fa,GCF_000000001.1,locus-01,query-3,WindowedKmer,0.01,500,2000,1000,1000,0,true,1,36,25,11,36",
                    "ref.fa,GCF_000000001.1,locus-01,query-4,WindowedKmer,0.01,500,2000,1000,1750,750,false,3,33,25,8,33",
                });

                var summaryPath = Path.Combine(directory, "divergence-summary.csv");
                BenchmarkSummary.WriteDivergenceSummary(rawPath, summaryPath);

                var lines = File.ReadAllLines(summaryPath);
                Assert.That(lines[0], Does.Contain("Method"));
                Assert.That(lines[0], Does.Contain("DivergenceRate"));
                Assert.That(lines[0], Does.Contain("Top1Accuracy"));
                Assert.That(lines.Skip(1).Any(line => line.Contains("ZCurve") && line.Contains("0.01") && line.Contains("0.5") || line.Contains("WindowedKmer") && line.Contains("0.01") && line.Contains("0.5")), Is.True);
            }
            finally
            {
                Directory.Delete(directory, recursive: true);
            }
        }
    }
}
