using SequenceAnalise.Benchmarks;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class WindowSizeBenchmarkTests
    {
        [Test]
        public void ExactQueryLocalizesAtTheCorrectWindowOffset()
        {
            var prefix = new string('A', 1_500);
            var unique = Generator.GetRandomDNA(8_000, seed: 4242);
            var suffix = new string('C', 1_500);
            var reference = prefix + unique + suffix;
            var trueStart = prefix.Length;
            var query = reference.Substring(trueStart, 2_000);

            var result = WindowSizeLocalizationBenchmark.EvaluateLocalization(reference, query, trueStart, 500, "ZCurve");

            Assert.That(result.ExactHit, Is.True);
            Assert.That(result.TrueHitRank, Is.EqualTo(1));
            Assert.That(result.BestOffsetBp, Is.EqualTo(trueStart));
            Assert.That(result.AbsoluteErrorBp, Is.EqualTo(0));
        }

        [Test]
        public void BenchmarkCsvIncludesMethodAndWindowColumns()
        {
            var directory = Path.Combine(Path.GetTempPath(), $"window_size_benchmark_{Guid.NewGuid():N}");
            Directory.CreateDirectory(directory);
            try
            {
                var manifestPath = Path.Combine(directory, "manifest.json");
                var outputPath = Path.Combine(directory, "window-size.csv");

                var sequence = string.Concat(Enumerable.Repeat("ACGT", 150_000));
                foreach (var fileName in PhaseOneDataset.ReferenceFileNames)
                {
                    using var file = File.Create(Path.Combine(directory, fileName));
                    using var gzip = new System.IO.Compression.GZipStream(file, System.IO.Compression.CompressionMode.Compress);
                    using var writer = new StreamWriter(gzip);
                    writer.WriteLine(">$fixture-record fixture genome");
                    writer.Write(sequence);
                }

                var manifest = PhaseOneDataset.Generate(directory, manifestPath, seed: 1234);
                WindowSizeLocalizationBenchmark.Run(directory, manifestPath, outputPath, new[] { 500 }, new[] { 2_000 });

                var lines = File.ReadAllLines(outputPath);
                Assert.That(lines, Is.Not.Empty);
                Assert.That(lines[0], Does.Contain("Method"));
                Assert.That(lines[0], Does.Contain("WindowSizeBp"));
                Assert.That(lines[0], Does.Contain("QueryLengthBp"));
                Assert.That(lines.Skip(1).Any(line => line.Contains("ZCurve") && line.Contains("500") && line.Contains("2000")), Is.True);
            }
            finally
            {
                Directory.Delete(directory, recursive: true);
            }
        }

        [Test]
        public void WindowSizeRobustnessBenchmarkWritesMultipleWindowSizes()
        {
            var directory = Path.Combine(Path.GetTempPath(), $"window_size_robustness_{Guid.NewGuid():N}");
            Directory.CreateDirectory(directory);
            try
            {
                var sequence = string.Concat(Enumerable.Repeat("ACGT", 20_000));
                var referenceFileName = "fixture_reference_GCF_000000001.1.fna.gz";
                var referencePath = Path.Combine(directory, referenceFileName);
                using (var file = File.Create(referencePath))
                using (var gzip = new System.IO.Compression.GZipStream(file, System.IO.Compression.CompressionMode.Compress))
                using (var writer = new StreamWriter(gzip))
                {
                    writer.WriteLine(">$fixture-record fixture genome");
                    writer.Write(sequence);
                }

                var queryLength = 2_000;
                var locusStart = 2_000;
                var manifest = new LocalizationDatasetManifest(
                    SchemaVersion: 1,
                    DatasetId: "phase1-exact-v1",
                    DatasetSeed: 4321,
                    WindowSizeBp: 500,
                    K: 4,
                    LociPerReference: 1,
                    CoordinateConvention: "0-based-half-open",
                    PrimaryQueryLengthsBp: new[] { queryLength },
                    StressQueryLengthsBp: Array.Empty<int>(),
                    References: new[]
                    {
                        new ReferenceManifest(
                            FileName: referenceFileName,
                            Accession: "GCF_000000001.1",
                            RecordId: "fixture-record",
                            RecordLengthBp: sequence.Length,
                            CompressedSha256: "unused",
                            SamplingSeed: 4321,
                            Loci: new[]
                            {
                                new LocusManifest(
                                    LocusId: "GCF_000000001.1:locus-01",
                                    StartBp: locusStart,
                                    Queries: new[]
                                    {
                                        new QueryManifest(
                                            QueryId: "GCF_000000001.1:locus-01:2000",
                                            Stratum: "Primary",
                                            LengthBp: queryLength,
                                            SequenceSha256: Convert.ToHexString(System.Security.Cryptography.SHA256.HashData(System.Text.Encoding.UTF8.GetBytes(sequence.Substring(locusStart, queryLength)))).ToLowerInvariant())
                                    })
                            })
                    });

                var manifestPath = Path.Combine(directory, "manifest.json");
                var outputPath = Path.Combine(directory, "window-size-sweep.csv");
                File.WriteAllText(manifestPath, System.Text.Json.JsonSerializer.Serialize(manifest));

                WindowSizeRobustnessBenchmark.Run(directory, manifestPath, outputPath, new[] { 250, 500 }, new[] { queryLength });

                var lines = File.ReadAllLines(outputPath);
                Assert.That(lines, Is.Not.Empty);
                Assert.That(lines[0], Does.Contain("WindowSizeBp"));
                Assert.That(lines.Skip(1).Any(line => line.Contains("ZCurve") && line.Contains("250") && line.Contains("2000")), Is.True);
            }
            finally
            {
                Directory.Delete(directory, recursive: true);
            }
        }
    }
}
