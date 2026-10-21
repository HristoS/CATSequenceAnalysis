using System.Text.Json;
using SequenceAnalise.Benchmarks;

namespace SequenceAnaliseTests
{
    public class DivergenceBenchmarkTests
    {
        [Test]
        public void DivergenceBenchmarkWritesMethodAndRateColumns()
        {
            var directory = Path.Combine(Path.GetTempPath(), $"divergence_benchmark_{Guid.NewGuid():N}");
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
                    PrimaryQueryLengthsBp: new[] { 2_000 },
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
                var outputPath = Path.Combine(directory, "divergence.csv");
                File.WriteAllText(manifestPath, JsonSerializer.Serialize(manifest));

                DivergenceLocalizationBenchmark.Run(directory, manifestPath, outputPath, windowSizeBp: 500, queryLengths: new[] { 2_000 }, divergenceRates: new[] { 0.00, 0.01, 0.05 });

                var lines = File.ReadAllLines(outputPath);
                Assert.That(lines, Is.Not.Empty);
                Assert.That(lines[0], Does.Contain("Method"));
                Assert.That(lines[0], Does.Contain("DivergenceRate"));
                Assert.That(lines.Skip(1).Any(line => line.Contains("ZCurve") && line.Contains("0.01")), Is.True);
            }
            finally
            {
                Directory.Delete(directory, recursive: true);
            }
        }
    }
}
