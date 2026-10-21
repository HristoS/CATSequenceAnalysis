using System.IO.Compression;
using System.Text;
using SequenceAnalise.Benchmarks;

namespace SequenceAnaliseTests
{
    public class PhaseOneDatasetTests
    {
        private const int fixtureLength = 600_000;

        [Test]
        public void GenerationIsDeterministicAndStoresOnlyQueryRecipes()
        {
            WithFixture((directory, fileName, sequence) =>
            {
                var firstPath = Path.Combine(directory, "first.json");
                var secondPath = Path.Combine(directory, "second.json");

                var first = PhaseOneDataset.Generate(directory, firstPath, seed: 1234);
                PhaseOneDataset.Generate(directory, secondPath, seed: 1234);

                Assert.That(File.ReadAllBytes(firstPath), Is.EqualTo(File.ReadAllBytes(secondPath)));
                Assert.That(first.References, Has.Length.EqualTo(30));
                Assert.That(first.References[0].Loci, Has.Length.EqualTo(20));
                Assert.That(first.References.SelectMany(r => r.Loci).SelectMany(l => l.Queries).Count(), Is.EqualTo(3_000));
                Assert.That(File.ReadAllText(firstPath), Does.Not.Contain(sequence));
                Assert.That(PhaseOneDataset.Validate(directory, firstPath), Is.EqualTo((30, 3_000)));
            });
        }

        [Test]
        public void ValidationRejectsChangedReferenceBytes()
        {
            WithFixture((directory, fileName, _) =>
            {
                var manifestPath = Path.Combine(directory, "manifest.json");
                PhaseOneDataset.Generate(directory, manifestPath, seed: 5678);
                WriteReference(Path.Combine(directory, fileName), new string('A', fixtureLength));

                Assert.Throws<InvalidDataException>(() => PhaseOneDataset.Validate(directory, manifestPath));
            });
        }

        [Test]
        public void DefaultReferenceSetIsExactlyTheBacterialCorpus()
        {
            Assert.That(PhaseOneDataset.ReferenceFileNames, Has.Count.EqualTo(30));
            Assert.That(PhaseOneDataset.ReferenceFileNames, Does.Not.Contain("phix174_GCF_000819615.1.fna.gz"));
            Assert.That(PhaseOneDataset.ReferenceFileNames, Does.Not.Contain("lambda_phage_GCF_000840245.1.fna.gz"));
        }

        private static void WithFixture(Action<string, string, string> action)
        {
            var directory = Path.Combine(Path.GetTempPath(), $"phase1_dataset_{Guid.NewGuid():N}");
            Directory.CreateDirectory(directory);
            try
            {
                var sequence = string.Concat(Enumerable.Repeat("ACGT", fixtureLength / 4));
                foreach (var referenceFileName in PhaseOneDataset.ReferenceFileNames)
                    WriteReference(Path.Combine(directory, referenceFileName), sequence);
                var fileName = PhaseOneDataset.ReferenceFileNames[0];
                action(directory, fileName, sequence);
            }
            finally
            {
                Directory.Delete(directory, recursive: true);
            }
        }

        private static void WriteReference(string path, string sequence)
        {
            using var file = File.Create(path);
            using var gzip = new GZipStream(file, CompressionMode.Compress);
            using var writer = new StreamWriter(gzip, new UTF8Encoding(false));
            writer.WriteLine(">fixture-record fixture genome");
            writer.Write(sequence);
        }
    }
}