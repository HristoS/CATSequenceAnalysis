using System.Security.Cryptography;
using System.Text;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class SequenceFeatureReaderTests
    {
        private string directory = null!;

        [SetUp]
        public void CreateDirectory() => directory = Directory.CreateTempSubdirectory("seqfeatures_").FullName;

        [TearDown]
        public void DeleteDirectory() => Directory.Delete(directory, recursive: true);

        [Test]
        public void LineWidthLetterCaseAndLineEndingsGiveTheSameFeatures()
        {
            var dna = Generator.GetRandomDNA(10_007, 20261016);
            var variants = new[]
            {
                WriteFasta("w60.fna", dna, 60, "\n"),
                WriteFasta("w80_crlf.fna", dna, 80, "\r\n"),
                WriteFasta("w70_lower.fna", dna.ToLowerInvariant(), 70, "\n"),
                WriteFasta("one_line.fna", dna, int.MaxValue, "\n"),
            };

            var features = variants.Select(Features).ToList();
            foreach (var f in features.Skip(1))
            {
                Assert.That(f.Sha256, Is.EqualTo(features[0].Sha256));
                Assert.That(f.DinucleotideCounts, Is.EqualTo(features[0].DinucleotideCounts));
            }
            Assert.That(features[0].Sha256, Is.EqualTo(SHA256.HashData(Encoding.ASCII.GetBytes(dna))));
            Assert.That(features[0].DinucleotideCounts.Sum(), Is.EqualTo(dna.Length - 1));
        }

        [Test]
        public void AmbiguityCodesAreHashedAndBreakDinucleotides()
        {
            var path = WriteFasta("n.fna", "acgtNacg", 4, "\n");
            var (sha256, counts) = Features(path);

            Assert.That(sha256, Is.EqualTo(SHA256.HashData(Encoding.ASCII.GetBytes("ACGTNACG"))));
            var expected = new long[16];
            expected[Index("AC")] = 2; // ACGT|N|ACG: AC, CG, GT, then AC, CG
            expected[Index("CG")] = 2;
            expected[Index("GT")] = 1;
            Assert.That(counts, Is.EqualTo(expected));
        }

        [Test]
        public void HexamersCountWithinRunsOfBases()
        {
            // AAAAAAC|N|ACGTACG: AAAAAA and AAAAAC, then ACGTAC and CGTACG; the N run can't start a 6-mer.
            using var reader = new SequenceFeatureReader(new StringReader("aaaaaac\nNacgtacg"));
            while (reader.Read(new char[3], 0, 3) > 0) { }
            reader.Finish();
            var counts = reader.HexamerCounts;

            var expected = new int[SequenceFeatureReader.HexamerCount];
            foreach (var word in new[] { "AAAAAA", "AAAAAC", "ACGTAC", "CGTACG" })
                expected[word.Aggregate(0, (code, c) => 4 * code + SequenceFeatureReader.BaseOrder.IndexOf(c))]++;
            Assert.That(counts, Is.EqualTo(expected));
        }

        [Test]
        public void HexamersOfALongSequenceSumToItsRuns()
        {
            var dna = Generator.GetRandomDNA(10_007, 20261104);
            using var reader = new SequenceFeatureReader(new StringReader(dna[..5000] + "NN" + dna[5000..]));
            ZCurveWindowedProfile.FromStream(reader);
            reader.Finish();

            Assert.That(reader.HexamerCounts.Sum(), Is.EqualTo((5000 - 5) + (5007 - 5)));
            Assert.That(() => new SequenceFeatureReader(new StringReader("ACGT")).HexamerCounts, Throws.InvalidOperationException);
        }

        [Test]
        public void WrappingTheReaderLeavesTheProfileUnchanged()
        {
            var path = WriteFasta("p.fna", Generator.GetRandomDNA(25_000, 7), 60, "\n");

            using var plain = FastaReader.OpenSequenceReader(path, "test");
            using var wrapped = new SequenceFeatureReader(FastaReader.OpenSequenceReader(path, "test"));
            var expected = ZCurveWindowedProfile.FromStream(plain);
            var actual = ZCurveWindowedProfile.FromStream(wrapped);

            Assert.That(actual.ToBinary(), Is.EqualTo(expected.ToBinary()));
            Assert.That(wrapped.CharactersRead, Is.EqualTo(25_000));
        }

        [Test]
        public void FinishBeforeTheEndThrows()
        {
            using var reader = new SequenceFeatureReader(new StringReader("ACGTACGT"));
            reader.Read(new char[4], 0, 4);
            Assert.That(() => reader.Finish(), Throws.InvalidOperationException);
        }

        private static int Index(string pair) =>
            4 * SequenceFeatureReader.BaseOrder.IndexOf(pair[0]) + SequenceFeatureReader.BaseOrder.IndexOf(pair[1]);

        private static (byte[] Sha256, long[] DinucleotideCounts) Features(string path)
        {
            using var reader = new SequenceFeatureReader(FastaReader.OpenSequenceReader(path, "test"));
            ZCurveWindowedProfile.FromStream(reader);
            return reader.Finish();
        }

        private string WriteFasta(string name, string sequence, int lineWidth, string newline)
        {
            var text = new StringBuilder(">test description").Append(newline);
            for (int i = 0; i < sequence.Length; i += lineWidth)
                text.Append(sequence, i, Math.Min(lineWidth, sequence.Length - i)).Append(newline);
            var path = Path.Combine(directory, name);
            File.WriteAllText(path, text.ToString());
            return path;
        }
    }
}
