using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class FastaReaderTests
    {
        // Klebsiella pneumoniae MGH 78578: chromosome + 5 plasmids. Expected id/length pairs below
        // were verified independently (zcat | awk over the raw file, outside this codebase) before
        // being hardcoded here, so this is a check against known-correct ground truth, not just
        // "whatever the parser happens to return".
        private static string KlebsiellaPath =>
            GenomeCatalog.PathTo("klebsiella_pneumoniae_mgh78578_GCF_000016305.1.fna.gz");

        [Test]
        public void ReadRecords_ParsesAllContigsWithCorrectIdsAndLengths()
        {
            var records = FastaReader.ReadRecords(KlebsiellaPath).ToList();

            var expected = new (string Id, int Length)[]
            {
                ("NC_009648.1", 5_315_120), // chromosome
                ("NC_009649.1", 175_879),   // pKPN3
                ("NC_009650.1", 107_576),   // pKPN4
                ("NC_009651.1", 88_582),    // pKPN5
                ("NC_009652.1", 4_259),     // pKPN6
                ("NC_009653.1", 3_478),     // pKPN7
            };

            Assert.That(records.Select(r => (r.Id, r.Length)), Is.EqualTo(expected));
        }

        [Test]
        public void ReadHeaders_MatchesReadRecordsIdsAndDescriptions()
        {
            var expected = FastaReader.ReadRecords(KlebsiellaPath).Select(r => (r.Id, r.Description)).ToList();

            Assert.That(FastaReader.ReadHeaders(KlebsiellaPath), Is.EqualTo(expected));
        }

        [Test]
        public void ReadRecords_SequenceContainsOnlyUppercaseLettersAndMatchesDeclaredLength()
        {
            foreach (var record in FastaReader.ReadRecords(KlebsiellaPath))
            {
                Assert.That(record.Sequence.Length, Is.EqualTo(record.Length));
                Assert.That(record.Sequence.All(c => c is >= 'A' and <= 'Z'), Is.True,
                    $"record {record.Id} contains a character outside A-Z after uppercasing");
            }
        }

        [Test]
        public void ReadPrimaryRecord_PicksTheChromosomeNotAPlasmid()
        {
            var primary = FastaReader.ReadPrimaryRecord(KlebsiellaPath);

            Assert.That(primary.Id, Is.EqualTo("NC_009648.1"));
            Assert.That(primary.Length, Is.EqualTo(5_315_120));
        }

        [Test]
        public void ReadRecords_HandlesManyContigsInOneFile()
        {
            // Borreliella burgdorferi B31: 1 chromosome + 21 plasmids, the most fragmented
            // assembly in data/genomes/ - a stress case for the header/body state machine.
            var path = GenomeCatalog.PathTo("borrelia_burgdorferi_b31_GCF_000008685.2.fna.gz");
            var records = FastaReader.ReadRecords(path).ToList();

            Assert.That(records.Count, Is.EqualTo(22));
            Assert.That(records.Sum(r => r.Length), Is.EqualTo(1_521_208),
                "total bp across every contig should match data/genomes/README.md's reported total");
        }

        [Test]
        public void OpenSequenceReader_MatchesReadRecordsCharacterForCharacter()
        {
            // Mirrors GenomeScaleStreamingTests.FromStreamMatchesStringConstructor's convention:
            // the streaming and materializing code paths must agree bit-for-bit, not just "look
            // similar". Mycoplasmoides genitalium (580 kbp, single record) keeps this test fast.
            var path = GenomeCatalog.PathTo("mycoplasma_genitalium_g37_GCF_000027325.1.fna.gz");
            var materialized = FastaReader.ReadPrimaryRecord(path);

            using var streamReader = FastaReader.OpenSequenceReader(path, materialized.Id);
            var streamed = streamReader.ReadToEnd();

            Assert.That(streamed, Is.EqualTo(materialized.Sequence));
        }

        [Test]
        public void OpenSequenceReader_StopsAtTheRequestedRecordEvenWhenOthersFollow()
        {
            // The chromosome (NC_009648.1) is listed first in this file, with 5 plasmid records
            // after it - OpenSequenceReader must stop at the first '>' following the target
            // record's body, not read past it into pKPN3's sequence.
            using var reader = FastaReader.OpenSequenceReader(KlebsiellaPath, "NC_009648.1");
            var streamed = reader.ReadToEnd();

            Assert.That(streamed.Length, Is.EqualTo(5_315_120));
        }

        [Test]
        public void ReadRecordStreams_MatchesReadRecordsInOnePass()
        {
            var expected = FastaReader.ReadRecords(KlebsiellaPath).Select(r => (r.Id, r.Description, r.Sequence)).ToList();

            using var file = new System.IO.Compression.GZipStream(File.OpenRead(KlebsiellaPath), System.IO.Compression.CompressionMode.Decompress);
            using var text = new StreamReader(file);
            var streamed = FastaReader.ReadRecordStreams(text).Select(r => (r.Id, r.Description, r.Sequence.ReadToEnd())).ToList();

            Assert.That(streamed, Is.EqualTo(expected));
        }

        [Test]
        public void ReadRecordStreams_DropsLineBreaksOfEveryKindAndUppercases()
        {
            var oneLine = new string('a', 100_000) + "cgt";
            var fasta = "\r\n>first some description\r\nacGT\r\nNNry\r\n\r\n>second\nGGGG\n" + $">third\n{oneLine}";

            var records = FastaReader.ReadRecordStreams(new StringReader(fasta)).Select(r => (r.Id, r.Description, r.Sequence.ReadToEnd())).ToList();

            Assert.That(records, Is.EqualTo(new[]
            {
                ("first", "some description", "ACGTNNRY"),
                ("second", "", "GGGG"),
                ("third", "", oneLine.ToUpperInvariant()),
            }));
        }

        [Test]
        public void ReadRecordStreams_SkipsAnUnreadBodyAndItsOldReaderReadsNothing()
        {
            var fasta = ">a\nAAAA\nCCCC\n>b\nGGGG\n>c\nTTTT\n";
            var readers = new List<TextReader>();
            var sequences = new List<string>();

            foreach (var (id, _, sequence) in FastaReader.ReadRecordStreams(new StringReader(fasta)))
            {
                readers.Add(sequence);
                if (id != "a") sequences.Add(sequence.ReadToEnd());
            }

            Assert.That(sequences, Is.EqualTo(new[] { "GGGG", "TTTT" }), "record a was never read, and b didn't start inside it");
            Assert.That(readers[0].ReadToEnd(), Is.Empty);
        }

        [Test]
        public void ReadRecordStreams_RejectsSequenceBeforeTheFirstHeader()
        {
            Assert.That(() => FastaReader.ReadRecordStreams(new StringReader("ACGT\n>a\nACGT\n")).ToList(), Throws.InstanceOf<FormatException>());
            Assert.That(FastaReader.ReadRecordStreams(new StringReader("")), Is.Empty);
        }
    }
}
