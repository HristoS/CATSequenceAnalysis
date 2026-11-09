using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// The gene layer for NCBI eukaryotes (S7.T): coding sequences from a CDS file and an assembly report, spliced genes on several
    /// sequences stored as ZCGP version 3, and a protein search that says which sequence its hit lies on.
    /// </summary>
    public class NcbiGeneLayerTests
    {
        private const string report =
            "# Assembly name:  Test1\n" +
            "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n" +
            "I\tassembled-molecule\tI\tChromosome\tCM000001.1\t=\tNC_000001.1\tPrimary Assembly\t500000\tna\n" +
            "II\tassembled-molecule\tII\tChromosome\tCM000002.1\t=\tNC_000002.1\tPrimary Assembly\t500000\tna\n" +
            "MT\tassembled-molecule\tMT\tMitochondrion\tna\t<>\tNC_000003.1\tnon-nuclear\t80000\tna\n";

        // A 600 aa gene spliced from two exons on the reverse strand of chromosome II, a shorter isoform of it, a 300 aa gene on
        // chromosome I, and one each of what the pass leaves out: partial, pseudo, ribosomal slippage, mitochondrial.
        private static readonly string spliced = ProteinSearchTests.RandomGene(600, seed: 1);
        private static readonly string isoform = ProteinSearchTests.RandomGene(400, seed: 2);
        private static readonly string plain = ProteinSearchTests.RandomGene(300, seed: 3);

        private static string CdsFile()
        {
            var other = ProteinSearchTests.RandomGene(200, seed: 4);
            string Record(string id, string attributes, string gene) => $">lcl|{id} {attributes} [gbkey=CDS]\n{gene}\n";
            return Record("NC_000002.1_cds_XP_1.1_1", "[gene=abc] [db_xref=GeneID:11] [protein_id=XP_1.1] [location=complement(join(1001..1900,3001..3903))]", spliced) +
                   Record("NC_000002.1_cds_XP_2.1_2", "[gene=abc] [db_xref=GeneID:11] [protein_id=XP_2.1] [location=1001..2203]", isoform) +
                   Record("NC_000001.1_cds_XP_3.1_3", "[locus_tag=T_3] [db_xref=GeneID:12] [protein_id=XP_3.1] [location=501..1403]", plain) +
                   Record("NC_000001.1_cds_XP_4.1_4", "[db_xref=GeneID:13] [partial=5'] [location=<1..603]", other) +
                   Record("NC_000001.1_cds_5", "[db_xref=GeneID:14] [pseudo=true] [location=2001..2603]", other) +
                   Record("NC_000001.1_cds_XP_6.1_6", "[db_xref=GeneID:15] [exception=ribosomal slippage] [location=4001..4603]", other) +
                   Record("NC_000003.1_cds_YP_7.1_7", "[db_xref=GeneID:16] [location=1..603]", other);
        }

        private static byte[] Gzip(string text)
        {
            using var output = new MemoryStream();
            using (var gzip = new GZipStream(output, CompressionMode.Compress))
                gzip.Write(Encoding.ASCII.GetBytes(text));
            return output.ToArray();
        }

        private static EukaryoteGenes Layer(out byte[] file)
        {
            file = Gzip(CdsFile());
            var nuclear = NcbiGeneLayer.NuclearSequences(new StringReader(report));
            var (kept, tally, sha256) = NcbiGeneLayer.Read(new MemoryStream(file), nuclear);
            return NcbiGeneLayer.Compute(kept, nuclear, ZCurveGeneProfile.standardCode, tally, sha256);
        }

        [Test]
        public void LocationsGiveTheSpanStrandAndExonLength()
        {
            Assert.That(NcbiGeneLayer.ParseLocation("complement(join(1001..1900,3001..3903))", out _), Is.EqualTo((1000, 3903, true, 1803L)));
            Assert.That(NcbiGeneLayer.ParseLocation("join(complement(5..10),complement(1..2))", out _), Is.EqualTo((0, 10, true, 8L)));
            Assert.That(NcbiGeneLayer.ParseLocation("501..1403", out _), Is.EqualTo((500, 1403, false, 903L)));
            Assert.That(NcbiGeneLayer.ParseLocation("join(1..5,complement(8..9))", out _), Is.Null, "two strands");
            Assert.That(NcbiGeneLayer.ParseLocation("join(NC_1.1:1..5,7..9)", out _), Is.Null, "another sequence");
            Assert.That(NcbiGeneLayer.ParseLocation("123", out _), Is.Null);
            Assert.That(NcbiGeneLayer.ParseLocation("<1..603", out bool partial), Is.Null);
            Assert.That(partial, Is.True);
        }

        [Test]
        public void TheReportListsNuclearSequencesByBothAccessions()
        {
            var nuclear = NcbiGeneLayer.NuclearSequences(new StringReader(report));
            Assert.That(nuclear, Is.EqualTo(new Dictionary<string, int>
            {
                ["CM000001.1"] = 0, ["NC_000001.1"] = 0, ["CM000002.1"] = 1, ["NC_000002.1"] = 1,
            }), "the mitochondrion is left out");
        }

        [Test]
        public void TheLongestIsoformOfEachNuclearGeneIsKept()
        {
            var layer = Layer(out var file);

            Assert.That(layer.CdsSha256, Is.EqualTo(SHA256.HashData(file)));
            var t = layer.Tally;
            Assert.That((t.Read, t.Isoforms, t.Partial, t.Pseudo, t.Exception, t.Organelle, t.Location, t.Frame), Is.EqualTo((7, 1, 1, 1, 1, 1, 0, 0)));
            Assert.That(layer.SequenceNames, Is.EqualTo(new[] { "NC_000001.1", "NC_000002.1" }));
            Assert.That(layer.Genes.Genes.Select(g => (g.Chromosome, g.Start, g.End, g.Reverse, g.Introns, g.Residues)),
                Is.EqualTo(new[] { (0, 500, 1403, false, 0, 300), (1, 1000, 3903, true, 1100, 600) }));
            Assert.That(layer.Proteins, Is.EqualTo(new[] { ZCurveGeneProfile.Translate(plain), ZCurveGeneProfile.Translate(spliced) }));
            Assert.That(layer.Genes.BinaryFormatVersion, Is.EqualTo(3));
        }

        [Test]
        public void OnlyTheCodesItCanReadAreAccepted()
        {
            Assert.That(NcbiGeneLayer.CodeTable(1), Is.EqualTo(ZCurveGeneProfile.standardCode));
            Assert.That(NcbiGeneLayer.CodeTable(12), Is.EqualTo(ZCurveGeneProfile.alternativeYeastCode));
            Assert.That(NcbiGeneLayer.CodeTable(6), Is.Null, "ciliates read TAA and TAG as Gln");
            Assert.That(NcbiGeneLayer.CodeTable(null), Is.Null);
        }

        [Test]
        public async Task StoredLayerIsSearchedAndNamesTheSequenceOfItsHit()
        {
            const string provider = "ncbi-gene-layer";
            var dataSource = ThrowawayDatabase.DataSource;
            long id = await ThrowawayDatabase.InsertRecordAsync(provider, "GCA_000000001.1", 1_000_000, kind: "eukaryote");
            var layer = Layer(out _);
            using var client = new NcbiGenomeClient(null);
            await new NcbiGeneLayer(dataSource, client).StoreAsync(id, layer,
                NcbiGeneLayer.LayerJson(layer, new NcbiAnnotation("Test annotation", "2026-01-01", "NCBI"), 1, "test://cds"));

            var stored = (await GeneProfiles.ReadAsync(dataSource, id))!;
            Assert.That(stored.Genes, Is.EqualTo(layer.Genes.Genes));
            // Each bucket row says its own version: the long gene is spliced on sequence 1, the short one plain on sequence 0.
            await using (var versions = dataSource.CreateCommand("select array_agg(format_version order by bucket) from gene_profile where sequence_id = $1"))
            {
                versions.Parameters.AddWithValue(id);
                Assert.That(await versions.ExecuteScalarAsync(), Is.EqualTo(new short[] { 3, 2 }));
            }
            Assert.That((await GeneProfiles.SequenceNamesAsync(dataSource, new[] { id }))[id], Is.EqualTo(new[] { "NC_000001.1", "NC_000002.1" }));

            var result = await ProteinSearch.SearchAsync(dataSource, ProteinQuery.FromProtein(layer.Proteins[1]), new CatalogueFilter(Provider: provider), top: 1);
            var top = result.Hits.Single();
            Assert.That((top.SequenceId, top.GeneSequence, top.GeneStart, top.GeneEnd, top.Reverse), Is.EqualTo((id, 1, 1000, 3903, true)));
        }
    }
}
