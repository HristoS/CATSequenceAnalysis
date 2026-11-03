using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    /// <summary>The eukaryote pilot: a genome becomes one record of nuclear counts, nothing of its DNA stays, and human assemblies are refused.</summary>
    public class NcbiIngestTests
    {
        private static readonly Lineage fungus = new("Eukaryota", null, "Fungi", "Ascomycota", "Saccharomycetes", "Saccharomycetales", "Saccharomycetaceae", "Saccharomyces", "Saccharomyces testus");

        private static NcbiAssembly Assembly(string accession, int taxid = 559292, long length = 12_000_000) =>
            new(accession, "R64", taxid, "Saccharomyces testus", length);

        private static byte[] Gzip(string text)
        {
            using var output = new MemoryStream();
            using (var gzip = new GZipStream(output, CompressionMode.Compress))
                gzip.Write(Encoding.ASCII.GetBytes(text));
            return output.ToArray();
        }

        [Test]
        public void CountsSumAcrossChromosomesAndLeaveOrganellesOut()
        {
            var fasta = ">NC_1 Saccharomyces chromosome I\nAACC\nGG\n>NC_2 Saccharomyces chromosome II\nTTNA\n" +
                        ">NC_3 Saccharomyces cerevisiae S288c mitochondrion, complete genome\nGGGGGGGG\n>NC_4 Arabidopsis chloroplast, complete genome\nCCCC\n";

            var counts = NcbiIngestor.CountNuclear(new StringReader(fasta));

            Assert.That(counts.Records, Is.EqualTo(2));
            Assert.That(counts.OrganellesLeftOut, Is.EqualTo(2));
            Assert.That(counts.LengthBp, Is.EqualTo(10));
            Assert.That(counts.Counts[0], Is.EqualTo(1), "AA");
            Assert.That(counts.Counts[1], Is.EqualTo(1), "AC");
            Assert.That(counts.Counts[5], Is.EqualTo(1), "CC");
            Assert.That(counts.Counts[6], Is.EqualTo(1), "CG");
            Assert.That(counts.Counts[10], Is.EqualTo(1), "GG, the organelle's seven left out");
            Assert.That(counts.Counts[15], Is.EqualTo(1), "TT; N breaks the pair");
            Assert.That(counts.Counts.Sum(), Is.EqualTo(6));
        }

        [Test]
        public void ReadingHashesTheFileAsDownloaded()
        {
            var file = Gzip(">NC_1 chromosome I\nACGTACGT\n");

            var (counts, sha256) = NcbiIngestor.Read(new MemoryStream(file));

            Assert.That(counts.LengthBp, Is.EqualTo(8));
            Assert.That(sha256, Is.EqualTo(SHA256.HashData(file)));
        }

        [Test]
        public void TheFastaUrlFollowsNcbisLayout()
        {
            Assert.That(NcbiGenomeClient.FastaUrl(Assembly("GCF_000146045.2")), Is.EqualTo(
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"));
            Assert.That(NcbiGenomeClient.FastaUrl(new NcbiAssembly("GCA_000001405.29", "GRCh38 p14/x", 1, "x", 1)), Does.EndWith(
                "GCA/000/001/405/GCA_000001405.29_GRCh38_p14_x/GCA_000001405.29_GRCh38_p14_x_genomic.fna.gz"));
        }

        [Test]
        public void PicksOneAssemblyPerFamilyWithinTheSizeLimitAndNeverHuman()
        {
            var lineages = new Dictionary<int, Lineage>
            {
                [1] = fungus,
                [2] = fungus with { Species = "Saccharomyces secundus" },
                [3] = fungus with { Family = "Debaryomycetaceae" },
                [4] = fungus with { Family = "Giganticaceae" },
                [9606] = fungus with { Family = "Hominidae", Species = "Homo sapiens" },
                [5] = fungus with { Family = "Otherhominidae", Species = "Homo sapiens" },
            };
            var assemblies = new[]
            {
                Assembly("GCF_1.1", 1), Assembly("GCF_2.1", 2), Assembly("GCF_3.1", 3),
                Assembly("GCF_4.1", 4, length: 5_000_000_000), Assembly("GCF_5.1", 9606), Assembly("GCF_6.1", 5),
            };

            var picked = NcbiIngestor.Pick(assemblies, lineages, 10, 2_000_000_000);

            Assert.That(picked.Select(p => p.Lineage.Family), Is.EquivalentTo(new[] { "Saccharomycetaceae", "Debaryomycetaceae" }));
            Assert.That(NcbiIngestor.Pick(assemblies, lineages, 10, 2_000_000_000).Select(p => p.Assembly.Accession), Is.EqualTo(picked.Select(p => p.Assembly.Accession)));
        }

        [Test]
        public void TheHumanGenomeIsRefusedEvenIfItIsAskedFor()
        {
            var ingestor = new NcbiIngestor(ThrowawayDatabase.DataSource, null);
            var human = new NuclearCounts(new long[16], 100, 1, 0);

            Assert.ThrowsAsync<ArgumentException>(() => ingestor.StoreAsync(Assembly("GCF_h.1", 9606), fungus, "test://h", new byte[32], 1, human));
            Assert.ThrowsAsync<ArgumentException>(() => ingestor.StoreAsync(Assembly("GCF_h.2", 1), fungus with { Species = "Homo sapiens" }, "test://h", new byte[32], 1, human));
        }

        [Test]
        public async Task AGenomeIsOneRecordOfCountsWithNoLevelsAndItsDnaIsNowhere()
        {
            var marker = Generator.GetRandomDNA(1000, seed: 77);
            var accession = $"GCF_{Guid.NewGuid():N}";
            await CrawlState.RegisterAsync(ThrowawayDatabase.DataSource, NcbiIngestor.Provider, "fungi", [accession]);
            var file = Gzip($">NC_1 chromosome I\n{marker}\n>NC_2 chromosome II\n{marker}\n>NC_3 mitochondrion, complete genome\n{marker}\n");
            var (counts, sha256) = NcbiIngestor.Read(new MemoryStream(file));

            await new NcbiIngestor(ThrowawayDatabase.DataSource, null).StoreAsync(Assembly(accession), fungus, "test://genome", sha256, file.Length, counts);

            await using var row = ThrowawayDatabase.DataSource.CreateCommand(
                "select r.id, r.kind, r.data_class, r.length_bp, r.dinucleotide_counts, (select count(*) from profile_level l where l.sequence_id = r.id), " +
                "position($2 in r::text) + position($2 in s::text), c.status " +
                "from sequence_record r join sequence_source s on s.id = r.source_id join crawl_item c on c.source_id = s.id where r.record_id = $1");
            row.Parameters.AddWithValue(accession);
            row.Parameters.AddWithValue(marker[..24]);
            await using var reader = await row.ExecuteReaderAsync();
            Assert.That(await reader.ReadAsync(), Is.True);
            Assert.That(reader.GetString(1), Is.EqualTo("eukaryote"));
            Assert.That(reader.GetString(2), Is.EqualTo("public_eukaryote"));
            Assert.That(reader.GetInt64(3), Is.EqualTo(2000), "two chromosomes, the organelle left out");
            Assert.That(reader.GetFieldValue<long[]>(4), Is.EqualTo(counts.Counts));
            Assert.That(reader.GetInt64(5), Is.EqualTo(0), "no profile levels");
            Assert.That(reader.GetInt32(6), Is.EqualTo(0), "none of the DNA in the record or its source");
            Assert.That(reader.GetString(7), Is.EqualTo(CrawlState.Done));
        }

        [Test]
        public async Task GenomesStoredInParallelAreAllKept()
        {
            var accessions = Enumerable.Range(0, 8).Select(_ => $"GCF_{Guid.NewGuid():N}").ToList();
            await CrawlState.RegisterAsync(ThrowawayDatabase.DataSource, NcbiIngestor.Provider, "fungi", accessions);
            var ingestor = new NcbiIngestor(ThrowawayDatabase.DataSource, null);
            var genome = new NuclearCounts(Enumerable.Repeat(5L, 16).ToArray(), 81, 1, 0);

            await Task.WhenAll(accessions.Select(a => ingestor.StoreAsync(Assembly(a), fungus, $"test://{a}", SHA256.HashData(Encoding.ASCII.GetBytes(a)), 1, genome)));

            await using var row = ThrowawayDatabase.DataSource.CreateCommand(
                "select count(distinct r.id), count(*) filter (where c.status = $2) " +
                "from sequence_record r join crawl_item c on c.source_id = r.source_id where r.record_id = any($1)");
            row.Parameters.AddWithValue(accessions.ToArray());
            row.Parameters.AddWithValue(CrawlState.Done);
            await using var reader = await row.ExecuteReaderAsync();
            await reader.ReadAsync();
            Assert.That(reader.GetInt64(0), Is.EqualTo(accessions.Count));
            Assert.That(reader.GetInt64(1), Is.EqualTo(accessions.Count));
        }

        [Test]
        public async Task AFailedGenomeIsDoneOnceARetryStoresIt()
        {
            var accession = $"GCF_{Guid.NewGuid():N}";
            await CrawlState.RegisterAsync(ThrowawayDatabase.DataSource, NcbiIngestor.Provider, "animals", [accession]);
            await CrawlState.MarkAsync(ThrowawayDatabase.DataSource, NcbiIngestor.Provider, accession, CrawlState.Failed, "IOException: EOF");
            Assert.That((await CrawlState.RegisterAsync(ThrowawayDatabase.DataSource, NcbiIngestor.Provider, "animals", [accession]))[accession], Is.EqualTo(CrawlState.Failed));

            await new NcbiIngestor(ThrowawayDatabase.DataSource, null).StoreAsync(Assembly(accession), fungus, $"test://{accession}",
                SHA256.HashData(Encoding.ASCII.GetBytes(accession)), 1, new NuclearCounts(Enumerable.Repeat(5L, 16).ToArray(), 81, 1, 0));

            Assert.That((await CrawlState.RegisterAsync(ThrowawayDatabase.DataSource, NcbiIngestor.Provider, "animals", [accession]))[accession], Is.EqualTo(CrawlState.Done));
        }

        [Test]
        public void ParallelCallersGetTheirOwnSlots()
        {
            using var http = new PoliteHttpClient("test", null, TimeSpan.FromSeconds(1));
            var now = new DateTime(2026, 11, 2, 0, 0, 0, DateTimeKind.Utc);
            var slots = new System.Collections.Concurrent.ConcurrentBag<DateTime>();

            Parallel.For(0, 100, _ => slots.Add(http.ReserveSlot(now)));

            Assert.That(slots.Order(), Is.EqualTo(Enumerable.Range(0, 100).Select(i => now.AddSeconds(i))));
        }

        [Test]
        public void EukaryotesAreTheirOwnGroup()
        {
            Assert.That(SignatureIndex.GroupOf("eukaryote"), Is.EqualTo(SignatureIndex.Eukaryotic));
        }

        [Test]
        public void ParsesAssembliesAcrossPages()
        {
            const string json = "{\"reports\":[{\"accession\":\"GCF_1.1\",\"assembly_info\":{\"assembly_name\":\"A1\"}," +
                "\"organism\":{\"tax_id\":7,\"organism_name\":\"Testus\"},\"assembly_stats\":{\"total_sequence_length\":\"1234\"}}],\"next_page_token\":\"more\"}";
            var assemblies = new List<NcbiAssembly>();

            var token = NcbiGenomeClient.ParseAssemblies(json, assemblies);

            Assert.That(token, Is.EqualTo("more"));
            Assert.That(assemblies.Single(), Is.EqualTo(new NcbiAssembly("GCF_1.1", "A1", 7, "Testus", 1234)));
            Assert.That(NcbiGenomeClient.ParseAssemblies("{\"reports\":[]}", assemblies), Is.Null);
        }
    }
}
