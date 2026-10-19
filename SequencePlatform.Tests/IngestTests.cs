using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    public class IngestTests
    {
        // Two small BV-BRC viruses (about 4 kbp each) with different contig accessions.
        private static readonly string[] genomeIds = { "1000373.10", "1000373.4" };

        [Test, Explicit("Downloads two small genomes from BV-BRC; run it by name.")]
        public async Task IngestLeavesNoFileBehindOnSuccessAndOnFailure()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var before = DownloadedFiles();
            using var client = new BvBrcClient(null);
            var genomes = await client.GetGenomesAsync(genomeIds);
            Assert.That(genomes.Select(g => g.GenomeId), Is.EquivalentTo(genomeIds));
            await CrawlState.RegisterAsync(dataSource, GenomeIngestor.Provider, "viruses", genomeIds);
            var ingestor = new GenomeIngestor(dataSource, client);

            var stored = await ingestor.IngestAsync(genomes[0]);
            Assert.That(stored.Status, Is.EqualTo(CrawlState.Done), stored.Error);
            Assert.That(DownloadedFiles(), Is.EquivalentTo(before));

            await using (var features = dataSource.CreateCommand(
                "select length(sequence_sha256), array_length(dinucleotide_counts, 1), " +
                "(select sum(c) from unnest(dinucleotide_counts) c) <= length_bp - 1 from sequence_record"))
            await using (var reader = await features.ExecuteReaderAsync())
            {
                Assert.That(await reader.ReadAsync(), Is.True);
                Assert.That((reader.GetInt32(0), reader.GetInt32(1), reader.GetBoolean(2)), Is.EqualTo((32, 16, true)));
            }

            // A failure halfway through the store transaction: nothing stored, no file left, outcome recorded.
            await Execute("create function test_fail() returns trigger language plpgsql as " +
                "$$ begin raise exception 'injected failure'; end $$; " +
                "create trigger test_fail before insert on profile_level for each row execute function test_fail()");
            try
            {
                var failed = await ingestor.IngestAsync(genomes[1]);
                Assert.That(failed.Status, Is.EqualTo(CrawlState.Failed));
                Assert.That(failed.Error, Does.Contain("injected failure"));
            }
            finally
            {
                await Execute("drop trigger test_fail on profile_level; drop function test_fail()");
            }
            Assert.That(DownloadedFiles(), Is.EquivalentTo(before));

            await using var sources = dataSource.CreateCommand("select count(*) from sequence_source where accession = $1");
            sources.Parameters.AddWithValue(genomeIds[1]);
            Assert.That((long)(await sources.ExecuteScalarAsync())!, Is.Zero);
        }

        // Echarate virus segment L listed twice by BV-BRC: GenBank HM119410 (1000646.3) and its RefSeq copy NC_055339 (1000646.8).
        [Test, Explicit("Downloads one small genome from BV-BRC; run it by name.")]
        public async Task RefSeqCopyOfAStoredGenBankRecordIsSkippedBeforeDownload()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var before = DownloadedFiles();
            using var client = new BvBrcClient(null);
            var genomes = await client.GetGenomesAsync(new[] { "1000646.3", "1000646.8" });
            Assert.That(genomes.Select(g => g.GenomeId), Is.EqualTo(new[] { "1000646.3", "1000646.8" }));
            await CrawlState.RegisterAsync(dataSource, GenomeIngestor.Provider, "viruses", genomes.Select(g => g.GenomeId).ToList());
            var ingestor = new GenomeIngestor(dataSource, client);

            var stored = await ingestor.IngestAsync(genomes[0]);
            Assert.That(stored.Status, Is.EqualTo(CrawlState.Done), stored.Error);
            var copy = await ingestor.IngestAsync(genomes[1]);
            Assert.That((copy.Status, copy.DuplicateBy), Is.EqualTo((CrawlState.Duplicate, "md5")), copy.Error);
            Assert.That(DownloadedFiles(), Is.EquivalentTo(before));

            await using var command = dataSource.CreateCommand(
                "select c.duplicate_of = s.id, c.source_id is null, (select count(*) from sequence_source where accession = $2) " +
                "from crawl_item c, sequence_source s where c.external_id = $2 and s.accession = $1");
            command.Parameters.AddWithValue("1000646.3");
            command.Parameters.AddWithValue("1000646.8");
            await using var reader = await command.ExecuteReaderAsync();
            Assert.That(await reader.ReadAsync(), Is.True);
            Assert.That((reader.GetBoolean(0), reader.GetBoolean(1), reader.GetInt64(2)), Is.EqualTo((true, true, 0L)));
        }

        private static string[] DownloadedFiles() => Directory.GetFiles(Path.GetTempPath(), "bvbrc_*");

        private static async Task Execute(string sql)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(sql);
            await command.ExecuteNonQueryAsync();
        }
    }
}
