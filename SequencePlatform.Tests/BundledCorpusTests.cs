using SequenceAnaliseTests;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class BundledCorpusTests
    {
        [Test]
        public void EveryBundledFileHasTaxonomyAndAnAccession()
        {
            var onDisk = Directory.EnumerateFiles(GenomeCatalog.DirectoryPath, "*.fna.gz", SearchOption.AllDirectories)
                .Select(f => Path.GetRelativePath(GenomeCatalog.DirectoryPath, f).Replace('\\', '/'))
                .ToList();

            Assert.That(BundledCorpus.Genomes.Select(g => g.FileName), Is.EquivalentTo(onDisk));
            Assert.That(BundledCorpus.Genomes.Select(g => BundledCorpus.Accession(g.FileName)), Is.Unique);
        }

        [Test]
        public async Task LoadStoresRecordsKindsLineageAndTheFileProfile()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            // phiX174 (one record), HPV16 (viral path) and Deinococcus (two chromosomes and two plasmids).
            var picked = new[] { "phix174_GCF_000819615.1.fna.gz", "viral/hpv16_NC_001526.4.fna.gz", "deinococcus_radiodurans_r1_GCF_020546685.1.fna.gz" }
                .Select(f => BundledCorpus.Genomes.Single(g => g.FileName == f))
                .ToList();
            foreach (var genome in picked)
                await BundledCorpus.LoadAsync(dataSource, genome);

            await using var command = dataSource.CreateCommand(
                "select s.accession, r.id, r.record_id, r.kind, r.tax_genus, r.tax_family, r.tax_phylum, r.ncbi_taxid is null " +
                "from sequence_record r join sequence_source s on s.id = r.source_id where s.provider = $1 order by r.id");
            command.Parameters.AddWithValue(BundledCorpus.Provider);
            var rows = new List<(string Accession, long Id, string RecordId, string Kind, string? Genus, string? Family, string? Phylum, bool NoTaxid)>();
            await using (var reader = await command.ExecuteReaderAsync())
                while (await reader.ReadAsync())
                    rows.Add((reader.GetString(0), reader.GetInt64(1), reader.GetString(2), reader.GetString(3),
                        reader.IsDBNull(4) ? null : reader.GetString(4), reader.IsDBNull(5) ? null : reader.GetString(5),
                        reader.IsDBNull(6) ? null : reader.GetString(6), reader.GetBoolean(7)));

            Assert.That(rows.Select(r => (r.Accession, r.RecordId, r.Kind)), Is.EqualTo(new[]
            {
                ("GCF_000819615.1", "NC_001422.1", "virus"),
                ("NC_001526.4", "NC_001526.4", "virus"),
                ("GCF_020546685.1", "NZ_CP038663.1", "bacterium"),
                ("GCF_020546685.1", "NZ_CP038664.1", "bacterium"),
                ("GCF_020546685.1", "NZ_CP038666.1", "plasmid"),
                ("GCF_020546685.1", "NZ_CP038665.1", "plasmid"),
            }));
            Assert.That(rows.Select(r => (r.Genus, r.Family, r.Phylum)).Distinct(), Is.EqualTo(new (string?, string?, string?)[]
            {
                ("Sinsheimervirus", "Eubullaviridae", null),
                ("Alphapapillomavirus", "Papillomaviridae", null),
                ("Deinococcus", null, "Deinococcota"),
            }));
            Assert.That(rows.All(r => r.NoTaxid), Is.True);

            // The stored 500 bp level is the profile of the file's record, built here independently.
            var hpv16 = rows.Single(r => r.RecordId == "NC_001526.4");
            var path = GenomeCatalog.PathTo("viral/hpv16_NC_001526.4.fna.gz");
            var expected = ZCurveWindowedProfile.FromStream(FastaReader.OpenSequenceReader(path, "NC_001526.4"));
            var stored = await ProfileLevels.ReadAsync(dataSource, hpv16.Id, RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
            Assert.That(stored!.Windows, Is.EqualTo(expected.ToHalfPrecision().Windows), "stored as float16");
        }

        [Test, Explicit("Profiles and stores all 61 bundled genomes (about 100 Mbp); run it by name.")]
        public async Task LoadsTheWholeBundledCorpus()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            await Execute(dataSource, "delete from profile_level where sequence_id in (select r.id from sequence_record r join sequence_source s on s.id = r.source_id where s.provider = 'NCBI RefSeq')");
            await Execute(dataSource, "delete from sequence_record where source_id in (select id from sequence_source where provider = 'NCBI RefSeq')");
            await Execute(dataSource, "delete from sequence_source where provider = 'NCBI RefSeq'");

            var started = DateTime.UtcNow;
            foreach (var genome in BundledCorpus.Genomes)
                await BundledCorpus.LoadAsync(dataSource, genome);
            TestContext.Out.WriteLine($"{BundledCorpus.Genomes.Count} genomes loaded in {(DateTime.UtcNow - started).TotalSeconds:F0} s");

            var expectedRecords = BundledCorpus.Genomes.Sum(g => FastaReader.ReadHeaders(GenomeCatalog.PathTo(g.FileName)).Count());
            await using var command = dataSource.CreateCommand(
                "select count(distinct s.id), count(r.id) from sequence_source s join sequence_record r on r.source_id = s.id where s.provider = $1");
            command.Parameters.AddWithValue(BundledCorpus.Provider);
            await using var reader = await command.ExecuteReaderAsync();
            await reader.ReadAsync();
            Assert.That((reader.GetInt64(0), reader.GetInt64(1)), Is.EqualTo(((long)BundledCorpus.Genomes.Count, (long)expectedRecords)));
        }

        private static async Task Execute(Npgsql.NpgsqlDataSource dataSource, string sql)
        {
            await using var command = dataSource.CreateCommand(sql);
            await command.ExecuteNonQueryAsync();
        }
    }
}
