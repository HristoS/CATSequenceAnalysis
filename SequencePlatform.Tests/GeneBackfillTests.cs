using System.Security.Cryptography;
using System.Text;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    /// <summary>The S7.R gene pass: what it computes per record from a download and BV-BRC's coding sequences, and how it's stored.</summary>
    public class GeneBackfillTests
    {
        private static readonly string gene = "ATG" + string.Concat(Enumerable.Repeat("GCTGAAAAACTG", 200)) + "TAA";

        [Test]
        public void MatchingRecordGetsGenesCountsAndHexamers()
        {
            var flank = Generator.GetRandomDNA(1_000, 21);
            var dna = flank + gene + flank + Reverse(gene) + flank;
            var stored = new[] { Stored(1, "REC1", dna, null) };
            var cds = new[]
            {
                new BvBrcCds("REC1", 1_000, 1_000 + gene.Length, false, "PGF_00421071", "Aspartokinase"),
                new BvBrcCds("REC1", 2_000 + gene.Length, 2_000 + 2 * gene.Length, true, null, null),
                new BvBrcCds("REC1", dna.Length - 10, dna.Length + 20, false, null, null),
            };

            var result = GeneBackfill.Compute(new StringReader($">accn|REC1\n{dna}\n"), stored, cds).Single();

            Assert.That(result.Status, Is.EqualTo(ProteinKmers.Done));
            Assert.That(result.Genes!.Genes.Select(g => (g.Start, g.Reverse, g.Family)),
                Is.EqualTo(new[] { (1_000, false, 421071), (2_000 + gene.Length, true, ZCurveGeneProfile.NoFamily) }));
            Assert.That((result.CodingSequences, result.Skipped, result.CleanTranslations), Is.EqualTo((3, 1, 2)));
            Assert.That(result.Proteins, Is.All.EqualTo("M" + string.Concat(Enumerable.Repeat("AEKL", 200))), "both strands read in coding orientation");
            Assert.That(result.ProteinCounts![4].Sum(), Is.EqualTo(2 * (801 - 3)));
            Assert.That(result.HexamerCounts!.Sum(), Is.EqualTo(dna.Length - 5));
            Assert.That(result.Families, Is.EqualTo(new[] { (421071, "Aspartokinase") }));
        }

        [Test]
        public void MismatchMissingAndNoCdsAreRecorded()
        {
            var kept = Generator.GetRandomDNA(3_000, 22);
            var changed = Generator.GetRandomDNA(3_000, 23);
            var stored = new[] { Stored(1, "NOCDS", kept, null), Stored(2, "CHANGED", changed, null), Stored(3, "GONE", "ACGT", null) };
            var download = $">accn|NOCDS\n{kept}\n>accn|CHANGED\n{changed[..2999]}{(changed[^1] == 'A' ? 'C' : 'A')}\n";

            var results = GeneBackfill.Compute(new StringReader(download), stored, Array.Empty<BvBrcCds>()).ToDictionary(r => r.SequenceId);

            Assert.That(results[1].Status, Is.EqualTo(ProteinKmers.NoCds));
            Assert.That(results[1].HexamerCounts, Is.Not.Null, "a record without coding sequences still gets its 6-mers");
            Assert.That((results[2].Status, results[2].Genes), Is.EqualTo((ProteinKmers.Mismatch, (ZCurveGeneProfile?)null)));
            Assert.That(results[3].Status, Is.EqualTo(ProteinKmers.Missing));
        }

        [Test]
        public void MycoplasmaOrdersReadTgaAsTryptophan()
        {
            var tgaGene = "ATG" + string.Concat(Enumerable.Repeat("TGAGCT", 300)) + "TAA";
            var stored = new[] { Stored(1, "MYCO", tgaGene, "Mycoplasmatales"), Stored(2, "ECOLI", tgaGene, "Enterobacterales") };
            var download = $">accn|MYCO\n{tgaGene}\n>accn|ECOLI\n{tgaGene}\n";
            var cds = new[] { new BvBrcCds("MYCO", 0, tgaGene.Length, false, null, null), new BvBrcCds("ECOLI", 0, tgaGene.Length, false, null, null) };

            var results = GeneBackfill.Compute(new StringReader(download), stored, cds).ToDictionary(r => r.SequenceId);

            Assert.That(results[1].CleanTranslations, Is.EqualTo(1));
            Assert.That(results[1].Proteins.Single(), Does.StartWith("MWAWA"));
            Assert.That(results[2].CleanTranslations, Is.EqualTo(0), "table 11 reads TGA as a stop");
        }

        [Test]
        public void ProteinCountsSurviveCompressionAndAreCapped()
        {
            var counts = new int[8000];
            counts[0] = 70_000;
            counts[7999] = 12;
            var decoded = ProteinKmers.Decode(ProteinKmers.Encode(counts), 3);
            Assert.That(decoded[0], Is.EqualTo(ushort.MaxValue));
            Assert.That(decoded[7999], Is.EqualTo(12));
            Assert.That(decoded.Sum(), Is.EqualTo(ushort.MaxValue + 12));
            Assert.That(ProteinFamilies.Parse("PGF_00421071"), Is.EqualTo(421071));
            Assert.That(ProteinFamilies.Parse("PLF_1"), Is.EqualTo(ZCurveGeneProfile.NoFamily));
        }

        [Test]
        public void ThreeMersComeFromTheFourMersShortOfEachGenesLastOne()
        {
            var proteins = new[] { "MAEKLWWAC", "MKKLL", "MAEKL" };
            var threes = AminoAcidKmers.Count(proteins, 3);

            var derived = ProteinKmers.Shorter(AminoAcidKmers.Count(proteins, 4));

            foreach (var protein in proteins)
                threes[AminoAcidKmers.Count(new[] { protein[^3..] }, 3).ToList().FindIndex(c => c == 1)]--;
            Assert.That(derived, Is.EqualTo(threes));
            Assert.That(ProteinKmers.Shorter(derived), Has.Length.EqualTo(400), "and the 2-mers from those");
        }

        [Test]
        public async Task GeneRowAndCountsRoundTripThroughTheDatabase()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var path = Path.Combine(Path.GetTempPath(), $"genes_{Guid.NewGuid():N}.fna");
            var dna = Generator.GetRandomDNA(40_000, 24) + gene;
            File.WriteAllText(path, ">accn|GB000001 test\n" + dna + "\n");
            List<ProfiledRecord> profiled;
            try { profiled = RecordStore.Profile(path).ToList(); }
            finally { File.Delete(path); }

            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            var stored = await RecordStore.StoreAsync(connection, transaction, new SourceRow("gene-roundtrip", "test://genes", "1.1", SHA256.HashData(Encoding.ASCII.GetBytes("genes")), 0, null),
                profiled.Select(p => (p, new RecordRow("GB000001", "test", "bacterium", null, new Lineage(null, null, null, null, null, null, null, null, null), "{}"))).ToList());
            long id = stored.SequenceIds[0];
            var genes = ZCurveGeneProfile.FromRecord(dna, new[] { (40_000, dna.Length, false, 7) });
            Assert.That(await GeneProfiles.WriteAsync(connection, transaction, id, genes), Is.True);
            Assert.That(await GeneProfiles.WriteAsync(connection, transaction, id, genes), Is.True, "a rerun replaces the row");
            var counts = AminoAcidKmers.Count(new[] { ZCurveGeneProfile.Translate(gene) }, 4);
            await ProteinKmers.WriteAsync(connection, transaction, id, 4, ProteinKmers.Done, 1, counts);
            await ProteinFamilies.AddAsync(connection, transaction, new[] { (7, "first"), (7, "duplicate") });
            await ProteinFamilies.AddAsync(connection, transaction, new[] { (7, "later") });
            await transaction.CommitAsync();

            var read = await GeneProfiles.ReadAsync(dataSource, id);
            Assert.That(read!.Genes, Is.EqualTo(genes.Genes));
            Assert.That(await ProteinKmers.ReadAsync(dataSource, id, 4), Is.EqualTo(counts));
            await using var product = dataSource.CreateCommand("select product from protein_family where id = 7");
            Assert.That(await product.ExecuteScalarAsync(), Is.EqualTo("first"), "a family's first product is kept");
            await using var levels = dataSource.CreateCommand("select count(*) from profile_level where sequence_id = $1 and method = $2");
            levels.Parameters.AddWithValue(id);
            levels.Parameters.AddWithValue(RecordStore.Method);
            Assert.That(await levels.ExecuteScalarAsync(), Is.EqualTo(4L), "the genome ladder is untouched");
        }

        [Test]
        public async Task GeneProfileIsSplitByLength()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var genes = new[] { 60, 300, 700 }.Select(aa => "ATG" + Generator.GetRandomDNA(3 * (aa - 1), 30 + aa).Replace("TAA", "TAC")
                .Replace("TAG", "TAC").Replace("TGA", "TGC") + "TAA").ToArray();
            var dna = string.Concat(genes);
            int a = genes[0].Length, b = a + genes[1].Length;
            var profile = ZCurveGeneProfile.FromRecord(dna, new[] { (0, a, false, 1), (a, b, false, 2), (b, dna.Length, false, 3) });
            var id = await ThrowawayDatabase.InsertRecordAsync("gene-buckets", "GB1", dna.Length, "bacterium");
            await using var connection = await dataSource.OpenConnectionAsync();
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                Assert.That(await GeneProfiles.WriteAsync(connection, transaction, id, profile), Is.True);
                await transaction.CommitAsync();
            }

            var buckets = new List<(short Bucket, int Genes)>();
            await using (var read = new Npgsql.NpgsqlCommand("select bucket, gene_count from gene_profile where sequence_id = $1 order by bucket", connection))
            {
                read.Parameters.AddWithValue(id);
                await using var reader = await read.ExecuteReaderAsync();
                while (await reader.ReadAsync()) buckets.Add((reader.GetInt16(0), reader.GetInt32(1)));
            }
            Assert.That(profile.Genes.Select(g => g.Residues), Is.EqualTo(new[] { 60, 300, 700 }), "the stop codon isn't a residue");
            Assert.That(buckets, Is.EqualTo(new[] { ((short)0, 1), ((short)GeneProfiles.BucketOf(60), 1), ((short)GeneProfiles.BucketOf(300), 1) }));
            Assert.That((GeneProfiles.BucketOf(50), GeneProfiles.BucketOf(499), GeneProfiles.BucketOf(500)), Is.EqualTo((1, GeneProfiles.MaxBucket, 0)));
            Assert.That(GeneProfiles.BucketsFor(300, 0.1), Is.EqualTo(Enumerable.Range(GeneProfiles.BucketOf(270), GeneProfiles.BucketOf(330) - GeneProfiles.BucketOf(270) + 1)));
            Assert.That(GeneProfiles.BucketsFor(800, 0.1), Is.EqualTo(new[] { 0 }));
            Assert.That((await GeneProfiles.ReadAsync(dataSource, id))!.Genes, Is.EqualTo(profile.Genes));
        }

        private static StoredGeneRecord Stored(long id, string recordId, string dna, string? order)
        {
            using var reader = new SequenceFeatureReader(new StringReader(dna));
            while (reader.Read(new char[4096], 0, 4096) > 0) { }
            var (_, dinucleotides) = reader.Finish();
            return new StoredGeneRecord(id, recordId, SHA256.HashData(Encoding.ASCII.GetBytes(dna)), dinucleotides, false, order);
        }

        private static string Reverse(string s) =>
            new string(s.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' }).ToArray());
    }
}
