using System.Security.Cryptography;
using System.Text;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    /// <summary>6-mer counts (migration 0013): stored with every new record, and backfilled only where a download is the stored sequence.</summary>
    public class RecordKmerTests
    {
        [Test]
        public async Task StoringARecordStoresItsHexamers()
        {
            var path = Path.Combine(Path.GetTempPath(), $"kmers_{Guid.NewGuid():N}.fna");
            File.WriteAllText(path, ">accn|KM000001 test\n" + Generator.GetRandomDNA(3_000, 11) + "\n");
            List<ProfiledRecord> profiled;
            try { profiled = RecordStore.Profile(path).ToList(); }
            finally { File.Delete(path); }

            var dataSource = ThrowawayDatabase.DataSource;
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            var stored = await RecordStore.StoreAsync(connection, transaction, new SourceRow("kmer-test", "test://kmer", "1.1", new byte[32], 0, null),
                profiled.Select(p => (p, new RecordRow("KM000001", "test", "virus", null, new Lineage(null, null, null, null, null, null, null, null, null), "{}"))).ToList());
            await transaction.CommitAsync();

            var counts = await RecordKmers.ReadAsync(dataSource, stored.SequenceIds[0]);
            Assert.That(counts, Is.EqualTo(profiled[0].HexamerCounts));
            Assert.That(counts!.Sum(), Is.EqualTo(3_000 - 5));
        }

        [Test]
        public void BackfillCountsOnlyTheSequenceThatWasProfiled()
        {
            var kept = Generator.GetRandomDNA(2_000, 1);
            var changed = Generator.GetRandomDNA(2_000, 2);
            var stored = new[]
            {
                Stored(1, "KEPT", kept),
                Stored(2, "CHANGED", changed),
                Stored(3, "GONE", Generator.GetRandomDNA(500, 3)),
            };
            var download = $">accn|KEPT\n{kept}\n>accn|CHANGED\n{changed[..1999]}A\n>accn|NEW\nACGTACGT\n";

            var results = KmerBackfill.Check(new StringReader(download), stored).ToDictionary(r => r.SequenceId);

            Assert.That(results.Keys, Is.EquivalentTo(new long[] { 1, 2, 3 }), "a record the catalogue doesn't hold is ignored");
            Assert.That(results[1].Status, Is.EqualTo(RecordKmers.Done));
            Assert.That(results[1].Counts!.Sum(), Is.EqualTo(2_000 - 5));
            Assert.That((results[2].Status, results[2].Counts), Is.EqualTo((RecordKmers.Mismatch, (int[]?)null)));
            Assert.That((results[3].Status, results[3].Counts), Is.EqualTo((RecordKmers.Missing, (int[]?)null)));
        }

        private static StoredRecordFeatures Stored(long id, string recordId, string dna)
        {
            using var reader = new SequenceFeatureReader(new StringReader(dna));
            while (reader.Read(new char[4096], 0, 4096) > 0) { }
            var (_, dinucleotides) = reader.Finish();
            return new StoredRecordFeatures(id, recordId, SHA256.HashData(Encoding.ASCII.GetBytes(dna)), dinucleotides);
        }
    }
}
