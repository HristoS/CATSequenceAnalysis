using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// Migrations 0014 and 0018 (S7.R, S7.S): every profile_level row has a window again once 0018 has dropped 0014's gene-row exception,
    /// and protein_kmer holds counts exactly when its status is done.
    /// </summary>
    public class GeneLayerMigrationTests
    {
        [Test]
        public async Task NoLevelRowWithoutAWindow()
        {
            var sequenceId = await StoreRecordAsync("GL000001");
            var gene = "ATG" + string.Concat(Enumerable.Repeat("GCT", 100)) + "TAA";
            var payload = ZCurveGeneProfile.FromRecord(gene, new[] { (0, gene.Length, false, 421071) }).ToBinary();

            foreach (var method in new[] { "zcurve-cds", "zcurve-windowed" })
            {
                var ex = Assert.ThrowsAsync<PostgresException>(() => InsertLevelAsync(sequenceId, method, 0, payload));
                Assert.That(ex!.ConstraintName, Is.EqualTo("profile_level_window_count_check"), method);
            }
        }

        [Test]
        public async Task ProteinKmerCountsExistExactlyWhenDone()
        {
            var sequenceId = await StoreRecordAsync("GL000002");
            await InsertProteinKmerAsync(sequenceId, 4, "done", new byte[] { 1, 2, 3 });
            var threeMers = Assert.ThrowsAsync<PostgresException>(() => InsertProteinKmerAsync(sequenceId, 3, "no_cds", null));
            Assert.That(threeMers!.SqlState, Is.EqualTo(PostgresErrorCodes.CheckViolation), "0016 keeps 4-mers only");
            var (other, another) = (await StoreRecordAsync("GL000003"), await StoreRecordAsync("GL000004"));
            var noCounts = Assert.ThrowsAsync<PostgresException>(() => InsertProteinKmerAsync(other, 4, "done", null));
            Assert.That(noCounts!.SqlState, Is.EqualTo(PostgresErrorCodes.CheckViolation));
            var countsWithoutDone = Assert.ThrowsAsync<PostgresException>(() => InsertProteinKmerAsync(another, 4, "mismatch", new byte[] { 1 }));
            Assert.That(countsWithoutDone!.SqlState, Is.EqualTo(PostgresErrorCodes.CheckViolation));
        }

        private static async Task<long> StoreRecordAsync(string recordId)
        {
            var path = Path.Combine(Path.GetTempPath(), $"genelayer_{Guid.NewGuid():N}.fna");
            File.WriteAllText(path, $">accn|{recordId} test\n" + Generator.GetRandomDNA(3_000, int.Parse(recordId[2..])) + "\n");
            List<ProfiledRecord> profiled;
            try { profiled = RecordStore.Profile(path).ToList(); }
            finally { File.Delete(path); }

            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            var stored = await RecordStore.StoreAsync(connection, transaction, new SourceRow($"gene-layer-{recordId}", $"test://{recordId}", "1.1", System.Security.Cryptography.SHA256.HashData(System.Text.Encoding.ASCII.GetBytes(recordId)), 0, null),
                profiled.Select(p => (p, new RecordRow(recordId, "test", "bacterium", null, new Lineage(null, null, null, null, null, null, null, null, null), "{}"))).ToList());
            await transaction.CommitAsync();
            return stored.SequenceIds[0];
        }

        private static async Task InsertLevelAsync(long sequenceId, string method, int windowCount, byte[] payload)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(
                "insert into profile_level (sequence_id, method, method_version, fine_window_size, factor, window_count, format_version, payload, byte_size, sha256) " +
                "values ($1, $2, '1', $3, 1, $4, 1, $5, $6, $7)");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(ZCurveGeneProfile.WindowBases);
            command.Parameters.AddWithValue(windowCount);
            command.Parameters.AddWithValue(payload);
            command.Parameters.AddWithValue((long)payload.Length);
            command.Parameters.AddWithValue(System.Security.Cryptography.SHA256.HashData(payload));
            await command.ExecuteNonQueryAsync();
        }

        private static async Task InsertProteinKmerAsync(long sequenceId, int k, string status, byte[]? payload)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(
                "insert into protein_kmer (sequence_id, k, format_version, gene_count, payload, status) values ($1, $2, 1, 1, $3, $4)");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue((short)k);
            command.Parameters.Add(new NpgsqlParameter { Value = (object?)payload ?? DBNull.Value, NpgsqlDbType = NpgsqlTypes.NpgsqlDbType.Bytea });
            command.Parameters.AddWithValue(status);
            await command.ExecuteNonQueryAsync();
        }
    }
}
