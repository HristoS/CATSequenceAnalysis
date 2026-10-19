using Npgsql;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    public class DuplicateDetectionTests
    {
        private static readonly byte[] hashA = Enumerable.Repeat((byte)0xA1, 32).ToArray();
        private static readonly byte[] hashB = Enumerable.Repeat((byte)0xB2, 32).ToArray();
        private static readonly byte[] hashUnknown = Enumerable.Repeat((byte)0xC3, 32).ToArray();

        [Test]
        public async Task StoredCopyIsFoundOnlyWhenEveryKeyMatches()
        {
            var sourceA = await InsertAsync(GenomeIngestor.Provider, "dup-test-a", "md5-dup-test-a", hashA);
            var sourceB = await InsertAsync(GenomeIngestor.Provider, "dup-test-b", "md5-dup-test-b", hashB);
            await InsertAsync("test-other-provider", "dup-test-c", "md5-dup-test-c", hashUnknown);
            var dataSource = ThrowawayDatabase.DataSource;

            Assert.That(await GenomeIngestor.FindStoredCopyAsync(dataSource, "provider_md5", new[] { "md5-dup-test-a" }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(1, 1, sourceA)));
            Assert.That(await GenomeIngestor.FindStoredCopyAsync(dataSource, "provider_md5", new[] { "md5-dup-test-a", "md5-dup-test-a" }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(1, 1, sourceA)));
            Assert.That(await GenomeIngestor.FindStoredCopyAsync(dataSource, "provider_md5", new[] { "md5-dup-test-b", "md5-dup-test-a" }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(2, 2, Math.Min(sourceA, sourceB))));
            Assert.That(await GenomeIngestor.FindStoredCopyAsync(dataSource, "provider_md5", new[] { "md5-dup-test-a", "md5-not-stored" }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(2, 1, null)), "a partial match isn't a duplicate");
            Assert.That(await GenomeIngestor.FindStoredCopyAsync(dataSource, "provider_md5", new[] { "md5-dup-test-c" }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(1, 0, null)), "another provider's MD5 isn't compared");
            Assert.That(await GenomeIngestor.FindStoredCopyAsync(dataSource, "record_id", Array.Empty<string>()),
                Is.EqualTo(new GenomeIngestor.StoredCopy(0, 0, null)));
            Assert.That(await GenomeIngestor.FindStoredCopyAsync(dataSource, "record_id", new[] { "dup-test-b" }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(1, 1, sourceB)));

            Assert.That(await GenomeIngestor.FindStoredCopyBySha256Async(dataSource, new[] { (byte[])hashB.Clone() }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(1, 1, sourceB)));
            Assert.That(await GenomeIngestor.FindStoredCopyBySha256Async(dataSource, new[] { hashA, hashUnknown }),
                Is.EqualTo(new GenomeIngestor.StoredCopy(2, 1, null)));
        }

        [Test]
        public async Task DuplicateOfIsOnlyAllowedOnDuplicateItems()
        {
            var source = await InsertAsync("test-duplicate-of", "dup-of-1", null, null);
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(
                "insert into crawl_item (provider, external_id, status, source_id, duplicate_of) values ('test-duplicate-of', $1, $2, $3, $4)");
            command.Parameters.AddWithValue("");
            command.Parameters.AddWithValue("");
            command.Parameters.AddWithValue(DBNull.Value);
            command.Parameters.AddWithValue(source);

            async Task Insert(string externalId, string status, long? sourceId)
            {
                command.Parameters[0].Value = externalId;
                command.Parameters[1].Value = status;
                command.Parameters[2].Value = sourceId is { } id ? id : DBNull.Value;
                await command.ExecuteNonQueryAsync();
            }

            await Insert("item-duplicate", CrawlState.Duplicate, null);
            var done = Assert.ThrowsAsync<PostgresException>(() => Insert("item-done", CrawlState.Done, source));
            Assert.That(done!.SqlState, Is.EqualTo(PostgresErrorCodes.CheckViolation));
            var pending = Assert.ThrowsAsync<PostgresException>(() => Insert("item-pending", "pending", null));
            Assert.That(pending!.SqlState, Is.EqualTo(PostgresErrorCodes.CheckViolation));
        }

        // Returns the source id.
        private static async Task<long> InsertAsync(string provider, string accession, string? providerMd5, byte[]? sha256)
        {
            var recordId = await ThrowawayDatabase.InsertRecordAsync(provider, accession, 1000);
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(
                "update sequence_record set provider_md5 = $2, sequence_sha256 = $3 where id = $1 returning source_id");
            command.Parameters.AddWithValue(recordId);
            command.Parameters.AddWithValue((object?)providerMd5 ?? DBNull.Value);
            command.Parameters.AddWithValue((object?)sha256 ?? DBNull.Value);
            return (long)(await command.ExecuteScalarAsync())!;
        }
    }
}
