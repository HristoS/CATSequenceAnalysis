global using NUnit.Framework;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// C5: every test run gets its own sequence_platform_test_&lt;random&gt; database on the server named by
    /// <see cref="ConnectionString"/>, migrated from scratch and dropped afterwards, so tests never touch
    /// the real catalogue and two runs can't collide.
    /// </summary>
    [SetUpFixture]
    public sealed class ThrowawayDatabase
    {
        public static NpgsqlDataSource DataSource { get; private set; } = null!;

        private static NpgsqlConnectionStringBuilder settings = null!;

        [OneTimeSetUp]
        public async Task CreateAsync()
        {
            settings = Database.ConnectionSettings(ConnectionString.Resolve());
            settings.Database = $"sequence_platform_test_{Guid.NewGuid():N}"[..32];
            await Database.EnsureDatabaseAsync(settings);
            DataSource = NpgsqlDataSource.Create(settings.ConnectionString);
            await Database.MigrateAsync(DataSource);
        }

        [OneTimeTearDown]
        public async Task DropAsync()
        {
            await DataSource.DisposeAsync();
            var maintenance = new NpgsqlConnectionStringBuilder(settings.ConnectionString) { Database = "postgres", Pooling = false };
            await using var connection = new NpgsqlConnection(maintenance.ConnectionString);
            await connection.OpenAsync();
            // The name was generated above from hex digits only, so it's safe to splice in.
            await using var drop = new NpgsqlCommand($"drop database if exists {settings.Database} with (force)", connection);
            await drop.ExecuteNonQueryAsync();
        }

        // Each test uses its own provider so tests sharing the database can filter to their own rows.
        public static async Task<long> InsertRecordAsync(string provider, string accession, long lengthBp, string kind = "virus")
        {
            await using var command = DataSource.CreateCommand(
                "with s as (insert into sequence_source (provider, url, accession, retrieved_at, raw_sha256, raw_byte_size) " +
                "values ($1, 'test://' || $2, $2, now(), '\\x00', 0) returning id) " +
                "insert into sequence_record (source_id, record_id, name, kind, data_class, length_bp) " +
                "select id, $2, 'test ' || $2, $3, 'public_microbe', $4 from s returning id");
            command.Parameters.AddWithValue(provider);
            command.Parameters.AddWithValue(accession);
            command.Parameters.AddWithValue(kind);
            command.Parameters.AddWithValue(lengthBp);
            return (long)(await command.ExecuteScalarAsync())!;
        }

        // A user's private upload (D3): provider 'user upload', no URL or accession, an owner.
        public static async Task<long> InsertUploadAsync(long ownerId, string recordId, long lengthBp, string kind = "virus")
        {
            await using var command = DataSource.CreateCommand(
                "with s as (insert into sequence_source (provider, retrieved_at, raw_sha256, raw_byte_size, owner_id) " +
                "values ('user upload', now(), '\\x00', 0, $1) returning id) " +
                "insert into sequence_record (source_id, record_id, name, kind, data_class, length_bp) " +
                "select id, $2, 'upload ' || $2, $3, 'restricted', $4 from s returning id");
            command.Parameters.AddWithValue(ownerId);
            command.Parameters.AddWithValue(recordId);
            command.Parameters.AddWithValue(kind);
            command.Parameters.AddWithValue(lengthBp);
            return (long)(await command.ExecuteScalarAsync())!;
        }
    }
}
