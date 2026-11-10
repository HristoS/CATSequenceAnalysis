using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    // S6.4: the role the hosted API and worker connect as (migration 0021). The API tests already run as it (ApiHost), which shows it has
    // every right they use; these show it has nothing beyond that.
    public class AppRoleTests
    {
        private const string insufficientPrivilege = "42501";

        [Test]
        public async Task TheApiDataSourceRunsAsTheAppRoleWithoutSuperuser()
        {
            await using var command = ThrowawayDatabase.AppDataSource.CreateCommand(
                "select current_user, (select rolsuper from pg_roles where rolname = current_user)");
            await using var reader = await command.ExecuteReaderAsync();
            Assert.That(await reader.ReadAsync(), Is.True);
            Assert.That(reader.GetString(0), Is.EqualTo(Database.AppRole));
            Assert.That(reader.GetBoolean(1), Is.False);
        }

        [TestCase("select count(*) from sequence_record")]
        [TestCase("select count(*) from gene_profile")]
        [TestCase("select count(*) from app_user")]
        public async Task ItReadsTheCatalogue(string sql)
        {
            await using var command = ThrowawayDatabase.AppDataSource.CreateCommand(sql);
            Assert.That(await command.ExecuteScalarAsync(), Is.GreaterThanOrEqualTo(0L));
        }

        // A sample of what the API never does: write the computed layers or users, change the schema, run a program on the server.
        [TestCase("delete from gene_profile")]
        [TestCase("update sequence_record set name = name")]
        [TestCase("delete from sequence_record")]
        [TestCase("insert into app_user (name, api_key_sha256) values ('intruder', sha256('x'))")]
        [TestCase("update app_user set disabled_at = null")]
        [TestCase("delete from score_null")]
        [TestCase("truncate search_run cascade")]
        [TestCase("drop table search_hit")]
        [TestCase("create table intruder (id int)")]
        [TestCase("alter table sequence_record add column intruder int")]
        [TestCase("copy (select 1) to program 'echo intruder'")]
        [TestCase("create role intruder")]
        public async Task ItCannotDoWhatTheApiNeverDoes(string sql)
        {
            await using var command = ThrowawayDatabase.AppDataSource.CreateCommand(sql);
            var error = Assert.ThrowsAsync<PostgresException>(async () => await command.ExecuteNonQueryAsync());
            Assert.That(error!.SqlState, Is.EqualTo(insufficientPrivilege), error.MessageText);
        }

        // migrate gives the role a login from SEQUENCE_APP_PASSWORD; the password is quoted on the server, so any character is safe.
        // The role belongs to the whole server, so the test takes the login away again.
        [Test]
        public async Task SetAppPasswordLetsTheRoleLogInWithThatPassword()
        {
            var password = $"it's \"{Guid.NewGuid():N}\" \\ %L";
            try
            {
                await Database.SetAppPasswordAsync(ThrowawayDatabase.DataSource, password);
                var settings = new NpgsqlConnectionStringBuilder(ThrowawayDatabase.DataSource.ConnectionString)
                {
                    Username = Database.AppRole,
                    Password = password,
                    Pooling = false,
                };
                await using var connection = new NpgsqlConnection(settings.ConnectionString);
                await connection.OpenAsync();
                await using var command = new NpgsqlCommand("select current_user", connection);
                Assert.That(await command.ExecuteScalarAsync(), Is.EqualTo(Database.AppRole));
            }
            finally
            {
                await using var reset = ThrowawayDatabase.DataSource.CreateCommand($"alter role {Database.AppRole} nologin password null");
                await reset.ExecuteNonQueryAsync();
            }
        }
    }
}
