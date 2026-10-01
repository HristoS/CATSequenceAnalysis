using System.Reflection;
using Npgsql;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// Connection handling and schema migrations. The connection string comes only from the
    /// SEQUENCE_PLATFORM_DB environment variable - never from source or config files in the repo.
    /// </summary>
    internal static class Database
    {
        public const string ConnectionVariable = "SEQUENCE_PLATFORM_DB";
        private const string defaultDatabaseName = "sequence_platform";

        public static NpgsqlConnectionStringBuilder ConnectionSettings()
        {
            var raw = Environment.GetEnvironmentVariable(ConnectionVariable);
            if (string.IsNullOrWhiteSpace(raw))
                throw new InvalidOperationException(
                    $"Set {ConnectionVariable} to a PostgreSQL connection string, e.g. " +
                    "\"Host=localhost;Port=5432;Username=...;Password=...;Database=sequence_platform\".");

            var settings = new NpgsqlConnectionStringBuilder(raw);
            if (string.IsNullOrEmpty(settings.Database)) settings.Database = defaultDatabaseName;
            return settings;
        }

        // Creates the platform database if it doesn't exist yet. Connects to the server's
        // maintenance database for that one statement; nothing else on the server is touched.
        public static async Task EnsureDatabaseAsync(NpgsqlConnectionStringBuilder settings)
        {
            var maintenance = new NpgsqlConnectionStringBuilder(settings.ConnectionString) { Database = "postgres" };
            await using var connection = new NpgsqlConnection(maintenance.ConnectionString);
            await connection.OpenAsync();

            await using var exists = new NpgsqlCommand("select 1 from pg_database where datname = $1", connection);
            exists.Parameters.AddWithValue(settings.Database!);
            if (await exists.ExecuteScalarAsync() is not null) return;

            // Identifiers can't be parameters; the name is validated to a plain identifier first.
            if (!System.Text.RegularExpressions.Regex.IsMatch(settings.Database!, "^[a-z_][a-z0-9_]*$"))
                throw new InvalidOperationException($"Refusing to create database with unusual name '{settings.Database}'.");
            await using var create = new NpgsqlCommand($"create database {settings.Database} encoding 'UTF8' template template0", connection);
            await create.ExecuteNonQueryAsync();
            Console.WriteLine($"Created database {settings.Database}.");
        }

        // Applies embedded Migrations/*.sql in name order, each in its own transaction, recording
        // applied versions so re-running is a no-op. Same files for local PG, Supabase and TencentDB.
        public static async Task MigrateAsync(NpgsqlDataSource dataSource)
        {
            await using var connection = await dataSource.OpenConnectionAsync();
            await using (var init = new NpgsqlCommand(
                "create table if not exists schema_migration (version text primary key, applied_at timestamptz not null default now())",
                connection))
                await init.ExecuteNonQueryAsync();

            var assembly = Assembly.GetExecutingAssembly();
            var scripts = assembly.GetManifestResourceNames()
                .Where(n => n.EndsWith(".sql", StringComparison.Ordinal) && n.Contains(".Migrations."))
                .OrderBy(n => n, StringComparer.Ordinal);

            foreach (var resource in scripts)
            {
                var version = resource.Substring(resource.IndexOf(".Migrations.", StringComparison.Ordinal) + ".Migrations.".Length);

                await using (var check = new NpgsqlCommand("select 1 from schema_migration where version = $1", connection))
                {
                    check.Parameters.AddWithValue(version);
                    if (await check.ExecuteScalarAsync() is not null) continue;
                }

                using var reader = new StreamReader(assembly.GetManifestResourceStream(resource)!);
                var sql = await reader.ReadToEndAsync();

                await using var transaction = await connection.BeginTransactionAsync();
                await using (var apply = new NpgsqlCommand(sql, connection, transaction))
                    await apply.ExecuteNonQueryAsync();
                await using (var record = new NpgsqlCommand("insert into schema_migration (version) values ($1)", connection, transaction))
                {
                    record.Parameters.AddWithValue(version);
                    await record.ExecuteNonQueryAsync();
                }
                await transaction.CommitAsync();
                Console.WriteLine($"Applied migration {version}.");
            }
        }
    }
}
