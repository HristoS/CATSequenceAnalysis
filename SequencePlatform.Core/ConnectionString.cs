using System.Text.Json;

namespace SequencePlatform.Core
{
    /// <summary>
    /// Where the Ingest and the platform tests find their connection string: the SEQUENCE_PLATFORM_DB
    /// environment variable, or else a git-ignored appsettings.json next to the executable - never a
    /// committed file. The API (P4) uses ASP.NET configuration instead and doesn't need this.
    /// </summary>
    public static class ConnectionString
    {
        public const string ConnectionVariable = "SEQUENCE_PLATFORM_DB";
        private const string settingsFileName = "appsettings.json";
        private const string connectionStringName = "SequencePlatform";

        public static string Resolve()
        {
            var raw = Environment.GetEnvironmentVariable(ConnectionVariable);
            if (string.IsNullOrWhiteSpace(raw)) raw = FromSettingsFile();
            if (string.IsNullOrWhiteSpace(raw))
                throw new InvalidOperationException(
                    $"Set {ConnectionVariable}, or ConnectionStrings:{connectionStringName} in {settingsFileName}, to a PostgreSQL " +
                    "connection string, e.g. \"Host=localhost;Port=5432;Username=...;Password=...;Database=sequence_platform\".");
            return raw;
        }

        // Plain System.Text.Json instead of Microsoft.Extensions.Configuration: one value doesn't justify the packages.
        private static string? FromSettingsFile()
        {
            var path = Path.Combine(AppContext.BaseDirectory, settingsFileName);
            if (!File.Exists(path)) return null;

            using var document = JsonDocument.Parse(File.ReadAllText(path));
            return document.RootElement.TryGetProperty("ConnectionStrings", out var connectionStrings)
                && connectionStrings.TryGetProperty(connectionStringName, out var value)
                ? value.GetString()
                : null;
        }
    }
}
