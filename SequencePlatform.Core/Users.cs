using System.Security.Cryptography;
using Npgsql;

namespace SequencePlatform.Core
{
    public sealed record AppUser(long Id, string Name);

    /// <summary>
    /// app_user (migration 0010, D3): who may call the API. A key is shown once, when the user is created or the key rotated,
    /// and only its SHA-256 is kept, so a lost key is replaced with <see cref="RotateKeyAsync"/>.
    /// </summary>
    public static class Users
    {
        public const string KeyPrefix = "sp_";

        private const int keyBytes = 32;

        /// <returns>The new user and their API key, which isn't stored anywhere.</returns>
        public static async Task<(AppUser User, string ApiKey)> CreateAsync(NpgsqlDataSource dataSource, string name)
        {
            if (string.IsNullOrWhiteSpace(name)) throw new ArgumentException("A user needs a name.", nameof(name));

            var apiKey = NewKey();
            await using var command = dataSource.CreateCommand("insert into app_user (name, api_key_sha256) values ($1, $2) returning id");
            command.Parameters.AddWithValue(name.Trim());
            command.Parameters.AddWithValue(Hash(apiKey));
            var id = (long)(await command.ExecuteScalarAsync())!;
            return (new AppUser(id, name.Trim()), apiKey);
        }

        /// <summary>A new key for an enabled user, keeping the user (and so their uploads and searches); the old key stops working at once.</summary>
        /// <returns>The new key, or null when no enabled user has this name.</returns>
        public static async Task<string?> RotateKeyAsync(NpgsqlDataSource dataSource, string name)
        {
            var apiKey = NewKey();
            await using var command = dataSource.CreateCommand("update app_user set api_key_sha256 = $2 where name = $1 and disabled_at is null");
            command.Parameters.AddWithValue(name);
            command.Parameters.AddWithValue(Hash(apiKey));
            return await command.ExecuteNonQueryAsync() == 1 ? apiKey : null;
        }

        private static string NewKey() => KeyPrefix + Base64Url(RandomNumberGenerator.GetBytes(keyBytes));

        private static byte[] Hash(string apiKey) => SHA256.HashData(System.Text.Encoding.ASCII.GetBytes(apiKey));

        /// <returns>null for an unknown key or a disabled user.</returns>
        public static async Task<AppUser?> AuthenticateAsync(NpgsqlDataSource dataSource, string apiKey, CancellationToken cancellationToken = default)
        {
            if (!apiKey.StartsWith(KeyPrefix, StringComparison.Ordinal) || apiKey.Length > 100) return null;

            // Lookup by hash: equality on a SHA-256 leaks nothing useful through timing.
            await using var command = dataSource.CreateCommand(
                "select id, name from app_user where api_key_sha256 = $1 and disabled_at is null");
            command.Parameters.AddWithValue(Hash(apiKey));
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            return await reader.ReadAsync(cancellationToken) ? new AppUser(reader.GetInt64(0), reader.GetString(1)) : null;
        }

        /// <returns>false when no enabled user has this name.</returns>
        public static async Task<bool> DisableAsync(NpgsqlDataSource dataSource, string name)
        {
            await using var command = dataSource.CreateCommand("update app_user set disabled_at = now() where name = $1 and disabled_at is null");
            command.Parameters.AddWithValue(name);
            return await command.ExecuteNonQueryAsync() == 1;
        }

        private static string Base64Url(byte[] bytes) =>
            Convert.ToBase64String(bytes).TrimEnd('=').Replace('+', '-').Replace('/', '_');
    }
}
