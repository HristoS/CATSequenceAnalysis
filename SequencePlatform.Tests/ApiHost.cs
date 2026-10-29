using System.Net.Http.Headers;
using Microsoft.AspNetCore.Mvc.Testing;
using Microsoft.AspNetCore.TestHost;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.DependencyInjection.Extensions;
using Npgsql;
using SequencePlatform.Api;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>The API in process over the throwaway database, and clients holding a new user's key.</summary>
    internal static class ApiHost
    {
        // Any public type of the API assembly names its entry point; Program would be ambiguous with the Ingest's.
        public static WebApplicationFactory<SequencePageResponse> Create() =>
            new WebApplicationFactory<SequencePageResponse>().WithWebHostBuilder(host => host.ConfigureTestServices(services =>
            {
                services.RemoveAll<NpgsqlDataSource>();
                services.AddSingleton(ThrowawayDatabase.DataSource);
            }));

        public static async Task<(AppUser User, HttpClient Client)> NewClientAsync(this WebApplicationFactory<SequencePageResponse> factory, string name)
        {
            var (user, apiKey) = await Users.CreateAsync(ThrowawayDatabase.DataSource, $"{name}-{Guid.NewGuid():N}");
            var client = factory.CreateClient();
            client.DefaultRequestHeaders.Authorization = new AuthenticationHeaderValue("Bearer", apiKey);
            return (user, client);
        }
    }
}
