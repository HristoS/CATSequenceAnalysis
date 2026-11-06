using System.Net;
using System.Net.Http.Headers;
using System.Net.Http.Json;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>The API in process over the throwaway database (S4.2a): authentication, visibility, and the catalogue endpoints.</summary>
    public class ApiTests
    {
        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        private Task<(AppUser User, HttpClient Client)> NewClientAsync(string name) => factory.NewClientAsync(name);

        private static async Task ExecuteAsync(string sql, params object[] values)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(sql);
            foreach (var value in values) command.Parameters.AddWithValue(value);
            await command.ExecuteNonQueryAsync();
        }

        [Test]
        public async Task RotatingAKeyKeepsTheUserAndRevokesTheOldKey()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var name = $"rotated-{Guid.NewGuid():N}";
            var (user, oldKey) = await Users.CreateAsync(dataSource, name);

            var newKey = await Users.RotateKeyAsync(dataSource, name);

            Assert.That(newKey, Is.Not.Null.And.Not.EqualTo(oldKey));
            Assert.That(await Users.AuthenticateAsync(dataSource, oldKey), Is.Null);
            Assert.That(await Users.AuthenticateAsync(dataSource, newKey!), Is.EqualTo(user));
            await Users.DisableAsync(dataSource, name);
            Assert.That(await Users.RotateKeyAsync(dataSource, name), Is.Null, "a disabled user stays disabled");
        }

        [Test]
        public async Task EveryEndpointNeedsAnEnabledKey()
        {
            var anonymous = factory.CreateClient();
            var wrong = factory.CreateClient();
            wrong.DefaultRequestHeaders.Authorization = new AuthenticationHeaderValue("Bearer", Users.KeyPrefix + "not-a-key");
            var name = $"revoked-{Guid.NewGuid():N}";
            var (_, apiKey) = await Users.CreateAsync(ThrowawayDatabase.DataSource, name);
            await Users.DisableAsync(ThrowawayDatabase.DataSource, name);
            var revoked = factory.CreateClient();
            revoked.DefaultRequestHeaders.Authorization = new AuthenticationHeaderValue("Bearer", apiKey);

            foreach (var client in new[] { anonymous, wrong, revoked })
                foreach (var route in new[] { "/api/sequences", "/api/sequences/1", "/api/sequences/1/profile" })
                {
                    var response = await client.GetAsync(route);
                    Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.Unauthorized), route);
                    Assert.That(response.Headers.WwwAuthenticate.ToString(), Is.EqualTo("Bearer"));
                }
        }

        [Test]
        public async Task ListPagesThePublicCatalogueAndYourOwnUploadsOnly()
        {
            var (alice, client) = await NewClientAsync("alice");
            var (bob, _) = await NewClientAsync("bob");
            var shared = new List<long>();
            for (int i = 0; i < 3; i++) shared.Add(await ThrowawayDatabase.InsertRecordAsync("api-list", $"apilist-{i}", 12_345));
            var mine = await ThrowawayDatabase.InsertUploadAsync(alice.Id, "apilist-alice", 12_345);
            await ThrowawayDatabase.InsertUploadAsync(bob.Id, "apilist-bob", 12_345);

            var ids = new List<long>();
            long? after = 0;
            while (after is not null)
            {
                var page = (await client.GetFromJsonAsync<JsonElement>($"/api/sequences?name=apilist&limit=2&after={after}"))!;
                ids.AddRange(page.GetProperty("items").EnumerateArray().Select(r => r.GetProperty("id").GetInt64()));
                after = page.TryGetProperty("nextAfter", out var next) ? next.GetInt64() : null;
            }

            Assert.That(ids, Is.EqualTo(shared.Append(mine).ToList()), "public records and Alice's upload, in id order, over two pages");
        }

        [Test]
        public async Task DetailCarriesLineagesLevelsAndProvenanceButNoOneElsesUpload()
        {
            var (alice, client) = await NewClientAsync("alice");
            var (bob, _) = await NewClientAsync("bob");
            var id = await ThrowawayDatabase.InsertRecordAsync("api-detail", "detail-1", 40_000, kind: "bacterium");
            await ExecuteAsync("update sequence_record set ncbi_taxid = 99999001, tax_genus = 'OldGenus', " +
                "annotations = '{\"bvbrc\": {\"assembly_accession\": \"GCF_999999001.1\"}}' where id = $1", id);
            await ExecuteAsync("insert into taxon (taxid, status, current_taxid, name, tax_superkingdom, tax_genus) " +
                "values (99999001, 'current', 99999001, 'Testus primus', 'Bacteria', 'Testus') on conflict do nothing");
            await ExecuteAsync("insert into gtdb_taxon (assembly_key, gtdb_accession, release, tax_domain, tax_genus) " +
                "values ('999999001.1', 'RS_GCF_999999001.1', 'v232', 'Bacteria', 'Testus_A') on conflict do nothing");
            var fine = new ZCurveWindowedProfile(Generator.GetRandomDNA(40_000, seed: 21));
            await using (var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync())
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                await ProfileLevels.InsertAsync(connection, transaction, id, RecordStore.Method, RecordStore.MethodVersion, 1, fine);
                await ProfileLevels.InsertAsync(connection, transaction, id, RecordStore.Method, RecordStore.MethodVersion, 4, fine.Downsample(4, 0));
                await transaction.CommitAsync();
            }
            var theirs = await ThrowawayDatabase.InsertUploadAsync(bob.Id, "detail-bob", 40_000);

            var detail = (await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}"))!;

            Assert.That(detail.GetProperty("sequence").GetProperty("providerLineage").GetProperty("genus").GetString(), Is.EqualTo("OldGenus"));
            Assert.That(detail.GetProperty("ncbi").GetProperty("lineage").GetProperty("genus").GetString(), Is.EqualTo("Testus"));
            Assert.That(detail.GetProperty("gtdb").GetProperty("lineage").GetProperty("genus").GetString(), Is.EqualTo("Testus_A"));
            Assert.That(detail.GetProperty("annotations").GetProperty("bvbrc").GetProperty("assembly_accession").GetString(), Is.EqualTo("GCF_999999001.1"));
            Assert.That(detail.GetProperty("levels").EnumerateArray().Select(l => (l.GetProperty("windowSize").GetInt32(), l.GetProperty("windowCount").GetInt32())),
                Is.EqualTo(new[] { (500, 80), (2_000, 20) }));
            Assert.That(detail.GetProperty("provenance").GetProperty("url").GetString(), Is.EqualTo("test://detail-1"));
            Assert.That(detail.GetProperty("owned").GetBoolean(), Is.False);
            Assert.That((await client.GetAsync($"/api/sequences/{theirs}")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await client.GetAsync("/api/sequences/999999999")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
        }

        [Test]
        public async Task ProfileIsTheStoredLevelAsZcwpBytes()
        {
            var (alice, client) = await NewClientAsync("alice");
            var (bob, _) = await NewClientAsync("bob");
            var mine = await ThrowawayDatabase.InsertUploadAsync(alice.Id, "profile-alice", 30_000);
            var theirs = await ThrowawayDatabase.InsertUploadAsync(bob.Id, "profile-bob", 30_000);
            var fine = new ZCurveWindowedProfile(Generator.GetRandomDNA(30_000, seed: 22));
            await using (var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync())
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                await ProfileLevels.InsertAsync(connection, transaction, mine, RecordStore.Method, RecordStore.MethodVersion, 1, fine);
                await ProfileLevels.InsertAsync(connection, transaction, theirs, RecordStore.Method, RecordStore.MethodVersion, 1, fine);
                await transaction.CommitAsync();
            }

            var response = await client.GetAsync($"/api/sequences/{mine}/profile");
            var missingLevel = await client.GetAsync($"/api/sequences/{mine}/profile?windowSize=32000");

            Assert.That(response.Content.Headers.ContentType!.MediaType, Is.EqualTo("application/x-zcwp"));
            // Stored as float16 (ZCWP version 2), served as version 1.
            Assert.That(await response.Content.ReadAsByteArrayAsync(), Is.EqualTo(fine.ToHalfPrecision().ToBinary()));
            Assert.That(missingLevel.StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That(missingLevel.Content.Headers.ContentType!.MediaType, Is.EqualTo("application/problem+json"));
            Assert.That((await client.GetAsync($"/api/sequences/{theirs}/profile")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
        }

        [Test]
        public async Task InvalidFiltersAreProblemDetails()
        {
            var (_, client) = await NewClientAsync("alice");

            var rankAlone = await client.GetAsync("/api/sequences?rank=family");
            var unknownRank = await client.GetAsync("/api/sequences?rank=tribe&taxon=X");
            var pageTooBig = await client.GetAsync("/api/sequences?limit=1001");

            foreach (var response in new[] { rankAlone, unknownRank, pageTooBig })
            {
                Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
                Assert.That(response.Content.Headers.ContentType!.MediaType, Is.EqualTo("application/problem+json"));
            }
            var problem = (await rankAlone.Content.ReadFromJsonAsync<JsonElement>())!;
            Assert.That(problem.GetProperty("detail").GetString(), Does.Contain("Rank and taxon filter go together"));
        }
    }
}
