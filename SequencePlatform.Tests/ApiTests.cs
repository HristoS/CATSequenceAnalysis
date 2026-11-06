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
        public async Task FacetsCountWhatTheCallerCanSee()
        {
            var (alice, client) = await NewClientAsync("alice");
            var (bob, _) = await NewClientAsync("bob");
            foreach (var (accession, kind) in new[] { ("facet-b", "bacterium"), ("facet-a", "archaeon"), ("facet-e", "eukaryote") })
                await ThrowawayDatabase.InsertRecordAsync("api-facets", accession, 1_000_000, kind);
            await ThrowawayDatabase.InsertUploadAsync(alice.Id, "facet-alice", 12_345);
            await ThrowawayDatabase.InsertUploadAsync(bob.Id, "facet-bob", 12_345);

            var facets = (await client.GetFromJsonAsync<JsonElement>("/api/sequences/facets"))!;
            Dictionary<string, long> Counts(string name) => facets.GetProperty(name).EnumerateArray()
                .ToDictionary(f => f.GetProperty("value").GetString()!, f => f.GetProperty("count").GetInt64());
            var (kinds, providers) = (Counts("kinds"), Counts("providers"));
            var group = facets.GetProperty("kindGroups").EnumerateArray().Single(g => g.GetProperty("value").GetString() == "prokaryote");

            Assert.That(providers["api-facets"], Is.EqualTo(3));
            Assert.That(providers["user upload"], Is.EqualTo(1), "Alice's upload, not Bob's");
            Assert.That(kinds["eukaryote"], Is.GreaterThanOrEqualTo(1));
            Assert.That(group.GetProperty("count").GetInt64(), Is.EqualTo(kinds["bacterium"] + kinds["archaeon"]));
            Assert.That(facets.GetProperty("total").GetInt64(), Is.EqualTo(kinds.Values.Sum()).And.EqualTo(providers.Values.Sum()));

            var prokaryotes = (await client.GetFromJsonAsync<JsonElement>("/api/sequences?provider=api-facets&kind=prokaryote"))!;
            Assert.That(prokaryotes.GetProperty("items").EnumerateArray().Select(r => r.GetProperty("kind").GetString()),
                Is.EquivalentTo(new[] { "bacterium", "archaeon" }));
        }

        [Test]
        public async Task GenesListTheGeneLayerAndOneGenesWindows()
        {
            var (alice, client) = await NewClientAsync("alice");
            var (bob, _) = await NewClientAsync("bob");
            var genes = new[] { 60, 300, 700 }.Select(aa => "ATG" + Generator.GetRandomDNA(3 * (aa - 1), 40 + aa).Replace("TAA", "TAC")
                .Replace("TAG", "TAC").Replace("TGA", "TGC") + "TAA").ToArray();
            // The second gene sits on the reverse strand, so the record holds its reverse complement.
            var dna = genes[0] + ZCurveGeneProfile.CodingSequence(genes[1], 0, genes[1].Length, reverse: true) + genes[2];
            int a = genes[0].Length, b = a + genes[1].Length;
            var profile = ZCurveGeneProfile.FromRecord(dna, new[] { (0, a, false, 990_000_101), (a, b, true, 990_000_102), (b, dna.Length, false, -1) });
            var id = await ThrowawayDatabase.InsertRecordAsync("api-genes", "genes-1", dna.Length, kind: "bacterium");
            var upload = await ThrowawayDatabase.InsertUploadAsync(bob.Id, "genes-bob", dna.Length);
            await using (var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync())
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                await ProteinFamilies.AddAsync(connection, transaction, [(990_000_101, "test gyrase")]);
                Assert.That(await GeneProfiles.WriteAsync(connection, transaction, id, profile), Is.True);
                Assert.That(await GeneProfiles.WriteAsync(connection, transaction, upload, profile), Is.True);
                await transaction.CommitAsync();
            }

            var list = (await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}/genes"))!;
            var listed = list.GetProperty("genes").EnumerateArray().ToList();
            Assert.That((list.GetProperty("geneCount").GetInt32(), list.GetProperty("windowCount").GetInt32()), Is.EqualTo((3, profile.WindowCount)));
            Assert.That(listed.Select(g => g.GetProperty("residues").GetInt32()), Is.EqualTo(new[] { 60, 300, 700 }));
            Assert.That(listed.Select(g => g.TryGetProperty("familyId", out var f) ? f.GetString() : null), Is.EqualTo(new[] { "PGF_990000101", "PGF_990000102", null }));
            Assert.That(listed[0].GetProperty("product").GetString(), Is.EqualTo("test gyrase"));
            Assert.That(listed[1].GetProperty("strand").GetString(), Is.EqualTo("Reverse"));
            Assert.That(list.GetProperty("meanWindow").EnumerateArray().Select(v => v.GetProperty("value").GetDouble()), Is.EqualTo(profile.MeanWindow()).Within(1e-9));

            var one = (await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}/genes/1"))!;
            Assert.That(one.GetProperty("windowPositions").GetInt32(), Is.EqualTo(120), "300 aa: 5 windows of 120 positions");
            Assert.That(one.GetProperty("windows").EnumerateArray().SelectMany(w => w.EnumerateArray().Select(v => v.GetDouble())),
                Is.EqualTo(profile.WindowsOf(1)).Within(1e-9));

            Assert.That((await client.GetAsync($"/api/sequences/{id}/genes/3")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await client.GetAsync($"/api/sequences/{upload}/genes")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound), "Bob's upload");
            var none = (await client.GetFromJsonAsync<JsonElement>(
                $"/api/sequences/{await ThrowawayDatabase.InsertRecordAsync("api-genes", "genes-none", 5_000)}/genes"))!;
            Assert.That(none.GetProperty("geneCount").GetInt32(), Is.Zero);
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
                await ProfileLevels.InsertAsync(connection, transaction, mine, RecordStore.Method, OriginRotation.RotatedVersion(RecordStore.MethodVersion), 4,
                    fine.Downsample(4, 1));
                await transaction.CommitAsync();
            }
            var rotated = await client.GetAsync($"/api/sequences/{mine}/profile?windowSize=2000&origin=true");
            Assert.That(await rotated.Content.ReadAsByteArrayAsync(), Is.EqualTo(fine.Downsample(4, 1).ToHalfPrecision().ToBinary()), "the origin-rotated level");
            Assert.That((await client.GetAsync($"/api/sequences/{mine}/profile?origin=true")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound),
                "no rotated 500 bp level");

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
