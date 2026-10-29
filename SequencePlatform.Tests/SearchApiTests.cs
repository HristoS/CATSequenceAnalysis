using System.Net;
using System.Net.Http.Headers;
using System.Net.Http.Json;
using System.Text;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// The search endpoints (S4.2c). Nothing runs queued jobs until the S4.3 worker, so <see cref="RunAsync"/> does what it will:
    /// rebuild the job's query, search with the job's own parameters, save the result.
    /// </summary>
    public class SearchApiTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "api-search";

        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        private static StringContent Fasta(string id, string dna) => new($">{id}\n{dna}\n", Encoding.ASCII, "text/x-fasta");

        private static ByteArrayContent Zcwp(params ZCurveWindowedProfile[] strands)
        {
            var content = new ByteArrayContent(strands.SelectMany(s => s.ToBinary()).ToArray());
            content.Headers.ContentType = new MediaTypeHeaderValue("application/x-zcwp");
            return content;
        }

        // The ladder and sequence_sha256 a real ingest stores, under this class's method version; an upload when an owner is given.
        private static async Task<long> StoreAsync(string provider, string accession, string dna, long? ownerId = null)
        {
            var id = ownerId is long owner
                ? await ThrowawayDatabase.InsertUploadAsync(owner, accession, dna.Length, "bacterium")
                : await ThrowawayDatabase.InsertRecordAsync(provider, accession, dna.Length, "bacterium");
            var query = QueryProfile.FromSequence(dna);
            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var update = new Npgsql.NpgsqlCommand("update sequence_record set sequence_sha256 = $1 where id = $2", connection, transaction))
            {
                update.Parameters.AddWithValue(query.Sha256!);
                update.Parameters.AddWithValue(id);
                await update.ExecuteNonQueryAsync();
            }
            foreach (var factor in query.Forward.Factors)
                if (query.Forward.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, id, method, methodVersion, factor, query.Forward.Level(factor));
            await transaction.CommitAsync();
            return id;
        }

        // What the S4.3 worker will do with a queued job.
        private static async Task RunAsync(Guid id)
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var run = (await SearchLog.GetAsync(dataSource, id))!;
            var query = (await SearchLog.QueryAsync(dataSource, id))!;
            if (run.Request.Kind == SearchKind.Similar)
            {
                var search = new SimilaritySearch(dataSource, useCache: false, method: method, methodVersion: methodVersion);
                await SearchLog.SaveAsync(dataSource, id, await search.SearchAsync(query, run.Request.Filter, run.Request.LengthRatio!.Value, run.Request.Top));
                return;
            }
            var localization = new Localization(dataSource, useCache: false, method, methodVersion);
            await SearchLog.SaveAsync(dataSource, id, run.Request.WithinSequenceId is long record
                ? await localization.LocalizeInRecordAsync(query, record)
                : await localization.LocalizeAsync(query, run.Request.Filter, run.Request.Top));
        }

        private static async Task<Guid> AcceptedAsync(HttpResponseMessage response)
        {
            Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.Accepted), await response.Content.ReadAsStringAsync());
            var body = await response.Content.ReadFromJsonAsync<JsonElement>();
            var id = body.GetProperty("id").GetGuid();
            Assert.That(response.Headers.Location!.ToString(), Is.EqualTo($"/api/searches/{id}"));
            Assert.That(body.GetProperty("status").GetString(), Is.EqualTo("Queued"));
            return id;
        }

        private static async Task<long> RunCountAsync(long ownerId)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand("select count(*) from search_run where owner_id = $1");
            command.Parameters.AddWithValue(ownerId);
            return (long)(await command.ExecuteScalarAsync())!;
        }

        [Test]
        public async Task FindSimilarIsQueuedThenPagesEachHitList()
        {
            const string provider = "api-similar";
            var (_, client) = await factory.NewClientAsync("alice");
            var dna = Enumerable.Range(0, 5).Select(i => Generator.GetRandomDNA(60_000, seed: 40 + i)).ToArray();
            var ids = new List<long>();
            for (int i = 0; i < dna.Length; i++) ids.Add(await StoreAsync(provider, $"s{i}", dna[i]));
            var noLevel = await ThrowawayDatabase.InsertRecordAsync(provider, "no-level", 60_000);

            var id = await AcceptedAsync(await client.PostAsync($"/api/searches/similar?provider={provider}&top=10", Fasta("query", dna[0])));
            var queued = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{id}"))!;
            Assert.That(queued.GetProperty("status").GetString(), Is.EqualTo("Queued"));
            Assert.That(queued.TryGetProperty("hits", out _) || queued.TryGetProperty("summary", out _), Is.False);
            Assert.That(queued.GetProperty("parameters").GetProperty("provider").GetString(), Is.EqualTo(provider));
            Assert.That(queued.GetProperty("queryLengthBp").GetInt64(), Is.EqualTo(60_000));

            await RunAsync(id);

            var done = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{id}"))!;
            Assert.That(done.GetProperty("status").GetString(), Is.EqualTo("Done"));
            var summary = done.GetProperty("summary");
            Assert.That((summary.GetProperty("bandMinBp").GetInt64(), summary.GetProperty("bandMaxBp").GetInt64(), summary.GetProperty("rankedTotal").GetInt32()),
                Is.EqualTo((30_000L, 120_000L, 5)));
            var ranked = done.GetProperty("hits").GetProperty("items").EnumerateArray().ToList();
            Assert.That(ranked.Select(h => h.GetProperty("sequenceId").GetInt64()), Is.EquivalentTo(ids));
            Assert.That((ranked[0].GetProperty("sequenceId").GetInt64(), ranked[0].GetProperty("score").GetDouble()), Is.EqualTo((ids[0], 1.0)));
            Assert.That(ranked[0].GetProperty("sequence").GetProperty("recordId").GetString(), Is.EqualTo("s0"));
            Assert.That(ranked.Select(h => h.GetProperty("position").GetInt32()), Is.EqualTo(new[] { 1, 2, 3, 4, 5 }));

            var first = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{id}?limit=2"))!.GetProperty("hits");
            var second = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{id}?limit=2&after={first.GetProperty("nextAfter").GetInt32()}"))!.GetProperty("hits");
            Assert.That(second.GetProperty("items").EnumerateArray().Select(h => h.GetProperty("position").GetInt32()), Is.EqualTo(new[] { 3, 4 }));
            var identical = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{id}?list=Identical"))!.GetProperty("hits");
            Assert.That(identical.GetProperty("items").EnumerateArray().Select(h => (h.GetProperty("sequenceId").GetInt64(), h.GetProperty("strand").GetString())),
                Is.EqualTo(new[] { (ids[0], "Forward") }));
            var notComparable = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{id}?list=NotComparable"))!.GetProperty("hits");
            Assert.That(notComparable.GetProperty("items").EnumerateArray().Select(h => h.GetProperty("sequenceId").GetInt64()), Is.EqualTo(new[] { noLevel }));
        }

        [Test]
        public async Task OnlyTheOwnerSeesOrDeletesASearchAndDeleteExpiresIt()
        {
            var (_, alice) = await factory.NewClientAsync("alice");
            var (_, bob) = await factory.NewClientAsync("bob");
            var id = await AcceptedAsync(await alice.PostAsync("/api/searches/similar?provider=nobody", Fasta("q", Generator.GetRandomDNA(20_000, seed: 50))));
            await RunAsync(id);

            Assert.That((await bob.GetAsync($"/api/searches/{id}")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await bob.DeleteAsync($"/api/searches/{id}")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await alice.GetAsync($"/api/searches/{Guid.NewGuid()}")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await alice.DeleteAsync($"/api/searches/{id}")).StatusCode, Is.EqualTo(HttpStatusCode.NoContent));
            Assert.That((await alice.DeleteAsync($"/api/searches/{id}")).StatusCode, Is.EqualTo(HttpStatusCode.NoContent), "deleting twice is harmless");

            var expired = (await alice.GetFromJsonAsync<JsonElement>($"/api/searches/{id}"))!;
            Assert.That(expired.GetProperty("status").GetString(), Is.EqualTo("Expired"));
            Assert.That(expired.TryGetProperty("hits", out _), Is.False);
            Assert.That(expired.GetProperty("summary").GetProperty("bandMinBp").GetInt64(), Is.EqualTo(10_000), "the summary outlives the results");
            Assert.That(await SearchLog.QueryAsync(ThrowawayDatabase.DataSource, id), Is.Null);
        }

        [Test]
        public async Task LocalizePlacesAFragmentOverTheCatalogueOrInsideAVisibleRecord()
        {
            const string provider = "api-localize";
            var (alice, client) = await factory.NewClientAsync("alice");
            var (bob, _) = await factory.NewClientAsync("bob");
            var genomes = Enumerable.Range(0, 3).Select(i => Generator.GetRandomDNA(200_000, seed: 60 + i)).ToArray();
            var ids = new List<long>();
            for (int i = 0; i < genomes.Length; i++) ids.Add(await StoreAsync(provider, $"l{i}", genomes[i]));
            var mine = await StoreAsync(provider, "mine", genomes[2], alice.Id);
            var theirs = await StoreAsync(provider, "theirs", genomes[2], bob.Id);

            var overCatalogue = await AcceptedAsync(await client.PostAsync($"/api/searches/localize?provider={provider}&top=3",
                Fasta("piece", genomes[1].Substring(50_000, 10_000))));
            var inRecord = await AcceptedAsync(await client.PostAsync($"/api/searches/localize?record={mine}",
                Fasta("short", genomes[2].Substring(120_000, 2_000))));
            var tooShort = await client.PostAsync($"/api/searches/localize?provider={provider}", Fasta("short", genomes[2].Substring(120_000, 2_000)));
            var hidden = await client.PostAsync($"/api/searches/localize?record={theirs}", Fasta("short", genomes[2].Substring(120_000, 2_000)));
            await RunAsync(overCatalogue);
            await RunAsync(inRecord);

            var top = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{overCatalogue}"))!.GetProperty("hits").GetProperty("items")[0];
            Assert.That((top.GetProperty("sequenceId").GetInt64(), top.GetProperty("strand").GetString()), Is.EqualTo((ids[1], "Forward")));
            Assert.That(top.GetProperty("offsetBp").GetInt64(), Is.EqualTo(50_000).Within(500));
            Assert.That(top.TryGetProperty("coarseScore", out _), Is.True);
            var placed = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{inRecord}"))!;
            Assert.That(placed.GetProperty("parameters").GetProperty("withinSequenceId").GetInt64(), Is.EqualTo(mine));
            Assert.That(placed.GetProperty("hits").GetProperty("items")[0].GetProperty("offsetBp").GetInt64(), Is.EqualTo(120_000).Within(500));
            Assert.That(tooShort.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
            Assert.That((await tooShort.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("detail").GetString(), Does.Contain("name a record"));
            Assert.That(hidden.StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
        }

        [Test]
        public async Task ProfileQueryCarriesOneStrandOrBothBackToBack()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var dna = Generator.GetRandomDNA(30_250, seed: 70);
            var forward = new ZCurveWindowedProfile(dna);
            var reverse = new ZCurveWindowedProfile(QueryProfile.ReverseComplement(dna));

            var one = await AcceptedAsync(await client.PostAsync("/api/searches/similar", Zcwp(forward)));
            var both = await AcceptedAsync(await client.PostAsync("/api/searches/similar?lengthBp=30250", Zcwp(forward, reverse)));
            var mismatched = await client.PostAsync("/api/searches/similar", Zcwp(forward, new ZCurveWindowedProfile(dna[..20_000])));
            var wrongLength = await client.PostAsync("/api/searches/similar?lengthBp=40000", Zcwp(forward));
            var coarse = await client.PostAsync("/api/searches/similar", Zcwp(new ZCurveWindowedProfile(dna, 2_000)));
            var three = await client.PostAsync("/api/searches/similar", Zcwp(forward, reverse, forward));

            var oneQuery = (await SearchLog.QueryAsync(ThrowawayDatabase.DataSource, one))!;
            var bothQuery = (await SearchLog.QueryAsync(ThrowawayDatabase.DataSource, both))!;
            Assert.That((oneQuery.Reverse, oneQuery.LengthBp), Is.EqualTo(((ZCurveProfileLadder?)null, 30_000L)), "windows × 500 without lengthBp");
            Assert.That(bothQuery.LengthBp, Is.EqualTo(30_250));
            Assert.That(bothQuery.Reverse!.Fine.Windows, Is.EqualTo(reverse.Windows));
            foreach (var rejected in new[] { mismatched, wrongLength, coarse, three })
                Assert.That(rejected.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest), await rejected.Content.ReadAsStringAsync());
        }

        [Test]
        public async Task InvalidSearchesAreRejectedBeforeAnyJobIsQueued()
        {
            var (alice, client) = await factory.NewClientAsync("alice");
            var dna = Generator.GetRandomDNA(10_000, seed: 80);

            var cases = new (string Route, HttpContent Body, HttpStatusCode Status)[]
            {
                ("similar?ratio=0.1", Fasta("q", dna), HttpStatusCode.BadRequest),
                ("similar?top=0", Fasta("q", dna), HttpStatusCode.BadRequest),
                ("similar?top=1001", Fasta("q", dna), HttpStatusCode.BadRequest),
                ("similar?rank=genus", Fasta("q", dna), HttpStatusCode.BadRequest),
                ("similar?rank=tribe&taxon=X", Fasta("q", dna), HttpStatusCode.BadRequest),
                ("similar", new StringContent($">a\n{dna}\n>b\n{dna}\n", Encoding.ASCII, "text/x-fasta"), HttpStatusCode.BadRequest),
                ("similar", new StringContent(">a\n", Encoding.ASCII, "text/x-fasta"), HttpStatusCode.BadRequest),
                ("similar", new StringContent(dna, Encoding.ASCII, "application/json"), HttpStatusCode.UnsupportedMediaType),
                ("localize?top=101", Fasta("q", dna), HttpStatusCode.BadRequest),
                ("localize?record=999999999", Fasta("q", dna), HttpStatusCode.NotFound),
            };
            foreach (var (route, body, status) in cases)
            {
                var response = await client.PostAsync("/api/searches/" + route, body);
                Assert.That(response.StatusCode, Is.EqualTo(status), route);
                Assert.That(response.Content.Headers.ContentType!.MediaType, Is.EqualTo("application/problem+json"), route);
            }
            Assert.That(await RunCountAsync(alice.Id), Is.Zero);
        }
    }
}
