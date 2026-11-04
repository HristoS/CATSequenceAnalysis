using System.Net;
using System.Net.Http.Json;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>P9 experiment results (U0): a publish replaces the run of its name, findings only about public records, served read-only.</summary>
    public class ExperimentTests
    {
        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        private static Task<(int Stored, int Skipped)> PublishAsync(string name, params ExperimentFinding[] findings) =>
            ExperimentStore.PublishAsync(ThrowawayDatabase.DataSource, name, "A test", "Does it work?", ExperimentStore.Pass, "snapshot_test",
                "test.cs", "{\"rule\": \"x < 1\", \"keyNumbers\": [{\"label\": \"x\", \"value\": 0.5}]}", findings);

        private static ExperimentFinding Finding(long id, double score, string kind = "cpg") => new(id, kind, score, $"{{\"class\": \"Testclass\", \"id\": {id}}}");

        [Test]
        public async Task APublishReplacesThePreviousRunAndSkipsRecordsThatArentPublic()
        {
            var name = $"test-{Guid.NewGuid():N}";
            var (user, _) = await factory.NewClientAsync("alice");
            var first = await ThrowawayDatabase.InsertRecordAsync("experiments", $"{name}-1", 20_000);
            var second = await ThrowawayDatabase.InsertRecordAsync("experiments", $"{name}-2", 20_000);
            var upload = await ThrowawayDatabase.InsertUploadAsync(user.Id, $"{name}-upload", 20_000);

            var published = await PublishAsync(name, Finding(first, 0.3), Finding(second, 0.9), Finding(upload, 0.1), Finding(long.MaxValue, 0.2));
            Assert.That(published, Is.EqualTo((2, 2)), "the upload and the missing record are skipped");

            await PublishAsync(name, Finding(second, 0.7));
            var run = await ExperimentStore.GetAsync(ThrowawayDatabase.DataSource, name);
            var (items, total) = await ExperimentStore.FindingsAsync(ThrowawayDatabase.DataSource, run!.Id, null, false, 0, 10);

            Assert.That(total, Is.EqualTo(1));
            Assert.That(items.Single().SequenceId, Is.EqualTo(second));
            Assert.That(items.Single().Score, Is.EqualTo(0.7));
            Assert.That(await ExperimentStore.ForRecordAsync(ThrowawayDatabase.DataSource, first), Is.Empty, "the old run's findings went with it");
        }

        [Test]
        public async Task ServesTheRunItsFindingsInScoreOrderAndARecordsFindings()
        {
            var name = $"test-{Guid.NewGuid():N}";
            var (_, client) = await factory.NewClientAsync("alice");
            var ids = new List<long>();
            for (int i = 0; i < 3; i++) ids.Add(await ThrowawayDatabase.InsertRecordAsync("experiments", $"{name}-{i}", 20_000));
            await PublishAsync(name, Finding(ids[0], 0.5), Finding(ids[1], 0.1), Finding(ids[2], 0.9), Finding(ids[2], 3.0, "other"));

            var list = await client.GetFromJsonAsync<JsonElement>("/api/experiments");
            Assert.That(list.GetProperty("items").EnumerateArray().Select(e => e.GetProperty("name").GetString()), Does.Contain(name));

            var page = await client.GetFromJsonAsync<JsonElement>($"/api/experiments/{name}?kind=cpg&order=-score&limit=2");
            Assert.That(page.GetProperty("total").GetInt64(), Is.EqualTo(3));
            Assert.That(page.GetProperty("items").EnumerateArray().Select(i => i.GetProperty("sequence").GetProperty("id").GetInt64()),
                Is.EqualTo(new[] { ids[2], ids[0] }));
            Assert.That(page.GetProperty("experiment").GetProperty("summary").GetProperty("rule").GetString(), Is.EqualTo("x < 1"));
            Assert.That(page.GetProperty("items")[0].GetProperty("payload").GetProperty("class").GetString(), Is.EqualTo("Testclass"));

            var own = await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{ids[2]}/findings");
            Assert.That(own.GetProperty("items").EnumerateArray().Select(i => i.GetProperty("kind").GetString()), Is.EquivalentTo(new[] { "cpg", "other" }));
        }

        [Test]
        public async Task NeedsAKeyAndAnswers404And400Plainly()
        {
            var (bob, bobClient) = await factory.NewClientAsync("bob");
            var (_, aliceClient) = await factory.NewClientAsync("alice");
            var bobsUpload = await ThrowawayDatabase.InsertUploadAsync(bob.Id, $"bob-{Guid.NewGuid():N}", 20_000);

            Assert.That((await factory.CreateClient().GetAsync("/api/experiments")).StatusCode, Is.EqualTo(HttpStatusCode.Unauthorized));
            Assert.That((await aliceClient.GetAsync("/api/experiments/no-such-experiment")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await aliceClient.GetAsync($"/api/sequences/{bobsUpload}/findings")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await bobClient.GetAsync($"/api/sequences/{bobsUpload}/findings")).StatusCode, Is.EqualTo(HttpStatusCode.OK));
            Assert.That((await aliceClient.GetAsync("/api/experiments/x?order=up")).StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
            Assert.That((await aliceClient.GetAsync("/api/experiments/x?limit=5000")).StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
        }
    }
}
