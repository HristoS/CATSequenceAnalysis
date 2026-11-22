using System.Diagnostics;
using System.Net;
using System.Net.Http.Json;
using System.Text;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;
using Npgsql;
using SequencePlatform.Core;
using SequencePlatform.Experiments;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// Experiment reruns (migration 0024): an admin queues one, the experiment worker runs it and publishes the server run beside the article
    /// run, as the least-privilege role, which can't touch article runs. The worker runs stand-in experiments here; the real ones need a
    /// catalogue (their port is checked against the scripts on the local database, results/README.md).
    /// </summary>
    public class ExperimentRerunTests
    {
        private const string insufficientPrivilege = "42501";

        private readonly BlockingExperiment blocking = new($"test-block-{Guid.NewGuid():N}");
        private readonly ThrowingExperiment explained = new($"test-explained-{Guid.NewGuid():N}", new ExperimentNotRunnableException("No query set is stored."));
        private readonly ThrowingExperiment unexplained = new($"test-unexplained-{Guid.NewGuid():N}", new FormatException("internal detail"));
        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        // Runs until the test releases it, so a test can see it running.
        private sealed class BlockingExperiment(string name) : IExperiment
        {
            public readonly TaskCompletionSource Release = new(TaskCreationOptions.RunContinuationsAsynchronously);

            public string Name => name;

            public async Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken)
            {
                await Release.Task.WaitAsync(cancellationToken);
                return new ExperimentResult(ExperimentStore.Pass, "live database, test", "test.cs", "{\"keyNumbers\": [{\"label\": \"x\", \"value\": 0.25}]}",
                    [new ExperimentFile("rerun.csv", Encoding.UTF8.GetBytes("a,b\n1,2\n"))], []);
            }
        }

        private sealed class ThrowingExperiment(string name, Exception error) : IExperiment
        {
            public string Name => name;

            public Task<ExperimentResult> RunAsync(NpgsqlDataSource dataSource, ExperimentContext context, CancellationToken cancellationToken) =>
                Task.FromException<ExperimentResult>(error);
        }

        [OneTimeSetUp]
        public async Task Start()
        {
            factory = ApiHost.Create(experiments: [blocking, explained, unexplained]);
            foreach (var name in new[] { blocking.Name, explained.Name, unexplained.Name }) await PublishArticleRunAsync(name, runnable: true);
        }

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        private static Task PublishArticleRunAsync(string name, bool runnable) =>
            ExperimentStore.PublishAsync(ThrowawayDatabase.DataSource, name, "A rerun test", "Does it rerun?", ExperimentStore.Fail, "snapshot_test", "test.cs",
                "{\"keyNumbers\": [{\"label\": \"x\", \"value\": 0.5}]}", [], new ExperimentDetails("B", Runnable: runnable),
                [new ExperimentFile("article.csv", Encoding.UTF8.GetBytes("a\n1\n"))]);

        private async Task<HttpClient> AdminAsync()
        {
            var (user, client) = await factory.NewClientAsync("admin");
            Assert.That(await Users.SetAdminAsync(ThrowawayDatabase.DataSource, user.Name, true), Is.True);
            return client;
        }

        private static async Task<JsonElement> WaitForAsync(HttpClient client, string href, params string[] statuses)
        {
            var watch = Stopwatch.StartNew();
            while (watch.Elapsed < TimeSpan.FromSeconds(60))
            {
                var job = await client.GetFromJsonAsync<JsonElement>(href);
                if (statuses.Contains(job.GetProperty("status").GetString())) return job;
                await Task.Delay(50);
            }
            Assert.Fail($"{href} didn't reach {string.Join(" or ", statuses)} within 60 s.");
            return default;
        }

        [Test]
        public async Task OnlyAnAdminQueuesARerunAndOnlyOfARunnableDatabaseExperiment()
        {
            var (_, member) = await factory.NewClientAsync("member");
            var admin = await AdminAsync();
            var notRunnable = $"test-offline-{Guid.NewGuid():N}";
            await PublishArticleRunAsync(notRunnable, runnable: false);

            Assert.That((await member.PostAsync($"/api/experiments/{blocking.Name}/runs", null)).StatusCode, Is.EqualTo(HttpStatusCode.Forbidden));
            var reruns = await member.GetFromJsonAsync<JsonElement>($"/api/experiments/{blocking.Name}/runs");
            Assert.That(reruns.GetProperty("rerunnable").GetBoolean(), Is.True);
            Assert.That(reruns.GetProperty("canRun").GetBoolean(), Is.False, "a member sees it reruns, but may not queue it");

            Assert.That((await admin.PostAsync($"/api/experiments/{notRunnable}/runs", null)).StatusCode, Is.EqualTo(HttpStatusCode.NotFound),
                "its article run needs DNA");
            Assert.That((await admin.PostAsync("/api/experiments/no-such-experiment/runs", null)).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That((await admin.GetFromJsonAsync<JsonElement>($"/api/experiments/{notRunnable}/runs")).GetProperty("canRun").GetBoolean(), Is.False);
        }

        [Test]
        public async Task ARerunRunsOnceAtATimeAndIsPublishedBesideTheArticleRun()
        {
            var admin = await AdminAsync();
            var response = await admin.PostAsync($"/api/experiments/{blocking.Name}/runs", null);
            Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.Accepted), await response.Content.ReadAsStringAsync());
            var href = response.Headers.Location!.ToString();
            await WaitForAsync(admin, href, "Running");

            var second = await admin.PostAsync($"/api/experiments/{blocking.Name}/runs", null);
            Assert.That(second.StatusCode, Is.EqualTo(HttpStatusCode.Conflict), "one queued or running rerun per experiment");

            blocking.Release.SetResult();
            var job = await WaitForAsync(admin, href, "Done", "Failed");
            Assert.That(job.GetProperty("status").GetString(), Is.EqualTo("Done"), job.ToString());
            Assert.That(job.GetProperty("runHref").GetString(), Is.EqualTo($"/api/experiments/{blocking.Name}?origin=server"));

            var server = (await admin.GetFromJsonAsync<JsonElement>(job.GetProperty("runHref").GetString()!)).GetProperty("experiment");
            var article = (await admin.GetFromJsonAsync<JsonElement>($"/api/experiments/{blocking.Name}")).GetProperty("experiment");
            Assert.That(server.GetProperty("origin").GetString(), Is.EqualTo(ExperimentStore.ServerOrigin));
            Assert.That(server.GetProperty("verdict").GetString(), Is.EqualTo(ExperimentStore.Pass));
            Assert.That(server.GetProperty("title").GetString(), Is.EqualTo("A rerun test"), "the article run's title and question");
            Assert.That(server.GetProperty("article").GetString(), Is.EqualTo("B"));
            Assert.That(server.GetProperty("runnable").GetBoolean(), Is.True);
            Assert.That(server.GetProperty("files")[0].GetProperty("name").GetString(), Is.EqualTo("rerun.csv"));
            Assert.That(article.GetProperty("verdict").GetString(), Is.EqualTo(ExperimentStore.Fail), "the article run is untouched");
            Assert.That(article.GetProperty("files")[0].GetProperty("name").GetString(), Is.EqualTo("article.csv"));

            var reruns = await admin.GetFromJsonAsync<JsonElement>($"/api/experiments/{blocking.Name}/runs");
            Assert.That(reruns.GetProperty("canRun").GetBoolean(), Is.True);
            Assert.That(reruns.GetProperty("jobs")[0].GetProperty("id").GetString(), Is.EqualTo(job.GetProperty("id").GetString()));
        }

        [Test]
        public async Task AFailedRerunSaysWhyOnlyWhenTheReasonIsWrittenForTheAdmin()
        {
            var admin = await AdminAsync();
            async Task<JsonElement> RunAsync(string name)
            {
                var response = await admin.PostAsync($"/api/experiments/{name}/runs", null);
                Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.Accepted));
                return await WaitForAsync(admin, response.Headers.Location!.ToString(), "Done", "Failed");
            }

            var withReason = await RunAsync(explained.Name);
            var withoutReason = await RunAsync(unexplained.Name);
            Assert.That(withReason.GetProperty("status").GetString(), Is.EqualTo("Failed"));
            Assert.That(withReason.GetProperty("error").GetString(), Is.EqualTo("No query set is stored."));
            Assert.That(withoutReason.GetProperty("error").GetString(), Does.Not.Contain("internal detail"));
            Assert.That((await admin.PostAsync($"/api/experiments/{explained.Name}/runs", null)).StatusCode, Is.EqualTo(HttpStatusCode.Accepted),
                "a finished job no longer blocks the next");
        }

        [Test]
        public async Task TheAppRolePublishesServerRunsButNeverTouchesArticleRuns()
        {
            var name = $"test-role-{Guid.NewGuid():N}";
            await PublishArticleRunAsync(name, runnable: true);
            var record = await ThrowawayDatabase.InsertRecordAsync("experiments", $"{name}-r", 20_000);
            var server = new ExperimentDetails(Origin: ExperimentStore.ServerOrigin);
            for (int i = 0; i < 2; i++)
                await ExperimentStore.PublishAsync(ThrowawayDatabase.AppDataSource, name, "t", "q", ExperimentStore.Pass, "s", "x.cs", "{}",
                    [new ExperimentFinding(record, "kind", 1.0, "{}")], server, [new ExperimentFile("f.csv", [1, 2])]);
            Assert.That((await ExperimentStore.GetAsync(ThrowawayDatabase.DataSource, name, ExperimentStore.ServerOrigin))!.Findings, Is.EqualTo(1),
                "a second publish replaced the first, findings and files with it");

            var error = Assert.ThrowsAsync<PostgresException>(() => ExperimentStore.PublishAsync(ThrowawayDatabase.AppDataSource, name, "t", "q",
                ExperimentStore.Pass, "s", "x.cs", "{}", []));
            Assert.That(error!.SqlState, Is.EqualTo(insufficientPrivilege), error.MessageText);

            var articleId = (await ExperimentStore.GetAsync(ThrowawayDatabase.DataSource, name))!.Id;
            foreach (var sql in new[]
            {
                $"insert into experiment_file (run_id, name, byte_size, content) values ({articleId}, 'x.csv', 0, ''::bytea)",
                $"insert into experiment_finding (run_id, sequence_id, kind, payload) values ({articleId}, {record}, 'k', '{{}}')",
            })
            {
                await using var insert = ThrowawayDatabase.AppDataSource.CreateCommand(sql);
                var denied = Assert.ThrowsAsync<PostgresException>(async () => await insert.ExecuteNonQueryAsync());
                Assert.That(denied!.SqlState, Is.EqualTo(insufficientPrivilege), sql);
            }
            await using var delete = ThrowawayDatabase.AppDataSource.CreateCommand($"delete from experiment_run where id = {articleId}");
            Assert.That(await delete.ExecuteNonQueryAsync(), Is.EqualTo(0), "an article run is invisible to the role's delete");
            Assert.That(await ExperimentStore.GetAsync(ThrowawayDatabase.DataSource, name), Is.Not.Null);
        }

        [Test]
        public async Task AStoredQuerySetRoundTripsInOrder()
        {
            var name = $"test-queries-{Guid.NewGuid():N}";
            var record = await ThrowawayDatabase.InsertRecordAsync("experiments", $"{name}-r", 20_000);
            var fields = new ProteinHomolog.QueryFields("bacterium", 120, 10, 370, true);
            await ExperimentQueries.ReplaceAsync(ThrowawayDatabase.DataSource, name,
                [ProteinHomolog.ToQuery(record, fields, [0.5, -0.25]), ProteinHomolog.ToQuery(record, fields with { Start = 400 }, [1.0])]);

            var read = await ExperimentQueries.ReadAsync(ThrowawayDatabase.DataSource, name);
            Assert.That(read.Select(q => q.Fields.Deserialize<ProteinHomolog.QueryFields>()!.Start), Is.EqualTo(new[] { 10, 400 }));
            Assert.That(read[0].Profile, Is.EqualTo(new[] { 0.5, -0.25 }));
            Assert.That(read[0].Fields.Deserialize<ProteinHomolog.QueryFields>(), Is.EqualTo(fields));
            Assert.That(await ExperimentQueries.ReadAsync(ThrowawayDatabase.AppDataSource, name), Has.Count.EqualTo(2), "the role reads them");
        }

        [Test]
        public void TheRankSumAucAveragesTies()
        {
            Assert.That(GcStratum.Auc([(0.1, false), (0.2, true), (0.3, false), (0.4, true)]).Auc, Is.EqualTo(0.75));
            Assert.That(GcStratum.Auc([(0.5, false), (0.5, true), (0.5, false), (0.5, true)]).Auc, Is.EqualTo(0.5));
            Assert.That(double.IsNaN(GcStratum.Auc([(0.1, true)]).Auc), Is.True, "no negatives");
        }
    }
}
