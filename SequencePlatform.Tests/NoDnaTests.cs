using System.Diagnostics;
using System.Net;
using System.Net.Http.Json;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// The platform never stores DNA (S4.4): a genome carrying a unique 1 kbp marker goes through every route that takes sequence
    /// (upload, find-similar, localize, rejected bodies), and afterwards no column of any table, no log entry and no temp file holds
    /// any piece of the marker. The API runs on Kestrel here, not the test server, so its request size limit is the real one. A
    /// marked protein goes through protein search (P7) the same way.
    /// </summary>
    public class NoDnaTests
    {
        private const int markerLength = 1_000;
        // Any 47 consecutive bases of the marker hold a whole aligned piece, so a copy wrapped at 60 or more per line still shows.
        private const int pieceLength = 24;

        private static readonly string[] tempVariables = ["TMP", "TEMP", "TMPDIR", "ASPNETCORE_TEMP"];

        private readonly LogCapture logs = new();
        private readonly List<string> tempFiles = [];
        private readonly Dictionary<string, string?> savedVariables = [];
        private string tempDirectory = null!;
        private FileSystemWatcher watcher = null!;
        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start()
        {
            // Everything that asks .NET or ASP.NET for a temp path while the API runs lands in this directory, which is watched,
            // so a file written and deleted within one request (a buffered body) is caught too.
            tempDirectory = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"no-dna-temp-{Guid.NewGuid():N}");
            Directory.CreateDirectory(tempDirectory);
            foreach (var variable in tempVariables)
            {
                savedVariables[variable] = Environment.GetEnvironmentVariable(variable);
                Environment.SetEnvironmentVariable(variable, tempDirectory);
            }
            watcher = new FileSystemWatcher(tempDirectory) { IncludeSubdirectories = true };
            watcher.Created += (_, e) => { lock (tempFiles) tempFiles.Add(e.FullPath); };
            watcher.Renamed += (_, e) => { lock (tempFiles) tempFiles.Add(e.FullPath); };
            watcher.EnableRaisingEvents = true;

            factory = ApiHost.Create(RecordStore.MethodVersion, logs);
            factory.UseKestrel(0);
            factory.StartServer();
        }

        [OneTimeTearDown]
        public void Stop()
        {
            factory.Dispose();
            watcher.Dispose();
            foreach (var (variable, value) in savedVariables) Environment.SetEnvironmentVariable(variable, value);
            Directory.Delete(tempDirectory, recursive: true);
        }

        private static StringContent Fasta(string text) => new(text, Encoding.ASCII, "text/x-fasta");

        private static string Wrap(string dna, int width) =>
            string.Join("\n", Enumerable.Range(0, (dna.Length + width - 1) / width).Select(i => dna.Substring(i * width, Math.Min(width, dna.Length - i * width))));

        private static string ReverseComplement(string dna) =>
            string.Concat(dna.Reverse().Select(b => b switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => b }));

        private static List<string> Pieces(string marker) =>
            new[] { marker, ReverseComplement(marker) }
                .SelectMany(strand => Enumerable.Range(0, strand.Length / pieceLength).Select(i => strand.Substring(i * pieceLength, pieceLength)))
                .ToList();

        private static bool Holds(string text, IEnumerable<string> pieces) => pieces.Any(p => text.Contains(p, StringComparison.OrdinalIgnoreCase));

        /// <summary>Every column of every table that holds one of the texts, ignoring case: bytea by its bytes, any other type as text.</summary>
        private static async Task<List<string>> ColumnsHoldingAsync(IEnumerable<string> texts)
        {
            var columns = new List<(string Table, string Column, string Type)>();
            await using (var command = ThrowawayDatabase.DataSource.CreateCommand(
                "select c.table_name, c.column_name, c.data_type from information_schema.columns c " +
                "join information_schema.tables t on t.table_schema = c.table_schema and t.table_name = c.table_name " +
                "where c.table_schema = 'public' and t.table_type = 'BASE TABLE' order by 1, 2"))
            await using (var reader = await command.ExecuteReaderAsync())
                while (await reader.ReadAsync()) columns.Add((reader.GetString(0), reader.GetString(1), reader.GetString(2)));

            var patterns = texts.Select(t => $"%{t}%").ToArray();
            var holding = new List<string>();
            foreach (var (table, column, type) in columns)
            {
                var value = type == "bytea" ? $"encode(\"{column}\", 'escape')" : $"\"{column}\"::text";
                await using var command = ThrowawayDatabase.DataSource.CreateCommand($"select exists (select 1 from \"{table}\" where {value} ilike any($1))");
                command.Parameters.AddWithValue(patterns);
                if ((bool)(await command.ExecuteScalarAsync())!) holding.Add($"{table}.{column}");
            }
            return holding;
        }

        private static async Task<Guid> PostSearchAsync(HttpClient client, string route, string fasta)
        {
            var response = await client.PostAsync(route, Fasta(fasta));
            Assert.That((int)response.StatusCode, Is.EqualTo(202), await response.Content.ReadAsStringAsync());
            return (await response.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("id").GetGuid();
        }

        private static async Task<SearchRun> DoneAsync(Guid id)
        {
            var watch = Stopwatch.StartNew();
            while (watch.Elapsed < TimeSpan.FromSeconds(60))
            {
                var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!;
                if (run.Status is SearchStatus.Done or SearchStatus.Failed)
                {
                    Assert.That(run.Status, Is.EqualTo(SearchStatus.Done), run.Error);
                    return run;
                }
                await Task.Delay(50);
            }
            Assert.Fail($"Search {id} didn't finish within 60 s.");
            return null!;
        }

        private void AssertNoTempFiles()
        {
            Assert.That(Path.GetTempPath().TrimEnd(Path.DirectorySeparatorChar), Is.EqualTo(tempDirectory), "the temp path is redirected");
            lock (tempFiles) Assert.That(tempFiles, Is.Empty, "files created in the temp directory");
            Assert.That(Directory.EnumerateFileSystemEntries(tempDirectory), Is.Empty);
        }

        [Test]
        public async Task AMarkedGenomeLeavesNoTraceInTheDatabaseTheLogsOrTheTempDirectory()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var marker = Generator.GetRandomDNA(markerLength, seed: 4401);
            // 121 kbp, past every buffering threshold (64 KB for a form file, 30 KB for request buffering).
            var genome = Generator.GetRandomDNA(60_000, seed: 4402) + marker + Generator.GetRandomDNA(60_000, seed: 4403);
            var recordId = $"marked-{Guid.NewGuid():N}";
            var fasta = $">{recordId} marked genome\n{Wrap(genome, 60)}\n";

            var upload = await client.PostAsync("/api/sequences?kind=bacterium", Fasta(fasta));
            Assert.That(upload.StatusCode, Is.EqualTo(HttpStatusCode.Created), await upload.Content.ReadAsStringAsync());
            var id = (await upload.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("records")[0].GetProperty("id").GetInt64();
            var similar = await PostSearchAsync(client, "/api/searches/similar?kind=bacterium", $">q\n{genome}\n");
            var inRecord = await PostSearchAsync(client, $"/api/searches/localize?record={id}", $">q\n{Wrap(marker, 80)}\n");
            var overCatalogue = await PostSearchAsync(client, "/api/searches/localize?kind=bacterium", $">q\n{genome.Substring(58_000, 5_000)}\n");
            // Rejected bodies whose error could echo what they carry.
            var noHeader = await client.PostAsync("/api/sequences?kind=virus", Fasta(marker + "\n"));
            var twoRecords = await client.PostAsync("/api/searches/similar", Fasta($">a\n{marker}\n>b\n{marker}\n"));
            await DoneAsync(similar);
            await DoneAsync(inRecord);
            await DoneAsync(overCatalogue);

            var placed = await SearchLog.HitsAsync(ThrowawayDatabase.DataSource, inRecord, HitList.Ranked);
            Assert.That((placed[0].SequenceId, placed[0].OffsetBp), Is.EqualTo((id, (long?)60_000)), "the searches did run on the marker");
            foreach (var rejected in new[] { noHeader, twoRecords })
            {
                Assert.That(rejected.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
                Assert.That(Holds(await rejected.Content.ReadAsStringAsync(), Pieces(marker)), Is.False, "the error repeats the body");
            }

            var pieces = Pieces(marker);
            Assert.That(await ColumnsHoldingAsync(pieces), Is.Empty, "columns holding a piece of the marker");
            Assert.That(await ColumnsHoldingAsync([recordId]), Does.Contain("sequence_record.record_id"), "the scan finds what is stored");

            var entries = logs.Entries;
            Assert.That(entries.Where(e => Holds(e, pieces)), Is.Empty, "log entries holding a piece of the marker");
            Assert.That(entries, Has.Some.Contains($"sha256 {Convert.ToHexStringLower(SHA256.HashData(Encoding.ASCII.GetBytes(fasta)))}")
                .And.Contains($"1 record(s), {genome.Length} bp"), "the upload is logged by hash and length");

            AssertNoTempFiles();
        }

        [Test]
        public async Task AMarkedProteinLeavesNoTraceInTheDatabaseTheLogsOrTheTempDirectory()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var marker = ZCurveGeneProfile.Translate(ProteinSearchTests.RandomGene(markerLength, seed: 4405));
            var wrapped = await PostSearchAsync(client, "/api/searches/protein", $">marked protein\n{Wrap(marker, 60)}\n");
            var bare = await client.PostAsync("/api/searches/protein", new StringContent(marker, Encoding.ASCII, "text/plain"));
            Assert.That((int)bare.StatusCode, Is.EqualTo(202), await bare.Content.ReadAsStringAsync());
            // Rejected bodies whose error could echo what they carry.
            var twoRecords = await client.PostAsync("/api/searches/protein", Fasta($">a\n{marker}\n>b\n{marker}\n"));
            var internalStop = await client.PostAsync("/api/searches/protein", Fasta($">a\n{marker[..500]}*{marker[500..]}\n"));
            await DoneAsync(wrapped);
            await DoneAsync((await bare.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("id").GetGuid());

            var pieces = Enumerable.Range(0, marker.Length / pieceLength).Select(i => marker.Substring(i * pieceLength, pieceLength)).ToList();
            foreach (var rejected in new[] { twoRecords, internalStop })
            {
                Assert.That(rejected.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
                Assert.That(Holds(await rejected.Content.ReadAsStringAsync(), pieces), Is.False, "the error repeats the body");
            }
            Assert.That(await ColumnsHoldingAsync(pieces), Is.Empty, "columns holding a piece of the marker");
            var entries = logs.Entries;
            Assert.That(entries.Where(e => Holds(e, pieces)), Is.Empty, "log entries holding a piece of the marker");
            Assert.That(entries, Has.Some.Contains($"{markerLength} aa, code 11, sha256 {Convert.ToHexStringLower(SHA256.HashData(Encoding.ASCII.GetBytes(marker)))}"),
                "the protein is logged by length and hash");
            AssertNoTempFiles();
        }

        [Test]
        public async Task ABodyOverTheSizeLimitIsRefusedAndStoresNothing()
        {
            var (alice, client) = await factory.NewClientAsync("alice");
            var lines = Wrap(Generator.GetRandomDNA(1_000_000, seed: 4404), 100);
            var body = new StringBuilder(">too-big\n");
            while (body.Length <= 32 * 1024 * 1024) body.Append(lines).Append('\n');

            var response = await client.PostAsync("/api/sequences?kind=eukaryote", Fasta(body.ToString()));

            Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.RequestEntityTooLarge));
            await using var command = ThrowawayDatabase.DataSource.CreateCommand("select count(*) from sequence_source where owner_id = $1");
            command.Parameters.AddWithValue(alice.Id);
            Assert.That((long)(await command.ExecuteScalarAsync())!, Is.Zero);
            AssertNoTempFiles();
        }
    }
}
