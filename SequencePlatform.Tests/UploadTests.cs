using System.Net;
using System.Net.Http.Headers;
using System.Net.Http.Json;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>POST /api/sequences (S4.2b): FASTA and ZCWP uploads become their owner's private records, profile only.</summary>
    public class UploadTests
    {
        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        private static StringContent Fasta(string text) => new(text, Encoding.ASCII, "text/x-fasta");

        private static ByteArrayContent Zcwp(byte[] bytes)
        {
            var content = new ByteArrayContent(bytes);
            content.Headers.ContentType = new MediaTypeHeaderValue("application/x-zcwp");
            return content;
        }

        private static string Wrap(string dna, int width) =>
            string.Join("\n", Enumerable.Range(0, (dna.Length + width - 1) / width).Select(i => dna.Substring(i * width, Math.Min(width, dna.Length - i * width))));

        private static async Task<long> UploadCountAsync(long ownerId)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand("select count(*) from sequence_source where owner_id = $1");
            command.Parameters.AddWithValue(ownerId);
            return (long)(await command.ExecuteScalarAsync())!;
        }

        [Test]
        public async Task FastaUploadStoresEachRecordAsTheIngestWouldVisibleToItsOwnerOnly()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var (_, bob) = await factory.NewClientAsync("bob");
            var dnaA = Generator.GetRandomDNA(12_345, seed: 31);
            var dnaB = Generator.GetRandomDNA(3_000, seed: 32);
            var text = $">contigA first contig\r\n{Wrap(dnaA, 70).Replace("\n", "\r\n")}\r\n>contigB\n{dnaB.ToLowerInvariant()}\n";
            // What the Ingest's path-based reader makes of the same file.
            var path = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"upload-{Guid.NewGuid():N}.fna");
            File.WriteAllText(path, text, Encoding.ASCII);
            List<ProfiledRecord> expected;
            try { expected = RecordStore.Profile(path).ToList(); }
            finally { File.Delete(path); }

            var response = await client.PostAsync("/api/sequences?kind=bacterium&taxid=562", Fasta(text));

            Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.Created));
            var records = (await response.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("records").EnumerateArray().ToList();
            Assert.That(records.Select(r => (r.GetProperty("recordId").GetString(), r.GetProperty("lengthBp").GetInt64(), r.GetProperty("windowCount").GetInt32())),
                Is.EqualTo(new[] { ("contigA", 12_345L, 24), ("contigB", 3_000L, 6) }));
            var ids = records.Select(r => r.GetProperty("id").GetInt64()).ToList();
            Assert.That(response.Headers.Location!.ToString(), Is.EqualTo($"/api/sequences/{ids[0]}"));

            var detail = (await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{ids[0]}"))!;
            Assert.That(detail.GetProperty("owned").GetBoolean(), Is.True);
            Assert.That(detail.GetProperty("dataClass").GetString(), Is.EqualTo("restricted"));
            Assert.That(detail.GetProperty("sequence").GetProperty("provider").GetString(), Is.EqualTo("user upload"));
            Assert.That(detail.GetProperty("sequence").GetProperty("name").GetString(), Is.EqualTo("first contig"));
            Assert.That(detail.GetProperty("sequence").GetProperty("ncbiTaxid").GetInt32(), Is.EqualTo(562));
            Assert.That(detail.GetProperty("sequenceSha256").GetString(), Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(Encoding.ASCII.GetBytes(dnaA)))));
            Assert.That(detail.GetProperty("sequenceSha256").GetString(), Is.EqualTo(Convert.ToHexStringLower(expected[0].SequenceSha256!)));
            var provenance = detail.GetProperty("provenance");
            Assert.That(provenance.TryGetProperty("url", out _), Is.False, "an upload has no URL to fetch again");
            Assert.That(provenance.GetProperty("rawSha256").GetString(), Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(Encoding.ASCII.GetBytes(text)))));
            Assert.That(provenance.GetProperty("rawByteSize").GetInt64(), Is.EqualTo(text.Length));
            Assert.That(detail.GetProperty("levels").EnumerateArray().Select(l => l.GetProperty("windowSize").GetInt32()), Is.EqualTo(new[] { 500, 2_000, 8_000 }));
            Assert.That(await (await client.GetAsync($"/api/sequences/{ids[0]}/profile")).Content.ReadAsByteArrayAsync(),
                Is.EqualTo(expected[0].Levels[0].Level.ToBinary()), "the same profile the Ingest's path-based reader gives");
            var second = (await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{ids[1]}"))!;
            Assert.That(second.GetProperty("sequence").GetProperty("name").GetString(), Is.EqualTo("contigB"));
            Assert.That(second.GetProperty("sequenceSha256").GetString(), Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(Encoding.ASCII.GetBytes(dnaB)))));
            foreach (var id in ids)
                Assert.That((await bob.GetAsync($"/api/sequences/{id}")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
        }

        [Test]
        public async Task ProfileUploadDerivesItsLadderAndMustMatchItsDeclaredLength()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var dna = Generator.GetRandomDNA(20_250, seed: 33);
            var fine = new ZCurveWindowedProfile(dna);

            var response = await client.PostAsync("/api/sequences?kind=virus&name=my%20virus&lengthBp=20250", Zcwp(fine.ToBinary()));
            var tooLong = await client.PostAsync("/api/sequences?kind=virus&lengthBp=21000", Zcwp(fine.ToBinary()));
            var noLength = await client.PostAsync("/api/sequences?kind=virus", Zcwp(fine.ToBinary()));
            var coarse = await client.PostAsync("/api/sequences?kind=virus&lengthBp=20250", Zcwp(new ZCurveWindowedProfile(dna, 2_000).ToBinary()));
            var garbage = await client.PostAsync("/api/sequences?kind=virus&lengthBp=20250", Zcwp(Encoding.ASCII.GetBytes("ACGTACGT")));

            Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.Created));
            var id = (await response.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("records")[0].GetProperty("id").GetInt64();
            var detail = (await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}"))!;
            Assert.That(detail.GetProperty("sequence").GetProperty("lengthBp").GetInt64(), Is.EqualTo(20_250));
            Assert.That(detail.GetProperty("sequence").GetProperty("name").GetString(), Is.EqualTo("my virus"));
            Assert.That(detail.TryGetProperty("sequenceSha256", out _), Is.False, "a profile carries no sequence to hash");
            Assert.That(detail.GetProperty("levels").EnumerateArray().Select(l => (l.GetProperty("windowSize").GetInt32(), l.GetProperty("windowCount").GetInt32())),
                Is.EqualTo(new[] { (500, 40), (2_000, 10), (8_000, 2) }));
            foreach (var rejected in new[] { tooLong, noLength, coarse, garbage })
            {
                Assert.That(rejected.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
                Assert.That(rejected.Content.Headers.ContentType!.MediaType, Is.EqualTo("application/problem+json"));
            }
            Assert.That((await tooLong.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("detail").GetString(),
                Is.EqualTo("A length of 21000 bp gives 42 windows of 500 bp, but the profile has 40."));
        }

        [Test]
        public async Task RejectedUploadsStoreNothing()
        {
            var (alice, client) = await factory.NewClientAsync("alice");
            var good = ">a\nACGTACGT\n";

            var cases = new (string Query, HttpContent Body, HttpStatusCode Status)[]
            {
                ("", Fasta(good), HttpStatusCode.BadRequest),                                  // no kind
                ("?kind=fungus", Fasta(good), HttpStatusCode.BadRequest),
                ("?kind=virus&dataClass=secret", Fasta(good), HttpStatusCode.BadRequest),
                ("?kind=virus", new StringContent("{}", Encoding.UTF8, "application/json"), HttpStatusCode.UnsupportedMediaType),
                ("?kind=virus", Fasta(""), HttpStatusCode.BadRequest),                         // no record
                ("?kind=virus", Fasta("ACGT\n>a\nACGT\n"), HttpStatusCode.BadRequest),         // bases before the first header
                ("?kind=virus", Fasta(">a\nACGT\n>a\nACGT\n"), HttpStatusCode.BadRequest),     // a repeated id
                ("?kind=virus", Fasta(">a\n>b\nACGT\n"), HttpStatusCode.BadRequest),           // a record without bases
            };
            foreach (var (query, body, status) in cases)
            {
                var response = await client.PostAsync("/api/sequences" + query, body);
                Assert.That(response.StatusCode, Is.EqualTo(status), query);
                Assert.That(response.Content.Headers.ContentType!.MediaType, Is.EqualTo("application/problem+json"), query);
            }
            Assert.That(await UploadCountAsync(alice.Id), Is.Zero);
        }
    }
}
