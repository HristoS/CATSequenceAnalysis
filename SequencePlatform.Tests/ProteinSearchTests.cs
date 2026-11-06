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
    /// Protein search (P7, migration 0015): which stored genes encode a protein, from the protein alone. Core's query checks and
    /// streamed search, the job's storage, and the API with the worker off (<see cref="RunAsync"/> does what it does).
    /// </summary>
    public class ProteinSearchTests
    {
        private static readonly string[] senseCodons = (from a in "TCAG" from b in "TCAG" from c in "TCAG" select $"{a}{b}{c}")
            .Where(c => c is not ("TAA" or "TAG" or "TGA")).ToArray();

        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        /// <summary>A start codon, codons - 1 random sense codons and a stop: a protein of <paramref name="codons"/> aa.</summary>
        internal static string RandomGene(int codons, int seed)
        {
            var random = new Random(seed);
            return "ATG" + string.Concat(Enumerable.Range(0, codons - 1).Select(_ => senseCodons[random.Next(senseCodons.Length)])) + "TAA";
        }

        /// <summary>A record holding these genes between random flanks, with its gene row as the S7.R pass stores it.</summary>
        internal static async Task<(long Id, List<(int Start, int End)> Genes)> StoreGenesAsync(string provider, string accession,
            IReadOnlyList<(string Gene, bool Reverse, int Family)> genes)
        {
            var flank = Generator.GetRandomDNA(300, seed: 91);
            var record = new StringBuilder(flank);
            var cds = new List<(int Start, int End, bool Reverse, int Family)>();
            foreach (var (gene, reverse, family) in genes)
            {
                int start = record.Length;
                record.Append(reverse ? QueryProfile.ReverseComplement(gene) : gene);
                cds.Add((start, record.Length, reverse, family));
                record.Append(flank);
            }
            var dna = record.ToString();
            var id = await ThrowawayDatabase.InsertRecordAsync(provider, accession, dna.Length, "bacterium");
            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            Assert.That(await GeneProfiles.WriteAsync(connection, transaction, id, ZCurveGeneProfile.FromRecord(dna, cds)), Is.True);
            await transaction.CommitAsync();
            return (id, cds.Select(c => (c.Start, c.End)).ToList());
        }

        private static async Task RunAsync(Guid id)
        {
            var run = (await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!;
            var query = (await SearchLog.ProteinQueryAsync(ThrowawayDatabase.DataSource, id))!;
            var result = await ProteinSearch.SearchAsync(ThrowawayDatabase.DataSource, query, run.Request.Filter, run.Request.Top);
            Assert.That(await SearchLog.SaveAsync(ThrowawayDatabase.DataSource, id, result), Is.True);
        }

        private static StringContent Body(string text, string contentType = "text/x-fasta") => new(text, Encoding.ASCII, contentType);

        [Test]
        public void AProteinQueryChecksItsLettersAndKeepsNoneOfThem()
        {
            var protein = ZCurveGeneProfile.Translate(RandomGene(600, seed: 1));
            var query = ProteinQuery.FromProtein(protein);
            var sloppy = ProteinQuery.FromProtein($"  {protein[..300].ToLowerInvariant()}\n{protein[300..]}*\n");

            Assert.That((query.LengthAa, query.WindowCount, query.GeneticCode), Is.EqualTo((600, 3, 11)));
            Assert.That(sloppy.Sha256, Is.EqualTo(query.Sha256), "case, white space and a final stop don't change the query");
            Assert.That(sloppy.Axes, Is.EqualTo(query.Axes));
            Assert.That(query.Sha256, Is.EqualTo(SHA256.HashData(Encoding.ASCII.GetBytes(protein))));
            Assert.That(ProteinQuery.AxesFromBinary(query.AxesToBinary()), Is.EqualTo(query.Axes));
            Assert.That((ProteinSearch.MinLengthAa, ProteinQuery.FromProtein(protein[..ProteinSearch.MinLengthAa]).WindowCount), Is.EqualTo((50, 5)), "50 aa: 5 windows of 20 (S7.S)");

            var rejected = new (string Protein, int Code, string Message)[]
            {
                (protein[..(ProteinSearch.MinLengthAa - 1)], 11, "can't be searched"),
                (protein[..200] + "*" + protein[200..], 11, "internal stop"),
                (protein[..200] + "1" + protein[200..], 11, "residue 201 is '1'"),
                (Generator.GetRandomDNA(900, seed: 2), 11, "looks like DNA"),
                (protein, 5, "genetic code"),
            };
            foreach (var (text, code, message) in rejected)
                Assert.That(() => ProteinQuery.FromProtein(text, code), Throws.ArgumentException.With.Message.Contains(message), message);
        }

        [Test]
        public async Task SearchFindsTheGeneOfAProteinAmongTheFilteredRecords()
        {
            const string provider = "protein-core";
            var genes = Enumerable.Range(0, 9).Select(i => RandomGene(600 + 50 * i, seed: 100 + i)).ToArray();
            var records = new List<(long Id, List<(int Start, int End)> Genes)>();
            for (int r = 0; r < 3; r++)
                records.Add(await StoreGenesAsync(provider, $"p{r}", Enumerable.Range(3 * r, 3).Select(g => (genes[g], g % 2 == 1, 40 + g)).ToList()));
            var filter = new CatalogueFilter(Provider: provider);
            var query = ProteinQuery.FromProtein(ZCurveGeneProfile.Translate(genes[5]));

            var result = await ProteinSearch.SearchAsync(ThrowawayDatabase.DataSource, query, filter);
            var top2 = await ProteinSearch.SearchAsync(ThrowawayDatabase.DataSource, query, filter, top: 2);
            var elsewhere = await ProteinSearch.SearchAsync(ThrowawayDatabase.DataSource, query, new CatalogueFilter(Provider: "nobody"));

            var own = records[1].Genes[2];
            Assert.That((result.Hits[0].SequenceId, result.Hits[0].GeneStart, result.Hits[0].GeneEnd, result.Hits[0].Reverse, result.Hits[0].Family),
                Is.EqualTo((records[1].Id, own.Start, own.End, true, 45)), "the gene the protein came from, on the reverse strand");
            Assert.That(result.Hits[0].Score, Is.GreaterThan(result.Hits[1].Score));
            Assert.That((result.Hits.Count, result.Candidates, result.WeakSignal), Is.EqualTo((9, 9L, false)));
            Assert.That(result.Hits.Select(h => h.Score), Is.Ordered.Descending);
            Assert.That(top2.Hits, Is.EqualTo(result.Hits.Take(2)), "the kept top is the same whichever worker scored what");
            Assert.That((elsewhere.Hits.Count, elsewhere.Candidates), Is.EqualTo((0, 0L)));
        }

        [Test]
        public async Task AProteinUnder500AaMeetsOnlyGenesWithinTenPercentOfItsLength()
        {
            const string provider = "protein-short";
            var gene = RandomGene(400, seed: 200);
            var (id, stored) = await StoreGenesAsync(provider, "short", new[]
            {
                (gene, false, 7), (RandomGene(420, seed: 203), false, 10), (RandomGene(300, seed: 201), false, 8), (RandomGene(700, seed: 202), false, 9),
            });
            var query = ProteinQuery.FromProtein(ZCurveGeneProfile.Translate(gene));

            var result = await ProteinSearch.SearchAsync(ThrowawayDatabase.DataSource, query, new CatalogueFilter(Provider: provider));

            Assert.That((query.WindowCount, result.WeakSignal), Is.EqualTo((5, false)));
            Assert.That(result.Candidates, Is.EqualTo(2), "the 420 aa gene is in the band; 300 aa is too short and 700 aa is long");
            Assert.That((result.Hits[0].SequenceId, result.Hits[0].GeneStart), Is.EqualTo((id, stored[0].Start)));
        }

        [Test]
        public async Task AProteinJobKeepsItsWindowsAndHashUntilItExpires()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var query = ProteinQuery.FromProtein(ZCurveGeneProfile.Translate(RandomGene(400, seed: 300)), ProteinQuery.MycoplasmaCode);
            var request = new SearchRequest(SearchKind.Protein, new CatalogueFilter(Provider: "nobody"), 5, GeneticCode: ProteinQuery.MycoplasmaCode);

            var id = await SearchLog.EnqueueAsync(dataSource, request, query);
            var run = (await SearchLog.GetAsync(dataSource, id))!;
            var kept = (await SearchLog.ProteinQueryAsync(dataSource, id))!;
            Assert.That((run.Request.Kind, run.Request.GeneticCode, run.QueryLengthAa, run.QueryLengthBp), Is.EqualTo((SearchKind.Protein, (int?)4, (int?)400, 1_200L)));
            Assert.That((kept.Axes, kept.LengthAa, kept.Sha256, kept.GeneticCode), Is.EqualTo((query.Axes, 400, query.Sha256, 4)));
            Assert.That(await SearchLog.QueryAsync(dataSource, id), Is.Null, "a protein job has no DNA query");
            Assert.That(() => SearchLog.EnqueueAsync(dataSource, request with { GeneticCode = 11 }, query), Throws.ArgumentException);
            Assert.That(() => SearchLog.EnqueueAsync(dataSource, request with { Kind = SearchKind.Similar }, query), Throws.ArgumentException);

            var hit = new ProteinHit(await ThrowawayDatabase.InsertRecordAsync("protein-job", "job", 2_000), 30, 1_233, true, ZCurveGeneProfile.NoFamily, 0.9);
            Assert.That(await SearchLog.SaveAsync(dataSource, id, new ProteinSearchResult(new[] { hit }, 12, true)), Is.True);
            var stored = (await SearchLog.HitsAsync(dataSource, id, HitList.Ranked)).Single();
            Assert.That((stored.SequenceId, stored.GeneStart, stored.GeneEnd, stored.Strand, stored.FamilyId, stored.Score),
                Is.EqualTo((hit.SequenceId, (int?)30, (int?)1_233, (Strand?)Strand.Reverse, (int?)null, (double?)0.9)));
            var done = (await SearchLog.GetAsync(dataSource, id))!;
            Assert.That((done.Candidates, done.WeakSignal), Is.EqualTo(((int?)12, (bool?)true)));

            Assert.That(await SearchLog.ExpireAsync(dataSource, id), Is.True);
            Assert.That(await SearchLog.ProteinQueryAsync(dataSource, id), Is.Null);
        }

        [Test]
        public async Task TheApiQueuesAProteinAndPagesItsGenesWithTheirFamilies()
        {
            const string provider = "protein-api";
            var (_, client) = await factory.NewClientAsync("alice");
            var genes = Enumerable.Range(0, 4).Select(i => RandomGene(620, seed: 400 + i)).ToArray();
            var (id, stored) = await StoreGenesAsync(provider, "api", genes.Select((g, i) => (g, i == 2, 9_000_000 + i)).ToList());
            await using (var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync())
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                await ProteinFamilies.AddAsync(connection, transaction, new[] { (9_000_002, "Test synthase") });
                await transaction.CommitAsync();
            }
            var protein = ZCurveGeneProfile.Translate(genes[2]);

            var withHeader = await client.PostAsync($"/api/searches/protein?provider={provider}&top=3", Body($">query some protein\n{protein[..300]}\n{protein[300..]}\n"));
            var bare = await client.PostAsync($"/api/searches/protein?provider={provider}", Body(protein, "text/plain"));
            Assert.That(withHeader.StatusCode, Is.EqualTo(HttpStatusCode.Accepted), await withHeader.Content.ReadAsStringAsync());
            Assert.That(bare.StatusCode, Is.EqualTo(HttpStatusCode.Accepted), await bare.Content.ReadAsStringAsync());
            var search = (await withHeader.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("id").GetGuid();
            await RunAsync(search);

            var done = (await client.GetFromJsonAsync<JsonElement>($"/api/searches/{search}"))!;
            Assert.That((done.GetProperty("kind").GetString(), done.GetProperty("status").GetString(), done.GetProperty("queryLengthAa").GetInt32(),
                done.GetProperty("parameters").GetProperty("geneticCode").GetInt32()), Is.EqualTo(("Protein", "Done", 620, 11)));
            Assert.That(done.GetProperty("summary").GetProperty("candidates").GetInt32(), Is.EqualTo(4));
            Assert.That(done.GetProperty("summary").TryGetProperty("weakSignal", out var weak) && weak.GetBoolean(), Is.False);
            var hits = done.GetProperty("hits").GetProperty("items").EnumerateArray().ToList();
            Assert.That(hits, Has.Count.EqualTo(3));
            var top = hits[0];
            Assert.That((top.GetProperty("sequenceId").GetInt64(), top.GetProperty("geneStartBp").GetInt32(), top.GetProperty("geneEndBp").GetInt32(),
                top.GetProperty("strand").GetString(), top.GetProperty("familyId").GetString(), top.GetProperty("product").GetString()),
                Is.EqualTo((id, stored[2].Start, stored[2].End, "Reverse", "PGF_09000002", "Test synthase")));
            Assert.That(hits[1].GetProperty("familyId").GetString(), Does.StartWith("PGF_0900000"));
            Assert.That(hits[1].TryGetProperty("product", out _), Is.False, "a family with no stored product has none");
        }

        [Test]
        public async Task InvalidProteinSearchesAreRejectedBeforeAnyJobIsQueued()
        {
            var (alice, client) = await factory.NewClientAsync("alice");
            var protein = ZCurveGeneProfile.Translate(RandomGene(600, seed: 500));

            var cases = new (string Route, HttpContent Body, HttpStatusCode Status)[]
            {
                ("protein", Body(protein[..40]), HttpStatusCode.BadRequest),
                ("protein", Body(">dna\n" + Generator.GetRandomDNA(3_000, seed: 501)), HttpStatusCode.BadRequest),
                ("protein", Body($">a\n{protein}\n>b\n{protein}\n"), HttpStatusCode.BadRequest),
                ("protein", Body($"{protein}\n>late header\n"), HttpStatusCode.BadRequest),
                ("protein", Body(protein[..300] + "*" + protein[300..]), HttpStatusCode.BadRequest),
                ("protein?code=5", Body(protein), HttpStatusCode.BadRequest),
                ("protein?top=0", Body(protein), HttpStatusCode.BadRequest),
                ("protein?rank=genus", Body(protein), HttpStatusCode.BadRequest),
                ("protein", Body(protein, "application/json"), HttpStatusCode.UnsupportedMediaType),
            };
            foreach (var (route, body, status) in cases)
            {
                var response = await client.PostAsync("/api/searches/" + route, body);
                Assert.That(response.StatusCode, Is.EqualTo(status), route);
                var problem = await response.Content.ReadAsStringAsync();
                Assert.That(problem, Does.Not.Contain(protein[..40]), "an error never repeats the protein");
            }
            await using var command = ThrowawayDatabase.DataSource.CreateCommand("select count(*) from search_run where owner_id = $1");
            command.Parameters.AddWithValue(alice.Id);
            Assert.That((long)(await command.ExecuteScalarAsync())!, Is.Zero);
        }
    }
}
