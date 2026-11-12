using System.Net;
using System.Net.Http.Json;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;

namespace SequencePlatform.Tests
{
    /// <summary>The map (GET /api/landscape): public genomes only, one per species, in parallel arrays, rebuilt when the catalogue changes.</summary>
    public class LandscapeTests
    {
        private const int taxid = 99999101;
        private const string species = "Testvirus primus";

        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        private static async Task ExecuteAsync(string sql, params object[] values)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(sql);
            foreach (var value in values) command.Parameters.AddWithValue(value);
            await command.ExecuteNonQueryAsync();
        }

        private static async Task WithCountsAsync(long id, int seed, int? ncbiTaxid = null)
        {
            var rng = new Random(seed);
            var counts = Enumerable.Range(0, 16).Select(_ => (long)rng.Next(1_000, 100_000)).ToArray();
            await ExecuteAsync("update sequence_record set dinucleotide_counts = $1 where id = $2", counts, id);
            if (ncbiTaxid is int value) await ExecuteAsync("update sequence_record set ncbi_taxid = $1 where id = $2", value, id);
        }

        private static long[] Ids(JsonElement map) => map.GetProperty("ids").EnumerateArray().Select(i => i.GetInt64()).ToArray();

        [Test]
        public async Task ShowsOnePublicGenomePerSpeciesAndNoUploads()
        {
            var (user, client) = await factory.NewClientAsync("alice");
            await ExecuteAsync("insert into taxon (taxid, status, current_taxid, name, tax_phylum, tax_family, tax_genus, tax_species) " +
                "values ($1, 'current', $1, $2, 'Testphylum', 'Testviridae', 'Testvirus', $2) on conflict do nothing", taxid, species);

            var first = await ThrowawayDatabase.InsertRecordAsync("landscape", "same-species-1", 20_000);
            var second = await ThrowawayDatabase.InsertRecordAsync("landscape", "same-species-2", 20_000);
            var unnamed = await ThrowawayDatabase.InsertRecordAsync("landscape", "no-species", 20_000);
            var tooShort = await ThrowawayDatabase.InsertRecordAsync("landscape", "too-short", 900);
            var upload = await ThrowawayDatabase.InsertUploadAsync(user.Id, "private", 20_000);
            await WithCountsAsync(first, 1, taxid);
            await WithCountsAsync(second, 2, taxid);
            await WithCountsAsync(unnamed, 3);
            await WithCountsAsync(tooShort, 4);
            await WithCountsAsync(upload, 5);

            var map = await client.GetFromJsonAsync<JsonElement>("/api/landscape?group=viral");
            var ids = Ids(map);

            Assert.That(ids.Count(id => id == first || id == second), Is.EqualTo(1), "one genome per species");
            Assert.That(ids, Does.Contain(unnamed), "a genome without a species stands for itself");
            Assert.That(ids, Does.Not.Contain(tooShort).And.Not.Contain(upload));

            foreach (var parallel in new[] { "names", "x", "y", "z", "gc", "lengthBp", "phylum", "family", "genus" })
                Assert.That(map.GetProperty(parallel).GetArrayLength(), Is.EqualTo(ids.Length), parallel);

            int at = Array.FindIndex(ids, id => id == first || id == second);
            Assert.That(map.GetProperty("genera")[map.GetProperty("genus")[at].GetInt32()].GetString(), Is.EqualTo("Testvirus"));
            Assert.That(map.GetProperty("phyla")[map.GetProperty("phylum")[at].GetInt32()].GetString(), Is.EqualTo("Testphylum"));
            Assert.That(map.GetProperty("genus")[Array.IndexOf(ids, unnamed)].GetInt32(), Is.EqualTo(-1), "no lineage: index -1");

            var axes = map.GetProperty("axes");
            Assert.That(axes.GetArrayLength(), Is.EqualTo(3), "three components, drawn two at a time");
            Assert.That(axes.EnumerateArray().Select(a => a.GetProperty("explainedVariance").GetDouble()), Is.Ordered.Descending);
            Assert.That(axes[0].GetProperty("topLoadings")[0].GetProperty("dinucleotide").GetString(), Has.Length.EqualTo(2));
        }

        [Test]
        public async Task CellularGenomesNeedAGtdbAssignmentAndTheirLargestRecordStands()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            await ExecuteAsync("insert into gtdb_taxon (assembly_key, gtdb_accession, release, tax_domain, tax_phylum, tax_genus, tax_species) " +
                "values ('999999102.1', 'RS_GCF_999999102.1', 'v232', 'Bacteria', 'Testota', 'Testus', 'Testus secundus') on conflict do nothing");
            var assigned = await ThrowawayDatabase.InsertRecordAsync("landscape", "gtdb", 2_000_000, kind: "bacterium");
            var unassigned = await ThrowawayDatabase.InsertRecordAsync("landscape", "no-gtdb", 2_000_000, kind: "bacterium");
            await ExecuteAsync("update sequence_record set annotations = '{\"bvbrc\": {\"assembly_accession\": \"GCF_999999102.1\"}}' where id = $1", assigned);
            await WithCountsAsync(assigned, 6);
            await WithCountsAsync(unassigned, 7);

            var ids = Ids(await client.GetFromJsonAsync<JsonElement>("/api/landscape"));

            Assert.That(ids, Does.Contain(assigned).And.Not.Contain(unassigned));
        }

        [Test]
        public async Task TheMapOfEverythingHoldsEveryPublicGenomeOnOneSetOfAxes()
        {
            var (user, client) = await factory.NewClientAsync("alice");
            var firstOfSpecies = await ThrowawayDatabase.InsertRecordAsync("landscape", "all-species-1", 20_000);
            var secondOfSpecies = await ThrowawayDatabase.InsertRecordAsync("landscape", "all-species-2", 20_000);
            var noAssignment = await ThrowawayDatabase.InsertRecordAsync("landscape", "all-bacterium", 2_000_000, kind: "bacterium");
            var archaeon = await ThrowawayDatabase.InsertRecordAsync("landscape", "all-archaeon", 2_000_000, kind: "archaeon");
            var plasmid = await ThrowawayDatabase.InsertRecordAsync("landscape", "all-plasmid", 50_000, kind: "plasmid");
            var upload = await ThrowawayDatabase.InsertUploadAsync(user.Id, "all-private", 20_000);
            await WithCountsAsync(firstOfSpecies, 11, taxid);
            await WithCountsAsync(secondOfSpecies, 12, taxid);
            await WithCountsAsync(noAssignment, 13);
            await WithCountsAsync(archaeon, 14);
            await WithCountsAsync(plasmid, 15);
            await WithCountsAsync(upload, 16);

            var map = await client.GetFromJsonAsync<JsonElement>("/api/landscape?group=all");
            var ids = Ids(map);

            Assert.That(ids, Does.Contain(firstOfSpecies).And.Contain(secondOfSpecies), "every genome, not one per species");
            Assert.That(ids, Does.Contain(noAssignment), "no GTDB assignment needed");
            Assert.That(ids, Does.Not.Contain(plasmid).And.Not.Contain(upload));

            string Domain(long id) => map.GetProperty("domains")[map.GetProperty("domain")[Array.IndexOf(ids, id)].GetInt32()].GetString()!;
            Assert.That(Domain(firstOfSpecies), Is.EqualTo("Viruses"));
            Assert.That(Domain(noAssignment), Is.EqualTo("Bacteria"));
            Assert.That(Domain(archaeon), Is.EqualTo("Archaea"));
            Assert.That(map.GetProperty("domain").GetArrayLength(), Is.EqualTo(ids.Length));
        }

        [Test]
        public async Task TheMapOfEverythingHoldsExactlyWhatTheSignatureIndexHolds()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            await WithCountsAsync(await ThrowawayDatabase.InsertRecordAsync("landscape", "parity-virus", 20_000), 21);
            await WithCountsAsync(await ThrowawayDatabase.InsertRecordAsync("landscape", "parity-yeast", 12_000_000, kind: "eukaryote"), 22);

            var index = new SequencePlatform.Core.SignatureIndex(ThrowawayDatabase.DataSource);
            await index.RefreshAsync();
            var map = await client.GetFromJsonAsync<JsonElement>("/api/landscape?group=all");

            Assert.That(Ids(map), Is.EquivalentTo(index.Entries.Values.Where(e => e.OwnerId is null).Select(e => e.Id)));
            Assert.That(map.GetProperty("domains").EnumerateArray().Select(d => d.GetString()), Does.Contain("Eukaryota"));
        }

        [Test]
        public async Task IsCompressedForAClientThatAsks()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            for (int i = 0; i < 40; i++)
                await WithCountsAsync(await ThrowawayDatabase.InsertRecordAsync("landscape", $"compress-{i}", 20_000), 100 + i);

            client.DefaultRequestHeaders.AcceptEncoding.ParseAdd("gzip");
            var response = await client.GetAsync("/api/landscape?group=all");

            Assert.That(response.Content.Headers.ContentEncoding, Does.Contain("gzip"));
        }

        [Test]
        public async Task IsStableUntilTheCatalogueChanges()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            await WithCountsAsync(await ThrowawayDatabase.InsertRecordAsync("landscape", "stable-1", 20_000), 8);
            await WithCountsAsync(await ThrowawayDatabase.InsertRecordAsync("landscape", "stable-2", 20_000), 9);

            var before = await client.GetStringAsync("/api/landscape?group=viral");
            Assert.That(await client.GetStringAsync("/api/landscape?group=viral"), Is.EqualTo(before));

            var added = await ThrowawayDatabase.InsertRecordAsync("landscape", "stable-3", 20_000);
            await WithCountsAsync(added, 10);

            Assert.That(Ids(await client.GetFromJsonAsync<JsonElement>("/api/landscape?group=viral")), Does.Contain(added));
        }

        [Test]
        public async Task NeedsAKeyAndAKnownGroup()
        {
            var (_, client) = await factory.NewClientAsync("alice");

            Assert.That((await factory.CreateClient().GetAsync("/api/landscape")).StatusCode, Is.EqualTo(HttpStatusCode.Unauthorized));
            Assert.That((await client.GetAsync("/api/landscape?group=fungal")).StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
        }
    }
}
