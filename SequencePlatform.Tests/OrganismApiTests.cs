using System.Diagnostics;
using System.Net;
using System.Net.Http.Json;
using System.Text;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// Nearest relatives, signature views and trees (P5 S5.5) over uploads of synthetic genomes whose neighbour structure is planted:
    /// three groups of four, drawn from first-order Markov chains that avoid CpG, avoid TpA, or favour A and T runs. Each genome's
    /// relatives and each tree must follow the groups, and nobody's upload may show up in anyone else's results.
    /// </summary>
    public class OrganismApiTests
    {
        private const int genomeLength = 60_000;
        private static readonly string[] groups = { "cpg", "tpa", "runs" };

        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        // P(next | current) over A, C, G, T: uniform, except the planted preference.
        private static string Markov(string group, int seed)
        {
            var transitions = new double[4, 4];
            for (int a = 0; a < 4; a++)
                for (int b = 0; b < 4; b++)
                    transitions[a, b] = (group, a, b) switch
                    {
                        ("cpg", 1, 2) => 0.05,      // C then G
                        ("tpa", 3, 0) => 0.05,      // T then A
                        ("runs", 0, 0) or ("runs", 3, 3) => 0.55,
                        _ => 1,
                    };
            var rng = new Random(seed);
            var chars = new char[genomeLength];
            int current = rng.Next(4);
            for (int i = 0; i < chars.Length; i++)
            {
                double total = 0;
                for (int b = 0; b < 4; b++) total += transitions[current, b];
                double r = rng.NextDouble() * total;
                int next = 0;
                while (next < 3 && (r -= transitions[current, next]) > 0) next++;
                chars[i] = "ACGT"[current = next];
            }
            return new string(chars);
        }

        private static StringContent Fasta(IEnumerable<(string Id, string Dna)> records) =>
            new(string.Concat(records.Select(r => $">{r.Id}\n{r.Dna}\n")), Encoding.ASCII, "text/x-fasta");

        // Uploads four genomes per group as one FASTA; returns each record's id and group.
        private static async Task<List<(long Id, string Group)>> UploadGroupsAsync(HttpClient client, string kind, int seed)
        {
            var genomes = groups.SelectMany((g, gi) => Enumerable.Range(0, 4).Select(i => (Id: $"{g}{i}", Group: g, Dna: Markov(g, seed + 10 * gi + i)))).ToList();
            var response = await client.PostAsync($"/api/sequences?kind={kind}", Fasta(genomes.Select(g => (g.Id, g.Dna))));
            Assert.That(response.StatusCode, Is.EqualTo(HttpStatusCode.Created));
            var records = (await response.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("records").EnumerateArray().ToList();
            return records.Select((r, i) => (r.GetProperty("id").GetInt64(), genomes[i].Group)).ToList();
        }

        private static async Task<List<long>> RelativeIdsAsync(HttpClient client, long id, int k)
        {
            var body = await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}/relatives?k={k}");
            return body.GetProperty("items").EnumerateArray().Select(r => r.GetProperty("sequence").GetProperty("id").GetInt64()).ToList();
        }

        [Test]
        public async Task RelativesAreTheGenomesWithTheSameNeighbourStructureAndStayPrivate()
        {
            var (_, alice) = await factory.NewClientAsync("alice");
            var (_, bob) = await factory.NewClientAsync("bob");
            var mine = await UploadGroupsAsync(alice, "bacterium", 100);
            var bobs = (await UploadGroupsAsync(bob, "bacterium", 500)).First(g => g.Group == "cpg").Id;

            foreach (var (id, group) in mine)
            {
                var body = await alice.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}/relatives?k=3");
                Assert.That(body.GetProperty("group").GetString(), Is.EqualTo("cellular"));
                var items = body.GetProperty("items").EnumerateArray().ToList();
                var ids = items.Select(r => r.GetProperty("sequence").GetProperty("id").GetInt64()).ToList();
                Assert.That(ids, Is.EquivalentTo(mine.Where(g => g.Group == group && g.Id != id).Select(g => g.Id)), $"{group} genome {id}");
                var distances = items.Select(r => r.GetProperty("distance").GetDouble()).ToList();
                Assert.That(distances, Is.Ordered);
            }

            // Bob's CpG-avoiding genome finds his own relatives and the public catalogue, never Alice's uploads; and the other way round.
            var bobsRelatives = await RelativeIdsAsync(bob, bobs, 20);
            Assert.That(bobsRelatives, Is.Not.Empty);
            Assert.That(bobsRelatives.Intersect(mine.Select(g => g.Id)), Is.Empty);
            Assert.That((await RelativeIdsAsync(alice, mine[0].Id, 20)), Does.Not.Contain(bobs));
            Assert.That((await alice.GetAsync($"/api/sequences/{bobs}/relatives")).StatusCode, Is.EqualTo(HttpStatusCode.NotFound));

            // Viruses are their own group: a viral upload's relatives are viruses only.
            var viral = await UploadGroupsAsync(alice, "virus", 900);
            var viralBody = await alice.GetFromJsonAsync<JsonElement>($"/api/sequences/{viral[0].Id}/relatives?k=20");
            Assert.That(viralBody.GetProperty("group").GetString(), Is.EqualTo("viral"));
            Assert.That(viralBody.GetProperty("items").EnumerateArray().Select(r => r.GetProperty("sequence").GetProperty("kind").GetString()).Distinct(),
                Is.EqualTo(new[] { "virus" }));

            // Warm, one relatives query reads the in-memory index: well under the 100 ms stop criterion even with HTTP around it.
            var watch = Stopwatch.StartNew();
            await RelativeIdsAsync(alice, mine[0].Id, 10);
            Assert.That(watch.ElapsedMilliseconds, Is.LessThan(100));
        }

        [Test]
        public async Task SignatureViewsAreOneVectorInTwoCoordinateSystems()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var mine = await UploadGroupsAsync(client, "bacterium", 200);
            var (cpg, tpa) = (mine.First(g => g.Group == "cpg").Id, mine.First(g => g.Group == "tpa").Id);

            async Task<List<JsonElement>> View(long id, string view) =>
                (await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}/signature?view={view}")).GetProperty("values").EnumerateArray().ToList();
            var (karlinA, karlinB) = (await View(cpg, "karlin"), await View(tpa, "karlin"));
            var (zcurveA, zcurveB) = (await View(cpg, "zcurve"), await View(tpa, "zcurve"));

            Assert.That(karlinA.Select(v => v.GetProperty("name").GetString()), Is.EqualTo(
                new[] { "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT" }));
            var cg = karlinA.Single(v => v.GetProperty("name").GetString() == "CG");
            Assert.That(cg.GetProperty("ratio").GetDouble(), Is.LessThan(0.2), "the planted CpG avoidance");
            Assert.That(cg.GetProperty("ratio").GetDouble(), Is.EqualTo(Math.Exp(cg.GetProperty("value").GetDouble())).Within(1e-9));
            Assert.That(zcurveA.Select(v => v.GetProperty("name").GetString()), Does.Contain("x(n)·x(n+1)").And.Contain("1(n)·1(n+1)"));

            double Distance(List<JsonElement> a, List<JsonElement> b) =>
                Math.Sqrt(a.Zip(b, (u, v) => Math.Pow(u.GetProperty("value").GetDouble() - v.GetProperty("value").GetDouble(), 2)).Sum());
            Assert.That(Distance(karlinA, karlinB), Is.EqualTo(4 * Distance(zcurveA, zcurveB)).Within(1e-9), "a rotation: distances scale by 1/4");
            Assert.That((await client.GetAsync($"/api/sequences/{cpg}/signature?view=other")).StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
        }

        [Test]
        public async Task TreesKeepEachPlantedGroupAsAClade()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var mine = await UploadGroupsAsync(client, "bacterium", 300);

            foreach (var method in new[] { "upgma", "nj" })
            {
                var body = await client.GetFromJsonAsync<JsonElement>($"/api/trees?ids={string.Join(",", mine.Select(g => g.Id))}&method={method}");
                Assert.That(body.GetProperty("leaves").GetArrayLength(), Is.EqualTo(12));
                var clades = Clades(body.GetProperty("newick").GetString()!);
                foreach (var group in groups)
                {
                    var members = mine.Where(g => g.Group == group).Select(g => g.Id.ToString()).ToHashSet();
                    var others = mine.Select(g => g.Id.ToString()).Except(members).ToHashSet();
                    Assert.That(clades.Any(c => c.SetEquals(members) || c.SetEquals(others)), Is.True, $"{method}: {group}");
                }
            }
        }

        // The leaf set below every parenthesis of a Newick string.
        private static List<HashSet<string>> Clades(string newick)
        {
            var clades = new List<HashSet<string>>();
            var open = new Stack<HashSet<string>>();
            for (int i = 0; i < newick.Length;)
            {
                char c = newick[i];
                if (c == '(') { open.Push(new HashSet<string>()); i++; }
                else if (c == ')')
                {
                    var clade = open.Pop();
                    clades.Add(clade);
                    if (open.Count > 0) open.Peek().UnionWith(clade);
                    i++;
                    while (i < newick.Length && ",);".IndexOf(newick[i]) < 0) i++;   // skip the length
                }
                else if (c == ',' || c == ';') i++;
                else
                {
                    int start = i;
                    while (i < newick.Length && ":,);".IndexOf(newick[i]) < 0) i++;
                    open.Peek().Add(newick[start..i]);
                    while (i < newick.Length && ",);".IndexOf(newick[i]) < 0) i++;
                }
            }
            return clades;
        }

        [Test]
        public async Task RequestsThatCantBeAnsweredSayWhy()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var mine = await UploadGroupsAsync(client, "bacterium", 400);
            var plasmid = (await (await client.PostAsync("/api/sequences?kind=plasmid", Fasta([("p", Markov("cpg", 1))])))
                .Content.ReadFromJsonAsync<JsonElement>()).GetProperty("records")[0].GetProperty("id").GetInt64();
            var virus = (await (await client.PostAsync("/api/sequences?kind=virus", Fasta([("v", Markov("cpg", 2))])))
                .Content.ReadFromJsonAsync<JsonElement>()).GetProperty("records")[0].GetProperty("id").GetInt64();
            var profileBody = new ByteArrayContent(new SequenceAnalyses.ZCurveWindowedProfile(Markov("tpa", 3)).ToBinary());
            profileBody.Headers.ContentType = new System.Net.Http.Headers.MediaTypeHeaderValue("application/x-zcwp");
            var profileOnly = (await (await client.PostAsync($"/api/sequences?kind=virus&lengthBp={genomeLength}", profileBody))
                .Content.ReadFromJsonAsync<JsonElement>()).GetProperty("records")[0].GetProperty("id").GetInt64();
            var ids = string.Join(",", mine.Select(g => g.Id));

            var cases = new (string Path, HttpStatusCode Status)[]
            {
                ($"/api/sequences/{plasmid}/relatives", HttpStatusCode.BadRequest),
                ($"/api/sequences/{profileOnly}/relatives", HttpStatusCode.BadRequest),
                ($"/api/sequences/{profileOnly}/signature", HttpStatusCode.BadRequest),
                ($"/api/sequences/{mine[0].Id}/relatives?k=0", HttpStatusCode.BadRequest),
                ($"/api/sequences/{mine[0].Id}/relatives?k=101", HttpStatusCode.BadRequest),
                ("/api/sequences/999999999/relatives", HttpStatusCode.NotFound),
                ($"/api/trees?ids={mine[0].Id},{mine[1].Id}", HttpStatusCode.BadRequest),
                ($"/api/trees?ids={string.Join(",", Enumerable.Range(1, 501))}", HttpStatusCode.BadRequest),
                ($"/api/trees?ids={ids},{virus}", HttpStatusCode.BadRequest),
                ($"/api/trees?ids={ids},999999999", HttpStatusCode.BadRequest),
                ($"/api/trees?ids={ids}&rank=genus&taxon=Escherichia", HttpStatusCode.BadRequest),
                ("/api/trees", HttpStatusCode.BadRequest),
                ($"/api/trees?ids={ids}&method=parsimony", HttpStatusCode.BadRequest),
                ("/api/trees?rank=tribe&taxon=x", HttpStatusCode.BadRequest),
            };
            foreach (var (path, status) in cases)
            {
                var response = await client.GetAsync(path);
                Assert.That(response.StatusCode, Is.EqualTo(status), path);
                Assert.That(response.Content.Headers.ContentType!.MediaType, Is.EqualTo("application/problem+json"), path);
            }
        }
    }
}
