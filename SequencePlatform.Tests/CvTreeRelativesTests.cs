using System.Net.Http.Json;
using System.Text.Json;
using Microsoft.AspNetCore.Mvc.Testing;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// Relatives and trees by amino-acid CVTree (2026-11-22): catalogue genomes whose groups differ only in their proteins' 4-mers
    /// (dinucleotide counts are random), so CVTree must find the groups and the dinucleotide signature can't.
    /// A genome without stored 4-mers, and a eukaryote, keep the dinucleotide signature.
    /// </summary>
    public class CvTreeRelativesTests
    {
        private const string provider = "cvtree";
        private const int groups = 3;

        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        // CVTree sees only what a word's 3-mers don't predict, so each group pairs the same prefixes and suffixes around a middle
        // differently: group g writes prefix i, middle, suffix (i + g) mod 3. Every group has the same 3-mers and 2-mers of these words,
        // and only their 4-mers differ.
        private static readonly string[] middles = { "CM", "HW", "QY", "FK", "NR" };
        private static readonly string[] prefixes = { "AED", "GIL", "PST", "VAE", "DGI" };
        private static readonly string[] suffixes = { "LPS", "TVA", "EDG", "ILP", "STV" };

        // 40 proteins of about 300 residues: the group's words in random order. (A random residue after each word buries the plant: CVTree
        // then measured about 0.28 between any two genomes, within a group or not.)
        private static int[] FourMers(int group, int seed)
        {
            var words = Enumerable.Range(0, middles.Length).SelectMany(m => Enumerable.Range(0, 3)
                .Select(i => $"{prefixes[m][i]}{middles[m]}{suffixes[m][(i + group) % 3]}")).ToArray();
            var rng = new Random(seed);
            var proteins = Enumerable.Range(0, 40).Select(_ =>
            {
                var protein = new System.Text.StringBuilder();
                while (protein.Length < 300) protein.Append(words[rng.Next(words.Length)]);
                return protein.ToString();
            });
            return AminoAcidKmers.Count(proteins, CvTreeSignature.K);
        }

        private static async Task<long> GenomeAsync(string name, string kind, int seed, int? group)
        {
            long id = await ThrowawayDatabase.InsertRecordAsync(provider, name, 1_000_000, kind);
            var rng = new Random(seed);
            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var command = new NpgsqlCommand("update sequence_record set dinucleotide_counts = $1 where id = $2", connection, transaction))
            {
                command.Parameters.AddWithValue(Enumerable.Range(0, 16).Select(_ => (long)rng.Next(1_000, 100_000)).ToArray());
                command.Parameters.AddWithValue(id);
                await command.ExecuteNonQueryAsync();
            }
            if (group is int g)
                await ProteinKmers.WriteAsync(connection, transaction, id, CvTreeSignature.K, ProteinKmers.Done, 40, FourMers(g, seed));
            await transaction.CommitAsync();
            return id;
        }

        // Three groups of four genomes of one kind; returns each genome's id and group.
        private static async Task<List<(long Id, int Group)>> PlantAsync(string prefix, string kind, int seed)
        {
            var genomes = new List<(long, int)>();
            for (int group = 0; group < groups; group++)
                for (int i = 0; i < 4; i++)
                    genomes.Add((await GenomeAsync($"{prefix}-{group}-{i}", kind, seed + 10 * group + i, group), group));
            return genomes;
        }

        // Asks until the API's index has built its vectors (in the background after the refresh a request triggers).
        private static async Task<JsonElement> WhenCvTreeAsync(HttpClient client, string path)
        {
            for (int attempt = 0; ; attempt++)
            {
                var body = await client.GetFromJsonAsync<JsonElement>(path);
                if (body.GetProperty("signature").GetString() == CvTreeSignature.Description || attempt == 100) return body;
                await Task.Delay(200);
            }
        }

        [Test]
        public void AVectorIsUnitLengthAndOneWithoutWordsIsNone()
        {
            var vector = CvTreeSignature.FromFourMers(FourMers(0, 1))!;
            Assert.That(vector, Has.Length.EqualTo(CvTreeSignature.Dimensions));
            // Within 16-bit rounding of unit length.
            Assert.That(CvTreeSignature.Distance(vector, vector), Is.EqualTo(0).Within(1e-3));
            Assert.That(CvTreeSignature.Distance(vector, vector.Select(x => (short)-x).ToArray()), Is.EqualTo(1).Within(1e-3));
            Assert.That(CvTreeSignature.FromFourMers(new int[160_000]), Is.Null);
        }

        [Test]
        public async Task RelativesWithAGeneLayerFollowTheProteinsAndTheRestKeepTheDinucleotideSignature()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var planted = await PlantAsync("relatives", "bacterium", 100);
            long withoutGenes = await GenomeAsync("relatives-no-genes", "bacterium", 900, null);
            long eukaryote = await GenomeAsync("relatives-eukaryote", "eukaryote", 901, 0);

            foreach (var (id, group) in planted)
            {
                var body = await WhenCvTreeAsync(client, $"/api/sequences/{id}/relatives?k=3");
                Assert.That(body.GetProperty("signature").GetString(), Is.EqualTo(CvTreeSignature.Description));
                var ids = body.GetProperty("items").EnumerateArray().Select(r => r.GetProperty("sequence").GetProperty("id").GetInt64()).ToList();
                Assert.That(ids, Is.EquivalentTo(planted.Where(g => g.Group == group && g.Id != id).Select(g => g.Id)), $"group {group} genome {id}");
                Assert.That(ids, Does.Not.Contain(withoutGenes), "only genomes with a vector are CVTree candidates");
            }

            foreach (var id in new[] { withoutGenes, eukaryote })
            {
                var body = await client.GetFromJsonAsync<JsonElement>($"/api/sequences/{id}/relatives?k=3");
                Assert.That(body.GetProperty("signature").GetString(), Is.EqualTo(SignatureIndex.Description), $"genome {id}");
            }
        }

        [Test]
        public async Task ATreeUsesCvTreeOnlyWhenEveryLeafHasAVector()
        {
            var (_, client) = await factory.NewClientAsync("alice");
            var planted = await PlantAsync("tree", "virus", 300);
            long withoutGenes = await GenomeAsync("tree-no-genes", "virus", 990, null);
            var ids = string.Join(",", planted.Select(g => g.Id));

            var body = await WhenCvTreeAsync(client, $"/api/trees?ids={ids}");
            Assert.That(body.GetProperty("signature").GetString(), Is.EqualTo(CvTreeSignature.Description));
            var newick = body.GetProperty("newick").GetString()!;
            foreach (var group in Enumerable.Range(0, groups))
            {
                // Each group's four leaves sit under one parenthesis (the root's other side holds the rest, which UPGMA also groups).
                var members = planted.Where(g => g.Group == group).Select(g => g.Id.ToString()).ToHashSet();
                Assert.That(Clades(newick).Any(c => c.SetEquals(members)), Is.True, $"group {group}: {newick}");
            }

            var mixed = await client.GetFromJsonAsync<JsonElement>($"/api/trees?ids={ids},{withoutGenes}");
            Assert.That(mixed.GetProperty("signature").GetString(), Is.EqualTo(SignatureIndex.Description));
        }

        // The leaf set below every parenthesis of a Newick string, as OrganismApiTests reads it.
        private static List<HashSet<string>> Clades(string newick)
        {
            var clades = new List<HashSet<string>>();
            var open = new Stack<HashSet<string>>();
            var label = new System.Text.StringBuilder();
            bool inLength = false;
            void Leaf()
            {
                if (label.Length > 0) open.Peek().Add(label.ToString());
                label.Clear();
                inLength = false;
            }
            foreach (char c in newick)
                switch (c)
                {
                    case '(': open.Push(new HashSet<string>()); break;
                    case ',': Leaf(); break;
                    case ')':
                        Leaf();
                        var clade = open.Pop();
                        clades.Add(clade);
                        if (open.Count > 0) open.Peek().UnionWith(clade);
                        inLength = true;   // a label after ')' is the clade's length, never a leaf
                        break;
                    case ':': inLength = true; break;
                    case ';': break;
                    default: if (!inLength) label.Append(c); break;
                }
            return clades;
        }
    }
}
