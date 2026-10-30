using System.Text.Json;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    /// <summary>--per-species: one genome for each species not yet stored, the best-curated one where there's a choice.</summary>
    public class SpeciesSelectionTests
    {
        private static BvBrcGenome Genome(string id, string? species, string? reference = null)
        {
            var fields = new Dictionary<string, string?> { ["species"] = species, ["reference_genome"] = reference };
            return new BvBrcGenome(id, $"genome {id}", JsonDocument.Parse(JsonSerializer.Serialize(fields)).RootElement.Clone());
        }

        [Test]
        public void TakesOneGenomePerNewSpeciesAndPrefersACuratedOne()
        {
            var genomes = new[]
            {
                Genome("1.1", "Alpha one"),
                Genome("1.2", "Alpha one"),
                Genome("2.3", "Beta two"),
                Genome("2.2", "Beta two", "Representative"),
                Genome("3.1", "Gamma stored"),
                Genome("4.1", null),
            };

            var picked = SpeciesSelection.PickNew(genomes, new HashSet<string> { "Gamma stored" }, 10);

            Assert.That(picked.Select(g => g.GenomeId), Is.EquivalentTo(new[] { "1.1", "2.2" }));
        }

        [Test]
        public void ALimitIsARepeatableSampleNotTheStartOfTheAlphabet()
        {
            var species = Enumerable.Range(0, 1000).Select(i => $"Species {i:D4}").ToList();

            var first = SpeciesSelection.Missing(species, new HashSet<string>(), 100);

            Assert.That(first, Has.Count.EqualTo(100));
            Assert.That(first, Is.EqualTo(SpeciesSelection.Missing(Enumerable.Reverse(species), new HashSet<string>(), 100)), "the same whatever order they arrive in");
            Assert.That(first.Max(), Is.GreaterThan("Species 0500"));
        }

        [Test]
        public async Task StoredSpeciesAreThePublicOnes()
        {
            var user = (await Users.CreateAsync(ThrowawayDatabase.DataSource, $"species-{Guid.NewGuid():N}")).User;
            var open = await ThrowawayDatabase.InsertRecordAsync("species", $"public-{Guid.NewGuid():N}", 20_000);
            var upload = await ThrowawayDatabase.InsertUploadAsync(user.Id, "private", 20_000);
            foreach (var (id, name) in new[] { (open, "Publicus speciei"), (upload, "Privatus speciei") })
            {
                await using var command = ThrowawayDatabase.DataSource.CreateCommand("update sequence_record set tax_species = $1 where id = $2");
                command.Parameters.AddWithValue(name);
                command.Parameters.AddWithValue(id);
                await command.ExecuteNonQueryAsync();
            }

            var stored = await Catalogue.StoredSpeciesAsync(ThrowawayDatabase.DataSource);

            Assert.That(stored, Does.Contain("Publicus speciei").And.Not.Contain("Privatus speciei"));
        }
    }
}
