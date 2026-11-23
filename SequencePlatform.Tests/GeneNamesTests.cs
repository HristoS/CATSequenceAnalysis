using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>Each stored gene's name (migration 0025): stored once per text, read back in the gene profile's order, found by where a gene lies.</summary>
    public class GeneNamesTests
    {
        private static string Gene(int aa, int seed) => "ATG" + Generator.GetRandomDNA(3 * (aa - 1), seed).Replace("TAA", "TAC")
            .Replace("TAG", "TAC").Replace("TGA", "TGC") + "TAA";

        [Test]
        public void NamesAreTrimmedCutAndEmptyMeansNone()
        {
            Assert.That(GeneNames.Clean("  holin "), Is.EqualTo("holin"));
            Assert.That(GeneNames.Clean("   "), Is.Null);
            Assert.That(GeneNames.Clean(null), Is.Null);
            Assert.That(GeneNames.Clean(new string('a', 1_500))!.Length, Is.EqualTo(GeneNames.MaxLength));
        }

        [Test]
        public async Task NamesRoundTripInProfileOrderAndShareTheirText()
        {
            var genes = new[] { Gene(60, 81), Gene(120, 82), Gene(80, 83) };
            var dna = string.Concat(genes);
            int a = genes[0].Length, b = a + genes[1].Length;
            var profile = ZCurveGeneProfile.FromRecord(dna, new[] { (0, a, false, -1), (a, b, false, -1), (b, dna.Length, true, -1) });
            long first = await ThrowawayDatabase.InsertRecordAsync("gene-names", "names-1", dna.Length);
            long second = await ThrowawayDatabase.InsertRecordAsync("gene-names", "names-2", dna.Length);
            await using (var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync())
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                Assert.That(await GeneProfiles.WriteAsync(connection, transaction, first, profile), Is.True);
                await GeneNames.WriteAsync(connection, transaction, first, GeneNames.BvBrc, new[] { " portal protein", null, "holin" }, unmatched: 1);
                await GeneNames.WriteAsync(connection, transaction, second, GeneNames.BvBrc, new[] { "holin", "" });
                await transaction.CommitAsync();
            }

            Assert.That(await GeneNames.ReadAsync(ThrowawayDatabase.DataSource, first), Is.EqualTo(new[] { "portal protein", null, "holin" }));
            Assert.That(await GeneNames.ReadAsync(ThrowawayDatabase.DataSource, second), Is.EqualTo(new[] { "holin", null }));
            Assert.That(await GeneNames.ReadAsync(ThrowawayDatabase.DataSource, -1), Is.Null);
            await using (var count = ThrowawayDatabase.DataSource.CreateCommand(
                "select (select count(*) from gene_name_text where name = 'holin'), named, unmatched from gene_name where sequence_id = $1"))
            {
                count.Parameters.AddWithValue(first);
                await using var reader = await count.ExecuteReaderAsync();
                Assert.That(await reader.ReadAsync(), Is.True);
                Assert.That((reader.GetInt64(0), reader.GetInt32(1), reader.GetInt32(2)), Is.EqualTo((1L, 2, 1)), "one shared text; two named; one unmatched");
            }

            // A protein hit keeps where the gene lies, not its index.
            var at = await GeneNames.AtAsync(ThrowawayDatabase.DataSource, new[] { (first, 0, b, dna.Length), (first, 0, a, b), (first, 0, 1, 2) });
            Assert.That(at, Is.EqualTo(new Dictionary<(long, int, int, int), string> { [(first, 0, b, dna.Length)] = "holin" }));

            // Rewriting replaces the record's names.
            await using (var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync())
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                await GeneNames.WriteAsync(connection, transaction, first, GeneNames.Ncbi, new string?[] { null, null, "terminase" });
                await transaction.CommitAsync();
            }
            Assert.That(await GeneNames.ReadAsync(ThrowawayDatabase.DataSource, first), Is.EqualTo(new[] { null, null, "terminase" }));
        }

        // Not names, but the per-gene rule the gene list shows beside them (2026-11-23): one rule per kind of organism.
        [Test]
        public void MembraneRuleDependsOnTheKindAndThePeakIsTheMostWaterRepellingWindow()
        {
            Assert.That(GeneProfiles.MembraneRule("bacterium"), Is.EqualTo(("mean", GeneProfiles.MembraneLikeCut)));
            Assert.That(GeneProfiles.MembraneRule("archaeon"), Is.EqualTo(("mean", GeneProfiles.MembraneLikeCut)));
            Assert.That(GeneProfiles.MembraneRule("virus"), Is.EqualTo(("peak", GeneProfiles.MembraneLikePeakCut)));
            Assert.That(GeneProfiles.MembraneRule("eukaryote"), Is.EqualTo(("model", 0.0)));
            Assert.That(GeneProfiles.MembraneRule("plasmid"), Is.Null, "not tested");

            var genes = new[] { Gene(60, 86), Gene(700, 87) };
            var dna = string.Concat(genes);
            var profile = ZCurveGeneProfile.FromRecord(dna, new[] { (0, genes[0].Length, false, -1), (genes[0].Length, dna.Length, false, -1) });
            for (int g = 0; g < 2; g++)
            {
                var w = profile.WindowsOf(g);
                var perWindow = Enumerable.Range(0, w.Length / 9).Select(k => -(w[9 * k] + w[9 * k + 1])).ToArray();
                Assert.That(GeneProfiles.HydropathyPeak(profile, g), Is.EqualTo(perWindow.Max()).Within(1e-12));
                Assert.That(GeneProfiles.MembraneScore(profile, g, "mean"), Is.EqualTo(perWindow.Average()).Within(1e-12));
                Assert.That(double.IsFinite(GeneProfiles.MembraneScore(profile, g, "model")), Is.True);
            }
            Assert.That(() => GeneProfiles.MembraneScore(profile, 0, "other"), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        [Test]
        public void MatchingByPlaceCountsGenesTheProviderLacks()
        {
            var genes = new[] { Gene(60, 84), Gene(70, 85) };
            var dna = string.Concat(genes);
            var profile = ZCurveGeneProfile.FromRecord(dna, new[] { (0, genes[0].Length, false, -1), (genes[0].Length, dna.Length, false, -1) });
            var (names, unmatched) = GeneNames.Match(profile, new Dictionary<(int, int, int, bool), string?> { [(0, 0, genes[0].Length, false)] = "kinase" });
            Assert.That(names, Is.EqualTo(new[] { "kinase", null }));
            Assert.That(unmatched, Is.EqualTo(1));
        }
    }
}
