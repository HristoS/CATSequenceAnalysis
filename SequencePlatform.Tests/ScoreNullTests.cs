using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class ScoreNullTests
    {
        private const string method = "zcurve-windowed";

        // A record with its ladder levels under this test's method version and a provider lineage the builder reads.
        private static async Task<long> StoreAsync(string methodVersion, string accession, string dna, string phylum)
        {
            var sequenceId = await ThrowawayDatabase.InsertRecordAsync($"null-{methodVersion}", accession, dna.Length, kind: "bacterium");
            var ladder = new ZCurveProfileLadder(new ZCurveWindowedProfile(dna));

            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var update = new Npgsql.NpgsqlCommand(
                "update sequence_record set tax_superkingdom = 'Bacteria', tax_phylum = $1 where id = $2", connection, transaction))
            {
                update.Parameters.AddWithValue(phylum);
                update.Parameters.AddWithValue(sequenceId);
                await update.ExecuteNonQueryAsync();
            }
            foreach (var factor in ladder.Factors)
                if (ladder.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, ladder.Level(factor));
            await transaction.CommitAsync();
            return sequenceId;
        }

        [Test]
        public async Task SaveReplacesTheCellsOfOneSourceAndLoadReadsThemBack()
        {
            var first = ScoreNullCell.FromScores(500, 6, 3, new[] { 0.5, 0.6, 0.7 });
            var second = ScoreNullCell.FromScores(2000, 7, 4, new[] { 0.9, 0.8 });

            await ScoreNullStore.SaveAsync(ThrowawayDatabase.DataSource, method, "store", ScoreNullStore.Genome, new[] { first });
            await ScoreNullStore.SaveAsync(ThrowawayDatabase.DataSource, method, "store", ScoreNullStore.Genome, new[] { second });
            await ScoreNullStore.SaveAsync(ThrowawayDatabase.DataSource, method, "store", ScoreNullStore.Fragment, new[] { first });

            var genome = await ScoreNullStore.LoadAsync(ThrowawayDatabase.DataSource, method, "store", ScoreNullStore.Genome);
            Assert.That(genome!.Cells, Has.Count.EqualTo(1));
            var cell = genome.Cells.Single();
            Assert.That((cell.WindowSize, cell.WindowsBucket, cell.OffsetsBucket, cell.SortedScores),
                Is.EqualTo((2000, 7, 4, new[] { 0.8f, 0.9f })));
            Assert.That(await ScoreNullStore.LoadAsync(ThrowawayDatabase.DataSource, method, "store", ScoreNullStore.Permuted), Is.Null);
        }

        // Every record of a phylum is the same sequence, so any related pair scores 1. The genome and fragment
        // nulls may only pair different phyla, so neither can hold a score of 1; the permuted control pairs anything.
        [Test]
        public async Task NullsPairOnlyUnrelatedRecords()
        {
            const string version = "build";
            var dnaA = Generator.GetRandomDNA(80_000, seed: 1);
            var dnaB = Generator.GetRandomDNA(80_000, seed: 2);
            for (int i = 0; i < 6; i++)
            {
                await StoreAsync(version, $"a{i}", dnaA, "PhylumA");
                await StoreAsync(version, $"b{i}", dnaB, "PhylumB");
            }

            var built = await ScoreNullBuilder.BuildAsync(ThrowawayDatabase.DataSource,
                new ScoreNullOptions(SamplesPerCell: 50, GenomeDraws: 2_000, FragmentDraws: 2_000, FragmentPool: 12, MinCellSamples: 5),
                _ => { }, method, version);

            Assert.That(built[ScoreNullStore.Genome], Is.Not.Empty);
            Assert.That(built[ScoreNullStore.Fragment], Is.Not.Empty);
            Assert.That(built[ScoreNullStore.Genome].SelectMany(c => c.SortedScores).Max(), Is.LessThan(0.999f));
            Assert.That(built[ScoreNullStore.Fragment].SelectMany(c => c.SortedScores).Max(), Is.LessThan(0.999f));
            Assert.That(built.Values.SelectMany(cells => cells).All(c => c.SortedScores.Length is >= 5 and <= 50), Is.True);
        }

        [Test]
        public async Task SearchReportsPValuesButRanksByScore()
        {
            const string version = "ranked";
            var dnaA = Generator.GetRandomDNA(80_000, seed: 3);
            for (int i = 0; i < 4; i++)
            {
                await StoreAsync(version, $"a{i}", Generator.GetRandomDNA(60_000 + 10_000 * i, seed: 10 + i), "PhylumA");
                await StoreAsync(version, $"b{i}", Generator.GetRandomDNA(60_000 + 10_000 * i, seed: 20 + i), "PhylumB");
            }
            var built = await ScoreNullBuilder.BuildAsync(ThrowawayDatabase.DataSource,
                new ScoreNullOptions(SamplesPerCell: 50, GenomeDraws: 2_000, FragmentDraws: 100, FragmentPool: 8, MinCellSamples: 5),
                _ => { }, method, version);
            await ScoreNullStore.SaveAsync(ThrowawayDatabase.DataSource, method, version, ScoreNullStore.Genome, built[ScoreNullStore.Genome].ToList());

            var search = new SimilaritySearch(ThrowawayDatabase.DataSource, method: method, methodVersion: version);
            var result = await search.SearchAsync(QueryProfile.FromSequence(dnaA), new CatalogueFilter(Provider: $"null-{version}"));

            Assert.That(result.Ranked, Is.Not.Empty);
            Assert.That(result.Ranked.All(h => h.P is not null && h.Log10E > h.P.Value.Log10), Is.True);
            Assert.That(result.Ranked.Select(h => h.Score), Is.Ordered.Descending);
        }
    }
}
