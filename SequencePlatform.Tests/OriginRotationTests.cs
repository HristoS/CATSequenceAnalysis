using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class OriginRotationTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "origin";

        // Random DNA with a replication-like GC skew: G-rich from 0 to the middle, C-rich after it, so the origin is at 0.
        private static string Skewed(int length, int seed)
        {
            var chars = Generator.GetRandomDNA(length, seed).ToCharArray();
            var rng = new Random(seed + 1);
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < 0.1)
                    chars[i] = i < length / 2 ? (chars[i] == 'C' ? 'G' : chars[i]) : (chars[i] == 'G' ? 'C' : chars[i]);
            return new string(chars);
        }

        // Stores the ladder levels a real ingest would, with BV-BRC's contig topology in the annotations.
        private static async Task<long> StoreAsync(string provider, string accession, string dna, string kind, string topology)
        {
            var sequenceId = await ThrowawayDatabase.InsertRecordAsync(provider, accession, dna.Length, kind);
            var ladder = RecordStore.Ladder(new ZCurveWindowedProfile(dna));

            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await using (var update = new NpgsqlCommand(
                "update sequence_record set annotations = jsonb_build_object('contig', jsonb_build_object('topology', $1::text)) where id = $2",
                connection, transaction))
            {
                update.Parameters.AddWithValue(topology);
                update.Parameters.AddWithValue(sequenceId);
                await update.ExecuteNonQueryAsync();
            }
            foreach (var factor in ladder.Factors)
                if (ladder.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, ladder.Level(factor));
            await transaction.CommitAsync();
            return sequenceId;
        }

        private static async Task<int[]> RotatedWindowSizesAsync(long sequenceId)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(
                "select window_size from profile_level where sequence_id = $1 and method_version = $2 order by window_size");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue(OriginRotation.RotatedVersion(methodVersion));
            var sizes = new List<int>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) sizes.Add(reader.GetInt32(0));
            return sizes.ToArray();
        }

        private static Task<(int Rotated, int Examined)> BuildAsync() =>
            OriginRotation.BuildAsync(ThrowawayDatabase.DataSource, _ => { }, method, methodVersion);

        [Test]
        public async Task OnlyCircularChromosomesWithAClearSkewAreRotated()
        {
            const string provider = "origin-build";
            var dna = Skewed(600_000, seed: 1);
            var chromosome = await StoreAsync(provider, "chromosome", dna, "bacterium", "circular");
            var noSkew = await StoreAsync(provider, "no-skew", Generator.GetRandomDNA(600_000, seed: 2), "bacterium", "circular");
            var linear = await StoreAsync(provider, "linear", Skewed(600_000, seed: 3), "bacterium", "linear");
            var tooShort = await StoreAsync(provider, "short", Skewed(400_000, seed: 4), "bacterium", "circular");
            var plasmid = await StoreAsync(provider, "plasmid", Skewed(600_000, seed: 5), "plasmid", "circular");

            await BuildAsync();

            Assert.That(await RotatedWindowSizesAsync(chromosome), Is.EqualTo(new[] { 2_000, 8_000, 32_000 }));
            foreach (var id in new[] { noSkew, linear, tooShort, plasmid })
                Assert.That(await RotatedWindowSizesAsync(id), Is.Empty, $"record {id}");

            var fine = new ZCurveWindowedProfile(dna);
            var stored = await ProfileLevels.ReadAsync(ThrowawayDatabase.DataSource, chromosome, method, OriginRotation.RotatedVersion(methodVersion), 32_000);
            Assert.That(stored!.Windows, Is.EqualTo(RecordStore.Ladder(OriginRotation.ToOrigin(fine.ToHalfPrecision(), dna.Length)!).Level(64).Windows));
        }

        [Test]
        public async Task ACopyCutElsewhereScoresAsItselfOnceRotated()
        {
            const string provider = "origin-search";
            var dna = Skewed(600_000, seed: 6);
            var original = await StoreAsync(provider, "original", dna, "bacterium", "circular");
            var copyDna = dna[300_000..] + dna[..300_000];   // the same ring cut half-way round, on the 500 bp grid
            var copy = await StoreAsync(provider, "copy", copyDna, "bacterium", "circular");
            await StoreAsync(provider, "other", Skewed(600_000, seed: 7), "bacterium", "circular");
            var filter = new CatalogueFilter(Provider: provider);

            var before = await new SimilaritySearch(ThrowawayDatabase.DataSource, method: method, methodVersion: methodVersion)
                .SearchRecordAsync(original, filter, useNeighbours: false);
            await BuildAsync();
            var search = new SimilaritySearch(ThrowawayDatabase.DataSource, method: method, methodVersion: methodVersion);
            var after = await search.SearchRecordAsync(original, filter, useNeighbours: false);
            var streamed = await new SimilaritySearch(ThrowawayDatabase.DataSource, useCache: false, method: method, methodVersion: methodVersion)
                .SearchRecordAsync(original, filter, useNeighbours: false);
            var fromDna = await search.SearchAsync(QueryProfile.FromSequence(copyDna), filter);

            var copyBefore = before.Ranked.Single(h => h.SequenceId == copy);
            var copyAfter = after.Ranked.Single(h => h.SequenceId == copy);
            Assert.That(copyBefore.Score, Is.LessThan(0.99));
            Assert.That(copyBefore.Rotated, Is.False);
            Assert.That(copyAfter.Score, Is.GreaterThan(0.999));
            Assert.That(copyAfter.Rotated, Is.True);
            Assert.That(after.Ranked[0].SequenceId, Is.EqualTo(copy));
            Assert.That(streamed.Ranked, Is.EqualTo(after.Ranked), "the cache and the database give the same result");
            Assert.That(fromDna.Ranked.Single(h => h.SequenceId == original).Score, Is.GreaterThan(0.999), "a FASTA query is rotated too");
            Assert.That(after.Ranked.Select(h => h.Score), Is.Ordered.Descending);
        }

        [Test]
        public void ShortOrUnskewedProfilesAreNotRotated()
        {
            Assert.That(OriginRotation.ToOrigin(new ZCurveWindowedProfile(Skewed(400_000, seed: 8)), 400_000), Is.Null);
            Assert.That(OriginRotation.ToOrigin(new ZCurveWindowedProfile(Generator.GetRandomDNA(600_000, seed: 9)), 600_000), Is.Null);
            Assert.That(OriginRotation.ToOrigin(new ZCurveWindowedProfile(Skewed(600_000, seed: 10)), 600_000), Is.Not.Null);
        }
    }
}
