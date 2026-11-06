using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class PlatformCoreTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "test";

        [Test]
        public async Task MigrationsAreIdempotent()
        {
            var scripts = typeof(Database).Assembly.GetManifestResourceNames().Count(n => n.Contains(".Migrations.") && n.EndsWith(".sql"));

            await Database.MigrateAsync(ThrowawayDatabase.DataSource);
            await Database.MigrateAsync(ThrowawayDatabase.DataSource);

            await using var command = ThrowawayDatabase.DataSource.CreateCommand("select count(*) from schema_migration");
            Assert.That((long)(await command.ExecuteScalarAsync())!, Is.EqualTo(scripts).And.GreaterThan(0));
        }

        [Test]
        public async Task ProfileLevelsRoundTripBitExact()
        {
            var dna = Generator.GetRandomDNA(300_000, 20261015);
            var ladder = RecordStore.Ladder(new ZCurveWindowedProfile(dna));
            var sequenceId = await ThrowawayDatabase.InsertRecordAsync("test-roundtrip", "rt-1", dna.Length);

            await using (var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync())
            await using (var transaction = await connection.BeginTransactionAsync())
            {
                foreach (var factor in ladder.Factors)
                    Assert.That(await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, ladder.Level(factor)), Is.True);
                await transaction.CommitAsync();
            }

            foreach (var factor in ladder.Factors)
            {
                var expected = ladder.Level(factor);
                var read = await ProfileLevels.ReadAsync(ThrowawayDatabase.DataSource, sequenceId, method, methodVersion, expected.WindowSize);
                Assert.That(read, Is.Not.Null);
                Assert.That(read!.ToBinary(), Is.EqualTo(expected.ToBinary()), $"{expected.WindowSize} bp level");

                var streamed = new List<long>();
                await foreach (var (id, level) in ProfileLevels.StreamAsync(ThrowawayDatabase.DataSource, method, methodVersion, expected.WindowSize))
                {
                    streamed.Add(id);
                    Assert.That(level.ToBinary(), Is.EqualTo(expected.ToBinary()));
                }
                Assert.That(streamed, Is.EqualTo(new[] { sequenceId }));
            }

            Assert.That(await ProfileLevels.ReadAsync(ThrowawayDatabase.DataSource, sequenceId, method, methodVersion, 123), Is.Null);
        }

        [Test]
        public async Task LevelOverInlineLimitIsRefusedAndNothingWritten()
        {
            // 16-byte header + 6 bytes per window (ZCWP version 2): one window past the 1 MB limit.
            int windows = (ProfileLevels.MaxInlineBytes - 16) / 6 + 1;
            var level = ZCurveWindowedProfile.FromPrecomputedWindows(500, Enumerable.Repeat((0.0, 0.0, 0.0), windows).ToArray());
            var sequenceId = await ThrowawayDatabase.InsertRecordAsync("test-oversize", "big-1", 500L * windows);

            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            Assert.That(await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, 1, level), Is.False);
            await transaction.CommitAsync();

            Assert.That(await ProfileLevels.ReadAsync(ThrowawayDatabase.DataSource, sequenceId, method, methodVersion, 500), Is.Null);
        }

        [Test]
        public async Task KeysetPagesStayStableWhileRowsAreInserted()
        {
            const string provider = "test-keyset";
            var filter = new CatalogueFilter(Provider: provider);
            var inserted = new List<long>();
            for (int i = 0; i < 25; i++)
                inserted.Add(await ThrowawayDatabase.InsertRecordAsync(provider, $"ks-{i}", 1000 + i));

            var seen = new List<long>();
            var page = await Catalogue.ListAsync(ThrowawayDatabase.DataSource, filter, after: 0, limit: 10);
            seen.AddRange(page.Records.Select(r => r.Id));

            // Rows arriving mid-walk get larger ids, so they land after the cursor, never inside a page already read.
            for (int i = 25; i < 30; i++)
                inserted.Add(await ThrowawayDatabase.InsertRecordAsync(provider, $"ks-{i}", 1000 + i));

            while (page.NextAfter is { } after)
            {
                page = await Catalogue.ListAsync(ThrowawayDatabase.DataSource, filter, after, limit: 10);
                seen.AddRange(page.Records.Select(r => r.Id));
            }

            Assert.That(seen, Is.EqualTo(inserted));
            Assert.That(seen, Is.Unique.And.Ordered);
        }

        [Test]
        public async Task CatalogueFiltersSelectOnlyMatchingRows()
        {
            const string provider = "test-filter";
            var shortVirus = await ThrowawayDatabase.InsertRecordAsync(provider, "f-short", 5_000);
            var longVirus = await ThrowawayDatabase.InsertRecordAsync(provider, "f-long", 50_000);
            var bacterium = await ThrowawayDatabase.InsertRecordAsync(provider, "f_bact%", 4_000_000, kind: "bacterium");

            async Task<long[]> Ids(CatalogueFilter filter) =>
                (await Catalogue.ListAsync(ThrowawayDatabase.DataSource, filter with { Provider = provider })).Records.Select(r => r.Id).ToArray();

            Assert.That(await Ids(new CatalogueFilter(Kind: "virus")), Is.EqualTo(new[] { shortVirus, longVirus }));
            var archaeon = await ThrowawayDatabase.InsertRecordAsync(provider, "f-arch", 2_000_000, kind: "archaeon");
            Assert.That(await Ids(new CatalogueFilter(Kind: "prokaryote")), Is.EqualTo(new[] { bacterium, archaeon }), "the group is bacterium + archaeon");
            Assert.That(await Ids(new CatalogueFilter(MinLengthBp: 10_000, MaxLengthBp: 100_000)), Is.EqualTo(new[] { longVirus }));
            // Wildcards in the search text are literal: "_bact%" must not match "f-short" or "f-long".
            Assert.That(await Ids(new CatalogueFilter(NameContains: "_BACT%")), Is.EqualTo(new[] { bacterium }));
            Assert.That(() => Catalogue.ListAsync(ThrowawayDatabase.DataSource, new CatalogueFilter(Rank: "tribe", Taxon: "x")),
                Throws.ArgumentException);
        }
    }
}
