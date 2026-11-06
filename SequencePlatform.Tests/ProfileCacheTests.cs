using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class ProfileCacheTests
    {
        private const string method = "zcurve-windowed";

        // Each test stores its levels under its own method version, so caches never see another test's rows.
        private static async Task<(long SequenceId, ZCurveProfileLadder Ladder)> StoreAsync(string methodVersion, string accession, int length, int seed)
        {
            var dna = Generator.GetRandomDNA(length, seed);
            var ladder = RecordStore.Ladder(new ZCurveWindowedProfile(dna));
            var sequenceId = await ThrowawayDatabase.InsertRecordAsync($"cache-{methodVersion}", accession, length);

            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var factor in ladder.Factors)
                if (ladder.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, ladder.Level(factor));
            await transaction.CommitAsync();
            return (sequenceId, ladder);
        }

        private static async Task DeleteRecordAsync(long sequenceId)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand("delete from sequence_record where id = $1");
            command.Parameters.AddWithValue(sequenceId);
            await command.ExecuteNonQueryAsync();
        }

        [Test]
        public async Task LoadsOneLevelOfEveryRecordBitExact()
        {
            var stored = new[]
            {
                await StoreAsync("load", "a", 40_000, seed: 1),
                await StoreAsync("load", "b", 300_000, seed: 2),
                await StoreAsync("load", "c", 1_500, seed: 3),   // no 2 kbp level
            };
            var cache = new ProfileCache(ThrowawayDatabase.DataSource, 2_000, method, "load");

            Assert.That(await cache.RefreshAsync(), Is.EqualTo((2, false)));
            Assert.That(cache.Profiles.Keys, Is.EquivalentTo(stored.Take(2).Select(s => s.SequenceId)));
            foreach (var (sequenceId, ladder) in stored.Take(2))
                Assert.That(cache.Profiles[sequenceId].ToBinary(), Is.EqualTo(ladder.Level(4).ToBinary()));
            Assert.That(cache.WindowCount, Is.EqualTo(20 + 150));
            Assert.That(cache.WindowBytes, Is.EqualTo(170 * 12));
        }

        [Test]
        public async Task RefreshAddsOnlyLevelsStoredSinceTheLastOne()
        {
            await StoreAsync("incremental", "a", 40_000, seed: 1);
            var cache = new ProfileCache(ThrowawayDatabase.DataSource, 500, method, "incremental");
            Assert.That(await cache.RefreshAsync(), Is.EqualTo((1, false)));

            var (added, _) = await StoreAsync("incremental", "b", 10_000, seed: 2);

            Assert.That(await cache.RefreshAsync(), Is.EqualTo((1, false)));
            Assert.That(await cache.RefreshAsync(), Is.EqualTo((0, false)));
            Assert.That(cache.Profiles, Has.Count.EqualTo(2).And.ContainKey(added));
            Assert.That(cache.WindowCount, Is.EqualTo(80 + 20));
        }

        [Test]
        public async Task DeletedRecordForcesAFullReload()
        {
            var (kept, _) = await StoreAsync("deleted", "a", 40_000, seed: 1);
            var (deleted, _) = await StoreAsync("deleted", "b", 10_000, seed: 2);
            var cache = new ProfileCache(ThrowawayDatabase.DataSource, 500, method, "deleted");
            await cache.RefreshAsync();

            await DeleteRecordAsync(deleted);

            Assert.That(await cache.RefreshAsync(), Is.EqualTo((1, true)), "an id watermark alone can't see a deletion");
            Assert.That(cache.Profiles.Keys, Is.EquivalentTo(new[] { kept }));
            Assert.That(cache.WindowCount, Is.EqualTo(80));
        }

        [Test]
        public async Task LengthRangeHoldsOnlyRecordsInIt()
        {
            await StoreAsync("range", "short", 10_000, seed: 1);
            var (inRange, _) = await StoreAsync("range", "middle", 40_000, seed: 2);
            await StoreAsync("range", "long", 300_000, seed: 3);
            var cache = new ProfileCache(ThrowawayDatabase.DataSource, 500, method, "range", new LengthRange(20_000, 300_000));

            Assert.That(await cache.RefreshAsync(), Is.EqualTo((1, false)), "the upper end is exclusive");
            Assert.That(cache.Profiles.Keys, Is.EqualTo(new[] { inRange }));
            Assert.That(await cache.RefreshAsync(), Is.EqualTo((0, false)), "records outside the range don't look like deletions");
        }

        [Test]
        public async Task SnapshotTakenBeforeARefreshDoesNotChange()
        {
            await StoreAsync("snapshot", "a", 40_000, seed: 1);
            var cache = new ProfileCache(ThrowawayDatabase.DataSource, 500, method, "snapshot");
            await cache.RefreshAsync();
            var before = cache.Profiles;

            await StoreAsync("snapshot", "b", 10_000, seed: 2);
            await cache.RefreshAsync();

            Assert.That(before, Has.Count.EqualTo(1));
            Assert.That(cache.Profiles, Has.Count.EqualTo(2));
        }
    }
}
