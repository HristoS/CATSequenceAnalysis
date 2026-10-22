using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    public class NeighbourGraphTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "neighbours";
        private const string provider = "neighbour-test";

        private static readonly List<long> records = new();
        private static long firstRecord, sameSourceRecord;

        private static async Task StoreLevelsAsync(long sequenceId, string dna)
        {
            var ladder = new ZCurveProfileLadder(new ZCurveWindowedProfile(dna));
            await using var connection = await ThrowawayDatabase.DataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var factor in ladder.Factors)
                if (ladder.Level(factor).WindowCount > 0)
                    await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, ladder.Level(factor));
            await transaction.CommitAsync();
        }

        // Ten records of 60-105 kbp, each its own source, plus a second record in the first one's source
        // (a chromosome and its plasmid, say), then the graph over all of them.
        [OneTimeSetUp]
        public async Task BuildCatalogueAndGraphAsync()
        {
            for (int i = 0; i < 10; i++)
            {
                var dna = Generator.GetRandomDNA(60_000 + 5_000 * i, seed: 300 + i);
                var id = await ThrowawayDatabase.InsertRecordAsync(provider, $"n{i}", dna.Length, kind: "bacterium");
                await StoreLevelsAsync(id, dna);
                records.Add(id);
            }
            firstRecord = records[0];

            var plasmid = Generator.GetRandomDNA(61_000, seed: 399);
            await using (var insert = ThrowawayDatabase.DataSource.CreateCommand(
                "insert into sequence_record (source_id, record_id, name, kind, data_class, length_bp) " +
                "select source_id, 'n0-plasmid', 'test n0 plasmid', 'plasmid', 'public_microbe', $2 from sequence_record where id = $1 returning id"))
            {
                insert.Parameters.AddWithValue(firstRecord);
                insert.Parameters.AddWithValue((long)plasmid.Length);
                sameSourceRecord = (long)(await insert.ExecuteScalarAsync())!;
            }
            await StoreLevelsAsync(sameSourceRecord, plasmid);
            records.Add(sameSourceRecord);

            Assert.That(await NeighbourGraph.BuildAsync(ThrowawayDatabase.DataSource, 0.5, SimilaritySearch.NeighbourCount, _ => { },
                method, methodVersion), Is.EqualTo(records.Count));
        }

        private static SimilaritySearch NewSearch() => new(ThrowawayDatabase.DataSource, method: method, methodVersion: methodVersion);

        [Test]
        public async Task SearchByRecordLeavesOutItsOwnSource()
        {
            var result = await NewSearch().SearchRecordAsync(firstRecord, useNeighbours: false, top: 100);
            var ids = result.Ranked.Select(h => h.SequenceId).ToList();

            Assert.That(ids, Does.Not.Contain(firstRecord).And.Not.Contain(sameSourceRecord));
            Assert.That(ids, Is.EquivalentTo(records.Except(new[] { firstRecord, sameSourceRecord })));
        }

        [Test]
        public async Task StoredListsAreEachRecordsOwnSearch()
        {
            var graph = await NeighbourGraph.LoadAsync(ThrowawayDatabase.DataSource, 0.5, method, methodVersion);
            var search = NewSearch();

            Assert.That(graph.Keys, Is.EquivalentTo(records));
            foreach (var id in records)
            {
                var live = await search.SearchRecordAsync(id, top: SimilaritySearch.NeighbourCount, useNeighbours: false);
                Assert.That(graph[id], Is.EqualTo(live.Ranked.Select(h => h.SequenceId)), $"record {id}");
            }
            Assert.That(await NeighbourGraph.LoadAsync(ThrowawayDatabase.DataSource, 0.75, method, methodVersion), Is.Empty);
        }

        [Test]
        public async Task TopHitsAreOrderedBySharedNeighboursThenScore()
        {
            var query = QueryProfile.FromSequence(Generator.GetRandomDNA(80_000, seed: 1));
            var filter = new CatalogueFilter(Provider: provider);

            var ranked = await NewSearch().SearchAsync(query, filter, top: 100);
            var unranked = await NewSearch().SearchAsync(query, filter, lengthRatio: 0.75, top: 100);   // no graph at 0.75

            Assert.That(ranked.NeighbourRanked, Is.True);
            Assert.That(ranked.Ranked.All(h => h.SharedNeighbours is not null), Is.True);
            Assert.That(ranked.Ranked, Is.Ordered.By(nameof(SimilarityHit.SharedNeighbours)).Descending
                .Then.By(nameof(SimilarityHit.Score)).Descending);
            Assert.That(unranked.NeighbourRanked, Is.False);
            Assert.That(unranked.Ranked.Select(h => h.Score), Is.Ordered.Descending);
            Assert.That(unranked.Ranked.All(h => h.SharedNeighbours is null), Is.True);
        }

        [Test]
        public async Task RecordMissingFromTheGraphIsSearchedLive()
        {
            var query = QueryProfile.FromSequence(Generator.GetRandomDNA(80_000, seed: 2));
            var filter = new CatalogueFilter(Provider: provider);
            var before = await NewSearch().SearchAsync(query, filter, top: 100);

            var removed = before.Ranked[0].SequenceId;
            var saved = (await NeighbourGraph.LoadAsync(ThrowawayDatabase.DataSource, 0.5, method, methodVersion))[removed];
            await Execute("delete from neighbour_list where method_version = $1 and sequence_id = $2", methodVersion, removed);
            try
            {
                var after = await NewSearch().SearchAsync(query, filter, top: 100);
                Assert.That(after.Ranked, Is.EqualTo(before.Ranked));
            }
            finally
            {
                await Execute("insert into neighbour_list (method, method_version, length_ratio, sequence_id, neighbour_ids) " +
                    "values ($1, $2, 0.5, $3, $4)", method, methodVersion, removed, saved);
            }
        }

        private static async Task Execute(string sql, params object[] values)
        {
            await using var command = ThrowawayDatabase.DataSource.CreateCommand(sql);
            foreach (var value in values) command.Parameters.AddWithValue(value);
            await command.ExecuteNonQueryAsync();
        }
    }
}
