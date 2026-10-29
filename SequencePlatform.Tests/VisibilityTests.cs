using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>D3: an upload is visible to its owner only, in listings, details and search candidates alike.</summary>
    public class VisibilityTests
    {
        private const string method = "zcurve-windowed";
        private const string methodVersion = "visibility";

        private static async Task<AppUser> NewUserAsync(string name) =>
            (await Users.CreateAsync(ThrowawayDatabase.DataSource, $"{name}-{Guid.NewGuid():N}")).User;

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

        [Test]
        public async Task OnlyTheKeyHashIsStoredAndADisabledKeyFails()
        {
            var name = $"keys-{Guid.NewGuid():N}";
            var (user, apiKey) = await Users.CreateAsync(ThrowawayDatabase.DataSource, name);

            Assert.That(apiKey, Does.StartWith(Users.KeyPrefix).And.Length.GreaterThan(40));
            Assert.That(await Users.AuthenticateAsync(ThrowawayDatabase.DataSource, apiKey), Is.EqualTo(user));
            Assert.That(await Users.AuthenticateAsync(ThrowawayDatabase.DataSource, apiKey + "x"), Is.Null);
            await using (var row = ThrowawayDatabase.DataSource.CreateCommand("select position($2 in t::text) from app_user t where name = $1"))
            {
                row.Parameters.AddWithValue(name);
                row.Parameters.AddWithValue(apiKey);
                Assert.That(await row.ExecuteScalarAsync(), Is.EqualTo(0), "the key itself is in no column");
            }

            Assert.That(await Users.DisableAsync(ThrowawayDatabase.DataSource, name), Is.True);
            Assert.That(await Users.AuthenticateAsync(ThrowawayDatabase.DataSource, apiKey), Is.Null);
        }

        [Test]
        public async Task CatalogueShowsAnUploadToItsOwnerOnly()
        {
            var alice = await NewUserAsync("alice");
            var bob = await NewUserAsync("bob");
            var shared = await ThrowawayDatabase.InsertRecordAsync("visibility", "public", 10_000);
            var mine = await ThrowawayDatabase.InsertUploadAsync(alice.Id, "alice-1", 10_000);
            var theirs = await ThrowawayDatabase.InsertUploadAsync(bob.Id, "bob-1", 10_000);

            var anonymous = await Catalogue.IdsAsync(ThrowawayDatabase.DataSource, new CatalogueFilter());
            var asAlice = await Catalogue.IdsAsync(ThrowawayDatabase.DataSource, new CatalogueFilter(ViewerId: alice.Id));

            Assert.That(anonymous, Does.Contain(shared).And.Not.Contain(mine).And.Not.Contain(theirs), "no viewer: public records only");
            Assert.That(asAlice, Does.Contain(shared).And.Contain(mine).And.Not.Contain(theirs));
            Assert.That((await Catalogue.DetailAsync(ThrowawayDatabase.DataSource, mine, alice.Id))!.Owned, Is.True);
            Assert.That((await Catalogue.DetailAsync(ThrowawayDatabase.DataSource, shared, alice.Id))!.Owned, Is.False);
            Assert.That(await Catalogue.DetailAsync(ThrowawayDatabase.DataSource, theirs, alice.Id), Is.Null);
            Assert.That(await Catalogue.IsVisibleAsync(ThrowawayDatabase.DataSource, theirs, alice.Id), Is.False);
            Assert.That(await Catalogue.IsVisibleAsync(ThrowawayDatabase.DataSource, theirs, bob.Id), Is.True);
        }

        // The search reaches its candidates through Catalogue.IdsAsync, so the same rule holds there.
        [Test]
        public async Task SearchCandidatesAreThePublicCatalogueAndTheViewersOwnUploads()
        {
            var alice = await NewUserAsync("alice");
            var bob = await NewUserAsync("bob");
            var dna = Generator.GetRandomDNA(60_000, seed: 11);
            var mine = await ThrowawayDatabase.InsertUploadAsync(alice.Id, "alice-search", dna.Length);
            var theirs = await ThrowawayDatabase.InsertUploadAsync(bob.Id, "bob-search", dna.Length);
            await StoreLevelsAsync(mine, dna);
            await StoreLevelsAsync(theirs, dna);
            var search = new SimilaritySearch(ThrowawayDatabase.DataSource, useCache: false, method: method, methodVersion: methodVersion);
            var query = QueryProfile.FromSequence(dna);

            var asAlice = await search.SearchAsync(query, new CatalogueFilter(ViewerId: alice.Id));
            var anonymous = await search.SearchAsync(query, new CatalogueFilter());

            Assert.That(asAlice.Ranked.Select(h => h.SequenceId), Is.EqualTo(new[] { mine }));
            Assert.That(asAlice.NotComparable, Does.Not.Contain(theirs));
            Assert.That(anonymous.Ranked, Is.Empty);
            Assert.That(anonymous.NotComparable, Does.Not.Contain(mine).And.Not.Contain(theirs));
        }

        [Test]
        public async Task SearchRunKeepsItsOwnerOutsideTheFilters()
        {
            var alice = await NewUserAsync("alice");
            var request = new SearchRequest(SearchKind.Similar, new CatalogueFilter(Kind: "virus", ViewerId: alice.Id), 10, LengthRatio: 0.5);

            var id = await SearchLog.EnqueueAsync(ThrowawayDatabase.DataSource, request, QueryProfile.FromSequence(Generator.GetRandomDNA(20_000, seed: 12)));

            Assert.That((await SearchLog.GetAsync(ThrowawayDatabase.DataSource, id))!.Request, Is.EqualTo(request));
            await using var row = ThrowawayDatabase.DataSource.CreateCommand("select owner_id, filters::text from search_run where public_id = $1");
            row.Parameters.AddWithValue(id);
            await using var reader = await row.ExecuteReaderAsync();
            await reader.ReadAsync();
            Assert.That((reader.GetInt64(0), reader.GetString(1)), Is.EqualTo((alice.Id, "{\"Kind\": \"virus\"}")));
        }
    }
}
