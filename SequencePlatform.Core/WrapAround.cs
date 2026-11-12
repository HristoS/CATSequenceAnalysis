using System.Runtime.CompilerServices;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// Wrap-around scoring for viruses. A virus genome is often deposited from an arbitrary start (circularly permuted phages), so a
    /// relative started elsewhere scored lower when the query only slid along it: Pantoea phage ROV85 against Escherichia phage Lidtsur
    /// 0.870 as deposited, 0.924 from a start 24.5 kbp on. A virus candidate is therefore scored as its windows followed by all but the
    /// last of themselves, so a query meets it from every start. Over 400 virus queries this raised the mean same-genus vs other-family
    /// AUC from 0.904 to 0.934 and same-genus records in the top 10 from 6.28 to 7.14, while 1.2% of other-family pairs gained more than
    /// 0.03 against 22.8% of same-genus pairs (results/CircularPermutation_Check.csv). Bacterial chromosomes get origin rotation instead
    /// (<see cref="OriginRotation"/>); eukaryote chromosomes are linear. The score null scores virus candidates the same way, so their
    /// p-values come from cells with the wrapped number of offsets.
    /// </summary>
    public static class WrapAround
    {
        // Each cached profile's wrapped form, made once and dropped with the profile when a cache refresh replaces it.
        private static readonly ConditionalWeakTable<ZCurveWindowedProfile, ZCurveWindowedProfile> wrapped = new();

        public static bool AppliesTo(string kind) => kind == "virus";

        /// <summary>The profile's windows followed by all but its last window again: every start, each once.</summary>
        public static ZCurveWindowedProfile Of(ZCurveWindowedProfile profile) => wrapped.GetValue(profile, Wrap);

        /// <summary>
        /// The ids of the records scored wrapped (viruses), read in full once and afterwards only those stored since: a record's kind
        /// never changes, and a deleted id is simply never a candidate. Each read publishes a new set, so a search keeps the one it took.
        /// </summary>
        public sealed class Ids(NpgsqlDataSource dataSource)
        {
            private readonly SemaphoreSlim reading = new(1, 1);
            private volatile HashSet<long> current = [];
            private long lastId;

            public async Task<IReadOnlySet<long>> CurrentAsync(CancellationToken cancellationToken = default)
            {
                await reading.WaitAsync(cancellationToken);
                try
                {
                    await using var command = dataSource.CreateCommand("select id from sequence_record where kind = 'virus' and id > $1");
                    command.Parameters.AddWithValue(lastId);
                    HashSet<long>? next = null;
                    await using var reader = await command.ExecuteReaderAsync(cancellationToken);
                    while (await reader.ReadAsync(cancellationToken))
                    {
                        long id = reader.GetInt64(0);
                        (next ??= new HashSet<long>(current)).Add(id);
                        lastId = Math.Max(lastId, id);
                    }
                    if (next is not null) current = next;
                    return current;
                }
                finally { reading.Release(); }
            }
        }

        private static ZCurveWindowedProfile Wrap(ZCurveWindowedProfile profile)
        {
            var windows = profile.Windows;
            return ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize,
                windows.Concat(windows.Take(Math.Max(0, windows.Count - 1))).Select(w => (w.X, w.Y, w.Z)).ToList());
        }
    }
}
