using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// One ladder level (method, version, window size) of every stored record, held in memory for the
    /// search scan (docs/database-decision.md section 5.2, P3 step S3.3). Each profile already keeps its
    /// windows as three per-axis float arrays, the layout <see cref="ZCurveWindowedProfile.BestAlignment"/>
    /// vectorizes over, so the cache is a map from sequence id to profile and adds no copy.
    ///
    /// <see cref="RefreshAsync"/> adds the levels inserted since the last refresh (by profile_level.id), then
    /// checks the stored count up to that id against the cache. A mismatch means levels were deleted
    /// (duplicates collapsed, a source reloaded), which an id watermark can't see, and the level is
    /// reloaded in full. Each refresh publishes a new snapshot, so a search that took
    /// <see cref="Profiles"/> keeps a consistent view while a refresh runs. Refreshes themselves must not overlap.
    /// </summary>
    public sealed class ProfileCache
    {
        private readonly NpgsqlDataSource dataSource;
        private readonly string method;
        private readonly string methodVersion;
        private volatile Snapshot current = Snapshot.Empty;

        public int WindowSize { get; }

        public IReadOnlyDictionary<long, ZCurveWindowedProfile> Profiles => current.Profiles;

        public long WindowCount => current.WindowCount;

        /// <summary>Bytes of window data held: 3 axes × 4 bytes per window, without object and dictionary overhead.</summary>
        public long WindowBytes => current.WindowCount * 12;

        public ProfileCache(NpgsqlDataSource dataSource, int windowSize,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion)
        {
            this.dataSource = dataSource;
            this.method = method;
            this.methodVersion = methodVersion;
            WindowSize = windowSize;
        }

        /// <returns>How many levels were added or replaced, and whether the whole level was reloaded.</returns>
        public async Task<(int Loaded, bool FullReload)> RefreshAsync(CancellationToken cancellationToken = default)
        {
            var snapshot = current;
            var (profiles, windowCount, lastLevelId, loaded) = await LoadAfterAsync(snapshot, cancellationToken);

            var stored = await ProfileLevels.CountUpToAsync(dataSource, method, methodVersion, WindowSize, lastLevelId, cancellationToken);
            if (stored == profiles.Count)
            {
                current = new Snapshot(profiles, windowCount, lastLevelId);
                return (loaded, false);
            }

            (profiles, windowCount, lastLevelId, loaded) = await LoadAfterAsync(Snapshot.Empty, cancellationToken);
            current = new Snapshot(profiles, windowCount, lastLevelId);
            return (loaded, true);
        }

        private async Task<(Dictionary<long, ZCurveWindowedProfile> Profiles, long WindowCount, long LastLevelId, int Loaded)> LoadAfterAsync(
            Snapshot from, CancellationToken cancellationToken)
        {
            var profiles = new Dictionary<long, ZCurveWindowedProfile>(from.Profiles);
            long windowCount = from.WindowCount;
            long lastLevelId = from.LastLevelId;
            int loaded = 0;

            await foreach (var (levelId, sequenceId, level) in ProfileLevels.StreamAfterAsync(
                dataSource, method, methodVersion, WindowSize, from.LastLevelId, cancellationToken))
            {
                if (level.WindowSize != WindowSize)
                    throw new InvalidDataException($"Level {levelId} has window size {level.WindowSize}, expected {WindowSize}.");

                // A sequence's level replaced under a new id: count only the newer one.
                if (profiles.TryGetValue(sequenceId, out var replaced)) windowCount -= replaced.WindowCount;
                profiles[sequenceId] = level;
                windowCount += level.WindowCount;
                lastLevelId = levelId;
                loaded++;
            }
            return (profiles, windowCount, lastLevelId, loaded);
        }

        private sealed record Snapshot(IReadOnlyDictionary<long, ZCurveWindowedProfile> Profiles, long WindowCount, long LastLevelId)
        {
            public static readonly Snapshot Empty = new(new Dictionary<long, ZCurveWindowedProfile>(), 0, 0);
        }
    }
}
