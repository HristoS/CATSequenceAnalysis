using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <param name="OffsetBp">Where the fragment starts on the record's stored strand, to about one 500 bp window
    /// (within ±0.6 windows over the local catalogue, results/S35_LocalizeCheckMin4.csv).</param>
    /// <param name="Strand">Reverse: the fragment's reverse complement is what matches the stored strand there.</param>
    /// <param name="CoarseScore">The score at the coarse level that brought this record into the refined top.</param>
    /// <param name="P">From score_null's fragment null (S3.6), built from full 500 bp scans; a refine tries fewer offsets, so this p is conservative.</param>
    /// <param name="Log10E">log10 of p × the records scanned.</param>
    public sealed record LocalizationHit(long SequenceId, long OffsetBp, Strand Strand, double Score, double CoarseScore,
        PValue? P = null, double? Log10E = null);

    /// <param name="CoarseWindowSize">The level the records were scanned at; 500 when the fragment is localized at the fine level only.</param>
    /// <param name="Candidates">Records scanned: at least as long as the fragment, with a stored level at that window size.</param>
    /// <param name="NotComparable">Records long enough but without a stored level at that window size.</param>
    /// <param name="WeakSignal">
    /// The fragment is under 8 kbp. At 4 kbp the source was among the catalogue's top 100 in 23 and 26 of 30 trials (bacteria,
    /// viruses) and first in 11 and 13 (results/S35_LocalizeCheckMin4_Summary.csv); inside a named record, 1.5 / 3 kbp fragments
    /// found their own position first in 38% / 77% of bacterial cases (results/R6b_GeneAligned_Summary.csv).
    /// </param>
    public sealed record LocalizationResult(int CoarseWindowSize, int Candidates, IReadOnlyList<LocalizationHit> Hits,
        IReadOnlyList<long> NotComparable, bool WeakSignal);

    /// <summary>
    /// Localize (docs/database-decision.md section 5.2, P3 step S3.5): where in the catalogue a fragment comes from. Candidates
    /// are every record at least as long as the fragment, with no length band: the score here is containment, which find-similar
    /// avoids. The scan runs at the coarsest stored level that leaves the fragment <see cref="ZCurveProfileLadder.DefaultMinFragmentWindows"/>
    /// windows (500 bp below 8 kbp, 2 kbp below 32 kbp, 8 kbp below 128 kbp, then 32 kbp), with every phase of the fragment at that
    /// level and both query strands; the best <see cref="RefineDepth"/> records are then refined at 500 bp within ±k fine windows
    /// of their coarse placement (<see cref="ZCurveWindowedProfile.RefineAlignment"/>), and ranked by the refined score.
    ///
    /// Over 540 fragments from the bundled genomes (results/S35_LocalizeCheckMin4_Summary.csv) the source was in the coarse top 100
    /// in every trial from 10 to 200 kbp, and the refined top hit was the source or, where that differed, a close relative (other
    /// strains of the species, other poxviruses: same genus in 28-30 of 30 at 10 kbp and above). Fewer coarse windows lost the
    /// source in about half the trials (results/S35_LocalizeCheck_Summary.csv), and so did scanning query phase 0 alone.
    ///
    /// Fragments under <see cref="MinCatalogueLengthBp"/> can't be localized over the catalogue: at 2 kbp the source was in the top 100
    /// in half the trials. They can be placed inside one named record (<see cref="LocalizeInRecordAsync"/>).
    /// Not handled: a fragment across the origin of a circular record, and records contained in the query (the other containment).
    /// </summary>
    public sealed class Localization
    {
        public const long MinCatalogueLengthBp = 4_000;
        public const long WeakSignalLengthBp = 8_000;
        public const int RefineDepth = 100;

        private readonly NpgsqlDataSource dataSource;
        private readonly string method;
        private readonly string methodVersion;
        private readonly int[] factors = ZCurveProfileLadder.DefaultFactors.ToArray();
        private readonly ProfileCache[]? caches;
        private readonly bool[] loaded;
        private ScoreNull? scoreNull;
        private bool nullLoaded;

        /// <summary>Bytes of window data held over every loaded level; 0 without a cache.</summary>
        public long WindowBytes => caches is null ? 0 : caches.Where((_, i) => loaded[i]).Sum(c => c.WindowBytes);

        /// <param name="useCache">false: every localization reads its level from the database (one query per process, as the CLI does).</param>
        public Localization(NpgsqlDataSource dataSource, bool useCache = true,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion)
        {
            this.dataSource = dataSource;
            this.method = method;
            this.methodVersion = methodVersion;
            loaded = new bool[factors.Length];
            if (useCache)
                caches = Enumerable.Range(0, factors.Length)
                    .Select(i => new ProfileCache(dataSource, WindowSizeAt(i), method, methodVersion, CachedLengths(i)))
                    .ToArray();
        }

        private static int WindowSizeAt(int level) => RecordStore.FineWindowSize * ZCurveProfileLadder.DefaultFactors[level];

        // A fragment is scanned at level i only when it keeps 4 windows there, so it is at least 4 of that level's windows long
        // (and at least MinCatalogueLengthBp), and its candidates are at least as long as it is.
        public static LengthRange CachedLengths(int level) =>
            new(Math.Max(MinCatalogueLengthBp, (long)ZCurveProfileLadder.DefaultMinFragmentWindows * WindowSizeAt(level)), null);

        /// <summary>Brings every level loaded so far up to date with the database. A level is first loaded when a localization needs it.</summary>
        public async Task RefreshAsync(CancellationToken cancellationToken = default)
        {
            nullLoaded = false;
            if (caches is null) return;
            for (int i = 0; i < caches.Length; i++)
                if (loaded[i]) await caches[i].RefreshAsync(cancellationToken);
        }

        public async Task<LocalizationResult> LocalizeAsync(QueryProfile fragment, CatalogueFilter? filter = null, int top = 20,
            CancellationToken cancellationToken = default)
        {
            if (top < 1 || top > RefineDepth) throw new ArgumentOutOfRangeException(nameof(top), $"top must be 1 to {RefineDepth}.");
            if (fragment.LengthBp < MinCatalogueLengthBp
                || ZCurveProfileLadder.SelectQueryFactor(fragment.Forward, ZCurveProfileLadder.DefaultMinFragmentWindows) is not int factor)
                throw new ArgumentException($"A fragment under {MinCatalogueLengthBp:N0} bp can't be localized over the catalogue; " +
                    "localize it inside a named record instead.", nameof(fragment));

            filter ??= new CatalogueFilter();
            var candidates = await Catalogue.IdsAsync(dataSource,
                filter with { MinLengthBp = Math.Max(filter.MinLengthBp ?? 0, fragment.LengthBp) }, cancellationToken);
            int level = Array.IndexOf(factors, factor);
            var profiles = await LevelAsync(level, new LengthRange(fragment.LengthBp, null), cancellationToken);
            var significance = await SignificanceAsync(cancellationToken);
            var comparable = candidates.Where(profiles.ContainsKey).ToArray();
            var notComparable = candidates.Where(id => !profiles.ContainsKey(id)).Order().ToList();

            var strands = Strands(fragment);
            var placements = new Placement[comparable.Length];
            var phased = strands.Select(s => Enumerable.Range(0, factor).Select(p => s.Downsample(factor, p)).ToArray()).ToArray();
            Parallel.For(0, comparable.Length, new ParallelOptions { CancellationToken = cancellationToken },
                i => placements[i] = Place(phased, profiles[comparable[i]]));

            var best = Enumerable.Range(0, comparable.Length)
                .Where(i => placements[i].Offset >= 0)
                .OrderByDescending(i => placements[i].Score).ThenBy(i => comparable[i])
                .Take(RefineDepth)
                .ToArray();
            var fine = factor == 1 ? profiles : await FineAsync(best.Select(i => comparable[i]).ToArray(), cancellationToken);

            var hits = new List<LocalizationHit>();
            foreach (var i in best)
            {
                var placement = placements[i];
                if (!fine.TryGetValue(comparable[i], out var reference)) continue;
                var (score, offset) = factor == 1
                    ? (placement.Score, placement.Offset)
                    : ZCurveWindowedProfile.RefineAlignment(strands[placement.Strand], reference,
                        placement.Offset * factor - placement.Phase, radius: factor);
                if (offset < 0) continue;

                var (shorter, offsets) = ScoreNull.Shape(strands[placement.Strand].WindowCount, reference.WindowCount);
                var p = significance?.PValue(RecordStore.FineWindowSize, shorter, offsets, score);
                hits.Add(new LocalizationHit(comparable[i], (long)offset * RecordStore.FineWindowSize, (Strand)placement.Strand, score,
                    placement.Score, p, p?.Log10 + Math.Log10(comparable.Length)));
            }

            return new LocalizationResult(WindowSizeAt(level), comparable.Length,
                hits.OrderByDescending(h => h.Score).ThenBy(h => h.SequenceId).Take(top).ToList(), notComparable,
                fragment.LengthBp < WeakSignalLengthBp);
        }

        // Index 0 is the forward strand, 1 the reverse when the query has one; the index is the Strand value.
        private static ZCurveWindowedProfile[] Strands(QueryProfile fragment) =>
            fragment.Reverse is null ? [fragment.Forward.Fine] : [fragment.Forward.Fine, fragment.Reverse.Fine];

        private readonly record struct Placement(double Score, int Strand, int Phase, int Offset);

        // The best coarse placement over both strands and every phase; Offset -1 when no phase fits in the record.
        private static Placement Place(ZCurveWindowedProfile[][] phased, ZCurveWindowedProfile candidate)
        {
            var best = new Placement(double.NegativeInfinity, 0, 0, -1);
            for (int strand = 0; strand < phased.Length; strand++)
                for (int phase = 0; phase < phased[strand].Length; phase++)
                {
                    var query = phased[strand][phase];
                    if (query.WindowCount == 0 || query.WindowCount > candidate.WindowCount) continue;
                    var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, candidate);
                    if (score > best.Score) best = new Placement(score, strand, phase, offset);
                }
            return best;
        }

        // The refine's 500 bp profiles: from the cache when the fine level is loaded, otherwise one database read for all of them.
        private async Task<IReadOnlyDictionary<long, ZCurveWindowedProfile>> FineAsync(long[] ids, CancellationToken cancellationToken)
        {
            if (caches is not null && loaded[0] && ids.All(caches[0].Profiles.ContainsKey)) return caches[0].Profiles;
            return await ProfileLevels.ReadManyAsync(dataSource, ids, method, methodVersion, RecordStore.FineWindowSize, cancellationToken);
        }

        /// <summary>
        /// Where in one record a fragment of any length sits: a full 500 bp scan of that record, both query strands. The way to
        /// place fragments under <see cref="MinCatalogueLengthBp"/>. No hit when the fragment has no full window or is longer than the record.
        /// </summary>
        public async Task<LocalizationResult> LocalizeInRecordAsync(QueryProfile fragment, long sequenceId, CancellationToken cancellationToken = default)
        {
            var reference = await ProfileLevels.ReadAsync(dataSource, sequenceId, method, methodVersion, RecordStore.FineWindowSize)
                ?? throw new ArgumentException($"Record {sequenceId} has no stored {RecordStore.FineWindowSize} bp level.", nameof(sequenceId));
            var significance = await SignificanceAsync(cancellationToken);
            bool weak = fragment.LengthBp < WeakSignalLengthBp;

            var strands = Strands(fragment);
            if (strands[0].WindowCount == 0 || strands[0].WindowCount > reference.WindowCount)
                return new LocalizationResult(RecordStore.FineWindowSize, 1, [], [], weak);

            var best = (Score: double.NegativeInfinity, Offset: -1, Strand: Strand.Forward);
            for (int strand = 0; strand < strands.Length; strand++)
            {
                var (score, offset) = ZCurveWindowedProfile.BestAlignment(strands[strand], reference);
                if (score > best.Score) best = (score, offset, (Strand)strand);
            }
            var (shorter, offsets) = ScoreNull.Shape(strands[0].WindowCount, reference.WindowCount);
            var p = significance?.PValue(RecordStore.FineWindowSize, shorter, offsets, best.Score);
            var hit = new LocalizationHit(sequenceId, (long)best.Offset * RecordStore.FineWindowSize, best.Strand, best.Score, best.Score, p, p?.Log10);
            return new LocalizationResult(RecordStore.FineWindowSize, 1, [hit], [], weak);
        }

        private async Task<ScoreNull?> SignificanceAsync(CancellationToken cancellationToken)
        {
            if (!nullLoaded)
            {
                scoreNull = await ScoreNullStore.LoadAsync(dataSource, method, methodVersion, ScoreNullStore.Fragment, cancellationToken);
                nullLoaded = true;
            }
            return scoreNull;
        }

        private async Task<IReadOnlyDictionary<long, ZCurveWindowedProfile>> LevelAsync(int level, LengthRange lengths,
            CancellationToken cancellationToken)
        {
            if (caches is not null)
            {
                if (!loaded[level])
                {
                    await caches[level].RefreshAsync(cancellationToken);
                    loaded[level] = true;
                }
                return caches[level].Profiles;
            }

            var profiles = new Dictionary<long, ZCurveWindowedProfile>();
            await foreach (var (_, sequenceId, profile) in ProfileLevels.StreamAfterAsync(
                dataSource, method, methodVersion, WindowSizeAt(level), 0, lengths, cancellationToken))
                profiles[sequenceId] = profile;
            return profiles;
        }
    }
}
