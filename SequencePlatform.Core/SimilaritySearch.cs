using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    public enum Strand { Forward, Reverse }

    /// <param name="P">The chance an unrelated record of this shape scores as high (score_null, S3.6), and how it was obtained; null without a null.</param>
    /// <param name="Log10E">log10 of the expected number of unrelated records scoring as high among this search's candidates: p × candidates.</param>
    /// <param name="SharedNeighbours">For the top hits of a neighbour-ranked search: how many of the query's top 20 are also among this hit's own top 20.</param>
    /// <param name="Rotated">The best score came from the query or the record laid out from its replication origin (<see cref="OriginRotation"/>).</param>
    public sealed record SimilarityHit(long SequenceId, double Score, Strand Strand, PValue? P = null, double? Log10E = null,
        int? SharedNeighbours = null, bool Rotated = false);

    /// <param name="QueryWindowSize">The one level every hit is scored at; null when the query has no full fine window.</param>
    /// <param name="Lengths">The candidate length band, [ratio × query length, query length / ratio].</param>
    /// <param name="Identical">Records whose sequence_sha256 equals the query's on either strand: sequence identity, not a score.</param>
    /// <param name="NotComparable">Candidates in the band with no stored level at the query's window size.</param>
    /// <param name="FewWindows">
    /// The query is too short for its shortest candidates to keep 8 windows even at 500 bp (under about 8 kbp at ratio 0.5).
    /// The ranking then reflects composition more than order: HPV16 (15 windows) ranks GC-matched strangers above its own
    /// family (results/S36_CatalogueRanking.csv).
    /// </param>
    /// <param name="NeighbourRanked">The top hits are ordered by shared neighbours, then score (S3.4b); false when no neighbour
    /// graph exists for this length ratio, and the order is by score.</param>
    public sealed record SimilarityResult(
        int? QueryWindowSize,
        LengthRange Lengths,
        IReadOnlyList<SimilarityHit> Ranked, int RankedTotal,
        IReadOnlyList<(long SequenceId, Strand Strand)> Identical,
        IReadOnlyList<long> NotComparable,
        bool FewWindows,
        bool NeighbourRanked = false);

    /// <summary>
    /// A query as the search sees it: its ladder on both strands, its length and, when it came from DNA, the
    /// SHA-256 of both strands normalized as <see cref="SequenceFeatureReader"/> does for stored records. The
    /// reverse strand is profiled from the reverse-complement DNA, which is exact; a stored profile is never
    /// flipped (P3 S3.2, results/C3_ProfileFlip.csv). The DNA itself isn't kept.
    ///
    /// A query long enough and with a clear GC skew also carries both strands laid out from its replication origin
    /// (<see cref="OriginRotation.ToOrigin"/>). Whether a FASTA query is circular isn't known, so it is searched both ways and
    /// keeps the better score, which costs nothing where the rotation is wrong (results/S32_RotationRanking_Summary.csv).
    /// </summary>
    /// <param name="LengthBp">As length_bp is for stored records; for a profile-only query, windows × window size (a lower bound).</param>
    /// <param name="Rotated">The forward strand rotated to its origin; null when the query isn't rotated.</param>
    /// <param name="RotatedReverse">The reverse strand rotated to the same origin, when there is a reverse strand.</param>
    public sealed record QueryProfile(ZCurveProfileLadder Forward, ZCurveProfileLadder? Reverse, long LengthBp, byte[]? Sha256, byte[]? ReverseSha256,
        ZCurveProfileLadder? Rotated = null, ZCurveProfileLadder? RotatedReverse = null)
    {
        public static QueryProfile FromSequence(string sequence)
        {
            var (forward, sha256, lengthBp) = Profile(sequence);
            var (reverse, reverseSha256, _) = Profile(ReverseComplement(sequence));
            return new QueryProfile(forward, reverse, lengthBp, sha256, reverseSha256).WithOrigin();
        }

        /// <summary>
        /// A stored record as a query (search by record, the neighbour graph). The DNA is gone, so the other strand is a
        /// profile-level flip, exact only when the length is a multiple of 500 bp (results/C3_ProfileFlip.csv); its length
        /// is the stored length_bp. It is rotated only if it is a circular chromosome, as <see cref="OriginRotation.BuildAsync"/> rotates the stored records.
        /// </summary>
        public static QueryProfile FromStoredRecord(ZCurveWindowedProfile fine, long lengthBp, bool circularChromosome = false)
        {
            var query = new QueryProfile(new ZCurveProfileLadder(fine), new ZCurveProfileLadder(Flip(fine)), lengthBp, null, null);
            return circularChromosome && OriginRotation.ToOrigin(fine, lengthBp) is { } turned
                ? query with { Rotated = new ZCurveProfileLadder(turned), RotatedReverse = new ZCurveProfileLadder(Flip(turned)) }
                : query;
        }

        /// <summary>A profile-only query (a ZCWP upload): no hash, and the reverse strand only if the uploader supplied it.</summary>
        public static QueryProfile FromProfiles(ZCurveWindowedProfile forward, ZCurveWindowedProfile? reverse) =>
            FromFineProfiles(forward, reverse, (long)forward.WindowCount * forward.WindowSize, null, null);

        /// <summary>
        /// A query from its fine profiles, length and hashes, as <see cref="SearchLog"/> keeps a job: everything else is derived
        /// from the fine level, so this equals the query <see cref="FromSequence"/> built from the DNA.
        /// </summary>
        public static QueryProfile FromFineProfiles(ZCurveWindowedProfile forward, ZCurveWindowedProfile? reverse, long lengthBp,
            byte[]? sha256, byte[]? reverseSha256) =>
            new QueryProfile(new ZCurveProfileLadder(forward), reverse is null ? null : new ZCurveProfileLadder(reverse),
                lengthBp, sha256, reverseSha256).WithOrigin();

        // Rotates both strands to the origin when the forward strand places one. The reverse strand's own skew minimum is the same
        // point of the ring.
        private QueryProfile WithOrigin()
        {
            if (OriginRotation.ToOrigin(Forward.Fine, LengthBp) is not { } turned) return this;
            var reverse = Reverse?.Fine;
            return this with
            {
                Rotated = new ZCurveProfileLadder(turned),
                RotatedReverse = reverse is null ? null : new ZCurveProfileLadder(reverse.RotateTo(reverse.GcSkewExtremes().OriginWindow)),
            };
        }

        private static ZCurveWindowedProfile Flip(ZCurveWindowedProfile fine) =>
            ZCurveWindowedProfile.FromPrecomputedWindows(fine.WindowSize, fine.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

        private static (ZCurveProfileLadder Ladder, byte[] Sha256, long LengthBp) Profile(string sequence)
        {
            using var features = new SequenceFeatureReader(new StringReader(sequence));
            var fine = ZCurveWindowedProfile.FromStream(features, RecordStore.FineWindowSize);
            return (new ZCurveProfileLadder(fine), features.Finish().Sha256, features.CharactersRead);
        }

        /// <summary>
        /// The other strand over the letters <see cref="SequenceFeatureReader"/> keeps (ASCII letters, uppercased),
        /// with IUPAC codes complemented too, so its hash matches a record deposited on that strand.
        /// </summary>
        public static string ReverseComplement(string sequence)
        {
            var letters = new char[sequence.Length];
            int count = 0;
            for (int i = sequence.Length - 1; i >= 0; i--)
            {
                if (!char.IsAsciiLetter(sequence[i])) continue;
                letters[count++] = char.ToUpperInvariant(sequence[i]) switch
                {
                    'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
                    'R' => 'Y', 'Y' => 'R', 'K' => 'M', 'M' => 'K',
                    'B' => 'V', 'V' => 'B', 'D' => 'H', 'H' => 'D',
                    var self => self, // S, W, N and anything else
                };
            }
            return new string(letters, 0, count);
        }
    }

    /// <summary>
    /// Find similar (docs/database-decision.md section 5.2, P3 step S3.4). Only candidates of comparable length
    /// are compared: across very different lengths the score measures containment, because the shorter sequence
    /// slides along the longer and keeps its best stretch (results/S34_LengthGate.csv: 270 kbp phages ranked
    /// above E. coli strains). With length ratio r, a query of length L meets records of length [rL, L/r],
    /// all scored at <see cref="ZCurveProfileLadder.SelectBandFactor"/>'s level, where each has at least
    /// 8 windows: one list on one scale, the better of the two query strands. There is no 500 bp re-rank: on the
    /// bundled corpora it made bacterial close-relative top hits worse (results/MultiResolution_Rerank.csv).
    ///
    /// Once score_null's genome null exists (S3.6), each hit also carries a p-value from the cell of its own shape and
    /// an E-value over this search's candidates. The ranking stays by score: ranking by p-value lost 0.008-0.043 bacterial
    /// AUC (as p or as z) and put an identical sequence at 416, because p is bounded by each cell's sample size
    /// (results/S36_CatalogueRanking.csv, MultiResolution_ScoreNullRanking.csv).
    ///
    /// The band also bounds what has to be read. A record of length X is only ever scored at the levels its
    /// band of queries uses, so each cached level holds one length range (<see cref="CachedLengths"/>): 19 MB
    /// at r = 0.5 on the local catalogue instead of 837 MB for every level (results/S34_BandLimitedCache.csv).
    /// Without a cache, a search streams just its band's rows of its one level from the database.
    ///
    /// Once neighbour_list holds a graph for the search's length ratio (S3.4b), the top 100 hits are re-ordered by how
    /// many of the query's top 20 each hit has among its own top 20, then by score. A GC look-alike turns up in many lists
    /// without its neighbours agreeing; real relatives share theirs. Over 300 random catalogue queries this raised
    /// same-genus precision@10 from 0.161 to 0.225 for bacteria and from 0.450 to 0.475 for viruses
    /// (results/Check_MutualNeighboursSample.csv; through this code 0.131 to 0.195 and 0.462 to 0.491,
    /// results/S34b_NeighbourRanking.csv). A hit with no stored list (ingested after the graph was built) is searched live.
    /// Shared neighbours can also lift a composition look-alike whose own neighbours are the query's: for SARS-CoV-2 a
    /// Vibrio fischeri plasmid with the coronavirus GC and keto skew comes fourth, with E about 32
    /// (results/S34b_PlasmidComposition.csv). Its E-value is what marks it.
    ///
    /// Once <see cref="OriginRotation.BuildAsync"/> has stored circular chromosomes laid out from their replication origin, each such
    /// record is scored both as stored and rotated, and so is a query that has a rotated form; the best score counts. Rotating only
    /// one side would cost the queries that can't be rotated (results/S32_RotationRanking_Summary.csv).
    /// </summary>
    public sealed class SimilaritySearch
    {
        public const double DefaultLengthRatio = 0.5;
        /// <summary>The smallest length ratio a search may use (S3.4: 0.25-1 allowed), and what the cache holds by default.</summary>
        public const double SmallestLengthRatio = 0.25;
        public const int NeighbourCount = 20;
        public const int NeighbourRankDepth = 100;

        private readonly NpgsqlDataSource dataSource;
        private readonly string method;
        private readonly string methodVersion;
        private readonly int[] factors = ZCurveProfileLadder.DefaultFactors.ToArray();
        private readonly ProfileCache[]? caches;
        private readonly ProfileCache[]? rotatedCaches;
        private readonly bool[] loaded;
        private ScoreNull? scoreNull;
        private bool nullLoaded;
        private readonly Dictionary<double, IReadOnlyDictionary<long, long[]>> graphs = new();

        /// <summary>The smallest length ratio a search may ask for; the cached length ranges are sized for it.</summary>
        public double MinLengthRatio { get; }

        /// <summary>Bytes of window data held over every loaded level; 0 without a cache.</summary>
        public long WindowBytes => caches is null ? 0 : caches.Concat(rotatedCaches!).Where((_, i) => loaded[i % loaded.Length]).Sum(c => c.WindowBytes);

        /// <param name="useCache">false: every search reads its band from the database (one query per process, as the CLI does).</param>
        public SimilaritySearch(NpgsqlDataSource dataSource, bool useCache = true, double minLengthRatio = SmallestLengthRatio,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion)
        {
            if (minLengthRatio <= 0 || minLengthRatio > 1) throw new ArgumentOutOfRangeException(nameof(minLengthRatio));
            this.dataSource = dataSource;
            this.method = method;
            this.methodVersion = methodVersion;
            MinLengthRatio = minLengthRatio;
            loaded = new bool[factors.Length];
            if (useCache)
            {
                caches = Enumerable.Range(0, factors.Length)
                    .Select(i => new ProfileCache(dataSource, WindowSizeAt(i), method, methodVersion, CachedLengths(i)))
                    .ToArray();
                rotatedCaches = Enumerable.Range(0, factors.Length)
                    .Select(i => new ProfileCache(dataSource, WindowSizeAt(i), method, OriginRotation.RotatedVersion(methodVersion), CachedLengths(i)))
                    .ToArray();
            }
        }

        private static int WindowSizeAt(int level) => RecordStore.FineWindowSize * ZCurveProfileLadder.DefaultFactors[level];

        public LengthRange CachedLengths(int level) => CachedLengths(level, MinLengthRatio);

        // Which record lengths level i can be asked for by any query with ratio r >= minLengthRatio. A query of
        // length L ranks at window w when r * floor(L / w) >= 8 and not at the next window w' (r * floor(L / w') < 8),
        // so L >= 8w / r and L < w' * (8 / r + 1). Its candidates X are in [rL, L / r]: X >= 8w, X < w' * (8 / r^2 + 1 / r).
        // The finest level is also every query's fallback, so it starts at 0; the coarsest has no upper end.
        public static LengthRange CachedLengths(int level, double minLengthRatio)
        {
            long min = level == 0 ? 0 : 8L * WindowSizeAt(level);
            if (level == ZCurveProfileLadder.DefaultFactors.Count - 1) return new LengthRange(min, null);
            double r = minLengthRatio;
            return new LengthRange(min, (long)Math.Ceiling(WindowSizeAt(level + 1) * (8 / (r * r) + 1 / r)));
        }

        /// <summary>
        /// The level index a query of this length ranks at for this ratio: <see cref="ZCurveProfileLadder.SelectBandFactor"/>
        /// from the window counts the length implies, for callers that have lengths but no query profile (the null builder).
        /// </summary>
        public static int BandLevel(long lengthBp, double lengthRatio)
        {
            long fineWindows = lengthBp / RecordStore.FineWindowSize;
            for (int i = ZCurveProfileLadder.DefaultFactors.Count - 1; i > 0; i--)
                if (lengthRatio * (fineWindows / ZCurveProfileLadder.DefaultFactors[i]) >= ZCurveProfileLadder.DefaultMinWindows)
                    return i;
            return 0;
        }

        /// <summary>The candidate lengths a query of this length meets: [ratio × L, L / ratio].</summary>
        public static LengthRange Band(long lengthBp, double lengthRatio) =>
            new((long)Math.Ceiling(lengthRatio * lengthBp), (long)Math.Floor(lengthBp / lengthRatio) + 1);

        /// <summary>Brings every level loaded so far up to date with the database. A level is first loaded when a search needs it.</summary>
        public async Task RefreshAsync(CancellationToken cancellationToken = default)
        {
            nullLoaded = false;
            graphs.Clear();
            if (caches is null) return;
            for (int i = 0; i < caches.Length; i++)
                if (loaded[i])
                {
                    await caches[i].RefreshAsync(cancellationToken);
                    await rotatedCaches![i].RefreshAsync(cancellationToken);
                }
        }

        /// <summary>Loads every cached level and the null now, so that searches running in parallel afterwards only read them.</summary>
        public async Task LoadAllAsync(CancellationToken cancellationToken = default)
        {
            for (int i = 0; i < factors.Length; i++)
                await LevelAsync(i, CachedLengths(i), cancellationToken);
            await SignificanceAsync(cancellationToken);
        }

        public Task<SimilarityResult> SearchAsync(QueryProfile query, CatalogueFilter? filter = null,
            double lengthRatio = DefaultLengthRatio, int top = 50, CancellationToken cancellationToken = default) =>
            SearchCoreAsync(query, null, filter, lengthRatio, top, useNeighbours: true, cancellationToken);

        /// <summary>
        /// Find similar to a stored record, from its stored profile (<see cref="QueryProfile.FromStoredRecord"/>), leaving out
        /// the records of its own source: its other chromosomes and plasmids aren't relatives, and it would find itself.
        /// </summary>
        public async Task<SimilarityResult> SearchRecordAsync(long sequenceId, CatalogueFilter? filter = null,
            double lengthRatio = DefaultLengthRatio, int top = 50, bool useNeighbours = true, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select source_id, length_bp, kind = 'bacterium' and coalesce(annotations->'contig'->>'topology' = 'circular', false) " +
                "from sequence_record where id = $1");
            command.Parameters.AddWithValue(sequenceId);
            long sourceId, lengthBp;
            bool circularChromosome;
            await using (var reader = await command.ExecuteReaderAsync(cancellationToken))
            {
                if (!await reader.ReadAsync(cancellationToken))
                    throw new ArgumentException($"No record {sequenceId}.", nameof(sequenceId));
                (sourceId, lengthBp, circularChromosome) = (reader.GetInt64(0), reader.GetInt64(1), reader.GetBoolean(2));
            }

            var fine = await ProfileLevels.ReadAsync(dataSource, sequenceId, method, methodVersion, RecordStore.FineWindowSize)
                ?? throw new ArgumentException($"Record {sequenceId} has no stored {RecordStore.FineWindowSize} bp level.", nameof(sequenceId));
            return await SearchCoreAsync(QueryProfile.FromStoredRecord(fine, lengthBp, circularChromosome), sourceId, filter, lengthRatio, top,
                useNeighbours, cancellationToken);
        }

        private async Task<SimilarityResult> SearchCoreAsync(QueryProfile query, long? excludedSourceId, CatalogueFilter? filter,
            double lengthRatio, int top, bool useNeighbours, CancellationToken cancellationToken)
        {
            if (top < 1) throw new ArgumentOutOfRangeException(nameof(top));
            if (lengthRatio < MinLengthRatio || lengthRatio > 1)
                throw new ArgumentOutOfRangeException(nameof(lengthRatio), $"The length ratio must be {MinLengthRatio} to 1.");

            var lengths = Band(query.LengthBp, lengthRatio);
            filter ??= new CatalogueFilter();
            var banded = filter with
            {
                MinLengthBp = Math.Max(filter.MinLengthBp ?? 0, lengths.Min),
                MaxLengthBp = Math.Min(filter.MaxLengthBp ?? long.MaxValue, lengths.MaxOrUnbounded - 1),
            };
            var candidates = await Catalogue.IdsAsync(dataSource, banded, cancellationToken);
            if (excludedSourceId is long excluded)
                candidates.ExceptWith(await SourceRecordsAsync(excluded, cancellationToken));
            var identical = await IdenticalAsync(query, candidates, cancellationToken);
            if (ZCurveProfileLadder.SelectBandFactor(query.Forward, lengthRatio) is not int queryFactor)
                return new SimilarityResult(null, lengths, [], 0, identical, [], FewWindows: true);

            int level = Array.IndexOf(factors, queryFactor);
            var (profiles, rotatedProfiles) = await LevelAsync(level, lengths, cancellationToken);
            var significance = await SignificanceAsync(cancellationToken);
            var forward = query.Forward.Level(queryFactor);
            var queryLevels = QueryLevels(query, queryFactor);

            var comparable = candidates.Where(profiles.ContainsKey).ToArray();
            var notComparable = candidates.Where(id => !profiles.ContainsKey(id)).Order().ToList();
            var hits = new SimilarityHit[comparable.Length];
            Parallel.For(0, comparable.Length, new ParallelOptions { CancellationToken = cancellationToken }, i =>
            {
                var candidate = profiles[comparable[i]];
                var (score, strand, rotated) = Best(queryLevels, candidate, false, (double.NegativeInfinity, Strand.Forward, false));
                if (rotatedProfiles.TryGetValue(comparable[i], out var turned))
                    (score, strand, rotated) = Best(queryLevels, turned, true, (score, strand, rotated));
                var (shorter, offsets) = ScoreNull.Shape(forward.WindowCount, candidate.WindowCount);
                var p = significance?.PValue(forward.WindowSize, shorter, offsets, score);
                hits[i] = new SimilarityHit(comparable[i], score, strand, p, p?.Log10 + Math.Log10(comparable.Length), Rotated: rotated);
            });

            var byScore = hits.OrderByDescending(h => h.Score).ThenBy(h => h.SequenceId).ToList();
            var graph = useNeighbours ? await GraphAsync(lengthRatio, cancellationToken) : null;
            bool neighbourRanked = graph is { Count: > 0 };
            if (neighbourRanked)
                byScore = await ByNeighboursAsync(byScore, graph!, lengthRatio, cancellationToken);

            bool fewWindows = lengthRatio * query.Forward.Fine.WindowCount < ZCurveProfileLadder.DefaultMinWindows;
            return new SimilarityResult(WindowSizeAt(level), lengths, byScore.Take(top).ToList(), hits.Length, identical, notComparable,
                fewWindows, neighbourRanked);
        }

        // The query's forms at one level, in the order ties are broken: forward, reverse, then the rotated pair.
        private static (ZCurveWindowedProfile Level, Strand Strand, bool Rotated)[] QueryLevels(QueryProfile query, int factor)
        {
            var levels = new List<(ZCurveWindowedProfile, Strand, bool)> { (query.Forward.Level(factor), Strand.Forward, false) };
            if (query.Reverse is not null) levels.Add((query.Reverse.Level(factor), Strand.Reverse, false));
            if (query.Rotated is not null) levels.Add((query.Rotated.Level(factor), Strand.Forward, true));
            if (query.RotatedReverse is not null) levels.Add((query.RotatedReverse.Level(factor), Strand.Reverse, true));
            return levels.ToArray();
        }

        // The best of the query's forms against one form of a candidate, kept only where it beats the score so far.
        private static (double Score, Strand Strand, bool Rotated) Best((ZCurveWindowedProfile Level, Strand Strand, bool Rotated)[] queryLevels,
            ZCurveWindowedProfile candidate, bool candidateRotated, (double Score, Strand Strand, bool Rotated) best)
        {
            foreach (var (level, strand, rotated) in queryLevels)
            {
                var score = ZCurveWindowedProfile.Compare(level, candidate);
                if (score > best.Score) best = (score, strand, rotated || candidateRotated);
            }
            return best;
        }

        // The top hits by shared neighbours, then score; the rest stay in score order below them. A hit missing from the
        // graph is searched live the way the graph was built: unfiltered, by its stored record.
        private async Task<List<SimilarityHit>> ByNeighboursAsync(List<SimilarityHit> byScore, IReadOnlyDictionary<long, long[]> graph,
            double lengthRatio, CancellationToken cancellationToken)
        {
            var queryNeighbours = byScore.Take(NeighbourCount).Select(h => h.SequenceId).ToHashSet();
            var head = new List<SimilarityHit>();
            foreach (var hit in byScore.Take(NeighbourRankDepth))
            {
                var own = graph.TryGetValue(hit.SequenceId, out var stored)
                    ? stored
                    : (await SearchRecordAsync(hit.SequenceId, null, lengthRatio, NeighbourCount, useNeighbours: false, cancellationToken))
                        .Ranked.Select(h => h.SequenceId).ToArray();
                head.Add(hit with { SharedNeighbours = own.Count(queryNeighbours.Contains) });
            }

            return head.OrderByDescending(h => h.SharedNeighbours).ThenByDescending(h => h.Score).ThenBy(h => h.SequenceId)
                .Concat(byScore.Skip(NeighbourRankDepth))
                .ToList();
        }

        private async Task<ScoreNull?> SignificanceAsync(CancellationToken cancellationToken)
        {
            if (!nullLoaded)
            {
                scoreNull = await ScoreNullStore.LoadAsync(dataSource, method, methodVersion, ScoreNullStore.Genome, cancellationToken);
                nullLoaded = true;
            }
            return scoreNull;
        }

        private async Task<IReadOnlyDictionary<long, long[]>> GraphAsync(double lengthRatio, CancellationToken cancellationToken)
        {
            if (!graphs.TryGetValue(lengthRatio, out var graph))
                graphs[lengthRatio] = graph = await NeighbourGraph.LoadAsync(dataSource, lengthRatio, method, methodVersion, cancellationToken);
            return graph;
        }

        private async Task<HashSet<long>> SourceRecordsAsync(long sourceId, CancellationToken cancellationToken)
        {
            await using var command = dataSource.CreateCommand("select id from sequence_record where source_id = $1");
            command.Parameters.AddWithValue(sourceId);
            var ids = new HashSet<long>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken)) ids.Add(reader.GetInt64(0));
            return ids;
        }

        // One level of the band as stored, and the rotated versions of the records that have one (empty before build-origins).
        private async Task<(IReadOnlyDictionary<long, ZCurveWindowedProfile> Stored, IReadOnlyDictionary<long, ZCurveWindowedProfile> Rotated)>
            LevelAsync(int level, LengthRange lengths, CancellationToken cancellationToken)
        {
            if (caches is not null)
            {
                if (!loaded[level])
                {
                    await caches[level].RefreshAsync(cancellationToken);
                    await rotatedCaches![level].RefreshAsync(cancellationToken);
                    loaded[level] = true;
                }
                return (caches[level].Profiles, rotatedCaches![level].Profiles);
            }

            async Task<Dictionary<long, ZCurveWindowedProfile>> Stream(string version)
            {
                var profiles = new Dictionary<long, ZCurveWindowedProfile>();
                await foreach (var (_, sequenceId, profile) in ProfileLevels.StreamAfterAsync(
                    dataSource, method, version, WindowSizeAt(level), 0, lengths, cancellationToken))
                    profiles[sequenceId] = profile;
                return profiles;
            }
            return (await Stream(methodVersion), await Stream(OriginRotation.RotatedVersion(methodVersion)));
        }

        private async Task<IReadOnlyList<(long SequenceId, Strand Strand)>> IdenticalAsync(QueryProfile query, HashSet<long> candidates,
            CancellationToken cancellationToken)
        {
            if (query.Sha256 is null) return [];

            await using var command = dataSource.CreateCommand(
                "select id, sequence_sha256 = $1 from sequence_record where sequence_sha256 = any($2) order by id");
            command.Parameters.AddWithValue(query.Sha256);
            command.Parameters.AddWithValue(new[] { query.Sha256, query.ReverseSha256 ?? query.Sha256 });

            var identical = new List<(long SequenceId, Strand Strand)>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                if (candidates.Contains(reader.GetInt64(0)))
                    identical.Add((reader.GetInt64(0), reader.GetBoolean(1) ? Strand.Forward : Strand.Reverse));
            return identical;
        }
    }
}
