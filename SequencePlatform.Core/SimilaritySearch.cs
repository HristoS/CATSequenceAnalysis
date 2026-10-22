using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    public enum Strand { Forward, Reverse }

    /// <param name="P">The chance an unrelated record of this shape scores as high (score_null, S3.6), and how it was obtained; null without a null.</param>
    /// <param name="Log10E">log10 of the expected number of unrelated records scoring as high among this search's candidates: p × candidates.</param>
    public sealed record SimilarityHit(long SequenceId, double Score, Strand Strand, PValue? P = null, double? Log10E = null);

    /// <param name="QueryWindowSize">The one level every hit is scored at; null when the query has no full fine window.</param>
    /// <param name="Lengths">The candidate length band, [ratio × query length, query length / ratio].</param>
    /// <param name="Identical">Records whose sequence_sha256 equals the query's on either strand: sequence identity, not a score.</param>
    /// <param name="NotComparable">Candidates in the band with no stored level at the query's window size.</param>
    /// <param name="FewWindows">
    /// The query is too short for its shortest candidates to keep 8 windows even at 500 bp (under about 8 kbp at ratio 0.5).
    /// The ranking then reflects composition more than order: HPV16 (15 windows) ranks GC-matched strangers above its own
    /// family (results/S36_CatalogueRanking.csv).
    /// </param>
    public sealed record SimilarityResult(
        int? QueryWindowSize,
        LengthRange Lengths,
        IReadOnlyList<SimilarityHit> Ranked, int RankedTotal,
        IReadOnlyList<(long SequenceId, Strand Strand)> Identical,
        IReadOnlyList<long> NotComparable,
        bool FewWindows);

    /// <summary>
    /// A query as the search sees it: its ladder on both strands, its length and, when it came from DNA, the
    /// SHA-256 of both strands normalized as <see cref="SequenceFeatureReader"/> does for stored records. The
    /// reverse strand is profiled from the reverse-complement DNA, which is exact; a stored profile is never
    /// flipped (P3 S3.2, results/C3_ProfileFlip.csv). The DNA itself isn't kept.
    /// </summary>
    /// <param name="LengthBp">As length_bp is for stored records; for a profile-only query, windows × window size (a lower bound).</param>
    public sealed record QueryProfile(ZCurveProfileLadder Forward, ZCurveProfileLadder? Reverse, long LengthBp, byte[]? Sha256, byte[]? ReverseSha256)
    {
        public static QueryProfile FromSequence(string sequence)
        {
            var (forward, sha256, lengthBp) = Profile(sequence);
            var (reverse, reverseSha256, _) = Profile(ReverseComplement(sequence));
            return new QueryProfile(forward, reverse, lengthBp, sha256, reverseSha256);
        }

        /// <summary>A profile-only query (a ZCWP upload): no hash, and the reverse strand only if the uploader supplied it.</summary>
        public static QueryProfile FromProfiles(ZCurveWindowedProfile forward, ZCurveWindowedProfile? reverse) =>
            new(new ZCurveProfileLadder(forward), reverse is null ? null : new ZCurveProfileLadder(reverse),
                (long)forward.WindowCount * forward.WindowSize, null, null);

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
    /// </summary>
    public sealed class SimilaritySearch
    {
        public const double DefaultLengthRatio = 0.5;

        private readonly NpgsqlDataSource dataSource;
        private readonly string method;
        private readonly string methodVersion;
        private readonly int[] factors = ZCurveProfileLadder.DefaultFactors.ToArray();
        private readonly ProfileCache[]? caches;
        private readonly bool[] loaded;
        private ScoreNull? scoreNull;
        private bool nullLoaded;

        /// <summary>The smallest length ratio a search may ask for; the cached length ranges are sized for it.</summary>
        public double MinLengthRatio { get; }

        /// <summary>Bytes of window data held over every loaded level; 0 without a cache.</summary>
        public long WindowBytes => caches is null ? 0 : caches.Where((_, i) => loaded[i]).Sum(c => c.WindowBytes);

        /// <param name="useCache">false: every search reads its band from the database (one query per process, as the CLI does).</param>
        public SimilaritySearch(NpgsqlDataSource dataSource, bool useCache = true, double minLengthRatio = 0.25,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion)
        {
            if (minLengthRatio <= 0 || minLengthRatio > 1) throw new ArgumentOutOfRangeException(nameof(minLengthRatio));
            this.dataSource = dataSource;
            this.method = method;
            this.methodVersion = methodVersion;
            MinLengthRatio = minLengthRatio;
            loaded = new bool[factors.Length];
            if (useCache)
                caches = Enumerable.Range(0, factors.Length)
                    .Select(i => new ProfileCache(dataSource, WindowSizeAt(i), method, methodVersion, CachedLengths(i)))
                    .ToArray();
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
            if (caches is null) return;
            for (int i = 0; i < caches.Length; i++)
                if (loaded[i]) await caches[i].RefreshAsync(cancellationToken);
        }

        public async Task<SimilarityResult> SearchAsync(QueryProfile query, CatalogueFilter? filter = null,
            double lengthRatio = DefaultLengthRatio, int top = 50, CancellationToken cancellationToken = default)
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
            var identical = await IdenticalAsync(query, candidates, cancellationToken);
            if (ZCurveProfileLadder.SelectBandFactor(query.Forward, lengthRatio) is not int queryFactor)
                return new SimilarityResult(null, lengths, [], 0, identical, [], FewWindows: true);

            int level = Array.IndexOf(factors, queryFactor);
            var profiles = await LevelAsync(level, lengths, cancellationToken);
            if (!nullLoaded)
            {
                scoreNull = await ScoreNullStore.LoadAsync(dataSource, method, methodVersion, ScoreNullStore.Genome, cancellationToken);
                nullLoaded = true;
            }
            var significance = scoreNull;
            var forward = query.Forward.Level(queryFactor);
            var reverse = query.Reverse?.Level(queryFactor);

            var comparable = candidates.Where(profiles.ContainsKey).ToArray();
            var notComparable = candidates.Where(id => !profiles.ContainsKey(id)).Order().ToList();
            var hits = new SimilarityHit[comparable.Length];
            Parallel.For(0, comparable.Length, new ParallelOptions { CancellationToken = cancellationToken }, i =>
            {
                var candidate = profiles[comparable[i]];
                var score = ZCurveWindowedProfile.Compare(forward, candidate);
                var strand = Strand.Forward;
                if (reverse is not null)
                {
                    var reversed = ZCurveWindowedProfile.Compare(reverse, candidate);
                    if (reversed > score) (score, strand) = (reversed, Strand.Reverse);
                }
                var (shorter, offsets) = ScoreNull.Shape(forward.WindowCount, candidate.WindowCount);
                var p = significance?.PValue(forward.WindowSize, shorter, offsets, score);
                hits[i] = new SimilarityHit(comparable[i], score, strand, p, p?.Log10 + Math.Log10(comparable.Length));
            });

            var ranked = hits.OrderByDescending(h => h.Score).ThenBy(h => h.SequenceId).Take(top).ToList();
            bool fewWindows = lengthRatio * query.Forward.Fine.WindowCount < ZCurveProfileLadder.DefaultMinWindows;
            return new SimilarityResult(WindowSizeAt(level), lengths, ranked, hits.Length, identical, notComparable, fewWindows);
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
