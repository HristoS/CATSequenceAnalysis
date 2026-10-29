using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <param name="SamplesPerCell">Scores kept per cell; a full cell takes no more.</param>
    /// <param name="GenomeDraws">Query-candidate pairs drawn for the genome null (and as many for the permuted control, divided by 4).</param>
    /// <param name="FragmentDraws">Fragment-reference pairs drawn for the fragment null.</param>
    /// <param name="FragmentPool">Records whose 500 bp level is read for the fragment null, as fragment sources and references.</param>
    public sealed record ScoreNullOptions(int SamplesPerCell = 10_000, int GenomeDraws = 4_000_000, int FragmentDraws = 2_000_000,
        int FragmentPool = 3_000, int MinCellSamples = 20, int Seed = 20261022);

    /// <summary>
    /// Builds the score_null cells (R5, P3 step S3.6) from stored profiles only; no DNA is needed or read. Every
    /// sample is scored the way a search would score it: the better of the query's two strands at the level the
    /// search would use. Here the reverse strand is a profile-level flip (windows reversed, x and y negated), which
    /// is only approximate when the length isn't a multiple of 500 bp (results/C3_ProfileFlip.csv); for a null
    /// that changes the window phase, not the statistics.
    ///
    /// Unrelated means: different superkingdom (or virus against cellular); for two cellular records, different
    /// phyla; for two viruses, different families. NCBI's current lineage (taxon) where known, else the record's own.
    /// Plasmids are left out (they carry their host's taxonomy and can cross phyla), and so are records with no
    /// name at the deciding rank. Where <see cref="OriginRotation.BuildAsync"/> has stored a record's rotated levels, the genome null
    /// scores both versions of either side and keeps the best, as find-similar does.
    /// </summary>
    public static class ScoreNullBuilder
    {
        private sealed record Candidate(long Id, string Kind, long LengthBp, string? Superkingdom, string? Phylum, string? Family);

        public static async Task<IReadOnlyDictionary<string, IReadOnlyList<ScoreNullCell>>> BuildAsync(NpgsqlDataSource dataSource,
            ScoreNullOptions options, Action<string> log, string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion)
        {
            var records = await RecordsAsync(dataSource, method, methodVersion);
            log($"{records.Length:N0} records (no plasmids) with a 500 bp level");

            var levels = new IReadOnlyDictionary<long, ZCurveWindowedProfile>[ZCurveProfileLadder.DefaultFactors.Count];
            var rotated = new IReadOnlyDictionary<long, ZCurveWindowedProfile>[levels.Length];
            for (int i = 0; i < levels.Length; i++)
            {
                int windowSize = RecordStore.FineWindowSize * ZCurveProfileLadder.DefaultFactors[i];
                var cache = new ProfileCache(dataSource, windowSize, method, methodVersion, SimilaritySearch.CachedLengths(i, 0.25));
                await cache.RefreshAsync();
                levels[i] = cache.Profiles;
                var rotatedCache = new ProfileCache(dataSource, windowSize, method, OriginRotation.RotatedVersion(methodVersion),
                    SimilaritySearch.CachedLengths(i, 0.25));
                await rotatedCache.RefreshAsync();
                rotated[i] = rotatedCache.Profiles;
            }

            var rng = new Random(options.Seed);
            var result = new Dictionary<string, IReadOnlyList<ScoreNullCell>>
            {
                [ScoreNullStore.Genome] = Genome(records, levels, rotated, options, rng, permuted: false, options.GenomeDraws, log),
                [ScoreNullStore.Permuted] = Genome(records, levels, rotated, options, rng, permuted: true, options.GenomeDraws / 4, log),
                [ScoreNullStore.Fragment] = await FragmentAsync(dataSource, records, options, rng, method, methodVersion, log),
            };
            return result;
        }

        private static async Task<Candidate[]> RecordsAsync(NpgsqlDataSource dataSource, string method, string methodVersion)
        {
            // Public records only (D3): no one's private upload shapes the p-values everyone gets.
            await using var command = dataSource.CreateCommand(
                "select r.id, r.kind, r.length_bp, coalesce(t.tax_superkingdom, r.tax_superkingdom), " +
                "coalesce(t.tax_phylum, r.tax_phylum), coalesce(t.tax_family, r.tax_family) " +
                "from sequence_record r join sequence_source s on s.id = r.source_id left join taxon t on t.taxid = r.ncbi_taxid " +
                "where s.owner_id is null and r.kind <> 'plasmid' and exists (select 1 from profile_level p where p.sequence_id = r.id " +
                "and p.method = $1 and p.method_version = $2 and p.window_size = $3) order by r.length_bp, r.id");
            command.Parameters.AddWithValue(method);
            command.Parameters.AddWithValue(methodVersion);
            command.Parameters.AddWithValue(RecordStore.FineWindowSize);

            var records = new List<Candidate>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync())
            {
                string? Text(int i) => reader.IsDBNull(i) ? null : reader.GetString(i);
                records.Add(new Candidate(reader.GetInt64(0), reader.GetString(1), reader.GetInt64(2), Text(3), Text(4), Text(5)));
            }
            return records.ToArray();
        }

        private static bool Unrelated(Candidate a, Candidate b)
        {
            bool aVirus = a.Kind == "virus", bVirus = b.Kind == "virus";
            if (aVirus != bVirus) return true;
            if (aVirus) return a.Family is not null && b.Family is not null && a.Family != b.Family;
            if (a.Superkingdom is not null && b.Superkingdom is not null && a.Superkingdom != b.Superkingdom) return true;
            return a.Phylum is not null && b.Phylum is not null && a.Phylum != b.Phylum;
        }

        // Records are sorted by length, so a band is one index range.
        private static (int From, int To) IndexRange(Candidate[] records, LengthRange lengths)
        {
            int From(long length)
            {
                int low = 0, high = records.Length;
                while (low < high)
                {
                    int middle = (low + high) / 2;
                    if (records[middle].LengthBp < length) low = middle + 1;
                    else high = middle;
                }
                return low;
            }
            return (From(lengths.Min), From(lengths.MaxOrUnbounded));
        }

        // Genome null: a random query, a random ratio in [0.25, 1], the level a search would rank at, and a random
        // unrelated candidate in the band. Permuted control: any other record in the band, its windows shuffled.
        private static IReadOnlyList<ScoreNullCell> Genome(Candidate[] records, IReadOnlyDictionary<long, ZCurveWindowedProfile>[] levels,
            IReadOnlyDictionary<long, ZCurveWindowedProfile>[] rotated, ScoreNullOptions options, Random rng, bool permuted, int draws,
            Action<string> log)
        {
            var planned = new Dictionary<(int, int, int), int>();
            var work = new List<Pair>();
            int noPartner = 0;
            for (int draw = 0; draw < draws; draw++)
            {
                var query = records[rng.Next(records.Length)];
                double ratio = 0.25 + 0.75 * rng.NextDouble();
                int level = SimilaritySearch.BandLevel(query.LengthBp, ratio);
                var (from, to) = IndexRange(records, SimilaritySearch.Band(query.LengthBp, ratio));

                Candidate? partner = null;
                for (int attempt = 0; attempt < 20 && to > from; attempt++)
                {
                    var pick = records[rng.Next(from, to)];
                    if (pick.Id != query.Id && (permuted || Unrelated(query, pick))) { partner = pick; break; }
                }
                if (partner is null || !levels[level].TryGetValue(query.Id, out var q) || !levels[level].TryGetValue(partner.Id, out var c))
                {
                    noPartner++;
                    continue;
                }

                var (shorter, offsets) = ScoreNull.Shape(q.WindowCount, c.WindowCount);
                var cell = (q.WindowSize, ScoreNull.Bucket(shorter), ScoreNull.Bucket(offsets));
                if (planned.GetValueOrDefault(cell) >= options.SamplesPerCell) continue;
                planned[cell] = planned.GetValueOrDefault(cell) + 1;
                work.Add(permuted
                    ? new Pair(cell, q, Permute(c, rng))
                    : new Pair(cell, q, c, rotated[level].GetValueOrDefault(query.Id), rotated[level].GetValueOrDefault(partner.Id)));
            }

            var cells = Score(work, options);
            log($"{(permuted ? "permuted" : "genome")} null: {draws:N0} draws, {work.Count:N0} scored, {noPartner:N0} without a partner, " +
                $"{cells.Count} cells");
            return cells;
        }

        // Fragment null: m consecutive 500 bp windows (m log-uniform in 1..400, i.e. 0.5 to 200 kbp) of one record,
        // against a longer unrelated record, over every offset. Localization (S3.5) refines a coarse hit, which can
        // only find a best at or below this full scan, so p-values from this null are conservative for it.
        private static async Task<IReadOnlyList<ScoreNullCell>> FragmentAsync(NpgsqlDataSource dataSource, Candidate[] records,
            ScoreNullOptions options, Random rng, string method, string methodVersion, Action<string> log)
        {
            var pool = records.OrderBy(_ => rng.Next()).Take(options.FragmentPool).ToArray();
            var windows = new Dictionary<long, IReadOnlyList<(double X, double Y, double Z)>>();
            var profiles = new Dictionary<long, ZCurveWindowedProfile>();
            foreach (var record in pool)
                if (await ProfileLevels.ReadAsync(dataSource, record.Id, method, methodVersion, RecordStore.FineWindowSize) is { } fine)
                {
                    profiles[record.Id] = fine;
                    windows[record.Id] = fine.Windows;
                }
            pool = pool.Where(r => profiles.ContainsKey(r.Id)).ToArray();

            var planned = new Dictionary<(int, int, int), int>();
            var work = new List<Pair>();
            for (int draw = 0; draw < options.FragmentDraws; draw++)
            {
                int length = (int)Math.Round(Math.Pow(2, rng.NextDouble() * Math.Log2(400)));
                var source = pool[rng.Next(pool.Length)];
                var reference = pool[rng.Next(pool.Length)];
                int sourceWindows = windows[source.Id].Count, referenceWindows = windows[reference.Id].Count;
                if (sourceWindows < length || referenceWindows <= length || !Unrelated(source, reference)) continue;

                var (shorter, offsets) = ScoreNull.Shape(length, referenceWindows);
                var cell = (RecordStore.FineWindowSize, ScoreNull.Bucket(shorter), ScoreNull.Bucket(offsets));
                if (planned.GetValueOrDefault(cell) >= options.SamplesPerCell) continue;
                planned[cell] = planned.GetValueOrDefault(cell) + 1;

                int start = rng.Next(sourceWindows - length + 1);
                var fragment = ZCurveWindowedProfile.FromPrecomputedWindows(RecordStore.FineWindowSize,
                    windows[source.Id].Skip(start).Take(length).ToList());
                work.Add(new Pair(cell, fragment, profiles[reference.Id]));
            }

            var cells = Score(work, options);
            log($"fragment null: {pool.Length:N0} records read at 500 bp, {options.FragmentDraws:N0} draws, {work.Count:N0} scored, {cells.Count} cells");
            return cells;
        }

        // One sampled pair; the rotated versions are set only where both a genome-null side and its stored rotation exist.
        private sealed record Pair((int WindowSize, int WindowsBucket, int OffsetsBucket) Cell, ZCurveWindowedProfile Query,
            ZCurveWindowedProfile Candidate, ZCurveWindowedProfile? RotatedQuery = null, ZCurveWindowedProfile? RotatedCandidate = null);

        private static IReadOnlyList<ScoreNullCell> Score(List<Pair> work, ScoreNullOptions options)
        {
            var scores = new double[work.Count];
            Parallel.For(0, work.Count, i =>
            {
                var pair = work[i];
                double best = double.NegativeInfinity;
                foreach (var query in new[] { pair.Query, pair.RotatedQuery })
                    foreach (var candidate in new[] { pair.Candidate, pair.RotatedCandidate })
                        if (query is not null && candidate is not null)
                            best = Math.Max(best, Math.Max(ZCurveWindowedProfile.Compare(query, candidate),
                                ZCurveWindowedProfile.Compare(Flip(query), candidate)));
                scores[i] = best;
            });

            return Enumerable.Range(0, work.Count)
                .GroupBy(i => work[i].Cell)
                .Where(g => g.Count() >= options.MinCellSamples)
                .Select(g => ScoreNullCell.FromScores(g.Key.WindowSize, g.Key.WindowsBucket, g.Key.OffsetsBucket, g.Select(i => scores[i])))
                .OrderBy(c => c.WindowSize).ThenBy(c => c.WindowsBucket).ThenBy(c => c.OffsetsBucket)
                .ToList();
        }

        private static ZCurveWindowedProfile Flip(ZCurveWindowedProfile profile) =>
            ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize,
                profile.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

        private static ZCurveWindowedProfile Permute(ZCurveWindowedProfile profile, Random rng) =>
            ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize, profile.Windows.OrderBy(_ => rng.Next()).ToList());
    }
}
