using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// Rotation for circular chromosomes (S3.2's open item): a file holds a ring cut at an arbitrary point, and two deposits cut
    /// at different points compare like different strains (results/C3_StrandRotation.csv). Every bacterial chromosome has a
    /// common reference point, its replication origin, and the stored 500 bp profile already places it: the minimum of the
    /// cumulative GC skew (<see cref="ZCurveWindowedProfile.GcSkewExtremes"/>; it finds E. coli K-12's oriC and the dnaA start
    /// of 20 bundled files, results/C3_OriginFromSkew.csv). 45% of the catalogue's clear-skew circular chromosomes start more
    /// than 5% of the ring from it (results/S32_OriginCatalogue_Summary.csv).
    ///
    /// <see cref="BuildAsync"/> stores, beside each qualifying record's levels, its 2, 8 and 32 kbp levels re-derived from the fine
    /// profile rotated to the origin, under <see cref="RotatedVersion"/>: the levels at which a record of at least
    /// <see cref="MinLengthBp"/> can be compared in find-similar (a band of ratio 0.25 or more reaches it only from queries scored at
    /// 2 kbp or coarser). No DNA is read, so this needs no re-download; rotating the DNA instead would only help copies of the same
    /// molecule cut elsewhere, and the catalogue holds none (results/S32_RotatedTwins.csv). Find-similar then scores both versions of
    /// a rotated record, and of the query, and keeps the best: over S3.4b's 150 bacterial queries that raised same-genus
    /// precision@10 from 0.131 to 0.169, while rotated versions alone lost on the queries that couldn't be rotated
    /// (results/S32_RotationRanking_Summary.csv). Localization keeps the stored orientation, so its offsets stay in the deposit's coordinates.
    /// </summary>
    public static class OriginRotation
    {
        public const long MinLengthBp = 500_000;
        // A clear skew: origin and terminus on opposite sides of the ring and at least 3 bp of G - C excess per kbp. In the bundled
        // corpus only Deinococcus radiodurans (1.2) falls below it, the one genome whose origin estimate looks wrong.
        public const double MinOriginToTerminusFraction = 0.4;
        public const double MinSkewBpPerKbp = 3;

        public static readonly IReadOnlyList<int> Factors = new[] { 4, 16, 64 };

        public static string RotatedVersion(string methodVersion) => methodVersion + "-origin";

        /// <summary>The fine profile rotated to start at its replication origin; null when it's too short or its skew too weak to place one.</summary>
        public static ZCurveWindowedProfile? ToOrigin(ZCurveWindowedProfile fine, long lengthBp) =>
            OriginWindow(fine, lengthBp) is int origin ? fine.RotateTo(origin) : null;

        /// <summary>
        /// Where a record's stored origin-rotated levels start, in bp of the deposit, recomputed from its stored fine level (exactly the
        /// rotation build-origins stored: results/GenePath_Check.csv); null when it has none.
        /// </summary>
        public static async Task<long?> OriginBpAsync(NpgsqlDataSource dataSource, long sequenceId, long lengthBp,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion)
        {
            await using (var exists = dataSource.CreateCommand("select 1 from profile_level where sequence_id = $1 and method = $2 and method_version = $3 limit 1"))
            {
                exists.Parameters.AddWithValue(sequenceId);
                exists.Parameters.AddWithValue(method);
                exists.Parameters.AddWithValue(RotatedVersion(methodVersion));
                if (await exists.ExecuteScalarAsync() is null) return null;
            }
            var fine = await ProfileLevels.ReadAsync(dataSource, sequenceId, method, methodVersion, RecordStore.FineWindowSize);
            return fine is not null && OriginWindow(fine, lengthBp) is int window ? (long)window * RecordStore.FineWindowSize : null;
        }

        /// <summary>The fine window the rotated levels start at; null when the record is too short or its skew too weak to place one.</summary>
        public static int? OriginWindow(ZCurveWindowedProfile fine, long lengthBp)
        {
            if (lengthBp < MinLengthBp || fine.WindowCount == 0) return null;
            var skew = fine.GcSkewExtremes();
            int count = fine.WindowCount, apart = Math.Abs(skew.OriginWindow - skew.TerminusWindow);
            if ((double)Math.Min(apart, count - apart) / count < MinOriginToTerminusFraction) return null;
            if (skew.AmplitudeBp / (lengthBp / 1000.0) < MinSkewBpPerKbp) return null;
            return skew.OriginWindow;
        }

        /// <summary>
        /// Replaces every <see cref="RotatedVersion"/> level: one per <see cref="Factors"/> for each circular bacterial chromosome
        /// (BV-BRC's contig topology) of at least <see cref="MinLengthBp"/> with a clear skew. Rebuild score_null and the neighbour
        /// graph afterwards, since both are scored the way searches are.
        /// </summary>
        /// <returns>Records rotated, and records examined.</returns>
        public static async Task<(int Rotated, int Examined)> BuildAsync(NpgsqlDataSource dataSource, Action<string> log,
            string method = RecordStore.Method, string methodVersion = RecordStore.MethodVersion, CancellationToken cancellationToken = default)
        {
            var lengths = new Dictionary<long, long>();
            await using (var command = dataSource.CreateCommand(
                "select id, length_bp from sequence_record where kind = 'bacterium' and length_bp >= $1 " +
                "and annotations->'contig'->>'topology' = 'circular'"))
            {
                command.Parameters.AddWithValue(MinLengthBp);
                await using var reader = await command.ExecuteReaderAsync(cancellationToken);
                while (await reader.ReadAsync(cancellationToken)) lengths[reader.GetInt64(0)] = reader.GetInt64(1);
            }

            var rotated = new List<(long SequenceId, ZCurveProfileLadder Ladder)>();
            await foreach (var (_, sequenceId, fine) in ProfileLevels.StreamAfterAsync(dataSource, method, methodVersion,
                RecordStore.FineWindowSize, 0, new LengthRange(MinLengthBp, null), cancellationToken))
                if (lengths.TryGetValue(sequenceId, out var lengthBp) && ToOrigin(fine, lengthBp) is { } turned)
                    rotated.Add((sequenceId, RecordStore.Ladder(turned)));
            log($"{lengths.Count:N0} circular bacterial chromosomes of at least {MinLengthBp:N0} bp, {rotated.Count:N0} with a clear skew");

            var version = RotatedVersion(methodVersion);
            await using var connection = await dataSource.OpenConnectionAsync(cancellationToken);
            await using var transaction = await connection.BeginTransactionAsync(cancellationToken);
            await using (var delete = new NpgsqlCommand("delete from profile_level where method = $1 and method_version = $2", connection, transaction))
            {
                delete.Parameters.AddWithValue(method);
                delete.Parameters.AddWithValue(version);
                await delete.ExecuteNonQueryAsync(cancellationToken);
            }
            foreach (var (sequenceId, ladder) in rotated)
                foreach (var factor in Factors)
                    if (ladder.Level(factor).WindowCount > 0)
                        await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, version, factor, ladder.Level(factor));
            await transaction.CommitAsync(cancellationToken);
            return (rotated.Count, lengths.Count);
        }
    }
}
