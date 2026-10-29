using System.Text.Json;
using System.Text.Json.Serialization;
using Npgsql;
using NpgsqlTypes;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    public enum SearchKind { Similar, Localize }

    public enum SearchStatus { Queued, Running, Done, Failed, Expired }

    public enum HitList { Ranked, Identical, NotComparable }

    /// <param name="Filter">Its <see cref="CatalogueFilter.ViewerId"/> is the run's owner, kept in search_run.owner_id rather than in the filters.</param>
    /// <param name="LengthRatio">Find-similar's band ratio; null for localize.</param>
    /// <param name="WithinSequenceId">Localize inside this one record (<see cref="Localization.LocalizeInRecordAsync"/>).</param>
    public sealed record SearchRequest(SearchKind Kind, CatalogueFilter Filter, int Top, double? LengthRatio = null, long? WithinSequenceId = null);

    /// <param name="WindowSize">Find-similar: the level every hit is scored at; localize: the coarse scan level. Null until done.</param>
    /// <param name="Lengths">Find-similar's candidate length band.</param>
    /// <param name="Candidates">Localize: the records scanned.</param>
    /// <param name="WeakSignal">Find-similar's <see cref="SimilarityResult.FewWindows"/> or localize's <see cref="LocalizationResult.WeakSignal"/>.</param>
    public sealed record SearchRun(Guid PublicId, SearchRequest Request, SearchStatus Status, string? Error, long QueryLengthBp,
        DateTimeOffset CreatedAt, DateTimeOffset? StartedAt, DateTimeOffset? CompletedAt, DateTimeOffset ExpiresAt,
        int? WindowSize, LengthRange? Lengths, int? RankedTotal, int? Candidates, bool? WeakSignal, bool? NeighbourRanked);

    /// <summary>One row of search_hit: a <see cref="SimilarityHit"/>, a <see cref="LocalizationHit"/>, an identical record or a not-comparable one.</summary>
    /// <param name="Position">1-based order within its list; for <see cref="HitList.Ranked"/>, the rank.</param>
    public sealed record StoredHit(HitList List, int Position, long SequenceId, Strand? Strand = null, double? Score = null,
        double? CoarseScore = null, long? OffsetBp = null, PValue? P = null, double? Log10E = null, int? SharedNeighbours = null,
        bool Rotated = false);

    /// <summary>
    /// search_run and search_hit I/O (migrations 0002 and 0009, P4 S4.1): a search as an asynchronous job. A job keeps what
    /// <see cref="QueryProfile.FromFineProfiles"/> needs to rebuild the query, never the DNA; its result is the hit lists of a
    /// <see cref="SimilarityResult"/> or <see cref="LocalizationResult"/> plus their summary. D7: <see cref="ExpireDueAsync"/>
    /// drops profiles and hits after 30 days and keeps the hashes, length, parameters and summary; <see cref="ExpireAsync"/>
    /// does the same early for one run.
    ///
    /// PostgreSQL is the queue (S4.3): <see cref="EnqueueAsync"/> notifies <see cref="Channel"/>, <see cref="ClaimAsync"/> takes the
    /// oldest queued job with <c>for update skip locked</c>, so two workers never take the same one, and <see cref="RequeueStuckAsync"/>
    /// returns jobs whose worker stopped mid-run.
    /// </summary>
    public static class SearchLog
    {
        /// <summary>The LISTEN/NOTIFY channel a new job is announced on; the payload is its public id.</summary>
        public const string Channel = "search_queued";
        public const int MaxAttempts = 3;

        private static readonly JsonSerializerOptions filterJson = new() { DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull };

        public static async Task<Guid> EnqueueAsync(NpgsqlDataSource dataSource, SearchRequest request, QueryProfile query,
            CancellationToken cancellationToken = default)
        {
            if (request.Top < 1) throw new ArgumentOutOfRangeException(nameof(request), "Top must be at least 1.");
            if (request.Kind == SearchKind.Similar && request.WithinSequenceId is not null)
                throw new ArgumentException("Only a localization runs inside one record.", nameof(request));
            if (request.Kind == SearchKind.Localize && request.LengthRatio is not null)
                throw new ArgumentException("A length ratio is find-similar's.", nameof(request));

            await using var command = dataSource.CreateCommand(
                "insert into search_run (search_kind, query_sha256, query_reverse_sha256, query_length_bp, query_profile, query_reverse_profile, " +
                "filters, top_k, length_ratio, within_sequence_id, owner_id) values ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11) returning public_id");
            var p = command.Parameters;
            p.AddWithValue(KindName(request.Kind));
            p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bytea, Value = (object?)query.Sha256 ?? DBNull.Value });
            p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bytea, Value = (object?)query.ReverseSha256 ?? DBNull.Value });
            p.AddWithValue(query.LengthBp);
            p.AddWithValue(query.Forward.Fine.ToBinary());
            p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bytea, Value = (object?)query.Reverse?.Fine.ToBinary() ?? DBNull.Value });
            p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Jsonb, Value = JsonSerializer.Serialize(request.Filter with { ViewerId = null }, filterJson) });
            p.AddWithValue(request.Top);
            p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Double, Value = (object?)request.LengthRatio ?? DBNull.Value });
            p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bigint, Value = (object?)request.WithinSequenceId ?? DBNull.Value });
            p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bigint, Value = (object?)request.Filter.ViewerId ?? DBNull.Value });
            var id = (Guid)(await command.ExecuteScalarAsync(cancellationToken))!;

            await using (var notify = dataSource.CreateCommand("select pg_notify($1, $2)"))
            {
                notify.Parameters.AddWithValue(Channel);
                notify.Parameters.AddWithValue(id.ToString());
                await notify.ExecuteNonQueryAsync(cancellationToken);
            }
            return id;
        }

        /// <summary>Marks the oldest queued job running and counts the attempt.</summary>
        /// <returns>Its public id; null when nothing is queued.</returns>
        public static async Task<Guid?> ClaimAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "update search_run set status = 'running', started_at = now(), attempts = attempts + 1 " +
                "where id = (select id from search_run where status = 'queued' order by created_at, id limit 1 for update skip locked) " +
                "returning public_id");
            return await command.ExecuteScalarAsync(cancellationToken) as Guid?;
        }

        /// <summary>
        /// Jobs running since before <paramref name="startedBefore"/>, whose worker stopped: back to queued, or failed once they have
        /// been claimed <see cref="MaxAttempts"/> times.
        /// </summary>
        public static async Task<(int Requeued, int Failed)> RequeueStuckAsync(NpgsqlDataSource dataSource, DateTimeOffset startedBefore,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "with stuck as (select id, attempts from search_run where status = 'running' and started_at < $1 for update skip locked), " +
                "failed as (update search_run r set status = 'failed', completed_at = now(), " +
                "error = 'The search stopped its worker ' || s.attempts || ' times; it was not retried.' " +
                "from stuck s where r.id = s.id and s.attempts >= $2 returning r.id), " +
                "requeued as (update search_run r set status = 'queued', started_at = null from stuck s where r.id = s.id and s.attempts < $2 returning r.id) " +
                "select (select count(*)::integer from requeued), (select count(*)::integer from failed)");
            command.Parameters.AddWithValue(startedBefore);
            command.Parameters.AddWithValue((short)MaxAttempts);
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            await reader.ReadAsync(cancellationToken);
            return (reader.GetInt32(0), reader.GetInt32(1));
        }

        /// <returns>null when no run has this id.</returns>
        public static async Task<SearchRun?> GetAsync(NpgsqlDataSource dataSource, Guid publicId, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select search_kind, filters::text, top_k, length_ratio, within_sequence_id, status, error, query_length_bp, " +
                "created_at, started_at, completed_at, expires_at, level_window_size, band_min_bp, band_max_bp, ranked_total, " +
                "candidates, weak_signal, neighbour_ranked, owner_id from search_run where public_id = $1");
            command.Parameters.AddWithValue(publicId);
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            if (!await reader.ReadAsync(cancellationToken)) return null;

            var request = new SearchRequest(ParseKind(reader.GetString(0)),
                JsonSerializer.Deserialize<CatalogueFilter>(reader.GetString(1), filterJson)! with { ViewerId = Nullable<long>(reader, 19) },
                reader.GetInt32(2),
                Nullable<double>(reader, 3), Nullable<long>(reader, 4));
            return new SearchRun(publicId, request, Enum.Parse<SearchStatus>(reader.GetString(5), ignoreCase: true),
                reader.IsDBNull(6) ? null : reader.GetString(6), reader.GetInt64(7),
                reader.GetFieldValue<DateTimeOffset>(8), Nullable<DateTimeOffset>(reader, 9), Nullable<DateTimeOffset>(reader, 10),
                reader.GetFieldValue<DateTimeOffset>(11), Nullable<int>(reader, 12),
                reader.IsDBNull(13) ? null : new LengthRange(reader.GetInt64(13), Nullable<long>(reader, 14)),
                Nullable<int>(reader, 15), Nullable<int>(reader, 16), Nullable<bool>(reader, 17), Nullable<bool>(reader, 18));
        }

        /// <returns>The job's query, rebuilt exactly; null when no run has this id or it has expired.</returns>
        public static async Task<QueryProfile?> QueryAsync(NpgsqlDataSource dataSource, Guid publicId, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select query_profile, query_reverse_profile, query_length_bp, query_sha256, query_reverse_sha256 " +
                "from search_run where public_id = $1 and query_profile is not null");
            command.Parameters.AddWithValue(publicId);
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            if (!await reader.ReadAsync(cancellationToken)) return null;

            return QueryProfile.FromFineProfiles(ZCurveWindowedProfile.FromBinary(reader.GetFieldValue<byte[]>(0)),
                reader.IsDBNull(1) ? null : ZCurveWindowedProfile.FromBinary(reader.GetFieldValue<byte[]>(1)), reader.GetInt64(2),
                reader.IsDBNull(3) ? null : reader.GetFieldValue<byte[]>(3), reader.IsDBNull(4) ? null : reader.GetFieldValue<byte[]>(4));
        }

        /// <returns>false when the run is no longer queued or running (expired or deleted meanwhile): nothing is written.</returns>
        public static Task<bool> SaveAsync(NpgsqlDataSource dataSource, Guid publicId, SimilarityResult result,
            CancellationToken cancellationToken = default)
        {
            var hits = new List<StoredHit>();
            hits.AddRange(result.Ranked.Select((h, i) => new StoredHit(HitList.Ranked, i + 1, h.SequenceId, h.Strand, h.Score,
                P: h.P, Log10E: h.Log10E, SharedNeighbours: h.SharedNeighbours, Rotated: h.Rotated)));
            hits.AddRange(result.Identical.Select((h, i) => new StoredHit(HitList.Identical, i + 1, h.SequenceId, h.Strand)));
            hits.AddRange(result.NotComparable.Select((id, i) => new StoredHit(HitList.NotComparable, i + 1, id)));
            return CompleteAsync(dataSource, publicId, result.QueryWindowSize, result.Lengths, result.RankedTotal, null,
                result.FewWindows, result.NeighbourRanked, hits, cancellationToken);
        }

        /// <returns>false when the run is no longer queued or running (expired or deleted meanwhile): nothing is written.</returns>
        public static Task<bool> SaveAsync(NpgsqlDataSource dataSource, Guid publicId, LocalizationResult result,
            CancellationToken cancellationToken = default)
        {
            var hits = new List<StoredHit>();
            hits.AddRange(result.Hits.Select((h, i) => new StoredHit(HitList.Ranked, i + 1, h.SequenceId, h.Strand, h.Score,
                h.CoarseScore, h.OffsetBp, h.P, h.Log10E)));
            hits.AddRange(result.NotComparable.Select((id, i) => new StoredHit(HitList.NotComparable, i + 1, id)));
            return CompleteAsync(dataSource, publicId, result.CoarseWindowSize, null, null, result.Candidates,
                result.WeakSignal, null, hits, cancellationToken);
        }

        private static async Task<bool> CompleteAsync(NpgsqlDataSource dataSource, Guid publicId, int? windowSize, LengthRange? lengths,
            int? rankedTotal, int? candidates, bool weakSignal, bool? neighbourRanked, List<StoredHit> hits, CancellationToken cancellationToken)
        {
            await using var connection = await dataSource.OpenConnectionAsync(cancellationToken);
            await using var transaction = await connection.BeginTransactionAsync(cancellationToken);
            long runId;
            await using (var update = new NpgsqlCommand(
                "update search_run set status = 'done', completed_at = now(), level_window_size = $2, band_min_bp = $3, band_max_bp = $4, " +
                "ranked_total = $5, candidates = $6, weak_signal = $7, neighbour_ranked = $8 " +
                "where public_id = $1 and status in ('queued','running') returning id", connection, transaction))
            {
                var p = update.Parameters;
                p.AddWithValue(publicId);
                p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Integer, Value = (object?)windowSize ?? DBNull.Value });
                p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bigint, Value = (object?)lengths?.Min ?? DBNull.Value });
                p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bigint, Value = (object?)lengths?.Max ?? DBNull.Value });
                p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Integer, Value = (object?)rankedTotal ?? DBNull.Value });
                p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Integer, Value = (object?)candidates ?? DBNull.Value });
                p.AddWithValue(weakSignal);
                p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Boolean, Value = (object?)neighbourRanked ?? DBNull.Value });
                if (await update.ExecuteScalarAsync(cancellationToken) is not long id) return false;
                runId = id;
            }

            // One statement for the whole result: unnest zips the column arrays into rows.
            await using (var insert = new NpgsqlCommand(
                "insert into search_hit (search_run_id, list_kind, position, sequence_id, strand, score, coarse_score, offset_bp, " +
                "log10_p, p_kind, log10_e, shared_neighbours, rotated) " +
                "select $1, * from unnest($2::text[], $3::integer[], $4::bigint[], $5::text[], $6::float8[], $7::float8[], $8::bigint[], " +
                "$9::float8[], $10::text[], $11::float8[], $12::smallint[], $13::boolean[])", connection, transaction))
            {
                var p = insert.Parameters;
                p.AddWithValue(runId);
                p.AddWithValue(hits.Select(h => ListName(h.List)).ToArray());
                p.AddWithValue(hits.Select(h => h.Position).ToArray());
                p.AddWithValue(hits.Select(h => h.SequenceId).ToArray());
                p.AddWithValue(hits.Select(h => h.Strand is { } s ? StrandName(s) : null).ToArray());
                p.AddWithValue(hits.Select(h => h.Score).ToArray());
                p.AddWithValue(hits.Select(h => h.CoarseScore).ToArray());
                p.AddWithValue(hits.Select(h => h.OffsetBp).ToArray());
                p.AddWithValue(hits.Select(h => h.P?.Log10).ToArray());
                p.AddWithValue(hits.Select(h => h.P is { } v ? v.Kind.ToString().ToLowerInvariant() : null).ToArray());
                p.AddWithValue(hits.Select(h => h.Log10E).ToArray());
                p.AddWithValue(hits.Select(h => (short?)h.SharedNeighbours).ToArray());
                p.AddWithValue(hits.Select(h => h.Rotated).ToArray());
                await insert.ExecuteNonQueryAsync(cancellationToken);
            }
            await transaction.CommitAsync(cancellationToken);
            return true;
        }

        /// <returns>false when the run is no longer queued or running.</returns>
        public static async Task<bool> FailAsync(NpgsqlDataSource dataSource, Guid publicId, string error, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "update search_run set status = 'failed', error = $2, completed_at = now() where public_id = $1 and status in ('queued','running')");
            command.Parameters.AddWithValue(publicId);
            command.Parameters.AddWithValue(error);
            return await command.ExecuteNonQueryAsync(cancellationToken) == 1;
        }

        /// <summary>One page of one hit list, keyset by position: the rows after <paramref name="afterPosition"/>.</summary>
        public static async Task<IReadOnlyList<StoredHit>> HitsAsync(NpgsqlDataSource dataSource, Guid publicId, HitList list,
            int afterPosition = 0, int limit = 50, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(
                "select position, sequence_id, strand, score, coarse_score, offset_bp, log10_p, p_kind, log10_e, shared_neighbours, rotated " +
                "from search_hit where search_run_id = (select id from search_run where public_id = $1) and list_kind = $2 and position > $3 " +
                "order by position limit $4");
            var p = command.Parameters;
            p.AddWithValue(publicId);
            p.AddWithValue(ListName(list));
            p.AddWithValue(afterPosition);
            p.AddWithValue(limit);

            var hits = new List<StoredHit>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                hits.Add(new StoredHit(list, reader.GetInt32(0), reader.GetInt64(1),
                    reader.IsDBNull(2) ? null : Enum.Parse<Strand>(reader.GetString(2), ignoreCase: true),
                    Nullable<double>(reader, 3), Nullable<double>(reader, 4), Nullable<long>(reader, 5),
                    reader.IsDBNull(6) ? null : new PValue(reader.GetDouble(6), Enum.Parse<PValueKind>(reader.GetString(7), ignoreCase: true)),
                    Nullable<double>(reader, 8), reader.IsDBNull(9) ? null : reader.GetInt16(9), reader.GetBoolean(10)));
            return hits;
        }

        /// <summary>DELETE of a search (D7 early): the run keeps its hashes, length, parameters and summary, nothing else.</summary>
        /// <returns>false when no run has this id or it had already expired.</returns>
        public static async Task<bool> ExpireAsync(NpgsqlDataSource dataSource, Guid publicId, CancellationToken cancellationToken = default) =>
            await ExpireWhereAsync(dataSource, "public_id = $1", publicId, cancellationToken) == 1;

        /// <returns>The runs expired: every one past its expires_at, whatever its status.</returns>
        public static Task<int> ExpireDueAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken = default) =>
            ExpireWhereAsync(dataSource, "expires_at <= now()", null, cancellationToken);

        // The hits are deleted in the same statement: a data-modifying CTE sees the run ids the update returned.
        private static async Task<int> ExpireWhereAsync(NpgsqlDataSource dataSource, string condition, Guid? publicId,
            CancellationToken cancellationToken)
        {
            await using var command = dataSource.CreateCommand(
                "with expired as (update search_run set status = 'expired', query_profile = null, query_reverse_profile = null " +
                $"where {condition} and status <> 'expired' returning id), " +
                "hits as (delete from search_hit h using expired e where h.search_run_id = e.id) " +
                "select count(*)::integer from expired");
            if (publicId is { } id) command.Parameters.AddWithValue(id);
            return (int)(await command.ExecuteScalarAsync(cancellationToken))!;
        }

        private static string KindName(SearchKind kind) => kind == SearchKind.Similar ? "similar" : "localize";

        private static SearchKind ParseKind(string name) => name == "similar" ? SearchKind.Similar : SearchKind.Localize;

        private static string ListName(HitList list) => list switch
        {
            HitList.Ranked => "ranked",
            HitList.Identical => "identical",
            _ => "not_comparable",
        };

        private static string StrandName(Strand strand) => strand == Strand.Forward ? "forward" : "reverse";

        private static T? Nullable<T>(NpgsqlDataReader reader, int ordinal) where T : struct =>
            reader.IsDBNull(ordinal) ? null : reader.GetFieldValue<T>(ordinal);
    }
}
