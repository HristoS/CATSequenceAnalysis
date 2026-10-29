using System.Text;
using Npgsql;

namespace SequencePlatform.Core
{
    /// <param name="ViewerId">Who is asking (D3): every filter sees the public catalogue, and a viewer also sees their own uploads.
    /// Null sees public records only, so a caller that forgets the viewer can't leak anyone's upload.</param>
    public sealed record CatalogueFilter(
        string? Kind = null,
        string? Rank = null,
        string? Taxon = null,
        long? MinLengthBp = null,
        long? MaxLengthBp = null,
        string? Provider = null,
        string? NameContains = null,
        long? ViewerId = null);

    public sealed record Lineage(string? Superkingdom, string? Realm, string? Kingdom, string? Phylum, string? Class,
        string? Order, string? Family, string? Genus, string? Species);

    public sealed record CatalogueRecord(long Id, string Provider, string? SourceAccession, string RecordId, string? Name,
        string Kind, long LengthBp, float? GcFraction, int? NcbiTaxid, Lineage Lineage);

    /// <param name="Release">The GTDB release the lineage comes from.</param>
    public sealed record GtdbLineage(string Release, string Accession, Lineage Lineage);

    public sealed record LevelInfo(string MethodVersion, int WindowSize, int WindowCount, long ByteSize);

    /// <summary>One record with everything known about it except its DNA, which isn't stored.</summary>
    /// <param name="Owned">The record is the viewer's own upload.</param>
    /// <param name="NcbiLineage">NCBI's current lineage of the record's taxid (taxon, S2.4); <see cref="CatalogueRecord.Lineage"/> is the provider's.</param>
    public sealed record CatalogueDetail(CatalogueRecord Record, string DataClass, string AnnotationsJson, byte[]? SequenceSha256,
        string? SourceUrl, DateTimeOffset RetrievedAt, byte[] RawSha256, long RawByteSize, string? License, bool Owned,
        string? NcbiName, Lineage? NcbiLineage, GtdbLineage? Gtdb, IReadOnlyList<LevelInfo> Levels);

    /// <param name="NextAfter">The <c>after</c> value for the next page; null on the last page.</param>
    public sealed record CataloguePage(IReadOnlyList<CatalogueRecord> Records, long? NextAfter);

    /// <summary>
    /// Catalogue listing shared by the Ingest's <c>list</c> command and the P4 API. Pages are keyset
    /// (<c>id &gt; after</c>, ordered by id), not OFFSET, so a page costs the same at any depth and
    /// rows inserted while paging can't shift or repeat the rows already returned.
    /// </summary>
    public static class Catalogue
    {
        public const int MaxPageSize = 1000;

        // The only way a rank name reaches SQL: a fixed map to its column.
        private static readonly Dictionary<string, string> rankColumns = new(StringComparer.OrdinalIgnoreCase)
        {
            ["superkingdom"] = "tax_superkingdom",
            ["realm"] = "tax_realm",
            ["kingdom"] = "tax_kingdom",
            ["phylum"] = "tax_phylum",
            ["class"] = "tax_class",
            ["order"] = "tax_order",
            ["family"] = "tax_family",
            ["genus"] = "tax_genus",
            ["species"] = "tax_species",
        };

        public static IReadOnlyCollection<string> Ranks => rankColumns.Keys;

        private const string recordSelect =
            "select r.id, s.provider, s.accession, r.record_id, r.name, r.kind, r.length_bp, r.gc_fraction, r.ncbi_taxid, " +
            "r.tax_superkingdom, r.tax_realm, r.tax_kingdom, r.tax_phylum, r.tax_class, r.tax_order, r.tax_family, r.tax_genus, r.tax_species";
        private const string recordFrom = " from sequence_record r join sequence_source s on s.id = r.source_id";
        private const string recordColumns = recordSelect + recordFrom;
        private const int recordColumnCount = 18;

        public static async Task<CataloguePage> ListAsync(NpgsqlDataSource dataSource, CatalogueFilter filter, long after = 0, int limit = 50,
            CancellationToken cancellationToken = default)
        {
            if (limit < 1 || limit > MaxPageSize)
                throw new ArgumentOutOfRangeException(nameof(limit), $"Page size must be 1 to {MaxPageSize}.");

            await using var command = dataSource.CreateCommand();
            var sql = new StringBuilder(recordColumns).Append(" where r.id > $1");
            command.Parameters.AddWithValue(after);
            AppendFilter(command, sql, filter);

            // One row past the page tells whether another page exists without a count query.
            sql.Append(" order by r.id limit ").Append(limit + 1);
            command.CommandText = sql.ToString();

            var records = await ReadRecordsAsync(command, limit + 1, cancellationToken);
            if (records.Count <= limit) return new CataloguePage(records, null);
            records.RemoveAt(limit);
            return new CataloguePage(records, records[^1].Id);
        }

        /// <summary>Every record id the filter matches, unpaged: a search's candidate set (P3 S3.4).</summary>
        public static async Task<HashSet<long>> IdsAsync(NpgsqlDataSource dataSource, CatalogueFilter filter, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand();
            var sql = new StringBuilder("select r.id from sequence_record r join sequence_source s on s.id = r.source_id where true");
            AppendFilter(command, sql, filter);
            command.CommandText = sql.ToString();

            var ids = new HashSet<long>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                ids.Add(reader.GetInt64(0));
            return ids;
        }

        /// <summary>
        /// The records with these ids, keyed by id; ids that don't exist are left out. No visibility check: this is for ids that a
        /// filtered listing or search already returned.
        /// </summary>
        public static async Task<IReadOnlyDictionary<long, CatalogueRecord>> ByIdsAsync(NpgsqlDataSource dataSource, IReadOnlyCollection<long> ids,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(recordColumns + " where r.id = any($1)");
            command.Parameters.AddWithValue(ids.ToArray());
            var records = await ReadRecordsAsync(command, ids.Count, cancellationToken);
            return records.ToDictionary(r => r.Id);
        }

        /// <returns>null when the record doesn't exist or the viewer may not see it, which a caller can't tell apart.</returns>
        public static async Task<CatalogueDetail?> DetailAsync(NpgsqlDataSource dataSource, long id, long? viewerId,
            CancellationToken cancellationToken = default)
        {
            CatalogueDetail? detail;
            await using (var command = dataSource.CreateCommand())
            {
                // GTDB is keyed by the assembly accession without its GCA_/GCF_ prefix (TaxonomyLoader.AssemblyKey).
                var sql = new StringBuilder(recordSelect).Append(
                    ", r.data_class, r.annotations::text, r.sequence_sha256, s.url, s.retrieved_at, s.raw_sha256, s.raw_byte_size, s.license, " +
                    "s.owner_id is not null, t.name, t.tax_superkingdom, t.tax_realm, t.tax_kingdom, t.tax_phylum, t.tax_class, t.tax_order, " +
                    "t.tax_family, t.tax_genus, t.tax_species, g.release, g.gtdb_accession, g.tax_domain, g.tax_phylum, g.tax_class, " +
                    "g.tax_order, g.tax_family, g.tax_genus, g.tax_species")
                    .Append(recordFrom)
                    .Append(" left join taxon t on t.taxid = r.ncbi_taxid" +
                        " left join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_'" +
                        " and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5)" +
                        " where r.id = $1");
                command.Parameters.AddWithValue(id);
                AppendFilter(command, sql, new CatalogueFilter(ViewerId: viewerId));
                command.CommandText = sql.ToString();

                await using var reader = await command.ExecuteReaderAsync(cancellationToken);
                if (!await reader.ReadAsync(cancellationToken)) return null;
                string? Text(int i) => reader.IsDBNull(i) ? null : reader.GetString(i);
                // NCBI has realm and kingdom; GTDB's domain stands in for the superkingdom and has neither.
                Lineage Ncbi(int i) => new(Text(i), Text(i + 1), Text(i + 2), Text(i + 3), Text(i + 4), Text(i + 5), Text(i + 6), Text(i + 7), Text(i + 8));
                Lineage Gtdb(int i) => new(Text(i), null, null, Text(i + 1), Text(i + 2), Text(i + 3), Text(i + 4), Text(i + 5), Text(i + 6));
                const int c = recordColumnCount;
                detail = new CatalogueDetail(ReadRecord(reader), reader.GetString(c), reader.GetString(c + 1),
                    reader.IsDBNull(c + 2) ? null : reader.GetFieldValue<byte[]>(c + 2), Text(c + 3), reader.GetFieldValue<DateTimeOffset>(c + 4),
                    reader.GetFieldValue<byte[]>(c + 5), reader.GetInt64(c + 6), Text(c + 7), reader.GetBoolean(c + 8),
                    Text(c + 9), reader.IsDBNull(c + 10) ? null : Ncbi(c + 10),
                    reader.IsDBNull(c + 19) ? null : new GtdbLineage(reader.GetString(c + 19), reader.GetString(c + 20), Gtdb(c + 21)),
                    []);
            }

            await using var levels = dataSource.CreateCommand(
                "select method_version, window_size, window_count, byte_size from profile_level where sequence_id = $1 " +
                "order by method_version, window_size");
            levels.Parameters.AddWithValue(id);
            var list = new List<LevelInfo>();
            await using var levelReader = await levels.ExecuteReaderAsync(cancellationToken);
            while (await levelReader.ReadAsync(cancellationToken))
                list.Add(new LevelInfo(levelReader.GetString(0), levelReader.GetInt32(1), levelReader.GetInt32(2), levelReader.GetInt64(3)));
            return detail with { Levels = list };
        }

        /// <summary>Whether the viewer may see this record (<see cref="CatalogueFilter.ViewerId"/>); false also when it doesn't exist.</summary>
        public static async Task<bool> IsVisibleAsync(NpgsqlDataSource dataSource, long id, long? viewerId, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand();
            var sql = new StringBuilder("select 1").Append(recordFrom).Append(" where r.id = $1");
            command.Parameters.AddWithValue(id);
            AppendFilter(command, sql, new CatalogueFilter(ViewerId: viewerId));
            command.CommandText = sql.ToString();
            return await command.ExecuteScalarAsync(cancellationToken) is not null;
        }

        /// <summary>The checks every query applies, for a caller that stores a filter to run later (a queued search).</summary>
        public static void Validate(CatalogueFilter filter)
        {
            if ((filter.Rank is null) != (filter.Taxon is null))
                throw new ArgumentException("Rank and taxon filter go together.", nameof(filter));
            if (filter.Rank is not null && !rankColumns.ContainsKey(filter.Rank))
                throw new ArgumentException($"Unknown rank '{filter.Rank}'. Known: {string.Join(", ", Ranks)}.", nameof(filter));
        }

        private static void AppendFilter(NpgsqlCommand command, StringBuilder sql, CatalogueFilter filter)
        {
            Validate(filter);

            void Where(string condition, object value)
            {
                command.Parameters.AddWithValue(value);
                sql.Append(" and ").Append(condition.Replace("@", "$" + command.Parameters.Count));
            }

            if (filter.Kind is not null) Where("r.kind = @", filter.Kind);
            if (filter.Rank is not null) Where($"r.{rankColumns[filter.Rank]} = @", filter.Taxon!);
            if (filter.MinLengthBp is not null) Where("r.length_bp >= @", filter.MinLengthBp.Value);
            if (filter.MaxLengthBp is not null) Where("r.length_bp <= @", filter.MaxLengthBp.Value);
            if (filter.Provider is not null) Where("s.provider = @", filter.Provider);
            if (filter.NameContains is not null) Where(@"r.name ilike '%' || @ || '%' escape '\'", EscapeLike(filter.NameContains));
            if (filter.ViewerId is not null) Where("(s.owner_id is null or s.owner_id = @)", filter.ViewerId.Value);
            else sql.Append(" and s.owner_id is null");
        }

        private static async Task<List<CatalogueRecord>> ReadRecordsAsync(NpgsqlCommand command, int capacity, CancellationToken cancellationToken = default)
        {
            var records = new List<CatalogueRecord>(capacity);
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                records.Add(ReadRecord(reader));
            return records;
        }

        // Columns in recordSelect order.
        private static CatalogueRecord ReadRecord(NpgsqlDataReader reader)
        {
            string? Text(int i) => reader.IsDBNull(i) ? null : reader.GetString(i);
            return new CatalogueRecord(
                reader.GetInt64(0), reader.GetString(1), Text(2), reader.GetString(3), Text(4), reader.GetString(5), reader.GetInt64(6),
                reader.IsDBNull(7) ? null : reader.GetFloat(7),
                reader.IsDBNull(8) ? null : reader.GetInt32(8),
                new Lineage(Text(9), Text(10), Text(11), Text(12), Text(13), Text(14), Text(15), Text(16), Text(17)));
        }

        private static string EscapeLike(string text) =>
            text.Replace(@"\", @"\\").Replace("%", @"\%").Replace("_", @"\_");
    }
}
