using System.Text;
using Npgsql;

namespace SequencePlatform.Core
{
    public sealed record CatalogueFilter(
        string? Kind = null,
        string? Rank = null,
        string? Taxon = null,
        long? MinLengthBp = null,
        long? MaxLengthBp = null,
        string? Provider = null,
        string? NameContains = null);

    public sealed record Lineage(string? Superkingdom, string? Realm, string? Kingdom, string? Phylum, string? Class,
        string? Order, string? Family, string? Genus, string? Species);

    public sealed record CatalogueRecord(long Id, string Provider, string? SourceAccession, string RecordId, string? Name,
        string Kind, long LengthBp, float? GcFraction, int? NcbiTaxid, Lineage Lineage);

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

        private const string recordColumns =
            "select r.id, s.provider, s.accession, r.record_id, r.name, r.kind, r.length_bp, r.gc_fraction, r.ncbi_taxid, " +
            "r.tax_superkingdom, r.tax_realm, r.tax_kingdom, r.tax_phylum, r.tax_class, r.tax_order, r.tax_family, r.tax_genus, r.tax_species " +
            "from sequence_record r join sequence_source s on s.id = r.source_id";

        public static async Task<CataloguePage> ListAsync(NpgsqlDataSource dataSource, CatalogueFilter filter, long after = 0, int limit = 50)
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

            var records = await ReadRecordsAsync(command, limit + 1);
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

        /// <summary>The records with these ids, keyed by id; ids that don't exist are left out.</summary>
        public static async Task<IReadOnlyDictionary<long, CatalogueRecord>> ByIdsAsync(NpgsqlDataSource dataSource, IReadOnlyCollection<long> ids,
            CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand(recordColumns + " where r.id = any($1)");
            command.Parameters.AddWithValue(ids.ToArray());
            var records = await ReadRecordsAsync(command, ids.Count, cancellationToken);
            return records.ToDictionary(r => r.Id);
        }

        private static void AppendFilter(NpgsqlCommand command, StringBuilder sql, CatalogueFilter filter)
        {
            if ((filter.Rank is null) != (filter.Taxon is null))
                throw new ArgumentException("Rank and taxon filter go together.", nameof(filter));

            void Where(string condition, object value)
            {
                command.Parameters.AddWithValue(value);
                sql.Append(" and ").Append(condition.Replace("@", "$" + command.Parameters.Count));
            }

            if (filter.Kind is not null) Where("r.kind = @", filter.Kind);
            if (filter.Rank is not null)
            {
                if (!rankColumns.TryGetValue(filter.Rank, out var column))
                    throw new ArgumentException($"Unknown rank '{filter.Rank}'. Known: {string.Join(", ", Ranks)}.", nameof(filter));
                Where($"r.{column} = @", filter.Taxon!);
            }
            if (filter.MinLengthBp is not null) Where("r.length_bp >= @", filter.MinLengthBp.Value);
            if (filter.MaxLengthBp is not null) Where("r.length_bp <= @", filter.MaxLengthBp.Value);
            if (filter.Provider is not null) Where("s.provider = @", filter.Provider);
            if (filter.NameContains is not null) Where(@"r.name ilike '%' || @ || '%' escape '\'", EscapeLike(filter.NameContains));
        }

        // Columns in recordColumns order.
        private static async Task<List<CatalogueRecord>> ReadRecordsAsync(NpgsqlCommand command, int capacity, CancellationToken cancellationToken = default)
        {
            var records = new List<CatalogueRecord>(capacity);
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
            {
                string? Text(int i) => reader.IsDBNull(i) ? null : reader.GetString(i);
                records.Add(new CatalogueRecord(
                    reader.GetInt64(0), reader.GetString(1), Text(2), reader.GetString(3), Text(4), reader.GetString(5), reader.GetInt64(6),
                    reader.IsDBNull(7) ? null : reader.GetFloat(7),
                    reader.IsDBNull(8) ? null : reader.GetInt32(8),
                    new Lineage(Text(9), Text(10), Text(11), Text(12), Text(13), Text(14), Text(15), Text(16), Text(17))));
            }
            return records;
        }

        private static string EscapeLike(string text) =>
            text.Replace(@"\", @"\\").Replace("%", @"\%").Replace("_", @"\_");
    }
}
