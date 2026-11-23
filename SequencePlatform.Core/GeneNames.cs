using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// Each stored gene's name as its provider gives it (migration 0025): BV-BRC's product, NCBI's [protein=...]. A record's names are one
    /// array in the order of its stored gene profile (<see cref="GeneProfiles.ReadAsync"/>), null where none was given. Annotation text only.
    /// </summary>
    public static class GeneNames
    {
        public const string BvBrc = "bvbrc";
        public const string Ncbi = "ncbi";
        public const int MaxLength = 1000;

        /// <summary>A provider's name as stored: trimmed, cut to <see cref="MaxLength"/>, null when empty.</summary>
        public static string? Clean(string? name)
        {
            var trimmed = name?.Trim();
            return string.IsNullOrEmpty(trimmed) ? null : trimmed.Length > MaxLength ? trimmed[..MaxLength] : trimmed;
        }

        /// <summary>
        /// Replaces the record's names inside the caller's transaction: one per gene of its stored profile, in that order. New names are
        /// added to gene_name_text in sorted order, so two writers adding the same names wait on each other rather than deadlock.
        /// </summary>
        /// <param name="unmatched">Stored genes the provider's list didn't have.</param>
        public static async Task WriteAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId, string source,
            IReadOnlyList<string?> names, int unmatched = 0)
        {
            var cleaned = names.Select(Clean).ToArray();
            var distinct = cleaned.OfType<string>().Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToArray();
            var ids = new Dictionary<string, int>(StringComparer.Ordinal);
            if (distinct.Length > 0)
            {
                await using (var insert = new NpgsqlCommand(
                    "insert into gene_name_text (name) select unnest($1::text[]) on conflict (name) do nothing", connection, transaction))
                {
                    insert.Parameters.AddWithValue(distinct);
                    await insert.ExecuteNonQueryAsync();
                }
                await using var select = new NpgsqlCommand("select id, name from gene_name_text where name = any($1)", connection, transaction);
                select.Parameters.AddWithValue(distinct);
                await using var reader = await select.ExecuteReaderAsync();
                while (await reader.ReadAsync()) ids[reader.GetString(1)] = reader.GetInt32(0);
            }

            await using var write = new NpgsqlCommand(
                "insert into gene_name (sequence_id, source, names, named, unmatched) values ($1, $2, $3, $4, $5) " +
                "on conflict (sequence_id) do update set source = excluded.source, names = excluded.names, named = excluded.named, " +
                "unmatched = excluded.unmatched, checked_at = now()", connection, transaction);
            write.Parameters.AddWithValue(sequenceId);
            write.Parameters.AddWithValue(source);
            write.Parameters.AddWithValue(cleaned.Select(n => n is null ? 0 : ids[n]).ToArray());
            write.Parameters.AddWithValue(cleaned.Count(n => n is not null));
            write.Parameters.AddWithValue(unmatched);
            await write.ExecuteNonQueryAsync();
        }

        /// <summary>The record's names in its gene profile's order; null when none are stored.</summary>
        public static async Task<string?[]?> ReadAsync(NpgsqlDataSource dataSource, long sequenceId, CancellationToken cancellationToken = default) =>
            (await ReadManyAsync(dataSource, [sequenceId], cancellationToken)).GetValueOrDefault(sequenceId);

        /// <summary>Several records' names, by record; records without stored names are left out.</summary>
        public static async Task<Dictionary<long, string?[]>> ReadManyAsync(NpgsqlDataSource dataSource, IReadOnlyCollection<long> ids,
            CancellationToken cancellationToken = default)
        {
            var result = new Dictionary<long, string?[]>();
            if (ids.Count == 0) return result;
            await using var command = dataSource.CreateCommand(
                "select n.sequence_id, array(select t.name from unnest(n.names) with ordinality u(id, i) left join gene_name_text t on t.id = u.id " +
                "order by u.i) from gene_name n where n.sequence_id = any($1)");
            command.Parameters.AddWithValue(ids.ToArray());
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                result[reader.GetInt64(0)] = reader.GetFieldValue<string?[]>(1);
            return result;
        }

        /// <summary>
        /// The names of genes given by where they lie (protein search hits keep coordinates, not a gene index): per record, the stored
        /// gene on that sequence with that start and end. A gene without a stored name, or not found, is left out.
        /// </summary>
        public static async Task<Dictionary<(long SequenceId, int Sequence, int Start, int End), string>> AtAsync(NpgsqlDataSource dataSource,
            IReadOnlyCollection<(long SequenceId, int Sequence, int Start, int End)> genes, CancellationToken cancellationToken = default)
        {
            var result = new Dictionary<(long, int, int, int), string>();
            var names = await ReadManyAsync(dataSource, genes.Select(g => g.SequenceId).ToHashSet(), cancellationToken);
            foreach (var group in genes.Where(g => names.ContainsKey(g.SequenceId)).GroupBy(g => g.SequenceId))
            {
                if (await GeneProfiles.ReadAsync(dataSource, group.Key) is not { } profile) continue;
                var recordNames = names[group.Key];
                var index = new Dictionary<(int, int, int), int>();
                for (int g = 0; g < profile.Genes.Count && g < recordNames.Length; g++)
                    index.TryAdd((profile.Genes[g].Chromosome, profile.Genes[g].Start, profile.Genes[g].End), g);
                foreach (var gene in group)
                    if (index.TryGetValue((gene.Sequence, gene.Start, gene.End), out int g) && recordNames[g] is { } name)
                        result[gene] = name;
            }
            return result;
        }

        /// <summary>Names for a profile's genes from a provider's list keyed by where each lies, and how many stored genes it lacked.</summary>
        public static (string?[] Names, int Unmatched) Match(ZCurveGeneProfile profile, IReadOnlyDictionary<(int Sequence, int Start, int End, bool Reverse), string?> byPlace)
        {
            var names = new string?[profile.Genes.Count];
            int unmatched = 0;
            for (int g = 0; g < names.Length; g++)
            {
                var gene = profile.Genes[g];
                if (byPlace.TryGetValue((gene.Chromosome, gene.Start, gene.End, gene.Reverse), out var name)) names[g] = name;
                else unmatched++;
            }
            return (names, unmatched);
        }
    }
}
