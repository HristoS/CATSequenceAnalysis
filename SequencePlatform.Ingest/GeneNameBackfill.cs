using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>A stored record with a gene profile and no stored names yet.</summary>
    internal sealed record PendingNames(long SequenceId, string Source, string RecordId, string Kind, string? FastaUrl, string? LayerSha256);

    /// <summary>
    /// backfill-gene-names (migration 0025): the names of genes stored before the gene passes kept them. BV-BRC records take their
    /// coding sequences' products from BV-BRC's feature list, which is metadata, so no genome is downloaded; NCBI eukaryotes take their
    /// [protein=] names from the CDS file streamed again, its bases read past and dropped. Each stored gene is matched by where it lies;
    /// a gene the provider's list no longer has keeps no name and is counted. Each genome commits alone, so a stopped run resumes.
    /// </summary>
    internal sealed class GeneNameBackfill
    {
        private const string pending = "exists (select 1 from gene_profile g where g.sequence_id = r.id) and " +
            "not exists (select 1 from gene_name n where n.sequence_id = r.id)";
        private readonly NpgsqlDataSource dataSource;

        public GeneNameBackfill(NpgsqlDataSource dataSource) => this.dataSource = dataSource;

        /// <summary>Public records still without names, by provider, optionally of one kind, in order of source.</summary>
        public async Task<List<PendingNames>> PendingAsync(string provider, string? kind, int limit)
        {
            await using var command = dataSource.CreateCommand(
                "select r.id, s.accession, r.record_id, r.kind, s.url, r.annotations #>> '{ncbi,gene_layer,cds_sha256}' " +
                "from sequence_record r join sequence_source s on s.id = r.source_id " +
                "where s.provider = $1 and s.owner_id is null and ($2::text is null or r.kind = $2) and " + pending + " order by s.id, r.id limit $3");
            command.Parameters.AddWithValue(provider);
            command.Parameters.Add(new NpgsqlParameter { Value = (object?)kind ?? DBNull.Value, NpgsqlDbType = NpgsqlTypes.NpgsqlDbType.Text });
            command.Parameters.AddWithValue(limit);
            var records = new List<PendingNames>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync())
                records.Add(new PendingNames(reader.GetInt64(0), reader.GetString(1), reader.GetString(2), reader.GetString(3),
                    reader.IsDBNull(4) ? null : reader.GetString(4), reader.IsDBNull(5) ? null : reader.GetString(5)));
            return records;
        }

        /// <summary>One BV-BRC genome's pending records named from its coding sequences (as listed for a batch), in one transaction.</summary>
        /// <returns>Per record: genes, named genes and genes BV-BRC's list lacked.</returns>
        public async Task<List<(long SequenceId, int Genes, int Named, int Unmatched)>> NameBvBrcAsync(IReadOnlyList<PendingNames> records,
            IReadOnlyList<BvBrcCds> codingSequences)
        {
            var done = new List<(long, int, int, int)>();
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var record in records)
            {
                if (await GeneProfiles.ReadAsync(dataSource, record.SequenceId) is not { } profile) continue;
                var onRecord = codingSequences.Where(c => string.Equals(c.Accession, record.RecordId, StringComparison.Ordinal));
                var (names, unmatched) = GeneNames.Match(profile, GeneBackfill.ProductsByPlace(onRecord));
                await GeneNames.WriteAsync(connection, transaction, record.SequenceId, GeneNames.BvBrc, names, unmatched);
                done.Add((record.SequenceId, names.Length, names.Count(n => GeneNames.Clean(n) is not null), unmatched));
            }
            await transaction.CommitAsync();
            return done;
        }

        /// <summary>One NCBI eukaryote named from its CDS file streamed again.</summary>
        /// <returns>Genes, named genes, genes the file lacked, and whether the file is the one the gene layer was built from.</returns>
        public async Task<(int Genes, int Named, int Unmatched, bool SameFile)> NameNcbiAsync(PendingNames record, NcbiGeneLayer layer)
        {
            if (record.FastaUrl is null) throw new InvalidOperationException($"Record {record.SequenceId} has no source URL.");
            var profile = await GeneProfiles.ReadAsync(dataSource, record.SequenceId)
                ?? throw new InvalidOperationException($"Record {record.SequenceId} has no gene profile.");
            var (names, unmatched, sha256) = await layer.NamesAsync(record.FastaUrl, profile);
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            await GeneNames.WriteAsync(connection, transaction, record.SequenceId, GeneNames.Ncbi, names, unmatched);
            await transaction.CommitAsync();
            return (names.Length, names.Count(n => GeneNames.Clean(n) is not null), unmatched,
                string.Equals(record.LayerSha256, Convert.ToHexStringLower(sha256), StringComparison.Ordinal));
        }
    }
}
