using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>A stored NCBI eukaryote whose chromosome levels aren't recorded yet.</summary>
    internal sealed record PendingGenome(long SequenceId, string Accession, string FastaUrl, long LengthBp, long[] DinucleotideCounts);

    /// <summary>
    /// chromosome_level for eukaryotes stored before migration 0022: each genome is downloaded again and profiled as it streams, as
    /// ingest-ncbi now does, and kept only if it is the genome the record was made from, i.e. its nuclear length and dinucleotide counts
    /// equal the stored ones (an accession.version's sequence never changes, but its gzip file can be rewritten). A genome that differs
    /// is recorded as a mismatch and not tried again. Each record commits alone, so a stopped run resumes where it was.
    /// </summary>
    internal sealed class ChromosomeBackfill
    {
        private readonly NpgsqlDataSource dataSource;
        private readonly NcbiGenomeClient client;

        public ChromosomeBackfill(NpgsqlDataSource dataSource, NcbiGenomeClient client)
        {
            this.dataSource = dataSource;
            this.client = client;
        }

        public async Task<List<PendingGenome>> PendingAsync(int limit, string? accession = null)
        {
            await using var command = dataSource.CreateCommand(
                "select r.id, s.accession, s.url, r.length_bp, r.dinucleotide_counts from sequence_record r join sequence_source s on s.id = r.source_id " +
                "where s.provider = $1 and s.owner_id is null and r.kind = $2 and r.dinucleotide_counts is not null and " +
                "($3::text is null or s.accession = $3) and r.annotations #> '{ncbi,chromosome_levels}' is null order by r.length_bp, r.id limit $4");
            command.Parameters.AddWithValue(NcbiIngestor.Provider);
            command.Parameters.AddWithValue(NcbiIngestor.Kind);
            command.Parameters.Add(new NpgsqlParameter { Value = (object?)accession ?? DBNull.Value, NpgsqlDbType = NpgsqlTypes.NpgsqlDbType.Text });
            command.Parameters.AddWithValue(limit);
            var pending = new List<PendingGenome>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync())
                pending.Add(new PendingGenome(reader.GetInt64(0), reader.GetString(1), reader.GetString(2), reader.GetInt64(3),
                    reader.GetFieldValue<long[]>(4)));
            return pending;
        }

        /// <returns>A line for the log: what was stored, or why not.</returns>
        public async Task<string> BackfillAsync(PendingGenome genome)
        {
            using var response = await client.OpenFastaAsync(genome.FastaUrl);
            await using var body = await response.Content.ReadAsStreamAsync();
            var (counts, _) = await Task.Run(() => NcbiIngestor.Read(body));
            return await StoreAsync(genome, counts);
        }

        internal async Task<string> StoreAsync(PendingGenome genome, NuclearCounts counts)
        {
            bool same = counts.LengthBp == genome.LengthBp && counts.Counts.AsSpan().SequenceEqual(genome.DinucleotideCounts);
            var summary = same ? counts.Summary() : counts.Summary() with { Status = ChromosomeLevels.Mismatch };
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            if (same) await ChromosomeLevels.InsertAsync(connection, transaction, genome.SequenceId, counts.Chromosomes!);
            await ChromosomeLevels.AnnotateAsync(connection, transaction, genome.SequenceId, summary);
            await transaction.CommitAsync();
            return same
                ? $"{summary.Kept,6:N0} sequences {summary.KeptBp,15:N0} bp kept, {summary.LeftOut:N0} shorter ones ({summary.LeftOutBp:N0} bp) left out"
                : $"MISMATCH: the file now gives {counts.LengthBp:N0} nuclear bp and other dinucleotide counts than the record ({genome.LengthBp:N0} bp); nothing stored";
        }
    }
}
