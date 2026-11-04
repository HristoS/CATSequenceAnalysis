using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>One stored record as the backfill checks it: the hash and dinucleotide counts it was ingested with.</summary>
    internal sealed record StoredRecordFeatures(long SequenceId, string RecordId, byte[]? Sha256, long[]? DinucleotideCounts);

    /// <summary>
    /// 6-mer counts for records stored before migration 0013 (D8 widened 2026-11-04). Each BV-BRC genome is downloaded again and
    /// read once through <see cref="SequenceFeatureReader"/>, which counts what's missing and recomputes what was stored. A record
    /// gets its counts only when the download is the sequence that was profiled (same SHA-256 and dinucleotide counts); otherwise
    /// its row says why, so a rerun skips it. Nothing else is recomputed, and the file is deleted as soon as it's read.
    /// </summary>
    internal sealed class KmerBackfill
    {
        private readonly NpgsqlDataSource dataSource;
        private readonly BvBrcClient client;

        public KmerBackfill(NpgsqlDataSource dataSource, BvBrcClient client)
        {
            this.dataSource = dataSource;
            this.client = client;
        }

        /// <summary>BV-BRC genomes with a public record that has no 6-mer row yet; with a seed, in a fixed random order (a trial sample).</summary>
        public async Task<List<string>> PendingGenomesAsync(int limit, int? seed)
        {
            await using var command = dataSource.CreateCommand(
                "select s.accession from sequence_source s where s.provider = $1 and s.owner_id is null and exists (" +
                "  select 1 from sequence_record r where r.source_id = s.id and not exists (" +
                "    select 1 from record_kmer k where k.sequence_id = r.id and k.k = $2)) " +
                (seed is null ? "order by s.id " : "order by md5(s.accession || $4) ") + "limit $3");
            command.Parameters.AddWithValue(GenomeIngestor.Provider);
            command.Parameters.AddWithValue((short)RecordKmers.K);
            command.Parameters.AddWithValue(limit);
            if (seed is not null) command.Parameters.AddWithValue(seed.Value.ToString());
            var genomes = new List<string>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) genomes.Add(reader.GetString(0));
            return genomes;
        }

        /// <returns>The status written for each of the genome's records that had none.</returns>
        public async Task<List<string>> BackfillAsync(string genomeId)
        {
            var stored = await StoredAsync(genomeId);
            var file = await client.DownloadFastaAsync(genomeId);
            List<(long SequenceId, string Status, int[]? Counts)> results;
            try
            {
                using var fasta = File.OpenText(file.Path);
                results = Check(fasta, stored);
            }
            finally
            {
                File.Delete(file.Path);
            }

            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var (sequenceId, status, counts) in results)
                await RecordKmers.WriteAsync(connection, transaction, sequenceId, status, counts);
            await transaction.CommitAsync();
            return results.Select(r => r.Status).ToList();
        }

        /// <summary>
        /// Reads every record of the FASTA once and matches it to a stored record by id: done with its counts when the hash and
        /// dinucleotide counts agree, mismatch when they don't, missing when the stored record isn't in the file.
        /// </summary>
        public static List<(long SequenceId, string Status, int[]? Counts)> Check(TextReader fasta, IReadOnlyList<StoredRecordFeatures> stored)
        {
            var byId = stored.ToDictionary(s => s.RecordId, StringComparer.Ordinal);
            var results = new List<(long, string, int[]?)>();
            var seen = new HashSet<string>(StringComparer.Ordinal);
            var buffer = new char[1 << 16];
            foreach (var (headerId, _, sequence) in FastaReader.ReadRecordStreams(fasta))
            {
                using var features = new SequenceFeatureReader(sequence);
                while (features.Read(buffer, 0, buffer.Length) > 0) { }
                var (sha256, dinucleotides) = features.Finish();
                var recordId = GenomeIngestor.RecordId(headerId);
                if (!byId.TryGetValue(recordId, out var record) || !seen.Add(recordId)) continue;
                bool same = record.Sha256 is not null && record.Sha256.AsSpan().SequenceEqual(sha256)
                    && record.DinucleotideCounts is not null && record.DinucleotideCounts.AsSpan().SequenceEqual(dinucleotides);
                results.Add(same ? (record.SequenceId, RecordKmers.Done, features.HexamerCounts) : (record.SequenceId, RecordKmers.Mismatch, null));
            }
            results.AddRange(stored.Where(s => !seen.Contains(s.RecordId)).Select(s => (s.SequenceId, RecordKmers.Missing, (int[]?)null)));
            return results;
        }

        private async Task<List<StoredRecordFeatures>> StoredAsync(string genomeId)
        {
            await using var command = dataSource.CreateCommand(
                "select r.id, r.record_id, r.sequence_sha256, r.dinucleotide_counts from sequence_record r " +
                "join sequence_source s on s.id = r.source_id where s.provider = $1 and s.accession = $2 and not exists (" +
                "  select 1 from record_kmer k where k.sequence_id = r.id and k.k = $3)");
            command.Parameters.AddWithValue(GenomeIngestor.Provider);
            command.Parameters.AddWithValue(genomeId);
            command.Parameters.AddWithValue((short)RecordKmers.K);
            var stored = new List<StoredRecordFeatures>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync())
                stored.Add(new StoredRecordFeatures(reader.GetInt64(0), reader.GetString(1),
                    reader.IsDBNull(2) ? null : reader.GetFieldValue<byte[]>(2), reader.IsDBNull(3) ? null : reader.GetFieldValue<long[]>(3)));
            return stored;
        }
    }
}
