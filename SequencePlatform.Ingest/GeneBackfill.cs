using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>One stored record as the gene pass checks it: its hash and dinucleotides, whether it has 6-mers, and its order (genetic code).</summary>
    internal sealed record StoredGeneRecord(long SequenceId, string RecordId, byte[]? Sha256, long[]? DinucleotideCounts, bool HasHexamers, string? TaxOrder);

    /// <summary>
    /// What the pass found for one record. Proteins, one per stored gene in the profile's order, are kept only in memory, for the trial's
    /// self-search check.
    /// </summary>
    internal sealed record RecordGeneResult(long SequenceId, string Status, ZCurveGeneProfile? Genes, Dictionary<int, int[]>? ProteinCounts,
        int[]? HexamerCounts, int CodingSequences, int CleanTranslations, int Skipped, List<string> Proteins, List<(int Id, string Product)> Families);

    /// <summary>
    /// The gene layer for records already stored (S7.R, decided 2026-11-04; S7.S, 2026-11-05): each BV-BRC genome is downloaded again,
    /// checked against what was profiled (same SHA-256 and dinucleotide counts, as <see cref="KmerBackfill"/>), and each record gets, in one
    /// transaction per genome, its gene profile (<see cref="GeneProfiles"/>, gene_profile rows by length bucket) from BV-BRC's coding
    /// sequences, its amino-acid k-mer counts (<see cref="ProteinKmers"/>) and, when it has none, its 6-mer counts. A record is pending
    /// when it has no k-mer row yet, or its k-mers are done but it has no gene_profile row (a record from before S7.S). The file is
    /// deleted as soon as it's read; no base, and no protein, is stored.
    /// </summary>
    internal sealed class GeneBackfill
    {
        // A record the pass still has to do, as SQL over sequence_record r, with $2 the stored k-mer length.
        private const string pending = "(not exists (select 1 from protein_kmer k where k.sequence_id = r.id and k.k = $2) or (exists (" +
            "select 1 from protein_kmer k where k.sequence_id = r.id and k.k = $2 and k.status = 'done') and not exists (" +
            "select 1 from gene_profile g where g.sequence_id = r.id)))";
        // NCBI orders whose members read TGA as Trp (translation table 4).
        private static readonly string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };
        private readonly NpgsqlDataSource dataSource;
        private readonly BvBrcClient client;

        public GeneBackfill(NpgsqlDataSource dataSource, BvBrcClient client)
        {
            this.dataSource = dataSource;
            this.client = client;
        }

        /// <summary>BV-BRC genomes with a public record still pending (see the class); with a seed, in a fixed random order (a trial sample).</summary>
        public async Task<List<string>> PendingGenomesAsync(int limit, int? seed)
        {
            await using var command = dataSource.CreateCommand(
                "select s.accession from sequence_source s where s.provider = $1 and s.owner_id is null and exists (" +
                "  select 1 from sequence_record r where r.source_id = s.id and " + pending + ") " +
                (seed is null ? "order by s.id " : "order by md5(s.accession || $4) ") + "limit $3");
            command.Parameters.AddWithValue(GenomeIngestor.Provider);
            command.Parameters.AddWithValue((short)ProteinKmers.Ks.Max());
            command.Parameters.AddWithValue(limit);
            if (seed is not null) command.Parameters.AddWithValue(seed.Value.ToString());
            var genomes = new List<string>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) genomes.Add(reader.GetString(0));
            return genomes;
        }

        public Task<Dictionary<string, List<BvBrcCds>>> CodingSequencesAsync(IReadOnlyCollection<string> genomeIds) =>
            client.ListCodingSequencesAsync(genomeIds);

        /// <returns>One result per stored record of the genome that had no amino-acid k-mer row.</returns>
        public async Task<List<RecordGeneResult>> BackfillAsync(string genomeId, IReadOnlyList<BvBrcCds> codingSequences)
        {
            var stored = await StoredAsync(genomeId);
            if (stored.Count == 0) return new List<RecordGeneResult>();
            var file = await client.DownloadFastaAsync(genomeId);
            List<RecordGeneResult> results;
            try
            {
                using var fasta = File.OpenText(file.Path);
                results = Compute(fasta, stored, codingSequences);
            }
            finally
            {
                File.Delete(file.Path);
            }

            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            var byId = stored.ToDictionary(s => s.SequenceId);
            foreach (var result in results)
            {
                if (result.Status == ProteinKmers.Done && !await GeneProfiles.WriteAsync(connection, transaction, result.SequenceId, result.Genes!))
                    throw new InvalidOperationException($"A gene profile bucket of record {result.SequenceId} is over the inline limit.");
                foreach (var k in ProteinKmers.Ks)
                    await ProteinKmers.WriteAsync(connection, transaction, result.SequenceId, k, result.Status,
                        result.Status == ProteinKmers.Done ? result.Genes!.Genes.Count : null, result.ProteinCounts?[k]);
                if (!byId[result.SequenceId].HasHexamers)
                    await RecordKmers.WriteAsync(connection, transaction, result.SequenceId,
                        result.HexamerCounts is null ? (result.Status == ProteinKmers.Missing ? RecordKmers.Missing : RecordKmers.Mismatch) : RecordKmers.Done,
                        result.HexamerCounts);
            }
            await ProteinFamilies.AddAsync(connection, transaction, results.SelectMany(r => r.Families).DistinctBy(f => f.Id).ToList());
            await transaction.CommitAsync();
            return results;
        }

        /// <summary>
        /// Reads every record of the FASTA once. A stored record whose download matches it (hash and dinucleotides) gets its 6-mers, its
        /// gene profile and its amino-acid counts: done when BV-BRC lists coding sequences on it, no_cds when not. One that doesn't match
        /// is mismatch, and one missing from the file is missing. A coding sequence outside the record or not a whole number of codons
        /// is skipped and counted.
        /// </summary>
        public static List<RecordGeneResult> Compute(TextReader fasta, IReadOnlyList<StoredGeneRecord> stored, IReadOnlyList<BvBrcCds> codingSequences)
        {
            var byId = stored.ToDictionary(s => s.RecordId, StringComparer.Ordinal);
            var onRecord = codingSequences.GroupBy(c => c.Accession, StringComparer.Ordinal).ToDictionary(g => g.Key, g => g.OrderBy(c => c.Start).ToList(), StringComparer.Ordinal);
            var results = new List<RecordGeneResult>();
            var seen = new HashSet<string>(StringComparer.Ordinal);
            var buffer = new char[1 << 16];
            foreach (var (headerId, _, stream) in FastaReader.ReadRecordStreams(fasta))
            {
                var sequence = stream.ReadToEnd();
                var recordId = GenomeIngestor.RecordId(headerId);
                if (!byId.TryGetValue(recordId, out var record) || !seen.Add(recordId)) continue;
                using var features = new SequenceFeatureReader(new StringReader(sequence));
                while (features.Read(buffer, 0, buffer.Length) > 0) { }
                var (sha256, dinucleotides) = features.Finish();
                bool same = record.Sha256 is not null && record.Sha256.AsSpan().SequenceEqual(sha256)
                    && record.DinucleotideCounts is not null && record.DinucleotideCounts.AsSpan().SequenceEqual(dinucleotides);
                if (!same)
                {
                    results.Add(Failed(record.SequenceId, ProteinKmers.Mismatch));
                    continue;
                }

                var code = record.TaxOrder is not null && mycoplasmaOrders.Contains(record.TaxOrder) ? ZCurveGeneProfile.mycoplasmaCode : ZCurveGeneProfile.standardCode;
                var cds = onRecord.GetValueOrDefault(recordId) ?? new List<BvBrcCds>();
                var usable = cds.Where(c => c.Start >= 0 && c.End <= sequence.Length && c.End > c.Start && (c.End - c.Start) % 3 == 0).ToList();
                var genes = ZCurveGeneProfile.FromRecord(sequence, usable.Select(c => (c.Start, c.End, c.Reverse, ProteinFamilies.Parse(c.FamilyId))), code);
                var proteins = usable.Select(c => ZCurveGeneProfile.Translate(ZCurveGeneProfile.CodingSequence(sequence, c.Start, c.End, c.Reverse), code)).ToList();
                var proteinAt = new Dictionary<(int, int, bool), string>();
                for (int c = 0; c < usable.Count; c++) proteinAt.TryAdd((usable[c].Start, usable[c].End, usable[c].Reverse), proteins[c]);
                var families = usable.Where(c => c.Product is not null && ProteinFamilies.Parse(c.FamilyId) != ZCurveGeneProfile.NoFamily)
                    .Select(c => (ProteinFamilies.Parse(c.FamilyId), c.Product!)).ToList();
                results.Add(new RecordGeneResult(record.SequenceId, usable.Count > 0 ? ProteinKmers.Done : ProteinKmers.NoCds, genes,
                    usable.Count > 0 ? ProteinKmers.Ks.ToDictionary(k => k, k => AminoAcidKmers.Count(proteins, k)) : null,
                    features.HexamerCounts, cds.Count, proteins.Count(p => !p.Contains('*')), cds.Count - usable.Count,
                    genes.Genes.Select(g => proteinAt[(g.Start, g.End, g.Reverse)]).ToList(), families));
            }
            results.AddRange(stored.Where(s => !seen.Contains(s.RecordId)).Select(s => Failed(s.SequenceId, ProteinKmers.Missing)));
            return results;
        }

        private static RecordGeneResult Failed(long sequenceId, string status) =>
            new(sequenceId, status, null, null, null, 0, 0, 0, new List<string>(), new List<(int, string)>());

        private async Task<List<StoredGeneRecord>> StoredAsync(string genomeId)
        {
            await using var command = dataSource.CreateCommand(
                "select r.id, r.record_id, r.sequence_sha256, r.dinucleotide_counts, " +
                "  exists (select 1 from record_kmer k where k.sequence_id = r.id and k.k = $3 and k.status = $4), r.tax_order " +
                "from sequence_record r join sequence_source s on s.id = r.source_id where s.provider = $1 and s.accession = $2 and " +
                pending.Replace("$2", "$5"));
            command.Parameters.AddWithValue(GenomeIngestor.Provider);
            command.Parameters.AddWithValue(genomeId);
            command.Parameters.AddWithValue((short)RecordKmers.K);
            command.Parameters.AddWithValue(RecordKmers.Done);
            command.Parameters.AddWithValue((short)ProteinKmers.Ks.Max());
            var stored = new List<StoredGeneRecord>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync())
                stored.Add(new StoredGeneRecord(reader.GetInt64(0), reader.GetString(1),
                    reader.IsDBNull(2) ? null : reader.GetFieldValue<byte[]>(2), reader.IsDBNull(3) ? null : reader.GetFieldValue<long[]>(3),
                    reader.GetBoolean(4), reader.IsDBNull(5) ? null : reader.GetString(5)));
            return stored;
        }
    }
}
