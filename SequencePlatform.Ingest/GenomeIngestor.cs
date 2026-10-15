using System.Text.Json;
using Npgsql;
using NpgsqlTypes;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// Download -> profile -> store -> delete, one BV-BRC genome at a time. The DNA is streamed from
    /// the temp file into ZCurveWindowedProfile.FromStream record by record (never held as a whole
    /// string), only the profile ladder, annotations and source are written, and the temp file is
    /// deleted whether or not ingestion succeeded (docs/database-decision.md section 4.1).
    /// </summary>
    internal sealed class GenomeIngestor
    {
        public const string Provider = "BV-BRC";
        private const string license = "BV-BRC public data; cite BV-BRC (Olson et al., NAR 2023) and the original submitters";
        private const string method = "zcurve-windowed";
        private const string methodVersion = "1";
        private const int fineWindowSize = ZCurveWindowedProfile.DefaultWindowSize;

        private static readonly string[] annotationFields =
        {
            "strain", "serovar", "host_name", "disease", "isolation_source", "isolation_country", "collection_year",
            "gram_stain", "genome_quality", "checkm_completeness", "checkm_contamination", "gc_content",
            "assembly_accession", "refseq_accessions", "genbank_accessions", "bioproject_accession", "biosample_accession",
            "genome_status", "reference_genome", "taxon_id", "genome_length",
        };

        private readonly NpgsqlDataSource dataSource;
        private readonly BvBrcClient client;

        public GenomeIngestor(NpgsqlDataSource dataSource, BvBrcClient client)
        {
            this.dataSource = dataSource;
            this.client = client;
        }

        // BV-BRC sometimes lists the same assembly under two genome ids (seen: 1002689.3 / .4, both
        // GCA_025264625.1). A genome whose every contig accession is already stored is a duplicate.
        private async Task<bool> AllContigsAlreadyStoredAsync(IReadOnlyCollection<string> accessions)
        {
            if (accessions.Count == 0) return false;
            await using var command = dataSource.CreateCommand(
                "select count(distinct r.record_id) from sequence_record r join sequence_source s on s.id = r.source_id " +
                "where s.provider = $1 and r.record_id = any($2)");
            command.Parameters.AddWithValue(Provider);
            command.Parameters.AddWithValue(accessions.ToArray());
            return (long)(await command.ExecuteScalarAsync())! == accessions.Count;
        }

        public sealed record Outcome(string Status, int Records = 0, int Levels = 0, long Bases = 0, string? Error = null);

        /// <summary>
        /// Ingests one registered crawl item and records the outcome in crawl_item. Never throws for a
        /// per-genome problem: a failure is recorded with its message and returned, so the batch goes on.
        /// </summary>
        public async Task<Outcome> IngestAsync(BvBrcGenome genome)
        {
            try
            {
                var contigs = (await client.ListContigsAsync(genome.GenomeId)).ToDictionary(c => c.Accession, StringComparer.Ordinal);
                if (await AllContigsAlreadyStoredAsync(contigs.Keys))
                {
                    await CrawlState.MarkAsync(dataSource, Provider, genome.GenomeId, CrawlState.Duplicate);
                    return new Outcome(CrawlState.Duplicate);
                }

                var file = await client.DownloadFastaAsync(genome.GenomeId);
                try
                {
                    var records = ProfileRecords(file.Path).ToList();
                    var tooLarge = await StoreAsync(genome, file, records, contigs);
                    return new Outcome(tooLarge ? CrawlState.TooLarge : CrawlState.Done,
                        records.Count, records.Sum(r => r.Levels.Count), records.Sum(r => r.LengthBp));
                }
                finally
                {
                    File.Delete(file.Path); // the no-sequence rule: nothing downloaded outlives ingestion
                }
            }
            catch (Exception e)
            {
                // The store transaction rolled back, so nothing partial is left; the item stays retryable.
                await CrawlState.MarkAsync(dataSource, Provider, genome.GenomeId, CrawlState.Failed, $"{e.GetType().Name}: {e.Message}");
                return new Outcome(CrawlState.Failed, Error: e.Message);
            }
        }

        private sealed record ProfiledRecord(string RecordId, string Description, long LengthBp, double? GcFraction,
            List<(int Factor, ZCurveWindowedProfile Level)> Levels);

        private static IEnumerable<ProfiledRecord> ProfileRecords(string fastaPath)
        {
            foreach (var (headerId, description) in ReadHeaders(fastaPath))
            {
                using var counting = new CountingReader(FastaReader.OpenSequenceReader(fastaPath, headerId));
                var fine = ZCurveWindowedProfile.FromStream(counting, fineWindowSize);

                var levels = new List<(int, ZCurveWindowedProfile)>();
                double? gcFraction = null;
                if (fine.WindowCount > 0)
                {
                    var ladder = new ZCurveProfileLadder(fine);
                    foreach (var factor in ladder.Factors)
                    {
                        var level = ladder.Level(factor);
                        if (level.WindowCount > 0) levels.Add((factor, level));
                    }
                    // z = +1 for G/C, -1 for A/T, so the mean z over all windows is 2*GC - 1.
                    gcFraction = (fine.Windows.Average(w => w.Z) + 1) / 2;
                }

                // BV-BRC headers look like ">accn|NC_003210   description"; store the bare accession.
                var recordId = headerId.StartsWith("accn|", StringComparison.Ordinal) ? headerId[5..] : headerId;
                yield return new ProfiledRecord(recordId, description, counting.CharactersRead, gcFraction, levels);
            }
        }

        // Header ids exactly as FastaReader parses them (text after '>' up to the first space),
        // so each can be handed back to FastaReader.OpenSequenceReader.
        private static IEnumerable<(string Id, string Description)> ReadHeaders(string fastaPath)
        {
            foreach (var line in File.ReadLines(fastaPath))
            {
                if (line.Length == 0 || line[0] != '>') continue;
                var header = line.Substring(1);
                var space = header.IndexOf(' ');
                yield return space < 0 ? (header, "") : (header.Substring(0, space), header.Substring(space + 1).Trim());
            }
        }

        // Returns true when some record's 500 bp level didn't fit inline (crawl status too_large: that record
        // can't be searched at full resolution until object storage exists). Coarser levels that fit are kept.
        private async Task<bool> StoreAsync(BvBrcGenome genome, DownloadedFile file, List<ProfiledRecord> records, Dictionary<string, BvBrcContig> contigs)
        {
            bool fineLevelSkipped = false;
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();

            long sourceId;
            await using (var insertSource = new NpgsqlCommand(
                "insert into sequence_source (provider, url, accession, retrieved_at, raw_sha256, raw_byte_size, license) " +
                "values ($1, $2, $3, now(), $4, $5, $6) returning id", connection, transaction))
            {
                insertSource.Parameters.AddWithValue(Provider);
                insertSource.Parameters.AddWithValue(file.Url);
                insertSource.Parameters.AddWithValue(genome.GenomeId);
                insertSource.Parameters.AddWithValue(file.Sha256);
                insertSource.Parameters.AddWithValue(file.ByteSize);
                insertSource.Parameters.AddWithValue(license);
                sourceId = (long)(await insertSource.ExecuteScalarAsync())!;
            }

            var superkingdom = genome.Field("superkingdom");
            foreach (var record in records)
            {
                contigs.TryGetValue(record.RecordId, out var contig);
                long sequenceId;
                await using (var insertRecord = new NpgsqlCommand(
                    "insert into sequence_record (source_id, record_id, name, kind, data_class, length_bp, gc_fraction, ncbi_taxid, " +
                    "tax_superkingdom, tax_realm, tax_kingdom, tax_phylum, tax_class, tax_order, tax_family, tax_genus, tax_species, annotations) " +
                    "values ($1, $2, $3, $4, 'public_microbe', $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17) returning id",
                    connection, transaction))
                {
                    var p = insertRecord.Parameters;
                    p.AddWithValue(sourceId);
                    p.AddWithValue(record.RecordId);
                    p.AddWithValue(genome.Name);
                    p.AddWithValue(Kind(superkingdom, contig?.SequenceType));
                    p.AddWithValue(record.LengthBp);
                    p.AddWithValue(record.GcFraction.HasValue ? (object)(float)record.GcFraction.Value : DBNull.Value);
                    p.AddWithValue(int.TryParse(genome.Field("taxon_id"), out var taxid) ? (object)taxid : DBNull.Value);
                    p.AddWithValue(OrNull(superkingdom));
                    p.AddWithValue(OrNull(Realm(superkingdom, genome.LineageNames)));
                    foreach (var rank in new[] { "kingdom", "phylum", "class", "order", "family", "genus", "species" })
                        p.AddWithValue(OrNull(genome.Field(rank)));
                    p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Jsonb, Value = Annotations(genome, contig, record.Description) });
                    sequenceId = (long)(await insertRecord.ExecuteScalarAsync())!;
                }

                foreach (var (factor, level) in record.Levels)
                    if (!await ProfileLevels.InsertAsync(connection, transaction, sequenceId, method, methodVersion, factor, level))
                    {
                        fineLevelSkipped |= factor == 1;
                        Console.WriteLine($"    skipped {level.WindowSize} bp level of {record.RecordId}: over {ProfileLevels.MaxInlineBytes:N0} B, needs object storage");
                    }
            }

            await CrawlState.MarkStoredAsync(connection, transaction, Provider, genome.GenomeId, sourceId, fineLevelSkipped);
            await transaction.CommitAsync();
            return fineLevelSkipped;
        }

        private static string Kind(string? superkingdom, string? sequenceType) => superkingdom switch
        {
            "Viruses" => "virus",
            "Archaea" => "archaeon",
            "Eukaryota" => "eukaryote",
            "Bacteria" when string.Equals(sequenceType, "plasmid", StringComparison.OrdinalIgnoreCase) => "plasmid",
            "Bacteria" => "bacterium",
            _ => "other",
        };

        // BV-BRC has no realm field; for viruses it's the lineage entry right after "Viruses" (ICTV realms end in "-viria").
        private static string? Realm(string? superkingdom, IReadOnlyList<string> lineage) =>
            superkingdom == "Viruses" && lineage.Count > 1 && lineage[1].EndsWith("viria", StringComparison.Ordinal) ? lineage[1] : null;

        private static object OrNull(string? value) => string.IsNullOrEmpty(value) ? DBNull.Value : value;

        private static string Annotations(BvBrcGenome genome, BvBrcContig? contig, string headerDescription)
        {
            var bvbrc = new Dictionary<string, object?>();
            foreach (var field in annotationFields)
                if (genome.Raw.TryGetProperty(field, out var value) && value.ValueKind != JsonValueKind.Null)
                    bvbrc[field] = value;
            bvbrc["genome_id"] = genome.GenomeId;

            var annotations = new Dictionary<string, object?>
            {
                ["bvbrc"] = bvbrc,
                ["contig"] = new Dictionary<string, object?>
                {
                    ["sequence_type"] = contig?.SequenceType,
                    ["topology"] = contig?.Topology,
                    ["description"] = contig?.Description ?? headerDescription,
                },
            };
            return JsonSerializer.Serialize(annotations);
        }

        // Counts bases handed to FromStream, so length_bp is exact without keeping the sequence.
        private sealed class CountingReader : TextReader
        {
            private readonly TextReader inner;

            public CountingReader(TextReader inner) => this.inner = inner;

            public long CharactersRead { get; private set; }

            public override int Read(char[] buffer, int index, int count)
            {
                int read = inner.Read(buffer, index, count);
                CharactersRead += read;
                return read;
            }

            protected override void Dispose(bool disposing)
            {
                if (disposing) inner.Dispose();
                base.Dispose(disposing);
            }
        }
    }
}
