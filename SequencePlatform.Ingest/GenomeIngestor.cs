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
    /// string), only the profile ladder, the D8 per-record hash and dinucleotide counts (same pass),
    /// annotations and source are written, and the temp file is
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

        // A genome is a duplicate when every record is already stored under this provider. BV-BRC lists the same
        // assembly under two genome ids (seen: 1002689.3 / .4, both GCA_025264625.1: same accessions) and a GenBank
        // record next to its RefSeq copy (HM119410 / NC_055339: different accessions, same sequence_md5). Both are
        // caught before download; a record without an MD5 is caught after profiling by its SHA-256. The assembly
        // accession isn't a key: the segments of one virus share it.
        internal static Task<StoredCopy> FindStoredCopyAsync(NpgsqlDataSource dataSource, string column, string[] keys) =>
            FindStoredCopyAsync(dataSource, column, keys, keys.Distinct(StringComparer.Ordinal).Count());

        internal static Task<StoredCopy> FindStoredCopyBySha256Async(NpgsqlDataSource dataSource, byte[][] hashes) =>
            FindStoredCopyAsync(dataSource, "sequence_sha256", hashes, hashes.Select(Convert.ToHexString).Distinct().Count());

        /// <param name="column">record_id, provider_md5 or sequence_sha256; spliced into the SQL, so never user input.</param>
        private static async Task<StoredCopy> FindStoredCopyAsync(NpgsqlDataSource dataSource, string column, object keys, int distinctKeys)
        {
            if (distinctKeys == 0) return new StoredCopy(0, 0, null);
            await using var command = dataSource.CreateCommand(
                $"select count(distinct r.{column}), min(s.id) from sequence_record r join sequence_source s on s.id = r.source_id " +
                $"where s.provider = $1 and r.{column} = any($2)");
            command.Parameters.AddWithValue(Provider);
            command.Parameters.AddWithValue(keys);
            await using var reader = await command.ExecuteReaderAsync();
            await reader.ReadAsync();
            var stored = (int)reader.GetInt64(0);
            return new StoredCopy(distinctKeys, stored, stored == distinctKeys ? reader.GetInt64(1) : null);
        }

        /// <param name="SourceId">The lowest stored source holding a match, set only when every key matched.</param>
        internal sealed record StoredCopy(int Keys, int Stored, long? SourceId);

        private async Task<(long SourceId, string By)?> FindStoredCopyBeforeDownloadAsync(IReadOnlyCollection<BvBrcContig> contigs)
        {
            if ((await FindStoredCopyAsync(dataSource, "record_id", contigs.Select(c => c.Accession).ToArray())).SourceId is { } byAccession)
                return (byAccession, "accession");
            if (contigs.Count == 0 || contigs.Any(c => c.SequenceMd5 is null)) return null;
            if ((await FindStoredCopyAsync(dataSource, "provider_md5", contigs.Select(c => c.SequenceMd5!).ToArray())).SourceId is { } byMd5)
                return (byMd5, "md5");
            return null;
        }

        /// <param name="DuplicateBy">How a duplicate was recognised: accession, md5 or sha256.</param>
        /// <param name="AlreadyStored">Records of a stored genome whose sequence was already stored (a partial duplicate).</param>
        public sealed record Outcome(string Status, int Records = 0, int Levels = 0, long Bases = 0, string? Error = null,
            long? DuplicateOf = null, string? DuplicateBy = null, int AlreadyStored = 0);

        /// <summary>
        /// Ingests one registered crawl item and records the outcome in crawl_item. Never throws for a
        /// per-genome problem: a failure is recorded with its message and returned, so the batch goes on.
        /// </summary>
        public async Task<Outcome> IngestAsync(BvBrcGenome genome)
        {
            try
            {
                var contigs = (await client.ListContigsAsync(genome.GenomeId)).ToDictionary(c => c.Accession, StringComparer.Ordinal);
                if (await FindStoredCopyBeforeDownloadAsync(contigs.Values) is { } before)
                    return await MarkDuplicateAsync(genome, before.SourceId, before.By);

                var file = await client.DownloadFastaAsync(genome.GenomeId);
                try
                {
                    var records = ProfileRecords(file.Path).ToList();
                    var stored = await FindStoredCopyBySha256Async(dataSource, records.Select(r => r.SequenceSha256).ToArray());
                    if (stored.SourceId is { } bySha256)
                        return await MarkDuplicateAsync(genome, bySha256, "sha256");

                    var tooLarge = await StoreAsync(genome, file, records, contigs);
                    return new Outcome(tooLarge ? CrawlState.TooLarge : CrawlState.Done,
                        records.Count, records.Sum(r => r.Levels.Count), records.Sum(r => r.LengthBp), AlreadyStored: stored.Stored);
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

        private async Task<Outcome> MarkDuplicateAsync(BvBrcGenome genome, long duplicateOf, string by)
        {
            await CrawlState.MarkAsync(dataSource, Provider, genome.GenomeId, CrawlState.Duplicate, duplicateOf: duplicateOf);
            return new Outcome(CrawlState.Duplicate, DuplicateOf: duplicateOf, DuplicateBy: by);
        }

        private sealed record ProfiledRecord(string RecordId, string Description, long LengthBp, double? GcFraction,
            List<(int Factor, ZCurveWindowedProfile Level)> Levels, byte[] SequenceSha256, long[] DinucleotideCounts);

        private static IEnumerable<ProfiledRecord> ProfileRecords(string fastaPath)
        {
            foreach (var (headerId, description) in ReadHeaders(fastaPath))
            {
                using var features = new SequenceFeatureReader(FastaReader.OpenSequenceReader(fastaPath, headerId));
                var fine = ZCurveWindowedProfile.FromStream(features, fineWindowSize);
                var (sequenceSha256, dinucleotideCounts) = features.Finish();

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
                yield return new ProfiledRecord(recordId, description, features.CharactersRead, gcFraction, levels,
                    sequenceSha256, dinucleotideCounts);
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
                    "tax_superkingdom, tax_realm, tax_kingdom, tax_phylum, tax_class, tax_order, tax_family, tax_genus, tax_species, annotations, " +
                    "sequence_sha256, dinucleotide_counts, provider_md5) " +
                    "values ($1, $2, $3, $4, 'public_microbe', $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20) returning id",
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
                    p.AddWithValue(record.SequenceSha256);
                    p.AddWithValue(record.DinucleotideCounts);
                    p.AddWithValue(OrNull(contig?.SequenceMd5));
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
    }
}
