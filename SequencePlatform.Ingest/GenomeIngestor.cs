using System.Text.Json;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// Download -> profile -> store -> delete, one BV-BRC genome at a time. Profiling and storing are
    /// <see cref="RecordStore"/>'s (streamed record by record, never a whole string; only the profile
    /// ladder, the D8 hash and dinucleotide counts, annotations and source are written). This class adds
    /// the BV-BRC side: metadata and annotations, duplicate checks, crawl state, and deleting the temp
    /// file whether or not ingestion succeeded (docs/database-decision.md section 4.1).
    /// </summary>
    internal sealed class GenomeIngestor
    {
        public const string Provider = "BV-BRC";
        private const string license = "BV-BRC public data; cite BV-BRC (Olson et al., NAR 2023) and the original submitters";

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
        /// <param name="listedContigs">The genome's contigs from a batched listing; fetched here when null.</param>
        public async Task<Outcome> IngestAsync(BvBrcGenome genome, IReadOnlyList<BvBrcContig>? listedContigs = null)
        {
            try
            {
                var contigs = (listedContigs ?? await client.ListContigsAsync(genome.GenomeId)).ToDictionary(c => c.Accession, StringComparer.Ordinal);
                if (await FindStoredCopyBeforeDownloadAsync(contigs.Values) is { } before)
                    return await MarkDuplicateAsync(genome, before.SourceId, before.By);

                var file = await client.DownloadFastaAsync(genome.GenomeId);
                try
                {
                    var records = RecordStore.Profile(file.Path).ToList();
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

        // BV-BRC headers look like ">accn|NC_003210   description"; the stored record id is the bare accession.
        private static string RecordId(ProfiledRecord record) =>
            record.HeaderId.StartsWith("accn|", StringComparison.Ordinal) ? record.HeaderId[5..] : record.HeaderId;

        // Returns true when some record's 500 bp level didn't fit inline (crawl status too_large: that record
        // can't be searched at full resolution until object storage exists). Coarser levels that fit are kept.
        private async Task<bool> StoreAsync(BvBrcGenome genome, DownloadedFile file, List<ProfiledRecord> records, Dictionary<string, BvBrcContig> contigs)
        {
            var superkingdom = genome.Field("superkingdom");
            var lineage = new Lineage(superkingdom, Realm(superkingdom, genome.LineageNames), genome.Field("kingdom"), genome.Field("phylum"),
                genome.Field("class"), genome.Field("order"), genome.Field("family"), genome.Field("genus"), genome.Field("species"));
            var taxid = int.TryParse(genome.Field("taxon_id"), out var parsed) ? parsed : (int?)null;
            var rows = records
                .Select(record =>
                {
                    var recordId = RecordId(record);
                    contigs.TryGetValue(recordId, out var contig);
                    return (record, new RecordRow(recordId, genome.Name, Kind(superkingdom, contig?.SequenceType), taxid, lineage,
                        Annotations(genome, contig, record.Description), contig?.SequenceMd5));
                })
                .ToList();

            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            var stored = await RecordStore.StoreAsync(connection, transaction,
                new SourceRow(Provider, file.Url, genome.GenomeId, file.Sha256, file.ByteSize, license), rows);
            foreach (var (recordId, windowSize) in stored.SkippedLevels)
                Console.WriteLine($"    skipped {windowSize} bp level of {recordId}: over {ProfileLevels.MaxInlineBytes:N0} B, needs object storage");

            await CrawlState.MarkStoredAsync(connection, transaction, Provider, genome.GenomeId, stored.SourceId, stored.FineLevelSkipped);
            await transaction.CommitAsync();
            return stored.FineLevelSkipped;
        }

        private static string Kind(string? superkingdom, string? sequenceType) => superkingdom switch
        {
            "Viruses" => "virus",
            "Bacteria" or "Archaea" when string.Equals(sequenceType, "plasmid", StringComparison.OrdinalIgnoreCase) => "plasmid",
            "Archaea" => "archaeon",
            "Eukaryota" => "eukaryote",
            "Bacteria" => "bacterium",
            _ => "other",
        };

        // BV-BRC has no realm field; for viruses it's the lineage entry right after "Viruses" (ICTV realms end in "-viria").
        private static string? Realm(string? superkingdom, IReadOnlyList<string> lineage) =>
            superkingdom == "Viruses" && lineage.Count > 1 && lineage[1].EndsWith("viria", StringComparison.Ordinal) ? lineage[1] : null;

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
