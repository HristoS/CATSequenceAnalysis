using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>A stored NCBI eukaryote record the gene pass hasn't done yet: no amino-acid k-mer row.</summary>
    internal sealed record PendingEukaryote(long SequenceId, string Accession, string FastaUrl, int? Taxid);

    /// <summary>One coding sequence of an assembly, spliced and in coding orientation, and where it lies on its sequence.</summary>
    internal sealed record NcbiCds(string SequenceName, string GeneKey, int Start, int End, bool Reverse, string Gene);

    /// <summary>How many coding sequences a CDS file had, and why those left out were.</summary>
    internal sealed class CdsTally
    {
        public int Read { get; set; }
        public int Partial { get; set; }
        public int Pseudo { get; set; }
        public int Exception { get; set; }
        public int Organelle { get; set; }
        public int Location { get; set; }
        public int Frame { get; set; }
        public int Isoforms { get; set; }
    }

    /// <summary>An assembly's gene layer before it's stored: the profile, its proteins for the 4-mers and checks, and its provenance.</summary>
    /// <param name="Proteins">One per gene of <paramref name="Genes"/>, in its order; kept in memory only.</param>
    /// <param name="AllProteins">Every kept coding sequence's protein, under 50 aa too, as the 4-mers count them.</param>
    internal sealed record EukaryoteGenes(ZCurveGeneProfile Genes, List<string> Proteins, List<string> AllProteins, List<string> SequenceNames,
        CdsTally Tally, byte[] CdsSha256);

    /// <summary>What the pass did with one record.</summary>
    internal sealed record EukaryoteGeneResult(PendingEukaryote Record, string Status, string Note, EukaryoteGenes? Layer = null, int MaxBucketBytes = 0);

    /// <summary>
    /// The gene layer for NCBI eukaryotes (P7 S7.T, 2026-11-06). The genome is never downloaded: NCBI's cds_from_genomic file already holds
    /// each annotated coding sequence spliced and in coding orientation, with its location, and the assembly report says which sequences
    /// are nuclear. Per gene (GeneID) only the longest isoform is kept, and partial, pseudo and exception-flagged coding sequences are left
    /// out, since the windows are anchored at the start codon. Genes lie on the assembly's sequences, numbered in assembly order (ZCGP version
    /// 3), and the record's annotations list those sequences' names. The record's identity is its accession.version, whose sequence NCBI
    /// never changes; the annotation can be re-released, so its name and the CDS file's SHA-256 are recorded. No base and no protein is stored.
    /// </summary>
    internal sealed class NcbiGeneLayer
    {
        public const string Skipped = "skipped";
        private const string reportSuffix = "_assembly_report.txt";
        private const string cdsSuffix = "_cds_from_genomic.fna.gz";
        private const string cdsMarker = "_cds_";
        private const int bufferChars = 1 << 20;

        private static readonly Regex attribute = new(@"\[([a-z_]+)=([^\]]*)\]", RegexOptions.Compiled);
        private static readonly Regex geneId = new(@"(?:^|,)GeneID:(\d+)", RegexOptions.Compiled);

        private readonly NpgsqlDataSource dataSource;
        private readonly NcbiGenomeClient client;

        public NcbiGeneLayer(NpgsqlDataSource dataSource, NcbiGenomeClient client)
        {
            this.dataSource = dataSource;
            this.client = client;
        }

        public async Task<List<PendingEukaryote>> PendingAsync(int limit, string? accession = null)
        {
            await using var command = dataSource.CreateCommand(
                "select r.id, s.accession, s.url, r.ncbi_taxid from sequence_record r join sequence_source s on s.id = r.source_id " +
                "where s.provider = $1 and s.owner_id is null and r.kind = $2 and ($3::text is null or s.accession = $3) and " +
                "not exists (select 1 from protein_kmer k where k.sequence_id = r.id and k.k = $4) order by s.id limit $5");
            command.Parameters.AddWithValue(NcbiIngestor.Provider);
            command.Parameters.AddWithValue(NcbiIngestor.Kind);
            command.Parameters.Add(new NpgsqlParameter { Value = (object?)accession ?? DBNull.Value, NpgsqlDbType = NpgsqlTypes.NpgsqlDbType.Text });
            command.Parameters.AddWithValue((short)ProteinKmers.Ks.Max());
            command.Parameters.AddWithValue(limit);
            var pending = new List<PendingEukaryote>();
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync())
                pending.Add(new PendingEukaryote(reader.GetInt64(0), reader.GetString(1), reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetInt32(3)));
            return pending;
        }

        /// <summary>The annotation and genetic code of each record, in as few Datasets calls as the batch allows.</summary>
        public async Task<(Dictionary<string, NcbiAnnotation?> Annotations, Dictionary<int, int> Codes)> LookUpAsync(IReadOnlyCollection<PendingEukaryote> batch) =>
            (await client.AnnotationsAsync(batch.Select(p => p.Accession).ToList()),
             await client.GeneticCodesAsync(batch.Where(p => p.Taxid is not null).Select(p => p.Taxid!.Value).ToList()));

        /// <summary>
        /// Builds one record's layer and, unless <paramref name="write"/> is false, stores it in one transaction. Never throws for a problem
        /// of the record: a record with nothing written is reported and stays pending.
        /// </summary>
        public async Task<EukaryoteGeneResult> RunAsync(PendingEukaryote record, NcbiAnnotation? annotation, int? geneticCode, bool write)
        {
            try
            {
                if (annotation is null)
                {
                    if (write) await StoreNoCdsAsync(record.SequenceId);
                    return new EukaryoteGeneResult(record, ProteinKmers.NoCds, "no annotation at NCBI");
                }
                if (CodeTable(geneticCode) is not { } table)
                    return new EukaryoteGeneResult(record, Skipped, $"genetic code {geneticCode?.ToString() ?? "unknown"} isn't supported");

                var nuclear = NuclearSequences(new StringReader(await client.ReadTextAsync(NcbiGenomeClient.SiblingUrl(record.FastaUrl, reportSuffix))));
                var cdsUrl = NcbiGenomeClient.SiblingUrl(record.FastaUrl, cdsSuffix);
                using var response = await client.OpenFastaAsync(cdsUrl);
                await using var body = await response.Content.ReadAsStreamAsync();
                var (kept, tally, sha256) = await Task.Run(() => Read(body, nuclear));
                if (kept.Count == 0)
                {
                    if (write) await StoreNoCdsAsync(record.SequenceId);
                    return new EukaryoteGeneResult(record, ProteinKmers.NoCds, $"no usable coding sequence of {tally.Read:N0}");
                }

                var layer = Compute(kept, nuclear, table, tally, sha256);
                int maxBucket = BucketPayloads(layer.Genes).Max(p => p.Length);
                if (write)
                    await StoreAsync(record.SequenceId, layer, LayerJson(layer, annotation, geneticCode!.Value, cdsUrl));
                return new EukaryoteGeneResult(record, ProteinKmers.Done, $"{layer.Genes.Genes.Count:N0} genes on {layer.SequenceNames.Count:N0} sequences", layer, maxBucket);
            }
            catch (Exception e) when (e is HttpRequestException or InvalidOperationException or InvalidDataException or IOException or ArgumentException)
            {
                return new EukaryoteGeneResult(record, CrawlState.Failed, $"{e.GetType().Name}: {e.Message}");
            }
        }

        // Table 1 and 11 translate alike; 12 is the CUG-Ser yeasts'. Any other nuclear code would read some codons wrongly.
        internal static string? CodeTable(int? geneticCode) => geneticCode switch
        {
            1 or 11 => ZCurveGeneProfile.standardCode,
            12 => ZCurveGeneProfile.alternativeYeastCode,
            _ => null,
        };

        /// <summary>
        /// The assembly report's nuclear sequences, by both their GenBank and RefSeq accessions, numbered in report order; organelles, as
        /// the genome counts leave them out, aren't in it.
        /// </summary>
        internal static Dictionary<string, int> NuclearSequences(TextReader report)
        {
            var order = new Dictionary<string, int>(StringComparer.Ordinal);
            int row = 0;
            for (var line = report.ReadLine(); line is not null; line = report.ReadLine())
            {
                if (line.Length == 0 || line[0] == '#') continue;
                var columns = line.Split('\t');
                if (columns.Length < 7) continue;
                if (!NcbiIngestor.organelle.IsMatch(columns[3]))
                    foreach (var accession in new[] { columns[4], columns[6] })
                        if (accession != "na") order.TryAdd(accession, row);
                row++;
            }
            return order;
        }

        /// <returns>The longest usable coding sequence of each gene, the tally and the SHA-256 of the gzip file as it came.</returns>
        internal static (List<NcbiCds> Kept, CdsTally Tally, byte[] Sha256) Read(Stream gzipFasta, IReadOnlyDictionary<string, int> nuclear)
        {
            var tally = new CdsTally();
            var best = new Dictionary<string, NcbiCds>(StringComparer.Ordinal);
            using var sha256 = SHA256.Create();
            using var hashing = new CryptoStream(gzipFasta, sha256, CryptoStreamMode.Read, leaveOpen: true);
            using (var gzip = new GZipStream(hashing, CompressionMode.Decompress, leaveOpen: true))
            using (var reader = new StreamReader(gzip, bufferSize: bufferChars))
                foreach (var (id, description, sequence) in FastaReader.ReadRecordStreams(reader))
                {
                    tally.Read++;
                    if (Parse(id, description, sequence.ReadToEnd(), nuclear, tally) is not { } cds) continue;
                    if (best.TryGetValue(cds.GeneKey, out var other))
                    {
                        tally.Isoforms++;
                        if (other.Gene.Length >= cds.Gene.Length) continue;
                    }
                    best[cds.GeneKey] = cds;
                }
            hashing.CopyTo(Stream.Null);
            return (best.Values.ToList(), tally, sha256.Hash!);
        }

        /// <summary>One cds_from_genomic record (>lcl|SEQ_cds_PROTEIN_N [key=value] ...), or null with the reason counted.</summary>
        internal static NcbiCds? Parse(string id, string description, string gene, IReadOnlyDictionary<string, int> nuclear, CdsTally tally)
        {
            int marker = id.IndexOf(cdsMarker, StringComparison.Ordinal);
            var name = marker < 0 ? null : id.StartsWith("lcl|", StringComparison.Ordinal) ? id[4..marker] : id[..marker];
            if (name is null)
            {
                tally.Location++;
                return null;
            }
            if (!nuclear.ContainsKey(name))
            {
                tally.Organelle++;
                return null;
            }

            var fields = new Dictionary<string, string>(StringComparer.Ordinal);
            foreach (Match match in attribute.Matches(description)) fields.TryAdd(match.Groups[1].Value, match.Groups[2].Value);
            if (fields.ContainsKey("pseudo"))
            {
                tally.Pseudo++;
                return null;
            }
            if (fields.ContainsKey("partial"))
            {
                tally.Partial++;
                return null;
            }
            if (fields.ContainsKey("exception"))
            {
                tally.Exception++;
                return null;
            }
            var placed = fields.TryGetValue("location", out var location) ? ParseLocation(location, out _) : null;
            if (placed is null || placed.Value.Bases != gene.Length)
            {
                if (location is not null && location.IndexOfAny(['<', '>']) >= 0) tally.Partial++;
                else tally.Location++;
                return null;
            }
            if (gene.Length % 3 != 0)
            {
                tally.Frame++;
                return null;
            }

            var key = fields.TryGetValue("db_xref", out var xref) && geneId.Match(xref) is { Success: true } g ? "GeneID:" + g.Groups[1].Value
                : fields.TryGetValue("locus_tag", out var tag) ? "locus:" + tag
                : fields.TryGetValue("gene", out var symbol) ? "gene:" + symbol
                : fields.TryGetValue("protein_id", out var protein) ? "protein:" + protein : "id:" + id;
            var (start, end, reverse, _) = placed.Value;
            return new NcbiCds(name, key, start, end, reverse, gene);
        }

        /// <summary>
        /// A GenBank location of exons on one sequence and one strand: a..b, join(a..b,c..d), complement(...) around either, or
        /// complement(a..b) parts inside a join. Null for a partial end (&lt;, &gt;), parts on other sequences or strands, or anything else.
        /// </summary>
        /// <returns>[Start, End) spanning every exon (0-based), the strand and the exons' total length.</returns>
        internal static (int Start, int End, bool Reverse, long Bases)? ParseLocation(string location, out bool partial)
        {
            partial = location.IndexOfAny(['<', '>']) >= 0;
            if (partial) return null;
            bool outer = location.StartsWith("complement(", StringComparison.Ordinal) && location.EndsWith(')');
            var body = outer ? location["complement(".Length..^1] : location;
            if (body.StartsWith("join(", StringComparison.Ordinal) && body.EndsWith(')')) body = body["join(".Length..^1];

            long min = long.MaxValue, max = -1, bases = 0;
            bool? strand = null;
            foreach (var raw in body.Split(','))
            {
                var part = raw.Trim();
                bool reverse = outer;
                if (part.StartsWith("complement(", StringComparison.Ordinal) && part.EndsWith(')'))
                {
                    if (outer) return null;
                    (part, reverse) = (part["complement(".Length..^1], true);
                }
                if (strand is not null && strand != reverse) return null;
                strand = reverse;
                int dots = part.IndexOf("..", StringComparison.Ordinal);
                if (dots < 0 || !long.TryParse(part[..dots], out var a) || !long.TryParse(part[(dots + 2)..], out var b) || a < 1 || b < a) return null;
                (min, max, bases) = (Math.Min(min, a), Math.Max(max, b), bases + b - a + 1);
            }
            return strand is null || max > int.MaxValue ? null : ((int)(min - 1), (int)max, strand.Value, bases);
        }

        /// <summary>The profile over the kept coding sequences, their sequences numbered in assembly order among those that carry one.</summary>
        internal static EukaryoteGenes Compute(List<NcbiCds> kept, IReadOnlyDictionary<string, int> nuclear, string table, CdsTally tally, byte[] sha256)
        {
            var names = kept.Select(c => c.SequenceName).Distinct(StringComparer.Ordinal).OrderBy(n => nuclear[n]).ThenBy(n => n, StringComparer.Ordinal).ToList();
            var number = names.Select((n, i) => (n, i)).ToDictionary(x => x.n, x => x.i, StringComparer.Ordinal);
            var genes = ZCurveGeneProfile.FromCodingSequences(
                kept.Select(c => (number[c.SequenceName], c.Start, c.End, c.Reverse, ZCurveGeneProfile.NoFamily, c.Gene)), table);
            var proteinAt = new Dictionary<(int, int, int, bool), string>();
            var all = new List<string>(kept.Count);
            foreach (var c in kept)
            {
                var protein = ZCurveGeneProfile.Translate(c.Gene, table);
                all.Add(protein);
                proteinAt.TryAdd((number[c.SequenceName], c.Start, c.End, c.Reverse), protein);
            }
            var proteins = genes.Genes.Select(g => proteinAt[(g.Chromosome, g.Start, g.End, g.Reverse)]).ToList();
            return new EukaryoteGenes(genes, proteins, all, names, tally, sha256);
        }

        // As GeneProfiles.WriteAsync splits a profile, to see the largest row before writing.
        private static IEnumerable<byte[]> BucketPayloads(ZCurveGeneProfile genes) =>
            Enumerable.Range(0, genes.Genes.Count).GroupBy(g => GeneProfiles.BucketOf(genes.Genes[g].Residues))
                .Select(b => genes.Select(b).ToBinary()).DefaultIfEmpty(genes.Select(Array.Empty<int>()).ToBinary());

        internal static string LayerJson(EukaryoteGenes layer, NcbiAnnotation annotation, int geneticCode, string cdsUrl) =>
            JsonSerializer.Serialize(new Dictionary<string, object?>
            {
                ["source"] = cdsUrl,
                ["cds_sha256"] = Convert.ToHexStringLower(layer.CdsSha256),
                ["annotation"] = annotation.Name,
                ["annotation_release"] = annotation.ReleaseDate,
                ["annotation_provider"] = annotation.Provider,
                ["genetic_code"] = geneticCode,
                ["coding_sequences"] = layer.Tally.Read,
                ["genes_kept"] = layer.AllProteins.Count,
                ["genes_profiled"] = layer.Genes.Genes.Count,
                ["left_out"] = new Dictionary<string, int>
                {
                    ["other_isoforms"] = layer.Tally.Isoforms, ["partial"] = layer.Tally.Partial, ["pseudo"] = layer.Tally.Pseudo,
                    ["exception"] = layer.Tally.Exception, ["organelle"] = layer.Tally.Organelle, ["location"] = layer.Tally.Location,
                    ["frame"] = layer.Tally.Frame,
                },
                ["sequences"] = layer.SequenceNames,
            });

        internal async Task StoreAsync(long sequenceId, EukaryoteGenes layer, string layerJson)
        {
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            if (!await GeneProfiles.WriteAsync(connection, transaction, sequenceId, layer.Genes))
                throw new InvalidOperationException($"A gene profile bucket of record {sequenceId} is over {GeneProfiles.MaxRowBytes:N0} bytes.");
            foreach (var k in ProteinKmers.Ks)
                await ProteinKmers.WriteAsync(connection, transaction, sequenceId, k, ProteinKmers.Done, layer.Genes.Genes.Count,
                    AminoAcidKmers.Count(layer.AllProteins, k));
            await GeneProfiles.AnnotateLayerAsync(connection, transaction, sequenceId, layerJson);
            await transaction.CommitAsync();
        }

        private async Task StoreNoCdsAsync(long sequenceId)
        {
            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            foreach (var k in ProteinKmers.Ks)
                await ProteinKmers.WriteAsync(connection, transaction, sequenceId, k, ProteinKmers.NoCds, null, null);
            await transaction.CommitAsync();
        }
    }
}
