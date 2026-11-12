using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    internal sealed record NcbiAssembly(string Accession, string AssemblyName, int Taxid, string Organism, long LengthBp);

    /// <summary>An assembly's gene annotation, as the Datasets genome report has it.</summary>
    internal sealed record NcbiAnnotation(string? Name, string? ReleaseDate, string? Provider);

    /// <summary>What a genome's nuclear records add up to.</summary>
    /// <param name="Counts">Dinucleotide counts summed over the nuclear records; a pair spanning two records isn't counted.</param>
    /// <param name="Chromosomes">The nuclear records of at least <see cref="ChromosomeLevels.MinLengthBp"/>, each with its 500 bp level, in
    /// file order; null when not asked for.</param>
    /// <param name="ShortLeftOut">Nuclear records too short to keep a level, and their bases.</param>
    internal sealed record NuclearCounts(long[] Counts, long LengthBp, int Records, int OrganellesLeftOut,
        IReadOnlyList<ChromosomeLevels.Chromosome>? Chromosomes = null, int ShortLeftOut = 0, long ShortLeftOutBp = 0)
    {
        public ChromosomeLevels.Summary Summary() => new(Chromosomes?.Count ?? 0, Chromosomes?.Sum(c => c.LengthBp) ?? 0, ShortLeftOut, ShortLeftOutBp);
    }

    /// <summary>
    /// Reference eukaryote assemblies from NCBI: the assembly list and lineages from the Datasets API, the genome from its FTP tree over
    /// HTTPS, which serves plain gzip, so a genome streams straight through without a temp file.
    /// </summary>
    internal sealed class NcbiGenomeClient : IDisposable
    {
        private const string datasetsBase = "https://api.ncbi.nlm.nih.gov/datasets/v2/";
        private const string fastaBase = "https://ftp.ncbi.nlm.nih.gov/genomes/all/";
        private const int pageSize = 1000;
        private const int taxaPerRequest = 200;
        private const int accessionsPerRequest = 50;
        private const string genomeSuffix = "_genomic.fna.gz";

        private readonly PoliteHttpClient http;

        public NcbiGenomeClient(string? contactEmail)
        {
            http = new PoliteHttpClient("NCBI", contactEmail, TimeSpan.FromSeconds(1));
        }

        public async Task<List<NcbiAssembly>> ListAssembliesAsync(int taxid)
        {
            var assemblies = new List<NcbiAssembly>();
            string? token = null;
            do
            {
                var url = $"{datasetsBase}genome/taxon/{taxid}/dataset_report?filters.reference_only=true&filters.assembly_level=chromosome" +
                          $"&filters.assembly_level=complete_genome&page_size={pageSize}" + (token is null ? "" : $"&page_token={Uri.EscapeDataString(token)}");
                using var response = await http.SendAsync(url, "application/json");
                token = ParseAssemblies(await response.Content.ReadAsStringAsync(), assemblies);
            } while (token is not null);
            return assemblies;
        }

        /// <returns>The next page's token, null on the last page.</returns>
        internal static string? ParseAssemblies(string json, List<NcbiAssembly> assemblies)
        {
            using var document = JsonDocument.Parse(json);
            if (document.RootElement.TryGetProperty("reports", out var reports))
                foreach (var report in reports.EnumerateArray())
                {
                    var organism = report.GetProperty("organism");
                    var length = report.GetProperty("assembly_stats").TryGetProperty("total_sequence_length", out var l)
                        ? long.Parse(l.GetString()!) : 0;
                    assemblies.Add(new NcbiAssembly(report.GetProperty("accession").GetString()!, report.GetProperty("assembly_info").GetProperty("assembly_name").GetString()!,
                        organism.GetProperty("tax_id").GetInt32(), organism.GetProperty("organism_name").GetString()!, length));
                }
            return document.RootElement.TryGetProperty("next_page_token", out var next) ? next.GetString() : null;
        }

        public async Task<Dictionary<int, Lineage>> LineagesAsync(IReadOnlyCollection<int> taxids)
        {
            var lineages = new Dictionary<int, Lineage>();
            foreach (var chunk in taxids.Distinct().Chunk(taxaPerRequest))
            {
                using var response = await http.SendAsync(
                    $"{datasetsBase}taxonomy/taxon/{string.Join(",", chunk)}/dataset_report?page_size={taxaPerRequest * 2}", "application/json");
                foreach (var taxon in TaxonomyLoader.ParseDatasetReport(await response.Content.ReadAsStringAsync()))
                    lineages[taxon.QueriedTaxid] = taxon.Lineage;
            }
            return lineages;
        }

        // GCF_000146045.2 with assembly R64 is GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz.
        internal static string FastaUrl(NcbiAssembly assembly)
        {
            var name = Regex.Replace(assembly.AssemblyName, "[ /,#()]", "_");
            var digits = assembly.Accession[4..13];
            var folder = $"{assembly.Accession}_{name}";
            return $"{fastaBase}{assembly.Accession[..3]}/{digits[..3]}/{digits[3..6]}/{digits[6..9]}/{folder}/{folder}_genomic.fna.gz";
        }

        public Task<HttpResponseMessage> OpenFastaAsync(string url) => http.SendAsync(url, "application/gzip", HttpCompletionOption.ResponseHeadersRead);

        /// <summary>Each assembly's gene annotation as the Datasets genome report names it; an assembly without one maps to null.</summary>
        public async Task<Dictionary<string, NcbiAnnotation?>> AnnotationsAsync(IReadOnlyCollection<string> accessions)
        {
            var annotations = accessions.ToDictionary(a => a, _ => (NcbiAnnotation?)null, StringComparer.Ordinal);
            foreach (var chunk in accessions.Chunk(accessionsPerRequest))
            {
                using var response = await http.SendAsync($"{datasetsBase}genome/accession/{string.Join(",", chunk)}/dataset_report?page_size={pageSize}",
                    "application/json");
                using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
                if (!document.RootElement.TryGetProperty("reports", out var reports)) continue;
                foreach (var report in reports.EnumerateArray())
                    if (report.TryGetProperty("annotation_info", out var info) && annotations.ContainsKey(report.GetProperty("accession").GetString()!))
                        annotations[report.GetProperty("accession").GetString()!] = new NcbiAnnotation(
                            info.TryGetProperty("name", out var name) ? name.GetString() : null,
                            info.TryGetProperty("release_date", out var date) ? date.GetString() : null,
                            info.TryGetProperty("provider", out var provider) ? provider.GetString() : null);
            }
            return annotations;
        }

        /// <summary>Each taxon's primary (nuclear) NCBI translation table, from the Datasets taxonomy report.</summary>
        public async Task<Dictionary<int, int>> GeneticCodesAsync(IReadOnlyCollection<int> taxids)
        {
            var codes = new Dictionary<int, int>();
            foreach (var chunk in taxids.Distinct().Chunk(taxaPerRequest))
            {
                using var response = await http.SendAsync(
                    $"{datasetsBase}taxonomy/taxon/{string.Join(",", chunk)}/dataset_report?page_size={taxaPerRequest * 2}", "application/json");
                using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
                if (!document.RootElement.TryGetProperty("reports", out var reports)) continue;
                // A taxid NCBI no longer recognizes comes back as an error report with no taxonomy (147272, 2026-11-09); its code stays
                // unknown, so the record is skipped and stays pending.
                foreach (var report in reports.EnumerateArray())
                {
                    if (!report.TryGetProperty("taxonomy", out var taxonomy)) continue;
                    if (taxonomy.TryGetProperty("genetic_code", out var code) && code.TryGetProperty("primary", out var primary))
                        codes[taxonomy.GetProperty("tax_id").GetInt32()] = primary.GetProperty("id").GetInt32();
                }
            }
            return codes;
        }

        // The assembly's other files sit next to its genome FASTA, with another suffix.
        internal static string SiblingUrl(string fastaUrl, string suffix) =>
            fastaUrl.EndsWith(genomeSuffix, StringComparison.Ordinal) ? fastaUrl[..^genomeSuffix.Length] + suffix
                : throw new ArgumentException($"Not an NCBI genome FASTA URL: {fastaUrl}");

        public async Task<string> ReadTextAsync(string url)
        {
            using var response = await http.SendAsync(url, "text/plain");
            return await response.Content.ReadAsStringAsync();
        }

        public void Dispose() => http.Dispose();
    }

    /// <summary>
    /// The eukaryote pilot: a genome becomes one record holding its nuclear dinucleotide counts and no profile levels, which is enough for
    /// relatives, the map and trees; whether the signature organizes eukaryotes at all is what the pilot is for. Since migration 0022 each
    /// nuclear sequence of at least 100 kbp also keeps its own 500 bp and 32 kbp levels in chromosome_level, for drawing it (12 MB per Gbp),
    /// built in the same pass; they are never compared as one genome. The DNA passes through memory only. Organelles are left out because
    /// their composition is bacterial in origin, and D9 already leaves plasmids out the same way. The human genome stays out until
    /// the consent and GDPR decision is made (decision log, 2026-10-30).
    /// </summary>
    internal sealed class NcbiIngestor
    {
        public const string Provider = "NCBI";
        public const string Kind = "eukaryote";
        private const string dataClass = "public_eukaryote";
        private const string license = "NCBI RefSeq and GenBank public data; cite the original submitters";
        private const int humanTaxid = 9606;
        private const string humanSpecies = "Homo sapiens";
        private const int seed = 20261030;
        private const int bufferChars = 1 << 20;

        internal static readonly Regex organelle = new(@"\b(mitochondri(on|al)|chloroplast|plastid|apicoplast|kinetoplast)\b",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private readonly NpgsqlDataSource dataSource;
        private readonly NcbiGenomeClient? client;

        public NcbiIngestor(NpgsqlDataSource dataSource, NcbiGenomeClient? client)
        {
            this.dataSource = dataSource;
            this.client = client;
        }

        internal static bool IsHuman(NcbiAssembly assembly, Lineage? lineage) => assembly.Taxid == humanTaxid || lineage?.Species == humanSpecies;

        /// <summary>One assembly per family, the families in a fixed random order, so a limit is a spread rather than one clade.</summary>
        internal static List<(NcbiAssembly Assembly, Lineage Lineage)> Pick(IEnumerable<NcbiAssembly> assemblies, IReadOnlyDictionary<int, Lineage> lineages,
            int limit, long maxLengthBp)
        {
            var rng = new Random(seed);
            return assemblies
                .Where(a => a.LengthBp <= maxLengthBp && lineages.ContainsKey(a.Taxid) && !IsHuman(a, lineages[a.Taxid]))
                .OrderBy(a => a.Accession, StringComparer.Ordinal)
                .GroupBy(a => lineages[a.Taxid].Family ?? lineages[a.Taxid].Genus ?? a.Organism, StringComparer.Ordinal)
                .OrderBy(family => family.Key, StringComparer.Ordinal)
                .Select(family => family.OrderBy(_ => rng.Next()).First())
                .OrderBy(_ => rng.Next())
                .Take(limit)
                .Select(a => (a, lineages[a.Taxid]))
                .ToList();
        }

        /// <param name="levels">Also each kept nuclear record's 500 bp level (chromosome_level), from the same pass.</param>
        internal static NuclearCounts CountNuclear(TextReader fasta, bool levels = true)
        {
            var counts = new long[16];
            var buffer = new char[bufferChars];
            var chromosomes = new List<ChromosomeLevels.Chromosome>();
            long lengthBp = 0, shortBp = 0;
            int records = 0, leftOut = 0, tooShort = 0;
            foreach (var (id, description, sequence) in FastaReader.ReadRecordStreams(fasta))
            {
                // Moving on skips what's left of a record, so an organelle costs nothing.
                if (organelle.IsMatch(description)) { leftOut++; continue; }

                using var features = new SequenceFeatureReader(sequence);
                ZCurveWindowedProfile? fine = null;
                if (levels) fine = ZCurveWindowedProfile.FromStream(features);
                else while (features.Read(buffer, 0, buffer.Length) > 0) { }
                var (_, recordCounts) = features.Finish();
                for (int i = 0; i < counts.Length; i++) counts[i] += recordCounts[i];
                lengthBp += features.CharactersRead;
                records++;
                if (fine is null) continue;
                if (features.CharactersRead >= ChromosomeLevels.MinLengthBp) chromosomes.Add(new(id, features.CharactersRead, fine));
                else { tooShort++; shortBp += features.CharactersRead; }
            }
            return new NuclearCounts(counts, lengthBp, records, leftOut, levels ? chromosomes : null, tooShort, shortBp);
        }

        /// <returns>The counts and the SHA-256 of the gzip file as it came.</returns>
        internal static (NuclearCounts Counts, byte[] Sha256) Read(Stream gzipFasta, bool levels = true)
        {
            using var sha256 = SHA256.Create();
            using var hashing = new CryptoStream(gzipFasta, sha256, CryptoStreamMode.Read, leaveOpen: true);
            NuclearCounts counts;
            using (var gzip = new GZipStream(hashing, CompressionMode.Decompress, leaveOpen: true))
            using (var reader = new StreamReader(gzip, bufferSize: bufferChars))
                counts = CountNuclear(reader, levels);
            hashing.CopyTo(Stream.Null);    // the hash covers the whole file, whatever the decompressor left unread
            return (counts, sha256.Hash!);
        }

        /// <summary>Never throws for a per-genome problem: it's recorded in crawl_item, and the batch goes on.</summary>
        public async Task<string> IngestAsync(NcbiAssembly assembly, Lineage lineage)
        {
            try
            {
                var url = NcbiGenomeClient.FastaUrl(assembly);
                using var response = await client!.OpenFastaAsync(url);
                await using var body = await response.Content.ReadAsStreamAsync();
                var (counts, sha256) = await Task.Run(() => Read(body));
                await StoreAsync(assembly, lineage, url, sha256, response.Content.Headers.ContentLength ?? 0, counts);
                return $"{counts.Records,6:N0} nuclear records {counts.LengthBp,15:N0} bp, {counts.OrganellesLeftOut} organelle(s) left out";
            }
            catch (Exception e) when (e is not ArgumentException)
            {
                await CrawlState.MarkAsync(dataSource, Provider, assembly.Accession, CrawlState.Failed, $"{e.GetType().Name}: {e.Message}");
                return $"FAILED: {e.Message}";
            }
        }

        internal async Task StoreAsync(NcbiAssembly assembly, Lineage lineage, string url, byte[] sha256, long byteSize, NuclearCounts counts)
        {
            if (IsHuman(assembly, lineage)) throw new ArgumentException("The human genome is out of the pilot until the consent and GDPR decision.");
            if (counts.Records == 0) throw new InvalidOperationException($"{assembly.Accession} has no nuclear records.");

            double pairs = counts.Counts.Sum();
            double gc = pairs == 0 ? 0 : Enumerable.Range(4, 8).Sum(i => (double)counts.Counts[i]) / pairs;   // first base C or G: indexes 4..11
            var annotations = JsonSerializer.Serialize(new Dictionary<string, object>
            {
                ["ncbi"] = new Dictionary<string, object>
                {
                    ["assembly_accession"] = assembly.Accession,
                    ["assembly_name"] = assembly.AssemblyName,
                    ["nuclear_records"] = counts.Records,
                    ["organelles_left_out"] = counts.OrganellesLeftOut,
                },
            });
            var profile = new ProfiledRecord(assembly.Accession, "", counts.LengthBp, gc, [], null, counts.Counts);
            var row = new RecordRow(assembly.Accession, assembly.Organism, Kind, assembly.Taxid, lineage, annotations, DataClass: dataClass);

            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            var stored = await RecordStore.StoreAsync(connection, transaction,
                new SourceRow(Provider, url, assembly.Accession, sha256, byteSize, license), [(profile, row)]);
            if (counts.Chromosomes is { } chromosomes)
            {
                await ChromosomeLevels.InsertAsync(connection, transaction, stored.SequenceIds[0], chromosomes);
                await ChromosomeLevels.AnnotateAsync(connection, transaction, stored.SequenceIds[0], counts.Summary());
            }
            await CrawlState.MarkStoredAsync(connection, transaction, Provider, assembly.Accession, stored.SourceId, false);
            await transaction.CommitAsync();
        }
    }
}
