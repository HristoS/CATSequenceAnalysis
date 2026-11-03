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

    /// <summary>What a genome's nuclear records add up to.</summary>
    /// <param name="Counts">Dinucleotide counts summed over the nuclear records; a pair spanning two records isn't counted.</param>
    internal sealed record NuclearCounts(long[] Counts, long LengthBp, int Records, int OrganellesLeftOut);

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

        public void Dispose() => http.Dispose();
    }

    /// <summary>
    /// The eukaryote pilot: a genome becomes one record holding its nuclear dinucleotide counts and no profile levels, which is enough for
    /// relatives, the map and trees. A genome would need 32 MB per Gbp of levels and object storage above 43 Mbp per chromosome, and whether the
    /// signature organizes eukaryotes at all is what the pilot is for. The DNA passes through memory only. Organelles are left out because
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

        private static readonly Regex organelle = new(@"\b(mitochondri(on|al)|chloroplast|plastid|apicoplast|kinetoplast)\b",
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

        internal static NuclearCounts CountNuclear(TextReader fasta)
        {
            var counts = new long[16];
            var buffer = new char[bufferChars];
            long lengthBp = 0;
            int records = 0, leftOut = 0;
            foreach (var (_, description, sequence) in FastaReader.ReadRecordStreams(fasta))
            {
                // Moving on skips what's left of a record, so an organelle costs nothing.
                if (organelle.IsMatch(description)) { leftOut++; continue; }

                using var features = new SequenceFeatureReader(sequence);
                while (features.Read(buffer, 0, buffer.Length) > 0) { }
                var (_, recordCounts) = features.Finish();
                for (int i = 0; i < counts.Length; i++) counts[i] += recordCounts[i];
                lengthBp += features.CharactersRead;
                records++;
            }
            return new NuclearCounts(counts, lengthBp, records, leftOut);
        }

        /// <returns>The counts and the SHA-256 of the gzip file as it came.</returns>
        internal static (NuclearCounts Counts, byte[] Sha256) Read(Stream gzipFasta)
        {
            using var sha256 = SHA256.Create();
            using var hashing = new CryptoStream(gzipFasta, sha256, CryptoStreamMode.Read, leaveOpen: true);
            NuclearCounts counts;
            using (var gzip = new GZipStream(hashing, CompressionMode.Decompress, leaveOpen: true))
            using (var reader = new StreamReader(gzip, bufferSize: bufferChars))
                counts = CountNuclear(reader);
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
            await CrawlState.MarkStoredAsync(connection, transaction, Provider, assembly.Accession, stored.SourceId, false);
            await transaction.CommitAsync();
        }
    }
}
