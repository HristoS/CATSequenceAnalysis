using System.Net;
using System.Security.Cryptography;
using System.Text.Json;

namespace SequencePlatform.Ingest
{
    internal sealed record BvBrcGenome(string GenomeId, string Name, JsonElement Raw)
    {
        public string? Field(string name) =>
            Raw.TryGetProperty(name, out var value) && value.ValueKind is not JsonValueKind.Null
                ? value.ValueKind == JsonValueKind.String ? value.GetString() : value.GetRawText()
                : null;

        public IReadOnlyList<string> LineageNames =>
            Raw.TryGetProperty("taxon_lineage_names", out var names) && names.ValueKind == JsonValueKind.Array
                ? names.EnumerateArray().Select(n => n.GetString() ?? "").ToList()
                : Array.Empty<string>();
    }

    internal sealed record BvBrcContig(string Accession, string? SequenceType, string? Topology, string? Description, long? Length,
        string? SequenceMd5);

    internal sealed record DownloadedFile(string Path, string Url, byte[] Sha256, long ByteSize);

    /// <summary>
    /// Minimal client for the BV-BRC data API (https://www.bv-brc.org/api/). Deliberately polite:
    /// requests are sequential and at least <see cref="minRequestInterval"/> apart, and 429/5xx
    /// responses are retried with exponential backoff that honours Retry-After. BV-BRC publishes
    /// no rate limit, so this stays well under anything a public service would object to.
    /// FASTA comes from the API (Accept: application/dna+fasta) because the FTP host isn't
    /// reachable from every network.
    /// </summary>
    internal sealed class BvBrcClient : IDisposable
    {
        private const string apiBase = "https://www.bv-brc.org/api/";
        private const int maxAttempts = 5;
        private static readonly TimeSpan minRequestInterval = TimeSpan.FromSeconds(1);

        private static readonly string[] genomeFields =
        {
            "genome_id", "genome_name", "taxon_id", "genome_length", "contigs", "genome_status", "reference_genome",
            "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "taxon_lineage_names",
            "strain", "serovar", "host_name", "disease", "isolation_source", "isolation_country", "collection_year",
            "gram_stain", "genome_quality", "checkm_completeness", "checkm_contamination", "gc_content",
            "assembly_accession", "refseq_accessions", "genbank_accessions", "bioproject_accession", "biosample_accession",
        };

        private readonly HttpClient http;
        private DateTime lastRequestUtc = DateTime.MinValue;

        public BvBrcClient(string? contactEmail)
        {
            http = new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
            var agent = "SequencePlatform.Ingest/0.1 (research test" + (contactEmail is null ? "" : $"; mailto:{contactEmail}") + ")";
            http.DefaultRequestHeaders.UserAgent.ParseAdd(agent);
        }

        public static string GroupFilter(string group) => group switch
        {
            "bacteria" => "and(eq(superkingdom,Bacteria),in(reference_genome,(Reference,Representative)),eq(genome_status,Complete))",
            "archaea" => "and(eq(superkingdom,Archaea),in(reference_genome,(Reference,Representative)),eq(genome_status,Complete))",
            "viruses" => "and(eq(superkingdom,Viruses),eq(reference_genome,Reference),eq(genome_status,Complete))",
            _ => throw new ArgumentException($"Unknown group '{group}'. Use 'bacteria', 'archaea' or 'viruses'.", nameof(group)),
        };

        public async Task<List<BvBrcGenome>> ListGenomesAsync(string group, int limit, int offset)
        {
            var url = $"{apiBase}genome/?{GroupFilter(group)}&select({string.Join(",", genomeFields)})&sort(+genome_id)&limit({limit},{offset})";
            return await QueryGenomesAsync(url);
        }

        // Metadata for known genome ids (retrying failed crawl items), in chunks that keep the URL short.
        public async Task<List<BvBrcGenome>> GetGenomesAsync(IReadOnlyList<string> genomeIds)
        {
            const int chunkSize = 50;
            var genomes = new List<BvBrcGenome>();
            for (int i = 0; i < genomeIds.Count; i += chunkSize)
            {
                var chunk = genomeIds.Skip(i).Take(chunkSize).ToList();
                var url = $"{apiBase}genome/?in(genome_id,({string.Join(",", chunk)}))&select({string.Join(",", genomeFields)})&sort(+genome_id)&limit({chunk.Count})";
                genomes.AddRange(await QueryGenomesAsync(url));
            }
            return genomes;
        }

        private async Task<List<BvBrcGenome>> QueryGenomesAsync(string url)
        {
            using var response = await SendAsync(url, "application/json");
            using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
            return document.RootElement.EnumerateArray()
                .Select(g => new BvBrcGenome(g.GetProperty("genome_id").GetString()!, g.GetProperty("genome_name").GetString() ?? "", g.Clone()))
                .ToList();
        }

        public async Task<List<BvBrcContig>> ListContigsAsync(string genomeId)
        {
            var url = $"{apiBase}genome_sequence/?eq(genome_id,{genomeId})&select(accession,sequence_type,topology,description,length,sequence_md5)&limit(10000)";
            using var response = await SendAsync(url, "application/json");
            using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
            string? Text(JsonElement e, string name) => e.TryGetProperty(name, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString() : null;
            return document.RootElement.EnumerateArray()
                .Select(c => new BvBrcContig(
                    Text(c, "accession") ?? "",
                    Text(c, "sequence_type"),
                    Text(c, "topology"),
                    Text(c, "description"),
                    c.TryGetProperty("length", out var l) && l.ValueKind == JsonValueKind.Number ? l.GetInt64() : null,
                    Text(c, "sequence_md5")))
                .ToList();
        }

        // Streams the genome's FASTA to a temp file, hashing as it goes. The caller must delete the file.
        public async Task<DownloadedFile> DownloadFastaAsync(string genomeId)
        {
            var url = $"{apiBase}genome_sequence/?eq(genome_id,{genomeId})&limit(10000)";
            using var response = await SendAsync(url, "application/dna+fasta", HttpCompletionOption.ResponseHeadersRead);

            var path = Path.Combine(Path.GetTempPath(), $"bvbrc_{genomeId}_{Guid.NewGuid():N}.fna");
            using var sha = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
            long size = 0;
            await using (var source = await response.Content.ReadAsStreamAsync())
            await using (var target = File.Create(path))
            {
                var buffer = new byte[81920];
                int read;
                while ((read = await source.ReadAsync(buffer)) > 0)
                {
                    sha.AppendData(buffer, 0, read);
                    await target.WriteAsync(buffer.AsMemory(0, read));
                    size += read;
                }
            }
            return new DownloadedFile(path, url, sha.GetHashAndReset(), size);
        }

        private async Task<HttpResponseMessage> SendAsync(string url, string accept, HttpCompletionOption completion = HttpCompletionOption.ResponseContentRead)
        {
            for (int attempt = 1; ; attempt++)
            {
                var wait = lastRequestUtc + minRequestInterval - DateTime.UtcNow;
                if (wait > TimeSpan.Zero) await Task.Delay(wait);
                lastRequestUtc = DateTime.UtcNow;

                using var request = new HttpRequestMessage(HttpMethod.Get, url);
                request.Headers.Accept.ParseAdd(accept);

                HttpResponseMessage? response = null;
                try
                {
                    response = await http.SendAsync(request, completion);
                    if (response.IsSuccessStatusCode) return response;
                    // Other 4xx won't succeed on retry; a distinct exception type keeps it out of the retry catch below.
                    if (response.StatusCode != HttpStatusCode.TooManyRequests && (int)response.StatusCode < 500)
                        throw new InvalidOperationException($"BV-BRC returned {(int)response.StatusCode} for {url}");
                }
                catch (Exception e) when (e is HttpRequestException or TaskCanceledException && attempt < maxAttempts)
                {
                    Console.WriteLine($"  request failed ({e.Message}), retrying");
                }

                if (attempt >= maxAttempts)
                    throw new HttpRequestException($"BV-BRC request failed after {maxAttempts} attempts: {url}");

                var backoff = response?.Headers.RetryAfter?.Delta ?? TimeSpan.FromSeconds(Math.Pow(2, attempt));
                response?.Dispose();
                await Task.Delay(backoff);
            }
        }

        public void Dispose() => http.Dispose();
    }
}
