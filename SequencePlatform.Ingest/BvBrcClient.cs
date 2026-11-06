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

    /// <summary>One coding sequence on a contig: [Start, End) 0-based, the stop codon included as BV-BRC gives it.</summary>
    internal sealed record BvBrcCds(string Accession, int Start, int End, bool Reverse, string? FamilyId, string? Product);

    internal sealed record DownloadedFile(string Path, string Url, byte[] Sha256, long ByteSize);

    /// <summary>
    /// Minimal client for the BV-BRC data API (https://www.bv-brc.org/api/). Deliberately polite
    /// (<see cref="PoliteHttpClient"/>): requests are sequential, one in flight at a time, and 429/5xx responses
    /// are retried with backoff. BV-BRC publishes no rate limit. There is no gap between requests (it was 1 s
    /// until 2026-11-04): with lookups at about 0.5 s and downloads at about 0.65 s, the gap had become half of
    /// each species' time, while one request at a time is no more load than a single browser tab.
    /// FASTA comes from the API (Accept: application/dna+fasta) because the FTP host isn't
    /// reachable from every network.
    /// </summary>
    internal sealed class BvBrcClient : IDisposable
    {
        private const string apiBase = "https://www.bv-brc.org/api/";
        private static readonly TimeSpan minRequestInterval = TimeSpan.Zero;

        private static readonly string[] genomeFields =
        {
            "genome_id", "genome_name", "taxon_id", "genome_length", "contigs", "genome_status", "reference_genome",
            "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "taxon_lineage_names",
            "strain", "serovar", "host_name", "disease", "isolation_source", "isolation_country", "collection_year",
            "gram_stain", "genome_quality", "checkm_completeness", "checkm_contamination", "gc_content",
            "assembly_accession", "refseq_accessions", "genbank_accessions", "bioproject_accession", "biosample_accession",
        };

        private readonly PoliteHttpClient http;

        public BvBrcClient(string? contactEmail)
        {
            http = new PoliteHttpClient("BV-BRC", contactEmail, minRequestInterval);
        }

        public static string GroupFilter(string group) => group switch
        {
            "bacteria" => "and(eq(superkingdom,Bacteria),in(reference_genome,(Reference,Representative)),eq(genome_status,Complete))",
            "archaea" => "and(eq(superkingdom,Archaea),in(reference_genome,(Reference,Representative)),eq(genome_status,Complete))",
            "viruses" => "and(eq(superkingdom,Viruses),eq(reference_genome,Reference),eq(genome_status,Complete))",
            _ => throw new ArgumentException($"Unknown group '{group}'. Use 'bacteria', 'archaea' or 'viruses'.", nameof(group)),
        };

        // Every complete genome of a group, reference or not: the pool --per-species picks from. Most are repeats of a few species
        // (10.3 M viruses, 71,364 bacteria, 2026-10-30), so nothing lists it whole for viruses.
        private static string CompleteFilter(string group) => group switch
        {
            "bacteria" => "and(eq(superkingdom,Bacteria),eq(genome_status,Complete))",
            "archaea" => "and(eq(superkingdom,Archaea),eq(genome_status,Complete))",
            "viruses" => "and(eq(superkingdom,Viruses),eq(genome_status,Complete))",
            _ => throw new ArgumentException($"Unknown group '{group}'. Use 'bacteria', 'archaea' or 'viruses'.", nameof(group)),
        };

        public async Task<List<BvBrcGenome>> ListGenomesAsync(string group, int limit, int offset)
        {
            var url = $"{apiBase}genome/?{GroupFilter(group)}&select({string.Join(",", genomeFields)})&sort(+genome_id)&limit({limit},{offset})";
            return await QueryGenomesAsync(url);
        }

        public async Task<List<BvBrcGenome>> ListCompleteGenomesAsync(string group, int limit, int offset)
        {
            var url = $"{apiBase}genome/?{CompleteFilter(group)}&select({string.Join(",", genomeFields)})&sort(+genome_id)&limit({limit},{offset})";
            return await QueryGenomesAsync(url);
        }

        // One request, however many genomes: the facet counts them on the server.
        public async Task<List<string>> ListSpeciesAsync(string group)
        {
            var url = $"{apiBase}genome/?{CompleteFilter(group)}&facet((field,species),(mincount,1),(limit,-1))&limit(1)&select(genome_id)";
            using var response = await SendAsync(url, "application/solr+json");
            using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
            return document.RootElement.GetProperty("facet_counts").GetProperty("facet_fields").GetProperty("species").EnumerateArray()
                .Where(v => v.ValueKind == JsonValueKind.String)    // name, count, name, count, ...
                .Select(v => v.GetString()!)
                .Where(name => name.Length > 0)
                .ToList();
        }

        public async Task<BvBrcGenome?> FirstGenomeOfSpeciesAsync(string group, string species)
        {
            var filter = CompleteFilter(group);
            var withSpecies = $"{filter[..^1]},eq(species,{Uri.EscapeDataString($"\"{species}\"")}))";
            var url = $"{apiBase}genome/?{withSpecies}&select({string.Join(",", genomeFields)})&sort(+genome_id)&limit(1)";
            return (await QueryGenomesAsync(url)).FirstOrDefault();
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

        public async Task<List<BvBrcContig>> ListContigsAsync(string genomeId) => (await ListContigsAsync(new[] { genomeId }))[genomeId];

        // One request per chunk of genomes instead of one per genome: at the 1 s request spacing that saves about 1 s per
        // genome, roughly 6 h over the D1 corpus. Every requested id gets an entry, empty when BV-BRC lists no contigs for it.
        public async Task<Dictionary<string, List<BvBrcContig>>> ListContigsAsync(IReadOnlyCollection<string> genomeIds)
        {
            const int chunkSize = 50;
            const int rowLimit = 25000;
            var contigs = genomeIds.Distinct(StringComparer.Ordinal).ToDictionary(id => id, _ => new List<BvBrcContig>(), StringComparer.Ordinal);
            var ids = contigs.Keys.ToList();
            for (int i = 0; i < ids.Count; i += chunkSize)
            {
                var chunk = ids.Skip(i).Take(chunkSize).ToList();
                var url = $"{apiBase}genome_sequence/?in(genome_id,({string.Join(",", chunk)}))" +
                          $"&select(genome_id,accession,sequence_type,topology,description,length,sequence_md5)&limit({rowLimit})";
                using var response = await SendAsync(url, "application/json");
                using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
                if (document.RootElement.GetArrayLength() >= rowLimit)
                    throw new InvalidOperationException($"Contig listing for {chunk.Count} genomes hit the {rowLimit}-row limit and may be truncated.");

                string? Text(JsonElement e, string name) => e.TryGetProperty(name, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString() : null;
                foreach (var c in document.RootElement.EnumerateArray())
                    if (Text(c, "genome_id") is { } genomeId && contigs.TryGetValue(genomeId, out var list))
                        list.Add(new BvBrcContig(
                            Text(c, "accession") ?? "",
                            Text(c, "sequence_type"),
                            Text(c, "topology"),
                            Text(c, "description"),
                            c.TryGetProperty("length", out var l) && l.ValueKind == JsonValueKind.Number ? l.GetInt64() : null,
                            Text(c, "sequence_md5")));
            }
            return contigs;
        }

        // Coding sequences (PATRIC annotation) of a chunk of genomes per request, paged by feature id so no genome is truncated. BV-BRC
        // coordinates are 1-based and inclusive, start below end on both strands, and include the stop codon. Every requested id gets
        // an entry, empty when BV-BRC lists no CDS for it.
        public async Task<Dictionary<string, List<BvBrcCds>>> ListCodingSequencesAsync(IReadOnlyCollection<string> genomeIds)
        {
            const int chunkSize = 10;
            const int rowLimit = 25000;
            var result = genomeIds.Distinct(StringComparer.Ordinal).ToDictionary(id => id, _ => new List<BvBrcCds>(), StringComparer.Ordinal);
            var ids = result.Keys.ToList();
            for (int i = 0; i < ids.Count; i += chunkSize)
            {
                var chunk = string.Join(",", ids.Skip(i).Take(chunkSize));
                for (int offset = 0; ; offset += rowLimit)
                {
                    var url = $"{apiBase}genome_feature/?and(in(genome_id,({chunk})),eq(feature_type,CDS),eq(annotation,PATRIC))" +
                              $"&select(genome_id,accession,start,end,strand,pgfam_id,product)&sort(+feature_id)&limit({rowLimit},{offset})";
                    using var response = await SendAsync(url, "application/json");
                    using var document = JsonDocument.Parse(await response.Content.ReadAsStringAsync());
                    string? Text(JsonElement e, string name) => e.TryGetProperty(name, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString() : null;
                    int rows = 0;
                    foreach (var f in document.RootElement.EnumerateArray())
                    {
                        rows++;
                        if (Text(f, "genome_id") is not { } genomeId || !result.TryGetValue(genomeId, out var list)) continue;
                        if (!f.TryGetProperty("start", out var start) || !f.TryGetProperty("end", out var end)) continue;
                        list.Add(new BvBrcCds(Text(f, "accession") ?? "", start.GetInt32() - 1, end.GetInt32(), Text(f, "strand") == "-",
                            Text(f, "pgfam_id"), Text(f, "product")));
                    }
                    if (rows < rowLimit) break;
                }
            }
            return result;
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

        private Task<HttpResponseMessage> SendAsync(string url, string accept, HttpCompletionOption completion = HttpCompletionOption.ResponseContentRead) =>
            http.SendAsync(url, accept, completion);

        public void Dispose() => http.Dispose();
    }
}
