using System.IO.Compression;
using System.Security.Cryptography;
using Npgsql;

// The eukaryote genomes the cheap checks read, kept on disk in data/genomes/eukaryotes/ (git-ignored) so a check doesn't download them
// again. Each file is fetched once from the URL the platform recorded, and the genome is kept only if its SHA-256 is the one the platform
// stored when it profiled it (the CDS file against the gene layer's cds_sha256). These are research files on this machine only: the
// platform itself still never stores DNA. Checks train on one genome of a pair and test on its relative, so they come in pairs.
internal sealed record EukaryoteGenome(long Id, string Label, string Fasta, string? Cds, string? Report);

internal static class EukaryoteFiles
{
    /// <summary>A compact fungus pair (35 Mbp each, genes ~1.5 kbp) and a plant pair (~820 Mbp each, long introns), as of 2026-11-11.</summary>
    public static readonly (long Id, string Label)[][] Pairs =
    [
        [(28677, "Podospora bellae-mahoneyi"), (31652, "Podospora pseudocomata")],
        [(61826, "Nelumbo nucifera"), (53150, "Nelumbo lutea")],
    ];

    private const string genomeSuffix = "_genomic.fna.gz";
    private static readonly HttpClient http = CreateClient();

    public static string Folder
    {
        get
        {
            var directory = new DirectoryInfo(Directory.GetCurrentDirectory());
            while (directory is not null && !Directory.Exists(Path.Combine(directory.FullName, ".git"))) directory = directory.Parent;
            return Path.Combine(directory?.FullName ?? Directory.GetCurrentDirectory(), "data", "genomes", "eukaryotes");
        }
    }

    /// <summary>The genome's local files, downloaded first if missing: genome FASTA, CDS FASTA and assembly report (null without a gene layer).</summary>
    public static async Task<EukaryoteGenome> GetAsync(NpgsqlDataSource dataSource, long id, string label, Action<string>? log = null)
    {
        string fastaUrl; byte[] fastaSha256; string? cdsUrl, cdsSha256;
        await using (var command = dataSource.CreateCommand(
            "select s.url, s.raw_sha256, r.annotations #>> '{ncbi,gene_layer,source}', r.annotations #>> '{ncbi,gene_layer,cds_sha256}' " +
            "from sequence_record r join sequence_source s on s.id = r.source_id where r.id = $1"))
        {
            command.Parameters.AddWithValue(id);
            await using var reader = await command.ExecuteReaderAsync();
            if (!await reader.ReadAsync()) throw new InvalidOperationException($"No record {id}.");
            fastaUrl = reader.GetString(0);
            fastaSha256 = reader.GetFieldValue<byte[]>(1);
            cdsUrl = reader.IsDBNull(2) ? null : reader.GetString(2);
            cdsSha256 = reader.IsDBNull(3) ? null : reader.GetString(3);
        }
        Directory.CreateDirectory(Folder);
        var fasta = await FetchAsync(fastaUrl, fastaSha256, log);
        string? cds = null, report = null;
        if (cdsUrl is not null)
        {
            cds = await FetchAsync(cdsUrl, cdsSha256 is null ? null : Convert.FromHexString(cdsSha256), log);
            report = await FetchAsync(fastaUrl[..^genomeSuffix.Length] + "_assembly_report.txt", null, log);
        }
        return new EukaryoteGenome(id, label, fasta, cds, report);
    }

    /// <summary>A local file as text, gunzipped when it ends in .gz.</summary>
    public static TextReader Open(string path)
    {
        Stream stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20);
        if (path.EndsWith(".gz", StringComparison.Ordinal)) stream = new GZipStream(stream, CompressionMode.Decompress);
        return new StreamReader(stream, bufferSize: 1 << 20);
    }

    private static async Task<string> FetchAsync(string url, byte[]? sha256, Action<string>? log)
    {
        var path = Path.Combine(Folder, url[(url.LastIndexOf('/') + 1)..]);
        if (File.Exists(path)) return path;
        var partial = path + ".part";
        log?.Invoke($"downloading {url}");
        using (var response = await http.GetAsync(url, HttpCompletionOption.ResponseHeadersRead))
        {
            response.EnsureSuccessStatusCode();
            await using var body = await response.Content.ReadAsStreamAsync();
            await using var file = File.Create(partial);
            await body.CopyToAsync(file);
        }
        if (sha256 is not null)
        {
            await using var file = File.OpenRead(partial);
            if (!SHA256.HashData(file).AsSpan().SequenceEqual(sha256))
            {
                file.Close();
                File.Delete(partial);
                throw new InvalidDataException($"{url} isn't the file the platform profiled (SHA-256 differs); not kept.");
            }
        }
        File.Move(partial, path);
        return path;
    }

    private static HttpClient CreateClient()
    {
        var client = new HttpClient { Timeout = TimeSpan.FromHours(2) };
        client.DefaultRequestHeaders.UserAgent.ParseAdd("sequence-platform-check/1.0");
        return client;
    }
}
