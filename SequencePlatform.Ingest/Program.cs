using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

// Small test app: populate the platform DB from BV-BRC. The connection string comes from the
// SEQUENCE_PLATFORM_DB environment variable, or else the git-ignored appsettings.json (never committed).
//
//   migrate                                   create the database if missing, apply migrations
//   ingest-bvbrc --group bacteria|archaea|viruses --limit N [--offset K] [--contact you@example.org] [--retry-failed]
//   load-taxonomy [--refresh] [--skip-gtdb]   NCBI's current lineage per referenced taxid, GTDB's per assembly (S2.4)
//   stats                                     what's in the DB, including crawl state
//   list [--kind K] [--rank R --taxon T] [--min-length N] [--max-length N] [--provider P] [--name TEXT] [--after ID] [--limit N]
//   profile --fasta X --out X.zcwp            the query-side profile of each record, both strands (X.zcwp and X.rc.zcwp); no database
//   search-similar --fasta X | --profile X.zcwp [list's filters] [--ratio R] [--top N]   records of comparable length (ratio 0.25-1, default 0.5) ranked by similarity (P3 S3.4)
const string usage =
    "usage: SequencePlatform.Ingest migrate | ingest-bvbrc --group bacteria|archaea|viruses --limit N [--offset K] [--contact EMAIL] [--retry-failed] | load-taxonomy [--refresh] [--skip-gtdb] | stats\n" +
    "       | list [--kind K] [--rank R --taxon T] [--min-length N] [--max-length N] [--provider P] [--name TEXT] [--after ID] [--limit N]\n" +
    "       | profile --fasta X --out X.zcwp | search-similar --fasta X | --profile X.zcwp [list's filters] [--ratio R] [--top N]";

if (args.Length == 0)
{
    Console.WriteLine(usage);
    return 1;
}

string? Option(string name)
{
    var index = Array.IndexOf(args, name);
    return index >= 0 && index + 1 < args.Length ? args[index + 1] : null;
}

// Profiling a query needs no database, so it runs before the connection is opened.
if (args[0] == "profile")
    return WriteQueryProfiles(Option("--fasta"), Option("--out"));

CatalogueFilter Filter() => new(
    Option("--kind"), Option("--rank"), Option("--taxon"),
    Option("--min-length") is { } min ? long.Parse(min) : null,
    Option("--max-length") is { } max ? long.Parse(max) : null,
    Option("--provider"), Option("--name"));

var settings = Database.ConnectionSettings(ConnectionString.Resolve());
await Database.EnsureDatabaseAsync(settings);
await using var dataSource = NpgsqlDataSource.Create(settings.ConnectionString);
await Database.MigrateAsync(dataSource);

switch (args[0])
{
    case "migrate":
        Console.WriteLine($"Database {settings.Database} is up to date.");
        return 0;

    case "ingest-bvbrc":
        var retryFailed = args.Contains("--retry-failed");
        return await IngestBvBrcAsync(
            dataSource,
            Option("--group") ?? (retryFailed ? null : "bacteria"),
            int.Parse(Option("--limit") ?? "10"),
            int.Parse(Option("--offset") ?? "0"),
            Option("--contact"),
            retryFailed);

    case "load-taxonomy":
        return await LoadTaxonomyAsync(dataSource, args.Contains("--refresh"), args.Contains("--skip-gtdb"), Option("--contact"));

    case "stats":
        await PrintStatsAsync(dataSource);
        return 0;

    case "list":
        await PrintListAsync(dataSource, Filter(), long.Parse(Option("--after") ?? "0"), int.Parse(Option("--limit") ?? "50"));
        return 0;

    case "search-similar":
        return await SearchSimilarAsync(dataSource, Option("--fasta"), Option("--profile"), Filter(),
            double.Parse(Option("--ratio") ?? SimilaritySearch.DefaultLengthRatio.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture),
            int.Parse(Option("--top") ?? "20"));

    default:
        Console.WriteLine(usage);
        return 1;
}

// Normal run: list a page of the group, register new items, process everything not in a final state
// except failed. --retry-failed: process failed items instead (of --group, or of every group if omitted).
static async Task<int> IngestBvBrcAsync(NpgsqlDataSource dataSource, string? group, int limit, int offset, string? contact, bool retryFailed)
{
    using var client = new BvBrcClient(contact);
    var ingestor = new GenomeIngestor(dataSource, client);
    var runStarted = DateTime.UtcNow;

    // A run killed mid-genome (Ctrl+C, reboot) skips the finally that deletes its download. The no-sequence rule says
    // nothing downloaded outlives ingestion, so sweep what an earlier run left. Only files over an hour old, so a
    // second run started alongside keeps its in-flight download.
    foreach (var leftover in Directory.EnumerateFiles(Path.GetTempPath(), "bvbrc_*.fna"))
        if (File.GetLastWriteTimeUtc(leftover) < runStarted.AddHours(-1))
        {
            File.Delete(leftover);
            Console.WriteLine($"  deleted leftover download {Path.GetFileName(leftover)}");
        }

    List<BvBrcGenome> genomes;
    if (retryFailed)
    {
        var ids = await CrawlState.FailedAsync(dataSource, GenomeIngestor.Provider, group, limit);
        genomes = ids.Count == 0 ? new() : await client.GetGenomesAsync(ids);
        Console.WriteLine($"BV-BRC {group ?? "all groups"}: retrying {genomes.Count} of {ids.Count} failed items.");
    }
    else
    {
        genomes = await client.ListGenomesAsync(group!, limit, offset);
        Console.WriteLine($"BV-BRC {group}: {genomes.Count} genomes listed (limit {limit}, offset {offset}).");
    }
    var statuses = await CrawlState.RegisterAsync(dataSource, GenomeIngestor.Provider, retryFailed ? null : group,
        genomes.Select(g => g.GenomeId).ToList());

    var pending = genomes.Where(g => !CrawlState.IsFinal(statuses[g.GenomeId]) && (retryFailed || statuses[g.GenomeId] != CrawlState.Failed)).ToList();
    Dictionary<string, List<BvBrcContig>>? listedContigs = null;
    try
    {
        listedContigs = await client.ListContigsAsync(pending.Select(g => g.GenomeId).ToList());
    }
    catch (Exception e) when (e is HttpRequestException or InvalidOperationException or TaskCanceledException)
    {
        Console.WriteLine($"  batched contig listing failed ({e.Message}); listing per genome instead");
    }

    var counts = new Dictionary<string, int>();
    int alreadyFinal = 0, failedSkipped = 0;
    foreach (var genome in genomes)
    {
        var status = statuses[genome.GenomeId];
        if (CrawlState.IsFinal(status)) { alreadyFinal++; continue; }
        if (status == CrawlState.Failed && !retryFailed) { failedSkipped++; continue; }

        var started = DateTime.UtcNow;
        var outcome = await ingestor.IngestAsync(genome, listedContigs?[genome.GenomeId]);
        counts[outcome.Status] = counts.GetValueOrDefault(outcome.Status) + 1;
        var detail = outcome.Status switch
        {
            CrawlState.Failed => $"FAILED: {outcome.Error}",
            CrawlState.Duplicate => $"duplicate of source {outcome.DuplicateOf} ({outcome.DuplicateBy}), skipped",
            _ => $"{outcome.Records,3} records {outcome.Levels,3} levels {outcome.Bases,12:N0} bp" +
                 (outcome.AlreadyStored > 0 ? $"  ({outcome.AlreadyStored} of {outcome.Records} records already stored)" : "") +
                 (outcome.Status == CrawlState.TooLarge ? "  (too large: a 500 bp level needs object storage)" : ""),
        };
        Console.WriteLine($"  {genome.GenomeId,-14} {genome.Name,-60} {detail}  {(DateTime.UtcNow - started).TotalSeconds,5:F1} s");
    }

    int processed = counts.Values.Sum();
    var elapsed = DateTime.UtcNow - runStarted;
    Console.WriteLine($"Done: {string.Join(", ", counts.Select(c => $"{c.Value} {c.Key}"))}{(processed == 0 ? "nothing to do" : "")}; " +
        $"{alreadyFinal} already handled, {failedSkipped} failed earlier (use --retry-failed). " +
        $"{elapsed.TotalSeconds:F0} s{(processed > 0 ? $", {elapsed.TotalSeconds / processed:F2} s per processed genome" : "")}.");
    return counts.ContainsKey(CrawlState.Failed) ? 2 : 0;
}

static async Task<int> LoadTaxonomyAsync(NpgsqlDataSource dataSource, bool refresh, bool skipGtdb, string? contact)
{
    using var loader = new TaxonomyLoader(dataSource, contact);
    var ncbi = await loader.LoadNcbiAsync(refresh);
    var unanswered = ncbi.Asked - ncbi.Current - ncbi.Merged - ncbi.Unknown;
    Console.WriteLine($"NCBI: {ncbi.Asked:N0} taxids asked: {ncbi.Current:N0} current, {ncbi.Merged:N0} merged into another, " +
        $"{ncbi.Unknown:N0} unknown, {unanswered:N0} unanswered (a re-run asks again).");
    // The S2.4 card's stop criterion: if the API can't resolve this many, the full dump with its merge tables is needed.
    if (ncbi.Asked > 0 && (ncbi.Unknown + unanswered) > 0.05 * ncbi.Asked)
        Console.WriteLine("  WARNING: over 5% of taxids unresolved; consider the full new_taxdump instead.");

    if (!skipGtdb)
    {
        var gtdb = await loader.LoadGtdbAsync();
        Console.WriteLine($"GTDB {gtdb.Release}: {gtdb.Matched:N0} of {gtdb.Referenced:N0} referenced bacterial and archaeal assemblies listed.");
    }

    await using var command = dataSource.CreateCommand(
        "select rank, outcome, count(*) from lineage_difference group by 1, 2 order by min(case rank when 'superkingdom' then 0 " +
        "when 'realm' then 1 when 'kingdom' then 2 when 'phylum' then 3 when 'class' then 4 when 'order' then 5 when 'family' then 6 " +
        "when 'genus' then 7 else 8 end), 2");
    await using var reader = await command.ExecuteReaderAsync();
    Console.WriteLine();
    Console.WriteLine($"{"rank",-13} {"outcome",-14} {"records",8}   (provider lineage vs NCBI's current one; view lineage_difference)");
    while (await reader.ReadAsync())
        Console.WriteLine($"{reader.GetString(0),-13} {reader.GetString(1),-14} {reader.GetInt64(2),8:N0}");
    return 0;
}

// One output pair per FASTA record: <out> and its .rc twin when the file has one record, <out minus .zcwp>.<record id>.zcwp
// otherwise. The reverse strand is profiled from the reverse-complement DNA (exact), so a profile-only search keeps both strands.
static int WriteQueryProfiles(string? fasta, string? output)
{
    if (fasta is null || output is null)
    {
        Console.WriteLine("usage: profile --fasta X --out X.zcwp");
        return 1;
    }

    var headers = FastaReader.ReadHeaders(fasta).ToList();
    foreach (var (recordId, _) in headers)
    {
        using var reader = FastaReader.OpenSequenceReader(fasta, recordId);
        var query = QueryProfile.FromSequence(reader.ReadToEnd());
        var path = headers.Count == 1 ? output : $"{Path.ChangeExtension(output, null)}.{recordId}.zcwp";
        File.WriteAllBytes(path, query.Forward.Fine.ToBinary());
        File.WriteAllBytes(ReverseStrandPath(path), query.Reverse!.Fine.ToBinary());
        Console.WriteLine($"{recordId}: {query.Forward.Fine.WindowCount:N0} windows -> {path} and {ReverseStrandPath(path)}");
    }
    return 0;
}

static string ReverseStrandPath(string path) => Path.ChangeExtension(path, null) + ".rc.zcwp";

static async Task<int> SearchSimilarAsync(NpgsqlDataSource dataSource, string? fasta, string? profile, CatalogueFilter filter, double ratio, int top)
{
    var queries = new List<(string Name, QueryProfile Query)>();
    if (fasta is not null)
        foreach (var (recordId, _) in FastaReader.ReadHeaders(fasta))
        {
            using var reader = FastaReader.OpenSequenceReader(fasta, recordId);
            queries.Add((recordId, QueryProfile.FromSequence(reader.ReadToEnd())));
        }
    else if (profile is not null)
    {
        var reversePath = ReverseStrandPath(profile);
        var reverse = File.Exists(reversePath) ? ZCurveWindowedProfile.FromBinary(File.ReadAllBytes(reversePath)) : null;
        if (reverse is null) Console.WriteLine($"No {reversePath}: searching the forward strand only.");
        queries.Add((Path.GetFileName(profile), QueryProfile.FromProfiles(ZCurveWindowedProfile.FromBinary(File.ReadAllBytes(profile)), reverse)));
    }
    else
    {
        Console.WriteLine("usage: search-similar --fasta X | --profile X.zcwp [list's filters] [--ratio R] [--top N]");
        return 1;
    }

    // One process, one or a few queries: reading each band from the database beats loading a cache.
    var search = new SimilaritySearch(dataSource, useCache: false);
    foreach (var (name, query) in queries)
    {
        var watch = System.Diagnostics.Stopwatch.StartNew();
        var result = await search.SearchAsync(query, filter, ratio, top);
        var elapsed = watch.Elapsed;
        var records = await Catalogue.ByIdsAsync(dataSource,
            result.Ranked.Select(h => h.SequenceId).Concat(result.Identical.Select(i => i.SequenceId)).ToHashSet());

        Console.WriteLine();
        Console.WriteLine($"Query {name}: {query.LengthBp:N0} bp, {(query.Reverse is null ? "forward strand only" : "both strands")}; " +
            (result.QueryWindowSize is { } w ? $"compared at {w:N0} bp" : "shorter than one window, comparable with nothing") +
            $" with records of {result.Lengths.Min:N0} to {result.Lengths.MaxOrUnbounded - 1:N0} bp (length ratio {ratio}). {elapsed.TotalMilliseconds:F0} ms.");
        foreach (var (id, strand) in result.Identical)
            Console.WriteLine($"  identical sequence ({strand} strand, by SHA-256): {Describe(records, id)}");

        PrintHits($"Ranked: top {result.Ranked.Count} of {result.RankedTotal:N0}", result.Ranked, records);
        if (result.NotComparable.Count > 0)
            Console.WriteLine($"Not comparable (no stored level at {result.QueryWindowSize:N0} bp): {result.NotComparable.Count:N0} records.");
    }
    return 0;
}

static void PrintHits(string title, IReadOnlyList<SimilarityHit> hits, IReadOnlyDictionary<long, CatalogueRecord> records)
{
    Console.WriteLine(title);
    if (hits.Count == 0) return;
    Console.WriteLine($"  {"#",3} {"score",7} {"strand",-7} {"id",7} {"kind",-10} {"length bp",12} {"genus",-20} name");
    for (int i = 0; i < hits.Count; i++)
    {
        var h = hits[i];
        var r = records.GetValueOrDefault(h.SequenceId);
        Console.WriteLine($"  {i + 1,3} {h.Score,7:F4} {h.Strand,-7} {h.SequenceId,7} {r?.Kind,-10} {r?.LengthBp,12:N0} {r?.Lineage.Genus,-20} {r?.Name}");
    }
}

static string Describe(IReadOnlyDictionary<long, CatalogueRecord> records, long id) =>
    records.TryGetValue(id, out var r) ? $"{id} {r.RecordId} {r.Name}" : id.ToString();

static async Task PrintListAsync(NpgsqlDataSource dataSource, CatalogueFilter filter, long after, int limit)
{
    var page = await Catalogue.ListAsync(dataSource, filter, after, limit);
    Console.WriteLine($"{"id",6} {"kind",-10} {"length bp",12} {"GC",6} {"genus",-20} {"record",-14} name");
    foreach (var r in page.Records)
        Console.WriteLine($"{r.Id,6} {r.Kind,-10} {r.LengthBp,12:N0} {r.GcFraction,6:F3} {r.Lineage.Genus,-20} {r.RecordId,-14} {r.Name}");
    Console.WriteLine(page.NextAfter is { } next ? $"{page.Records.Count} shown; next page: --after {next}" : $"{page.Records.Count} shown; last page");
}

static async Task PrintStatsAsync(NpgsqlDataSource dataSource)
{
    await using (var command = dataSource.CreateCommand(
        "select kind, count(*), sum(length_bp), count(distinct tax_family), count(distinct tax_genus) from sequence_record group by kind order by kind"))
    await using (var reader = await command.ExecuteReaderAsync())
    {
        Console.WriteLine($"{"kind",-10} {"records",8} {"bases",16} {"families",9} {"genera",7}");
        while (await reader.ReadAsync())
            Console.WriteLine($"{reader.GetString(0),-10} {reader.GetInt64(1),8:N0} {reader.GetDecimal(2),16:N0} {reader.GetInt64(3),9} {reader.GetInt64(4),7}");
    }

    await using (var command = dataSource.CreateCommand(
        "select window_size, count(*), sum(byte_size), sum(window_count) from profile_level group by window_size order by window_size"))
    await using (var reader = await command.ExecuteReaderAsync())
    {
        Console.WriteLine();
        Console.WriteLine($"{"level",-10} {"profiles",8} {"bytes",14} {"windows",12}");
        while (await reader.ReadAsync())
            Console.WriteLine($"{reader.GetInt32(0) + " bp",-10} {reader.GetInt64(1),8:N0} {reader.GetDecimal(2),14:N0} {reader.GetInt64(3),12:N0}");
    }

    await using (var command = dataSource.CreateCommand(
        "select provider, coalesce(item_group, '-'), status, count(*), max(attempts) from crawl_item group by 1, 2, 3 order by 1, 2, 3"))
    await using (var reader = await command.ExecuteReaderAsync())
    {
        Console.WriteLine();
        Console.WriteLine($"{"crawl",-10} {"group",-10} {"status",-10} {"items",8} {"max attempts",13}");
        while (await reader.ReadAsync())
            Console.WriteLine($"{reader.GetString(0),-10} {reader.GetString(1),-10} {reader.GetString(2),-10} {reader.GetInt64(3),8:N0} {reader.GetInt32(4),13}");
    }
}
