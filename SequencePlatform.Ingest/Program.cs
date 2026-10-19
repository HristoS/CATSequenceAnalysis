using Npgsql;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

// Small test app: populate the platform DB from BV-BRC. The connection string comes from the
// SEQUENCE_PLATFORM_DB environment variable, or else the git-ignored appsettings.json (never committed).
//
//   migrate                                   create the database if missing, apply migrations
//   ingest-bvbrc --group bacteria|viruses --limit N [--offset K] [--contact you@example.org] [--retry-failed]
//   stats                                     what's in the DB, including crawl state
//   list [--kind K] [--rank R --taxon T] [--min-length N] [--max-length N] [--provider P] [--name TEXT] [--after ID] [--limit N]
const string usage =
    "usage: SequencePlatform.Ingest migrate | ingest-bvbrc --group bacteria|viruses --limit N [--offset K] [--contact EMAIL] [--retry-failed] | stats\n" +
    "       | list [--kind K] [--rank R --taxon T] [--min-length N] [--max-length N] [--provider P] [--name TEXT] [--after ID] [--limit N]";

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

    case "stats":
        await PrintStatsAsync(dataSource);
        return 0;

    case "list":
        var filter = new CatalogueFilter(
            Option("--kind"), Option("--rank"), Option("--taxon"),
            Option("--min-length") is { } min ? long.Parse(min) : null,
            Option("--max-length") is { } max ? long.Parse(max) : null,
            Option("--provider"), Option("--name"));
        await PrintListAsync(dataSource, filter, long.Parse(Option("--after") ?? "0"), int.Parse(Option("--limit") ?? "50"));
        return 0;

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

    var counts = new Dictionary<string, int>();
    int alreadyFinal = 0, failedSkipped = 0;
    foreach (var genome in genomes)
    {
        var status = statuses[genome.GenomeId];
        if (CrawlState.IsFinal(status)) { alreadyFinal++; continue; }
        if (status == CrawlState.Failed && !retryFailed) { failedSkipped++; continue; }

        var started = DateTime.UtcNow;
        var outcome = await ingestor.IngestAsync(genome);
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
