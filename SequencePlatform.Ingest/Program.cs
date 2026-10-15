using Npgsql;
using SequencePlatform.Ingest;

// Small test app: populate the platform DB from BV-BRC. The connection string comes from the
// SEQUENCE_PLATFORM_DB environment variable, or else the git-ignored appsettings.json (never committed).
//
//   migrate                                   create the database if missing, apply migrations
//   ingest-bvbrc --group bacteria|viruses --limit N [--offset K] [--contact you@example.org]
//   stats                                     what's in the DB
const string usage =
    "usage: SequencePlatform.Ingest migrate | ingest-bvbrc --group bacteria|viruses --limit N [--offset K] [--contact EMAIL] | stats";

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

var settings = Database.ConnectionSettings();
await Database.EnsureDatabaseAsync(settings);
await using var dataSource = NpgsqlDataSource.Create(settings.ConnectionString);
await Database.MigrateAsync(dataSource);

switch (args[0])
{
    case "migrate":
        Console.WriteLine($"Database {settings.Database} is up to date.");
        return 0;

    case "ingest-bvbrc":
        return await IngestBvBrcAsync(
            dataSource,
            Option("--group") ?? "bacteria",
            int.Parse(Option("--limit") ?? "10"),
            int.Parse(Option("--offset") ?? "0"),
            Option("--contact"));

    case "stats":
        await PrintStatsAsync(dataSource);
        return 0;

    default:
        Console.WriteLine(usage);
        return 1;
}

static async Task<int> IngestBvBrcAsync(NpgsqlDataSource dataSource, string group, int limit, int offset, string? contact)
{
    using var client = new BvBrcClient(contact);
    var ingestor = new GenomeIngestor(dataSource, client);

    var genomes = await client.ListGenomesAsync(group, limit, offset);
    Console.WriteLine($"BV-BRC {group}: {genomes.Count} genomes listed (limit {limit}, offset {offset}).");

    int ingested = 0, skipped = 0, duplicates = 0, failed = 0;
    foreach (var genome in genomes)
    {
        if (await ingestor.AlreadyIngestedAsync(genome.GenomeId))
        {
            skipped++;
            continue;
        }

        try
        {
            var started = DateTime.UtcNow;
            var result = await ingestor.IngestAsync(genome);
            if (result is null)
            {
                duplicates++;
                Console.WriteLine($"  {genome.GenomeId,-14} {genome.Name,-60} duplicate of an assembly already stored, skipped");
                continue;
            }

            var (records, levels, bases) = result.Value;
            ingested++;
            Console.WriteLine($"  {genome.GenomeId,-14} {genome.Name,-60} {records,3} records {levels,3} levels {bases,12:N0} bp  {(DateTime.UtcNow - started).TotalSeconds,5:F1} s");
        }
        catch (Exception e)
        {
            // One bad genome shouldn't stop the batch; it stays un-ingested, so a re-run retries it.
            failed++;
            Console.WriteLine($"  {genome.GenomeId,-14} FAILED: {e.Message}");
        }
    }

    Console.WriteLine($"Done: {ingested} ingested, {skipped} already present, {duplicates} duplicates skipped, {failed} failed.");
    return failed == 0 ? 0 : 2;
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
}
