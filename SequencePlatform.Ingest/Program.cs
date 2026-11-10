using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;
using SequencePlatform.Ingest;

// Small test app: populate the platform DB from BV-BRC. The connection string comes from the
// SEQUENCE_PLATFORM_DB environment variable, or else the git-ignored appsettings.json (never committed).
//
//   migrate                                   create the database if missing, apply migrations; with SEQUENCE_APP_PASSWORD set, let sequence_app log in with it
//   ingest-bvbrc --group bacteria|archaea|viruses --limit N [--offset K] [--contact you@example.org] [--retry-failed] [--per-species]
//   ingest-ncbi --taxon fungi|plants|animals --limit N [--max-bp B] [--parallel P] [--retry-failed] [--contact you@example.org]   the eukaryote pilot: reference genomes, dinucleotide counts only
//   load-taxonomy [--refresh] [--skip-gtdb]   NCBI's current lineage per referenced taxid, GTDB's per assembly (S2.4)
//   stats                                     what's in the DB, including crawl state
//   list [--kind K] [--rank R --taxon T] [--min-length N] [--max-length N] [--provider P] [--name TEXT] [--after ID] [--limit N]
//   build-null [--per-cell N] [--genome-draws N] [--fragment-draws N] [--pool N] [--seed S] [--csv PATH]   score_null from stored profiles (P3 S3.6)
//   build-neighbours [--ratio R]             every record's own top 20 (other assemblies), for re-ranking by shared neighbours (P3 S3.4b)
//   build-origins                             circular bacterial chromosomes' levels laid out from their replication origin, for find-similar
//   convert-levels-half                       rewrite stored ZCWP version 1 levels as version 2 (float16), then rebuild the origin rows; no download
//   backfill-kmers [--limit N] [--sample SEED] [--contact you@example.org]   6-mer counts for BV-BRC records stored before migration 0013, by re-downloading
//   backfill-genes [--limit N] [--sample SEED] [--contact EMAIL] [--check-search N] [--parallel P]   the gene layer (S7.R, S7.S): gene_profile rows, amino-acid 4-mers, missing 6-mers
//   backfill-genes-ncbi [--limit N] [--accession A] [--dry-run] [--check-search N] [--contact EMAIL]   the gene layer for NCBI eukaryotes (S7.T), from each annotated assembly's CDS file; no genome download
//   build-protein-null [--per-bucket N] [--seed S] [--contact EMAIL]   p and E for protein hits: real BV-BRC proteins (not stored) scored against other families
//   verify-genes                              every gene row against its stored SHA-256 (protein search doesn't hash rows); exit 2 on a mismatch
//   cluster-genes                             gene_profile rewritten in bucket order, so protein search reads less scattered rows; rerun after a large ingest
//   profile --fasta X --out X.zcwp            the query-side profile of each record, both strands (X.zcwp and X.rc.zcwp); no database
//   search-similar --fasta X | --profile X.zcwp [list's filters] [--ratio R] [--top N]   records of comparable length (ratio 0.25-1, default 0.5) ranked by similarity (P3 S3.4)
//   localize --fasta X | --profile X.zcwp [--record ID] [list's filters] [--top N]   where each fragment sits: over records at least as long (>= 4 kbp), or inside one record (P3 S3.5)
//   add-user NAME                             an API user (D3); prints their key once, only its hash is stored
//   rotate-key NAME                           a new key for an enabled user, printed once; the old key stops working
//   disable-user NAME                         revokes that user's key; their uploads stay
const string usage =
    "usage: SequencePlatform.Ingest migrate | ingest-bvbrc --group bacteria|archaea|viruses --limit N [--offset K] [--contact EMAIL] [--retry-failed] [--per-species]\n" +
    "       | ingest-ncbi --taxon fungi|plants|animals --limit N [--max-bp B] [--parallel P] [--retry-failed] [--contact EMAIL] | load-taxonomy [--refresh] [--skip-gtdb] | stats\n" +
    "       | list [--kind K] [--rank R --taxon T] [--min-length N] [--max-length N] [--provider P] [--name TEXT] [--after ID] [--limit N]\n" +
    "       | build-null [--per-cell N] [--genome-draws N] [--fragment-draws N] [--pool N] [--seed S] [--csv PATH] | build-neighbours [--ratio R] | build-origins | convert-levels-half\n" +
    "       | backfill-kmers [--limit N] [--sample SEED] [--contact EMAIL] | backfill-genes [--limit N] [--sample SEED] [--contact EMAIL] [--check-search N] [--parallel P]\n" +
    "       | backfill-genes-ncbi [--limit N] [--accession A] [--dry-run] [--check-search N] [--contact EMAIL]\n" +
    "       | build-protein-null [--per-bucket N] [--seed S] [--contact EMAIL] | verify-genes | cluster-genes\n" +
    "       | profile --fasta X --out X.zcwp | search-similar --fasta X | --profile X.zcwp [list's filters] [--ratio R] [--top N]\n" +
    "       | localize --fasta X | --profile X.zcwp [--record ID] [list's filters] [--top N] | add-user NAME | rotate-key NAME | disable-user NAME";

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
        if (Environment.GetEnvironmentVariable(Database.AppPasswordVariable) is { Length: > 0 } appPassword)
        {
            await Database.SetAppPasswordAsync(dataSource, appPassword);
            Console.WriteLine($"Role {Database.AppRole} can log in with the password from {Database.AppPasswordVariable}.");
        }
        return 0;

    case "ingest-bvbrc":
        var retryFailed = args.Contains("--retry-failed");
        return await IngestBvBrcAsync(
            dataSource,
            Option("--group") ?? (retryFailed ? null : "bacteria"),
            int.Parse(Option("--limit") ?? "10"),
            int.Parse(Option("--offset") ?? "0"),
            Option("--contact"),
            retryFailed,
            args.Contains("--per-species"));

    case "ingest-ncbi":
        return await IngestNcbiAsync(dataSource, Option("--taxon") ?? "fungi", int.Parse(Option("--limit") ?? "100"),
            long.Parse(Option("--max-bp") ?? "2000000000"), int.Parse(Option("--parallel") ?? "1"), args.Contains("--retry-failed"), Option("--contact"), args.Contains("--dry-run"));

    case "load-taxonomy":
        return await LoadTaxonomyAsync(dataSource, args.Contains("--refresh"), args.Contains("--skip-gtdb"), Option("--contact"));

    case "stats":
        await PrintStatsAsync(dataSource);
        return 0;

    case "list":
        await PrintListAsync(dataSource, Filter(), long.Parse(Option("--after") ?? "0"), int.Parse(Option("--limit") ?? "50"));
        return 0;

    case "build-null":
        var defaults = new ScoreNullOptions();
        return await BuildNullAsync(dataSource, defaults with
        {
            SamplesPerCell = int.Parse(Option("--per-cell") ?? defaults.SamplesPerCell.ToString()),
            GenomeDraws = int.Parse(Option("--genome-draws") ?? defaults.GenomeDraws.ToString()),
            FragmentDraws = int.Parse(Option("--fragment-draws") ?? defaults.FragmentDraws.ToString()),
            FragmentPool = int.Parse(Option("--pool") ?? defaults.FragmentPool.ToString()),
            Seed = int.Parse(Option("--seed") ?? defaults.Seed.ToString()),
        }, Option("--csv"));

    case "build-neighbours":
        var neighbourRatio = double.Parse(Option("--ratio") ?? SimilaritySearch.DefaultLengthRatio.ToString(CultureInfo.InvariantCulture),
            CultureInfo.InvariantCulture);
        var neighbourWatch = System.Diagnostics.Stopwatch.StartNew();
        var listed = await NeighbourGraph.BuildAsync(dataSource, neighbourRatio, SimilaritySearch.NeighbourCount,
            message => Console.WriteLine($"  {neighbourWatch.Elapsed.TotalSeconds,6:F0} s  {message}"));
        Console.WriteLine($"{listed:N0} records listed at length ratio {neighbourRatio} in {neighbourWatch.Elapsed.TotalSeconds:F0} s.");
        return 0;

    case "build-origins":
        var originWatch = System.Diagnostics.Stopwatch.StartNew();
        var (rotatedCount, examined) = await OriginRotation.BuildAsync(dataSource,
            message => Console.WriteLine($"  {originWatch.Elapsed.TotalSeconds,6:F0} s  {message}"));
        Console.WriteLine($"{rotatedCount:N0} of {examined:N0} circular chromosomes stored rotated to their origin in {originWatch.Elapsed.TotalSeconds:F0} s. " +
            "Rebuild the null (build-null) and the neighbour graph (build-neighbours): both score the way searches do.");
        return 0;

    case "convert-levels-half":
        var convertWatch = System.Diagnostics.Stopwatch.StartNew();
        var convertedCount = await ProfileLevels.ConvertToHalfAsync(dataSource,
            message => Console.WriteLine($"  {convertWatch.Elapsed.TotalSeconds,6:F0} s  {message}"));
        Console.WriteLine($"{convertedCount:N0} records' levels rewritten as ZCWP version 2 in {convertWatch.Elapsed.TotalSeconds:F0} s; rebuilding the origin rows.");
        var (turnedCount, circular) = await OriginRotation.BuildAsync(dataSource,
            message => Console.WriteLine($"  {convertWatch.Elapsed.TotalSeconds,6:F0} s  {message}"));
        Console.WriteLine($"{turnedCount:N0} of {circular:N0} circular chromosomes rotated again in {convertWatch.Elapsed.TotalSeconds:F0} s. " +
            "Rebuild the null (build-null) and the neighbour graph (build-neighbours), then VACUUM FULL profile_level to give the space back.");
        return 0;

    case "backfill-kmers":
        return await BackfillKmersAsync(dataSource, int.Parse(Option("--limit") ?? "100"),
            Option("--sample") is { } sampleSeed ? int.Parse(sampleSeed) : null, Option("--contact"));

    case "backfill-genes":
        return await BackfillGenesAsync(dataSource, int.Parse(Option("--limit") ?? "100"),
            Option("--sample") is { } geneSeed ? int.Parse(geneSeed) : null, Option("--contact"), int.Parse(Option("--check-search") ?? "0"),
            int.Parse(Option("--parallel") ?? "1"));

    case "backfill-genes-ncbi":
        using (var ncbi = new NcbiGenomeClient(Option("--contact")))
        {
            var geneLayer = new NcbiGeneLayer(dataSource, ncbi);
            return await NcbiGenePassAsync(dataSource, geneLayer, await geneLayer.PendingAsync(int.Parse(Option("--limit") ?? "1000"), Option("--accession")),
                args.Contains("--dry-run"), int.Parse(Option("--check-search") ?? "0"));
        }

    case "verify-genes":
    {
        var watch = System.Diagnostics.Stopwatch.StartNew();
        var (checkedRows, mismatches) = await GeneProfiles.VerifyAsync(dataSource, rows => Console.WriteLine($"  {rows:N0} rows checked"));
        foreach (var (sequenceId, bucket) in mismatches) Console.WriteLine($"  MISMATCH: record {sequenceId}, bucket {bucket}");
        Console.WriteLine($"{checkedRows:N0} gene rows checked, {mismatches.Count} mismatches, {watch.Elapsed.TotalSeconds:F0} s.");
        return mismatches.Count == 0 ? 0 : 2;
    }
    case "cluster-genes":
    {
        var watch = System.Diagnostics.Stopwatch.StartNew();
        Console.WriteLine("Rewriting gene_profile in bucket order (exclusive lock: no protein search or gene write until it's done)...");
        await GeneProfiles.ClusterAsync(dataSource);
        Console.WriteLine($"gene_profile clustered by bucket, {watch.Elapsed.TotalSeconds:F0} s.");
        return 0;
    }
    case "build-protein-null":
        return await BuildProteinNullAsync(dataSource, int.Parse(Option("--per-bucket") ?? "30"), int.Parse(Option("--seed") ?? "20261106"), Option("--contact"));

    case "search-similar":
        return await SearchSimilarAsync(dataSource, Option("--fasta"), Option("--profile"), Filter(),
            double.Parse(Option("--ratio") ?? SimilaritySearch.DefaultLengthRatio.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture),
            int.Parse(Option("--top") ?? "20"));

    case "localize":
        return await LocalizeAsync(dataSource, Option("--fasta"), Option("--profile"), Option("--record") is { } record ? long.Parse(record) : null,
            Filter(), int.Parse(Option("--top") ?? "20"));

    case "add-user" when args.Length > 1:
        var (user, apiKey) = await Users.CreateAsync(dataSource, args[1]);
        Console.WriteLine($"User {user.Name} (id {user.Id}). API key, shown only now: {apiKey}");
        return 0;

    case "rotate-key" when args.Length > 1:
        if (await Users.RotateKeyAsync(dataSource, args[1]) is { } rotatedKey)
        {
            Console.WriteLine($"User {args[1]}: new API key, shown only now: {rotatedKey}. The old key no longer works.");
            return 0;
        }
        Console.WriteLine($"No enabled user named {args[1]}.");
        return 1;

    case "disable-user" when args.Length > 1:
        if (await Users.DisableAsync(dataSource, args[1])) return 0;
        Console.WriteLine($"No enabled user named {args[1]}.");
        return 1;

    default:
        Console.WriteLine(usage);
        return 1;
}

// Normal run: list a page of the group, register new items, process everything not in a final state
// except failed. --retry-failed: process failed items instead (of --group, or of every group if omitted).
// --per-species: one genome for each species of the group's complete genomes not yet stored, up to --limit; failed items of such a run
// are retried with --group <group>-species.
static async Task<int> IngestBvBrcAsync(NpgsqlDataSource dataSource, string? group, int limit, int offset, string? contact, bool retryFailed, bool perSpecies)
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

    var counts = new Dictionary<string, int>();
    int alreadyFinal = 0, failedSkipped = 0;

    async Task ProcessAsync(List<BvBrcGenome> genomes, string? crawlGroup)
    {
        var statuses = await CrawlState.RegisterAsync(dataSource, GenomeIngestor.Provider, crawlGroup, genomes.Select(g => g.GenomeId).ToList());

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
    }

    if (retryFailed)
    {
        var ids = await CrawlState.FailedAsync(dataSource, GenomeIngestor.Provider, group, limit);
        var genomes = ids.Count == 0 ? new List<BvBrcGenome>() : await client.GetGenomesAsync(ids);
        Console.WriteLine($"BV-BRC {group ?? "all groups"}: retrying {genomes.Count} of {ids.Count} failed items.");
        await ProcessAsync(genomes, null);
    }
    else if (perSpecies)
    {
        const int batchSize = 50;
        var stored = await Catalogue.StoredSpeciesAsync(dataSource);
        var crawlGroup = $"{group}-species";
        if (group == "viruses")
        {
            // Millions of genomes can't be listed to pick from, so each missing species is one request, made as the run reaches it: an
            // interrupted run loses one batch, and the next run sees those species as stored.
            var missing = SpeciesSelection.Missing(await client.ListSpeciesAsync(group), stored, limit);
            Console.WriteLine($"BV-BRC {group}: {missing.Count} new species to fetch (limit {limit}).");
            foreach (var batch in missing.Chunk(batchSize))
            {
                var genomes = new List<BvBrcGenome>();
                foreach (var species in batch)
                    try
                    {
                        if (await client.FirstGenomeOfSpeciesAsync(group, species) is { } genome) genomes.Add(genome);
                    }
                    catch (InvalidOperationException e)
                    {
                        // BV-BRC rejects some names outright, such as any with an apostrophe; one species shouldn't end the run.
                        Console.WriteLine($"  skipped species {species}: {e.Message}");
                    }
                await ProcessAsync(genomes, crawlGroup);
            }
        }
        else
        {
            const int pageSize = 25000;
            var listed = new List<BvBrcGenome>();
            for (int page = 0; ; page += pageSize)
            {
                var genomes = await client.ListCompleteGenomesAsync(group!, pageSize, page);
                listed.AddRange(genomes);
                if (genomes.Count < pageSize) break;
            }
            var picked = SpeciesSelection.PickNew(listed, stored, limit);
            Console.WriteLine($"BV-BRC {group}: {listed.Count} complete genomes listed, {picked.Count} new species to fetch (limit {limit}).");
            foreach (var batch in picked.Chunk(batchSize)) await ProcessAsync(batch.ToList(), crawlGroup);
        }
    }
    else
    {
        var genomes = await client.ListGenomesAsync(group!, limit, offset);
        Console.WriteLine($"BV-BRC {group}: {genomes.Count} genomes listed (limit {limit}, offset {offset}).");
        await ProcessAsync(genomes, group);
    }

    int processed = counts.Values.Sum();
    var elapsed = DateTime.UtcNow - runStarted;
    Console.WriteLine($"Done: {string.Join(", ", counts.Select(c => $"{c.Value} {c.Key}"))}{(processed == 0 ? "nothing to do" : "")}; " +
        $"{alreadyFinal} already handled, {failedSkipped} failed earlier (use --retry-failed). " +
        $"{elapsed.TotalSeconds:F0} s{(processed > 0 ? $", {elapsed.TotalSeconds / processed:F2} s per processed genome" : "")}.");
    return counts.ContainsKey(CrawlState.Failed) ? 2 : 0;
}

// The eukaryote pilot: one reference genome per family, counted and discarded; see NcbiIngestor.
static async Task<int> IngestNcbiAsync(NpgsqlDataSource dataSource, string taxon, int limit, long maxLengthBp, int parallel, bool retryFailed, string? contact, bool dryRun)
{
    int taxid = taxon switch
    {
        "fungi" => 4751,
        "plants" => 33090,
        "animals" => 33208,
        _ => throw new ArgumentException($"Unknown taxon '{taxon}'. Use fungi, plants or animals."),
    };
    using var client = new NcbiGenomeClient(contact);
    var ingestor = new NcbiIngestor(dataSource, client);

    var assemblies = await client.ListAssembliesAsync(taxid);
    var picked = NcbiIngestor.Pick(assemblies, await client.LineagesAsync(assemblies.Select(a => a.Taxid).ToList()), limit, maxLengthBp);
    Console.WriteLine($"NCBI {taxon}: {assemblies.Count} reference assemblies, {picked.Count} picked, {picked.Sum(p => p.Assembly.LengthBp) / 1e9:F1} Gbp in all.");
    if (dryRun)
    {
        foreach (var classes in picked.GroupBy(p => p.Lineage.Class ?? p.Lineage.Phylum ?? "unknown").OrderByDescending(c => c.Count()))
            Console.WriteLine($"  {classes.Key,-28} {classes.Count(),4} genomes {classes.Sum(p => p.Assembly.LengthBp) / 1e9,7:F1} Gbp");
        return 0;
    }
    var statuses = await CrawlState.RegisterAsync(dataSource, NcbiIngestor.Provider, taxon, picked.Select(p => p.Assembly.Accession).ToList());

    int failed = 0;
    var pending = picked.Where(p => statuses[p.Assembly.Accession] == "pending" || (retryFailed && statuses[p.Assembly.Accession] == CrawlState.Failed));
    var stored = new System.Collections.Concurrent.ConcurrentBag<string>();
    await Parallel.ForEachAsync(pending, new ParallelOptions { MaxDegreeOfParallelism = parallel }, async (item, _) =>
    {
        var (assembly, lineage) = item;
        var started = DateTime.UtcNow;
        var outcome = await ingestor.IngestAsync(assembly, lineage);
        if (outcome.StartsWith("FAILED", StringComparison.Ordinal)) Interlocked.Increment(ref failed);
        else stored.Add(assembly.Accession);
        Console.WriteLine($"  {assembly.Accession,-18} {assembly.Organism,-44} {outcome}  {(DateTime.UtcNow - started).TotalSeconds,6:F0} s");
    });
    Console.WriteLine($"Done: {failed} failed.");

    // The gene layer comes from the CDS file, not the genome just read, so it follows as its own pass; a failure there leaves the
    // record pending for backfill-genes-ncbi.
    var geneLayer = new NcbiGeneLayer(dataSource, client);
    var forGenes = new List<PendingEukaryote>();
    foreach (var accession in stored) forGenes.AddRange(await geneLayer.PendingAsync(1, accession));
    int genePass = forGenes.Count == 0 ? 0 : await NcbiGenePassAsync(dataSource, geneLayer, forGenes, dryRun: false, checkSearch: 0);
    return failed > 0 || genePass != 0 ? 2 : 0;
}

// The gene layer for NCBI eukaryotes (S7.T): annotated assemblies' coding sequences from their CDS files, never the genome; see
// NcbiGeneLayer. --check-search N: N proteins of each assembly searched through ProteinSearch over the whole catalogue afterwards.
static async Task<int> NcbiGenePassAsync(NpgsqlDataSource dataSource, NcbiGeneLayer layer, List<PendingEukaryote> pending, bool dryRun, int checkSearch)
{
    const int chunk = 50;
    string[] checkBands = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
    static int CheckBand(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;

    Console.WriteLine($"{pending.Count} NCBI eukaryote records for the gene pass{(dryRun ? " (dry run: nothing is written)" : "")}.");
    var statuses = new Dictionary<string, int>();
    var checks = new List<(PendingEukaryote Record, EukaryoteGenes Layer, int Gene)>();
    var picker = new Random(20261106);
    var watch = System.Diagnostics.Stopwatch.StartNew();
    long genes = 0;
    int largest = 0, done = 0;
    foreach (var batch in pending.Chunk(chunk))
    {
        var (annotations, codes) = await layer.LookUpAsync(batch);
        foreach (var record in batch)
        {
            var started = watch.Elapsed;
            int? code = record.Taxid is int taxid && codes.TryGetValue(taxid, out var known) ? known : null;
            var result = await layer.RunAsync(record, annotations.GetValueOrDefault(record.Accession), code, !dryRun);
            statuses[result.Status] = statuses.GetValueOrDefault(result.Status) + 1;
            done++;
            if (result.Layer is { } built)
            {
                var t = built.Tally;
                genes += built.Genes.Genes.Count;
                largest = Math.Max(largest, result.MaxBucketBytes);
                Console.WriteLine($"  {record.Accession,-18} {result.Note}; {t.Read:N0} CDS read, left out: isoforms {t.Isoforms:N0}, partial {t.Partial:N0}, " +
                    $"pseudo {t.Pseudo:N0}, exception {t.Exception:N0}, organelle {t.Organelle:N0}, location {t.Location:N0}, frame {t.Frame:N0}; " +
                    $"largest row {result.MaxBucketBytes:N0} bytes  {(watch.Elapsed - started).TotalSeconds:F0} s");
                if (checkSearch > 0)
                    checks.AddRange(Enumerable.Range(0, built.Genes.Genes.Count)
                        .Where(g => built.Proteins[g].IndexOfAny(['*', 'X']) < 0).OrderBy(_ => picker.Next()).Take(checkSearch).Select(g => (record, built, g)));
            }
            else if (result.Status != ProteinKmers.NoCds || pending.Count <= chunk)
                Console.WriteLine($"  {record.Accession,-18} {result.Status}: {result.Note}");
        }
        Console.WriteLine($"  {done}/{pending.Count} records, {watch.Elapsed.TotalMinutes:F1} min; " +
            string.Join(", ", statuses.OrderBy(s => s.Key).Select(s => $"{s.Key} {s.Value}")));
    }
    Console.WriteLine($"Gene pass: {genes:N0} genes profiled; largest row {largest:N0} bytes (limit {GeneProfiles.MaxRowBytes:N0}); " +
        $"{watch.Elapsed.TotalMinutes:F1} min.");

    if (checks.Count > 0 && !dryRun)
    {
        // Through the real search, against every stored gene: the right gene at rank 1, or tied with it (identical paralogs score alike).
        var outcomes = new List<(string Accession, int Band, bool Exact, bool Tied, double Seconds)>();
        foreach (var (record, built, g) in checks)
        {
            var gene = built.Genes.Genes[g];
            var query = ProteinQuery.FromProtein(built.Proteins[g]);
            var searchWatch = System.Diagnostics.Stopwatch.StartNew();
            var top = (await ProteinSearch.SearchAsync(dataSource, query, top: 1)).Hits[0];
            bool exact = top.SequenceId == record.SequenceId && top.GeneSequence == gene.Chromosome && top.GeneStart == gene.Start && top.GeneEnd == gene.End;
            bool tied = exact || built.Genes.Score(query.Axes, g) >= top.Score;
            outcomes.Add((record.Accession, CheckBand(gene.Residues), exact, tied, searchWatch.Elapsed.TotalSeconds));
        }
        Console.WriteLine($"Self-search through ProteinSearch over the catalogue ({outcomes.Count} proteins, median {outcomes.Select(o => o.Seconds).Order().ElementAt(outcomes.Count / 2):F1} s each):");
        foreach (var group in outcomes.GroupBy(o => o.Accession).Concat(outcomes.GroupBy(o => checkBands[o.Band]).OrderBy(b => Array.IndexOf(checkBands, b.Key))))
            Console.WriteLine($"  {group.Key,-18} {group.Count(),4} proteins: right gene first {100.0 * group.Count(o => o.Exact) / group.Count():F1}%, " +
                $"first or tied {100.0 * group.Count(o => o.Tied) / group.Count():F1}%");
    }
    return statuses.ContainsKey(CrawlState.Failed) ? 2 : 0;
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

static async Task<int> BuildNullAsync(NpgsqlDataSource dataSource, ScoreNullOptions options, string? csvPath)
{
    var watch = System.Diagnostics.Stopwatch.StartNew();
    var built = await ScoreNullBuilder.BuildAsync(dataSource, options, message => Console.WriteLine($"  {watch.Elapsed.TotalSeconds,6:F1} s  {message}"));
    foreach (var (source, cells) in built)
    {
        await ScoreNullStore.SaveAsync(dataSource, RecordStore.Method, RecordStore.MethodVersion, source, cells.ToList());
        Console.WriteLine($"{source,-9} {cells.Count,4} cells, {cells.Sum(c => c.SortedScores.Length),8:N0} scores, " +
            $"{cells.Count(c => c.SortedScores.Length >= options.SamplesPerCell)} full");
    }
    if (csvPath is not null)
    {
        using var writer = File.CreateText(csvPath);
        ScoreNullStore.WriteCsv(writer, built.SelectMany(b => b.Value.Select(c => (b.Key, c))));
        Console.WriteLine($"Written to {csvPath}.");
    }
    Console.WriteLine($"Done in {watch.Elapsed.TotalSeconds:F0} s.");
    return 0;
}

static string ReverseStrandPath(string path) => Path.ChangeExtension(path, null) + ".rc.zcwp";

static async Task<int> SearchSimilarAsync(NpgsqlDataSource dataSource, string? fasta, string? profile, CatalogueFilter filter, double ratio, int top)
{
    if (ReadQueries(fasta, profile) is not { } queries)
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
        if (result.FewWindows)
            Console.WriteLine($"  WARNING: under {ZCurveProfileLadder.DefaultMinWindows} windows for the shortest candidates; " +
                "this ranking mostly reflects base composition (GC), not relatedness.");
        foreach (var (id, strand) in result.Identical)
            Console.WriteLine($"  identical sequence ({strand} strand, by SHA-256): {Describe(records, id)}");
        if (result.Ranked.Count > 0 && result.Ranked[0].P is null)
            Console.WriteLine("  no score_null yet (build-null): p and E aren't shown.");

        PrintHits($"Ranked: top {result.Ranked.Count} of {result.RankedTotal:N0}" + (result.NeighbourRanked
            ? $", the first {SimilaritySearch.NeighbourRankDepth} by shared neighbours, then score"
            : " by score (no neighbour graph for this ratio: build-neighbours)"), result.Ranked, records);
        if (result.NotComparable.Count > 0)
            Console.WriteLine($"Not comparable (no stored level at {result.QueryWindowSize:N0} bp): {result.NotComparable.Count:N0} records.");
    }
    return 0;
}

// Each FASTA record, or one profile (with its X.rc.zcwp reverse strand if present); null when neither was given.
static List<(string Name, QueryProfile Query)>? ReadQueries(string? fasta, string? profile)
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
    else return null;
    return queries;
}

static async Task<int> LocalizeAsync(NpgsqlDataSource dataSource, string? fasta, string? profile, long? record, CatalogueFilter filter, int top)
{
    if (ReadQueries(fasta, profile) is not { } queries)
    {
        Console.WriteLine("usage: localize --fasta X | --profile X.zcwp [--record ID] [list's filters] [--top N]");
        return 1;
    }

    var localization = new Localization(dataSource, useCache: false);
    foreach (var (name, query) in queries)
    {
        Console.WriteLine();
        if (record is null && query.LengthBp < Localization.MinCatalogueLengthBp)
        {
            Console.WriteLine($"Fragment {name}: {query.LengthBp:N0} bp is under {Localization.MinCatalogueLengthBp:N0} bp, too short to localize " +
                "over the catalogue; name a record with --record ID.");
            continue;
        }

        var watch = System.Diagnostics.Stopwatch.StartNew();
        var result = record is long id
            ? await localization.LocalizeInRecordAsync(query, id)
            : await localization.LocalizeAsync(query, filter, top);
        var elapsed = watch.Elapsed;
        var records = await Catalogue.ByIdsAsync(dataSource, result.Hits.Select(h => h.SequenceId).ToHashSet());

        Console.WriteLine($"Fragment {name}: {query.LengthBp:N0} bp, {(query.Reverse is null ? "forward strand only" : "both strands")}; " +
            (record is long r
                ? $"placed inside record {r} at 500 bp"
                : $"scanned at {result.CoarseWindowSize:N0} bp over {result.Candidates:N0} records at least as long, the best {Localization.RefineDepth} refined at 500 bp") +
            $". {elapsed.TotalMilliseconds:F0} ms.");
        if (result.WeakSignal)
            Console.WriteLine($"  WARNING: under {Localization.WeakSignalLengthBp:N0} bp the placement is weak; check it against an aligner.");
        if (result.Hits.Count > 0 && result.Hits[0].P is null)
            Console.WriteLine("  no score_null yet (build-null): p and E aren't shown.");
        if (result.Hits.Count == 0)
            Console.WriteLine("  no placement: the fragment has no full 500 bp window or is longer than the record.");
        else
        {
            Console.WriteLine($"  {"#",3} {"score",7} {"log10 E",9} {"strand",-7} {"offset bp",12} {"id",7} {"kind",-10} {"length bp",12} {"genus",-20} name");
            Console.WriteLine("      (offset: where the fragment starts on the record's stored strand, to about 500 bp; Reverse: its reverse complement matches)");
            for (int i = 0; i < result.Hits.Count; i++)
            {
                var h = result.Hits[i];
                var c = records.GetValueOrDefault(h.SequenceId);
                var e = h.P is not { } p ? "" : p.Kind == PValueKind.Bound ? $"<={h.Log10E:F1}" : $"{h.Log10E:F1}";
                Console.WriteLine($"  {i + 1,3} {h.Score,7:F4} {e,9} {h.Strand,-7} {h.OffsetBp,12:N0} {h.SequenceId,7} {c?.Kind,-10} {c?.LengthBp,12:N0} {c?.Lineage.Genus,-20} {c?.Name}");
            }
        }
        if (result.NotComparable.Count > 0)
            Console.WriteLine($"Not comparable (no stored level at {result.CoarseWindowSize:N0} bp): {result.NotComparable.Count:N0} records.");
    }
    return 0;
}

static void PrintHits(string title, IReadOnlyList<SimilarityHit> hits, IReadOnlyDictionary<long, CatalogueRecord> records)
{
    Console.WriteLine(title);
    if (hits.Count == 0) return;
    Console.WriteLine($"  {"#",3} {"score",7} {"log10 E",9} {"shared",6} {"strand",-8} {"id",7} {"kind",-10} {"length bp",12} {"genus",-20} name");
    Console.WriteLine("      (log10 E: expected unrelated records scoring as high; '<=' only an upper bound: the null sample has too few scores that high." +
        " strand*: best laid out from the replication origin)");
    for (int i = 0; i < hits.Count; i++)
    {
        var h = hits[i];
        var r = records.GetValueOrDefault(h.SequenceId);
        Console.WriteLine($"  {i + 1,3} {h.Score,7:F4} {Evalue(h),9} {h.SharedNeighbours,6} {h.Strand + (h.Rotated ? "*" : ""),-8} {h.SequenceId,7} {r?.Kind,-10} {r?.LengthBp,12:N0} {r?.Lineage.Genus,-20} {r?.Name}");
    }
}

static string Evalue(SimilarityHit hit) => hit.P is not { } p ? "" : p.Kind == PValueKind.Bound ? $"<={hit.Log10E:F1}" : $"{hit.Log10E:F1}";

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

static async Task<int> BackfillKmersAsync(NpgsqlDataSource dataSource, int limit, int? seed, string? contact)
{
    using var client = new BvBrcClient(contact);
    var backfill = new KmerBackfill(dataSource, client);
    var genomes = await backfill.PendingGenomesAsync(limit, seed);
    Console.WriteLine($"{genomes.Count} genomes to backfill{(seed is null ? "" : $" (random sample, seed {seed})")}.");
    var statuses = new Dictionary<string, int>();
    int failed = 0;
    var watch = System.Diagnostics.Stopwatch.StartNew();
    for (int i = 0; i < genomes.Count; i++)
    {
        try
        {
            foreach (var status in await backfill.BackfillAsync(genomes[i])) statuses[status] = statuses.GetValueOrDefault(status) + 1;
        }
        catch (Exception e)
        {
            // Nothing was written for this genome, so the next run tries it again.
            failed++;
            Console.WriteLine($"  {genomes[i]}: {e.GetType().Name}: {e.Message}");
        }
        if ((i + 1) % 25 == 0 || i + 1 == genomes.Count)
            Console.WriteLine($"  {i + 1}/{genomes.Count} genomes, {watch.Elapsed.TotalMinutes:F1} min, {watch.Elapsed.TotalSeconds / (i + 1):F2} s each; " +
                string.Join(", ", statuses.OrderBy(s => s.Key).Select(s => $"{s.Key} {s.Value}")) + (failed > 0 ? $", failed {failed}" : ""));
    }
    return failed == 0 ? 0 : 1;
}

// The S7.R pass. Genomes go in chunks of 10 so one BV-BRC request lists their coding sequences. With --check-search N, N genes of
// at least 1.5 kbp per genome keep their protein in memory, and at the end each is searched from the protein alone (uniform codon
// weights, first 1,500 bp) among every gene of the pass with three windows: the trial's self-search check. Nothing is written for it.
// p and E for protein hits (migration 0019): a seeded sample of real BV-BRC proteins with a family, per length bucket, scored as protein
// search scores them against the genes of other families. The proteins aren't stored; only the score counts are.
static async Task<int> BuildProteinNullAsync(NpgsqlDataSource dataSource, int perBucket, int seed, string? contact)
{
    var watch = System.Diagnostics.Stopwatch.StartNew();
    using var client = new BvBrcClient(contact);
    var queries = await ProteinNullSample.CollectAsync(dataSource, client, perBucket, seed);
    Console.WriteLine($"{queries.Count:N0} proteins in {watch.Elapsed.TotalMinutes:F1} min; scoring them against the gene layer.");
    var cells = await ProteinScoreNull.BuildAsync(dataSource, queries, Console.WriteLine);
    Console.WriteLine($"{"bucket",6} {"proteins",8} {"pairs",14} {"score with 1 pair in 1,000",26}");
    foreach (var cell in cells)
    {
        int bin = Enumerable.Range(0, ProteinScoreNull.Bins + 1).First(b => cell.Exceed[b] * 1000 < cell.Pairs);
        Console.WriteLine($"{cell.Bucket,6} {cell.Queries,8} {cell.Pairs,14:N0} {bin / (double)ProteinScoreNull.Bins,26:F3}");
    }
    Console.WriteLine($"Protein score null: {cells.Count} buckets in {watch.Elapsed.TotalMinutes:F1} min.");
    return 0;
}

static async Task<int> BackfillGenesAsync(NpgsqlDataSource dataSource, int limit, int? seed, string? contact, int checkSearch, int parallel)
{
    string[] checkBands = ["50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500+ aa"];
    static int CheckBand(int residues) => residues < 100 ? 0 : residues < 167 ? 1 : residues < 333 ? 2 : residues < 500 ? 3 : 4;
    // Levina-Bickel tail below the 100 nearest competitors: the share of a pool this size expected closer than `distance`.
    static double TailShare(double[] sorted, double distance)
    {
        int n = sorted.Length, k = Math.Min(100, n);
        double kth = Math.Max(sorted[k - 1], 1e-12);
        if (distance >= kth) return sorted.Count(d => d < distance) / (double)n;
        double logSum = 0;
        for (int i = 0; i < k - 1; i++) logSum += Math.Log(kth / Math.Max(sorted[i], 1e-12));
        return k / (double)n * Math.Pow(Math.Max(distance, 0) / kth, (k - 1) / logSum);
    }

    const int chunk = 10;
    using var client = new BvBrcClient(contact);
    var backfill = new GeneBackfill(dataSource, client);
    var genomes = await backfill.PendingGenomesAsync(limit, seed);
    Console.WriteLine($"{genomes.Count} genomes for the gene pass{(seed is null ? "" : $" (random sample, seed {seed})")}, {parallel} at a time.");
    var statuses = new Dictionary<string, int>();
    long cds = 0, clean = 0, skipped = 0, records = 0;
    int failed = 0, done = 0;
    var checks = new List<(long SequenceId, int Gene, string Protein)>();
    var picker = new Random(seed ?? 0);
    var progress = new object();
    var watch = System.Diagnostics.Stopwatch.StartNew();
    // The BV-BRC client spaces and retries its requests across callers, so several genomes can download at once; each one still
    // commits in its own transaction, which keeps a rerun safe.
    for (int i = 0; i < genomes.Count; i += chunk)
    {
        var batch = genomes.Skip(i).Take(chunk).ToList();
        Dictionary<string, List<BvBrcCds>> codingSequences;
        try { codingSequences = await backfill.CodingSequencesAsync(batch); }
        catch (Exception e)
        {
            failed += batch.Count;
            Console.WriteLine($"  CDS listing for {batch[0]}..: {e.GetType().Name}: {e.Message}");
            continue;
        }
        await Parallel.ForEachAsync(batch, new ParallelOptions { MaxDegreeOfParallelism = parallel }, async (genome, _) =>
        {
            List<RecordGeneResult>? results = null;
            string? error = null;
            try { results = await backfill.BackfillAsync(genome, codingSequences[genome]); }
            catch (Exception e) { error = $"  {genome}: {e.GetType().Name}: {e.Message}"; }
            lock (progress)
            {
                if (error is not null)
                {
                    // Nothing was written for this genome, so the next run tries it again.
                    failed++;
                    Console.WriteLine(error);
                }
                foreach (var result in results ?? new List<RecordGeneResult>())
                {
                    statuses[result.Status] = statuses.GetValueOrDefault(result.Status) + 1;
                    (cds, clean, skipped, records) = (cds + result.CodingSequences, clean + result.CleanTranslations, skipped + result.Skipped, records + 1);
                    // checkSearch whole proteins per record and length band, as a user would send them.
                    if (checkSearch > 0 && result.Genes is not null)
                        foreach (var band in Enumerable.Range(0, checkBands.Length))
                            checks.AddRange(Enumerable.Range(0, result.Genes.Genes.Count)
                                .Where(g => CheckBand(result.Genes.Genes[g].Residues) == band && !result.Proteins[g].Contains('*'))
                                .OrderBy(_ => picker.Next()).Take(checkSearch).Select(g => (result.SequenceId, g, result.Proteins[g])));
                }
                done++;
                if (done % 25 == 0 || done == genomes.Count)
                    Console.WriteLine($"  {done}/{genomes.Count} genomes, {watch.Elapsed.TotalMinutes:F1} min, {watch.Elapsed.TotalSeconds / done:F2} s each; " +
                        string.Join(", ", statuses.OrderBy(s => s.Key).Select(s => $"{s.Key} {s.Value}")) + (failed > 0 ? $", failed {failed}" : "") +
                        $"; CDS {cds:N0}, clean translation {(cds == 0 ? 0 : 100.0 * clean / (cds - skipped)):F2}%, skipped {skipped}");
            }
        });
    }

    if (checks.Count > 0)
    {
        // The pool: every record the pass wrote, not only those a check came from.
        var profiles = new Dictionary<long, ZCurveGeneProfile>();
        var ids = new List<long>();
        long catalogueRecords;
        await using (var command = dataSource.CreateCommand(
            "select k.sequence_id from protein_kmer k where k.k = $1 and k.status = 'done' and k.checked_at >= now() - make_interval(secs => $2)"))
        {
            command.Parameters.AddWithValue((short)ProteinKmers.Ks.Max());
            command.Parameters.AddWithValue(watch.Elapsed.TotalSeconds + 60);
            await using var reader = await command.ExecuteReaderAsync();
            while (await reader.ReadAsync()) ids.Add(reader.GetInt64(0));
        }
        await using (var command = dataSource.CreateCommand("select count(*) from protein_kmer where k = $1 and status = 'done'"))
        {
            command.Parameters.AddWithValue((short)ProteinKmers.Ks.Max());
            catalogueRecords = (long)(await command.ExecuteScalarAsync())!;
        }
        foreach (var id in ids)
            if (await GeneProfiles.ReadAsync(dataSource, id) is { } profile) profiles[id] = profile;
        var pool = profiles.SelectMany(p => Enumerable.Range(0, p.Value.Genes.Count).Select(g => (Record: p.Key, Gene: g))).ToList();
        double scale = catalogueRecords / (double)Math.Max(1, profiles.Count);

        // As ProteinSearch compares: from 500 aa against genes that long, below it against genes within 10% of the length. Measured
        // top-1 within the pool, and extrapolated to the catalogue as T3 was (results/T3_Pool300.csv): a query some other gene already
        // beats keeps that share, else the share comes from the tail below its nearest competitors.
        var outcomes = new (int Band, bool Top, double Catalogue)[checks.Count];
        Parallel.For(0, checks.Count, j =>
        {
            var check = checks[j];
            var own = profiles[check.SequenceId];
            int residues = own.Genes[check.Gene].Residues;
            var query = ZCurveGeneProfile.QueryFromProtein(check.Protein);
            double trueDistance = 1 - own.Score(query, check.Gene);
            int closerOwn = 0;
            var others = new List<double>();
            foreach (var (record, gene) in pool)
            {
                if (record == check.SequenceId && gene == check.Gene) continue;
                var candidate = profiles[record].Genes[gene];
                bool comparable = residues >= ZCurveGeneProfile.LongResidues ? candidate.Residues >= ZCurveGeneProfile.LongResidues
                    : candidate.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(candidate.Residues - residues) <= ProteinSearch.LengthBand * residues;
                if (!comparable) continue;
                double distance = 1 - profiles[record].Score(query, gene);
                if (record == check.SequenceId) { if (distance < trueDistance) closerOwn++; }
                else others.Add(distance);
            }
            var sorted = others.Order().ToArray();
            int closer = sorted.Count(d => d < trueDistance);
            double share = closerOwn > 0 ? 1 : closer > 0 || sorted.Length < 2 ? closer / (double)Math.Max(1, sorted.Length) : TailShare(sorted, trueDistance);
            outcomes[j] = (CheckBand(residues), closerOwn == 0 && closer == 0, closerOwn > 0 ? 0 : Math.Exp(-scale * sorted.Length * share));
        });
        Console.WriteLine($"Self-search among {pool.Count:N0} genes of {profiles.Count:N0} records, from the protein alone (catalogue: {catalogueRecords:N0} records):");
        for (int band = 0; band < checkBands.Length; band++)
        {
            var inBand = outcomes.Where(o => o.Band == band).ToList();
            if (inBand.Count > 0)
                Console.WriteLine($"  {checkBands[band],-12} {inBand.Count,5} proteins: top-1 {100.0 * inBand.Count(o => o.Top) / inBand.Count:F1}% here, " +
                    $"{100 * inBand.Average(o => o.Catalogue):F1}% extrapolated to the catalogue");
        }
    }

    await using (var command = dataSource.CreateCommand(
        "select k, count(*), avg(octet_length(payload)) from protein_kmer where status = 'done' group by k order by k"))
    await using (var reader = await command.ExecuteReaderAsync())
        while (await reader.ReadAsync())
            Console.WriteLine($"Amino-acid {reader.GetInt16(0)}-mers: {reader.GetInt64(1):N0} records, mean {reader.GetDecimal(2):N0} bytes.");
    await using (var command = dataSource.CreateCommand(
        "select count(distinct sequence_id), count(*), sum(byte_size), sum(gene_count), sum(window_count) from gene_profile"))
    {
        await using var reader = await command.ExecuteReaderAsync();
        if (await reader.ReadAsync() && reader.GetInt64(1) > 0)
            Console.WriteLine($"Gene profiles: {reader.GetInt64(0):N0} records in {reader.GetInt64(1):N0} bucket rows, {reader.GetInt64(2):N0} bytes, " +
                $"{reader.GetInt64(3):N0} genes, {reader.GetInt64(4):N0} windows.");
    }
    return failed == 0 ? 0 : 1;
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
        "select method, window_size, count(*), sum(byte_size), sum(window_count) from profile_level group by 1, 2 order by 1 desc, 2"))
    await using (var reader = await command.ExecuteReaderAsync())
    {
        Console.WriteLine();
        Console.WriteLine($"{"level",-10} {"profiles",8} {"bytes",14} {"windows",12}");
        while (await reader.ReadAsync())
            Console.WriteLine($"{reader.GetInt32(1) + " bp",-10} {reader.GetInt64(2),8:N0} " +
                $"{reader.GetDecimal(3),14:N0} {reader.GetInt64(4),12:N0}");
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
