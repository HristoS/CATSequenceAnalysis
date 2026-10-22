#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

var root = "data/genomes/";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var cached = new SimilaritySearch(dataSource);
var streamed = new SimilaritySearch(dataSource, useCache: false);

GC.Collect(); long heapBefore = GC.GetTotalMemory(true);
Console.WriteLine("Query, LengthBp, WindowSize, Ranked, CachedFirstMs, CachedWarmMs, StreamedMs");
foreach (var (name, file) in new[] { ("E. coli K-12", "ecoli_k12_mg1655_GCF_000005845.2.fna.gz"), ("SARS-CoV-2", "viral/sars_cov_2_NC_045512.2.fna.gz"), ("HPV16", "viral/hpv16_NC_001526.4.fna.gz") })
{
    var query = QueryProfile.FromSequence(FastaReader.ReadPrimaryRecord(root + file).Sequence);
    var watch = Stopwatch.StartNew();
    var first = await cached.SearchAsync(query);
    var firstMs = watch.Elapsed.TotalMilliseconds;
    var warm = new List<double>();
    for (int i = 0; i < 5; i++) { watch.Restart(); await cached.SearchAsync(query); warm.Add(watch.Elapsed.TotalMilliseconds); }
    await streamed.SearchAsync(query);
    var stream = new List<double>();
    for (int i = 0; i < 3; i++) { watch.Restart(); await streamed.SearchAsync(query); stream.Add(watch.Elapsed.TotalMilliseconds); }
    Console.WriteLine($"{name}, {query.LengthBp}, {first.QueryWindowSize}, {first.RankedTotal}, {firstMs:F0}, {warm.Order().ElementAt(2):F0}, {stream.Order().ElementAt(1):F0}");
}
long heapAfter = GC.GetTotalMemory(true);
Console.WriteLine($"# cached window data {cached.WindowBytes / 1e6:F1} MB; heap delta {(heapAfter - heapBefore) / 1e6:F1} MB; " +
    $"levels held for lengths {string.Join("; ", Enumerable.Range(0, 4).Select(i => $"{500 << (2 * i)} bp: {cached.CachedLengths(i).Min}-{cached.CachedLengths(i).Max}"))}");
