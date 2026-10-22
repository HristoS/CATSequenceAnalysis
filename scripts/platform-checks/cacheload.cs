#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using Npgsql;
using SequencePlatform.Core;

// One level per process, so the heap delta is that level's alone.
int windowSize = int.Parse(args[0]);
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var cache = new ProfileCache(dataSource, windowSize);

GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();
long heapBefore = GC.GetTotalMemory(true);
var watch = Stopwatch.StartNew();
var (loaded, _) = await cache.RefreshAsync();
double loadSeconds = watch.Elapsed.TotalSeconds;

watch.Restart();
var (again, full) = await cache.RefreshAsync();
double refreshMs = watch.Elapsed.TotalMilliseconds;

long heapAfter = GC.GetTotalMemory(true);
var process = Process.GetCurrentProcess();
Console.WriteLine($"{windowSize}, {loaded}, {cache.WindowCount}, {cache.WindowBytes / 1e6:F1}, {(heapAfter - heapBefore) / 1e6:F1}, " +
    $"{process.WorkingSet64 / 1e6:F0}, {process.PeakWorkingSet64 / 1e6:F0}, {loadSeconds:F2}, {refreshMs:F1}, {again}, {full}");
