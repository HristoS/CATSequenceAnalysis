#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Diagnostics;
using SequenceAnalyses;
using SequencePlatform.Core;
using Npgsql;

// Storage check (2026-11-04): can the 2 / 8 / 32 kbp rows and their origin-rotated versions be dropped and derived from the 500 bp row
// when a level is loaded? Times the load the worker does today (stream the stored 2 kbp level and its "-origin" version) against the
// derived path (stream the 500 bp level, Downsample(4, 0), and for records that have an origin row, rotate to the origin first), and
// checks that the derived levels equal the stored ones (largest absolute difference over all windows). Each path runs twice, alternating,
// so the database cache is warm for both. Rule fixed before the run: keep the stored coarse rows if the derived load takes more than
// twice the stored one (warm). Prints one line per run and the verdict.
const string method = "zcurve-windowed", version = "1";
const int fine = 500, factor = 4;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var stored = new Dictionary<long, ZCurveWindowedProfile>();
var storedOrigin = new Dictionary<long, ZCurveWindowedProfile>();
var timesStored = new List<double>();
var timesDerived = new List<double>();
for (int run = 0; run < 2; run++)
{
    var watch = Stopwatch.StartNew();
    stored.Clear();
    storedOrigin.Clear();
    await foreach (var (id, level) in ProfileLevels.StreamAsync(dataSource, method, version, fine * factor)) stored[id] = level;
    await foreach (var (id, level) in ProfileLevels.StreamAsync(dataSource, method, OriginRotation.RotatedVersion(version), fine * factor)) storedOrigin[id] = level;
    timesStored.Add(watch.Elapsed.TotalSeconds);
    Console.WriteLine($"stored, run {run + 1}: {stored.Count} + {storedOrigin.Count} levels in {timesStored[^1]:F1} s");

    watch.Restart();
    var derived = new Dictionary<long, ZCurveWindowedProfile>();
    var derivedOrigin = new Dictionary<long, ZCurveWindowedProfile>();
    await foreach (var (id, level) in ProfileLevels.StreamAsync(dataSource, method, version, fine))
    {
        if (level.WindowCount >= factor) derived[id] = level.Downsample(factor, 0);
        if (storedOrigin.ContainsKey(id) && OriginRotation.ToOrigin(level, (long)level.WindowCount * fine) is { } rotated)
            derivedOrigin[id] = rotated.Downsample(factor, 0);
    }
    timesDerived.Add(watch.Elapsed.TotalSeconds);
    Console.WriteLine($"derived, run {run + 1}: {derived.Count} + {derivedOrigin.Count} levels in {timesDerived[^1]:F1} s");

    if (run == 1)
    {
        double worst = 0;
        int compared = 0, missing = 0;
        foreach (var (pairs, label) in new[] { (stored, "plain"), (storedOrigin, "origin") })
        {
            var source = label == "plain" ? derived : derivedOrigin;
            foreach (var (id, level) in pairs)
            {
                if (!source.TryGetValue(id, out var other) || other.WindowCount != level.WindowCount) { missing++; continue; }
                var a = level.Windows;
                var b = other.Windows;
                for (int w = 0; w < a.Count; w++)
                    worst = Math.Max(worst, Math.Max(Math.Abs(a[w].X - b[w].X), Math.Max(Math.Abs(a[w].Y - b[w].Y), Math.Abs(a[w].Z - b[w].Z))));
                compared++;
            }
        }
        Console.WriteLine($"compared {compared} levels, {missing} without a matching derived level, largest difference {worst:E2}");
    }
}
double ratio = timesDerived[^1] / timesStored[^1];
Console.WriteLine($"warm ratio derived / stored: {ratio:F2} -> {(ratio > 2 ? "keep the stored coarse rows" : "the coarse rows can be derived")}");
