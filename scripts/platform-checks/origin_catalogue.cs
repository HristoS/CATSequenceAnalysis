#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using Npgsql;
using SequencePlatform.Core;

// Rotation check (S3.2's open item): where do stored bacterial chromosomes start relative to their replication origin? The
// origin is estimated from the stored 500 bp profile alone, as the minimum of the cumulative GC skew, sum of (x - y)
// (StrandAndRotationBenchmark.OriginFromGcSkew; on the 30 bundled bacteria it found E. coli K-12's oriC at 3.93 Mbp and
// position 0 in the 20 files that start at dnaA). Records of kind 'bacterium' of at least 500 kbp, with BV-BRC's topology.
// Pass an output path (default S32_OriginCatalogue.csv); a summary goes to the console.
var output = args.Length > 0 ? args[0] : "S32_OriginCatalogue.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var meta = new Dictionary<long, (string? Topology, long Length)>();
await using (var cmd = dataSource.CreateCommand(
    "select id, annotations->'contig'->>'topology', length_bp from sequence_record where kind = 'bacterium' and length_bp >= 500000"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync()) meta[r.GetInt64(0)] = (r.IsDBNull(1) ? null : r.GetString(1), r.GetInt64(2));

var rows = new List<(string? Topology, double FromStart, double OriginToTerminus, double SkewPerKbp)>();
using var report = File.CreateText(output);
report.WriteLine("SequenceId, Topology, LengthBp, OriginBp, TerminusBp, OriginFromStartFraction, OriginToTerminusFraction, SkewAmplitudePerKbp");
await foreach (var (_, id, profile) in ProfileLevels.StreamAfterAsync(dataSource, RecordStore.Method, RecordStore.MethodVersion,
    RecordStore.FineWindowSize, 0, new LengthRange(500_000, null)))
{
    if (!meta.TryGetValue(id, out var m)) continue;
    var windows = profile.Windows;
    double sum = 0, min = 0, max = 0;
    int origin = 0, terminus = 0;
    for (int i = 0; i < windows.Count; i++)
    {
        sum += windows[i].X - windows[i].Y;
        if (sum < min) (min, origin) = (sum, (i + 1) % windows.Count);
        if (sum > max) (max, terminus) = (sum, (i + 1) % windows.Count);
    }
    long n = windows.Count;
    double Circular(long a, long b) { long d = Math.Abs(a - b) % n; return (double)Math.Min(d, n - d) / n; }
    // G - C excess between origin and terminus, in bases per kbp of genome: (x - y) / 2 × 500 bp per window.
    double skewPerKbp = (max - min) * 250 / (m.Length / 1000.0);
    rows.Add((m.Topology, Circular(origin, 0), Circular(origin, terminus), skewPerKbp));
    report.WriteLine($"{id}, {m.Topology}, {m.Length}, {origin * 500L}, {terminus * 500L}, {Circular(origin, 0):F4}, {Circular(origin, terminus):F3}, {skewPerKbp:F2}");
}

// A clear skew: the origin and terminus sit on opposite sides (0.4-0.5 of the ring apart) with at least 3 bp of G - C excess
// per kbp. In the bundled corpus only Deinococcus radiodurans (1.2) falls below it, the one genome whose estimate looks wrong;
// M. pneumoniae (3.8) and M. smegmatis (4.2) are found at their dnaA start. The split hardly moves between 0 and 7
// (results/S32_OriginCatalogue_Summary.csv).
Console.WriteLine("Topology, Records, ClearSkew, OriginWithin1PctOfStart, Within5Pct, Beyond5Pct, Beyond25Pct");
foreach (var g in rows.GroupBy(r => r.Topology ?? "(none)").Append(rows.GroupBy(_ => "all").Single()))
{
    var clear = g.Where(r => r.OriginToTerminus >= 0.4 && r.SkewPerKbp >= 3).ToList();
    Console.WriteLine($"{g.Key}, {g.Count()}, {clear.Count}, {clear.Count(r => r.FromStart <= 0.01)}, {clear.Count(r => r.FromStart <= 0.05)}, " +
        $"{clear.Count(r => r.FromStart > 0.05)}, {clear.Count(r => r.FromStart > 0.25)}");
}
