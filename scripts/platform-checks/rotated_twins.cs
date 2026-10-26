#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Rotation check A: how many stored records are the same molecule as another one, deposited from a different start (or strand)?
// Only those gain from rotating the DNA itself before profiling; for relatives, rotating the stored profile does as well. A
// rotation keeps the length and changes the dinucleotide counts (migration 0004) in two pairs at most, the old junction pair
// gained and the new one lost; the other strand maps each pair (a, b) to (b', a'). So a twin is a pair of records with the same
// length_bp, a different sequence_sha256 and dinucleotide counts within 2 in total, on either strand. "Near" (≤ 20) shows
// how many are a few changes away. Each twin is scored as stored and with both profiles rotated to their skew-minimum origin.
// Pass an output path (default S32_RotatedTwins.csv); a summary goes to the console.
var output = args.Length > 0 ? args[0] : "S32_RotatedTwins.csv";

// Sanity check of the rule on the bundled E. coli K-12 and its half rotation.
{
    var dna = FastaReader.ReadPrimaryRecord("data/genomes/ecoli_k12_mg1655_GCF_000005845.2.fna.gz").Sequence;
    var rotated = dna[(dna.Length / 2)..] + dna[..(dna.Length / 2)];
    var (shaA, countsA) = Features(dna);
    var (shaB, countsB) = Features(rotated);
    var (_, countsRc) = Features(QueryProfile.ReverseComplement(rotated));
    Console.Error.WriteLine($"K-12 vs half rotation: SHA-256 equal {shaA.SequenceEqual(shaB)}, dinucleotide L1 {L1(countsA, countsB)}; " +
        $"vs its reverse complement, flipped back: L1 {L1(countsA, Flip(countsRc))}");
}

await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var records = new List<(long Id, long SourceId, string Kind, string? Topology, long Length, byte[] Sha, long[] Counts)>();
await using (var cmd = dataSource.CreateCommand(
    "select id, source_id, kind, annotations->'contig'->>'topology', length_bp, sequence_sha256, dinucleotide_counts from sequence_record " +
    "where kind in ('bacterium', 'plasmid', 'archaeon') and length_bp >= 10000 and sequence_sha256 is not null and dinucleotide_counts is not null"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
        records.Add((r.GetInt64(0), r.GetInt64(1), r.GetString(2), r.IsDBNull(3) ? null : r.GetString(3), r.GetInt64(4),
            r.GetFieldValue<byte[]>(5), r.GetFieldValue<long[]>(6)));

var pairs = new List<(int A, int B, long Distance, bool OtherStrand)>();
foreach (var group in records.Select((rec, i) => (rec, i)).GroupBy(x => x.rec.Length).Where(g => g.Count() > 1))
{
    var members = group.ToArray();
    for (int a = 0; a < members.Length; a++)
        for (int b = a + 1; b < members.Length; b++)
        {
            var (x, y) = (members[a].rec, members[b].rec);
            if (x.Sha.SequenceEqual(y.Sha)) continue;
            long same = L1(x.Counts, y.Counts), other = L1(x.Counts, Flip(y.Counts));
            if (Math.Min(same, other) <= 20) pairs.Add((members[a].i, members[b].i, Math.Min(same, other), other < same));
        }
}

using var report = File.CreateText(output);
report.WriteLine("IdA, IdB, Kind, TopologyA, TopologyB, LengthBp, SameSource, DinucleotideL1, OtherStrand, " +
    "Stored500bp, Stored32kbp, Rotated500bp, Rotated32kbp");
foreach (var (a, b, distance, otherStrand) in pairs)
{
    var (x, y) = (records[a], records[b]);
    var fineX = await ProfileLevels.ReadAsync(dataSource, x.Id, RecordStore.Method, RecordStore.MethodVersion, 500);
    var fineY = await ProfileLevels.ReadAsync(dataSource, y.Id, RecordStore.Method, RecordStore.MethodVersion, 500);
    if (fineX is null || fineY is null) continue;
    if (otherStrand) fineY = FlipProfile(fineY);
    var stored = Score(fineX, fineY);
    var rotated = Score(ToOrigin(fineX), ToOrigin(fineY));
    report.WriteLine($"{x.Id}, {y.Id}, {x.Kind}, {x.Topology}, {y.Topology}, {x.Length}, {x.SourceId == y.SourceId}, {distance}, {otherStrand}, " +
        $"{stored.Fine:F4}, {stored.Coarse:F4}, {rotated.Fine:F4}, {rotated.Coarse:F4}");
}

var circularChromosomes = records.Count(r => r.Kind == "bacterium" && r.Topology == "circular" && r.Length >= 500_000);
Console.WriteLine($"records {records.Count:N0} (circular bacterial chromosomes >= 500 kbp: {circularChromosomes:N0}); same-length pairs with " +
    "different SHA-256 and dinucleotide L1 <= 20, by kind:");
Console.WriteLine("Kind, TwinsL1Max2, OfWhichOtherStrand, RecordsInTwins, NearL1Max20");
foreach (var g in pairs.GroupBy(p => records[p.A].Kind))
{
    var twins = g.Where(p => p.Distance <= 2).ToList();
    Console.WriteLine($"{g.Key}, {twins.Count}, {twins.Count(p => p.OtherStrand)}, {twins.SelectMany(p => new[] { p.A, p.B }).Distinct().Count()}, {g.Count()}");
}

static (byte[] Sha, long[] Counts) Features(string dna)
{
    using var reader = new SequenceFeatureReader(new StringReader(dna));
    ZCurveWindowedProfile.FromStream(reader);
    return reader.Finish();
}

static long L1(long[] a, long[] b) => a.Zip(b).Sum(p => Math.Abs(p.First - p.Second));

// Counts of the other strand: pair (a, b) becomes (3 - b, 3 - a) in A, C, G, T order.
static long[] Flip(long[] counts)
{
    var flipped = new long[16];
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++) flipped[4 * (3 - b) + (3 - a)] = counts[4 * a + b];
    return flipped;
}

static ZCurveWindowedProfile FlipProfile(ZCurveWindowedProfile p) =>
    ZCurveWindowedProfile.FromPrecomputedWindows(p.WindowSize, p.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

static (double Fine, double Coarse) Score(ZCurveWindowedProfile x, ZCurveWindowedProfile y) =>
    (ZCurveWindowedProfile.Compare(x, y), ZCurveWindowedProfile.Compare(x.Downsample(64), y.Downsample(64)));

// As StrandAndRotationBenchmark.ToOrigin: start at the window after the cumulative GC skew's minimum, sum of (x - y).
static ZCurveWindowedProfile ToOrigin(ZCurveWindowedProfile profile)
{
    var windows = profile.Windows;
    double sum = 0, min = 0;
    int origin = 0;
    for (int i = 0; i < windows.Count; i++)
    {
        sum += windows[i].X - windows[i].Y;
        if (sum < min) (min, origin) = (sum, (i + 1) % windows.Count);
    }
    return ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize, windows.Skip(origin).Concat(windows.Take(origin)).ToList());
}
