#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// Origin rotation at scale (P3, for Paper B): build-origins rotates every circular bacterial chromosome >= 500 kbp with a clear skew
// to the minimum of its cumulative GC skew, a rule checked on K-12 and the 20 bundled genomes only (C3_OriginFromSkew.csv). Here every
// record the database holds rotated levels for is checked against an independent key, the start codon of its dnaA (BV-BRC product
// "Chromosomal replication initiator protein DnaA", any PGfam), read from the gene headers of gene_profile; no DNA. oriC sits next to
// dnaA in most bacteria but not all (E. coli: about 43 kbp away), so this measures agreement with dnaA, a lower bound on agreement with
// oriC. The unit is the 32 kbp window: rotated chromosomes of 512 kbp and up are scored at that level (SimilaritySearch.BandLevel at
// r = 0.5). Baselines: where the deposit starts (position 0, often dnaA by convention) and chance (2k windows of a ring of L).
// Stop rule (plan, 2026-11-12): under 80% within +-5 windows restricts rotation to the confident subset. dnaA gave 74.3%, but in
// Bacteroidota and Cyanobacteriota dnaA lies far from oriC; origin_doric_check.py takes this script's CSV against DoriC 12 instead.
// Args: output path (default S32_OriginDnaA.csv). Summary on the console.
const string dnaAProduct = "Chromosomal replication initiator protein DnaA";
const int coarseWindow = 32_000;
var output = args.Length > 0 ? args[0] : "S32_OriginDnaA.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var watch = System.Diagnostics.Stopwatch.StartNew();

var dnaAFamilies = new HashSet<int>();
await using (var cmd = dataSource.CreateCommand("select id from protein_family where lower(product) = lower($1)"))
{
    cmd.Parameters.AddWithValue(dnaAProduct);
    await using var r = await cmd.ExecuteReaderAsync();
    while (await r.ReadAsync()) dnaAFamilies.Add(r.GetInt32(0));
}

var records = new Dictionary<long, Record>();
await using (var cmd = dataSource.CreateCommand(
    "select r.id, r.length_bp, coalesce(g.tax_phylum, t.tax_phylum, '?'), coalesce(t.tax_genus, '?'), r.record_id from sequence_record r " +
    "left join taxon t on t.taxid = r.ncbi_taxid " +
    "left join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where exists (select 1 from profile_level p where p.sequence_id = r.id and p.method = $1 and p.method_version = $2)"))
{
    cmd.Parameters.AddWithValue(RecordStore.Method);
    cmd.Parameters.AddWithValue(OriginRotation.RotatedVersion(RecordStore.MethodVersion));
    await using var r = await cmd.ExecuteReaderAsync();
    while (await r.ReadAsync())
        records[r.GetInt64(0)] = new Record(r.GetInt64(0), r.GetInt64(1), r.GetString(2), r.GetString(3), r.GetString(4));
}
Console.Error.WriteLine($"{dnaAFamilies.Count} DnaA families; {records.Count:N0} records with rotated levels ({watch.Elapsed.TotalSeconds:F0} s)");

// The origin exactly as build-origins placed it, from the stored fine level.
await foreach (var (_, sequenceId, fine) in ProfileLevels.StreamAfterAsync(dataSource, RecordStore.Method, RecordStore.MethodVersion,
    RecordStore.FineWindowSize, 0, new LengthRange(OriginRotation.MinLengthBp, null)))
    if (records.TryGetValue(sequenceId, out var record) && OriginRotation.OriginWindow(fine, record.LengthBp) is int window)
    {
        record.OriginBp = (long)window * RecordStore.FineWindowSize;
        record.SkewBpPerKbp = fine.GcSkewExtremes().AmplitudeBp / (record.LengthBp / 1000.0);
    }
Console.Error.WriteLine($"origins placed for {records.Values.Count(r => r.OriginBp >= 0):N0} ({watch.Elapsed.TotalSeconds:F0} s)");

// dnaA starts from the gene headers: ScoreRow's candidate filter sees every gene and decodes no window.
var scratch = new ZCurveGeneProfile.Gene[1 << 16];
await using (var cmd = dataSource.CreateCommand("select sequence_id, payload from gene_profile where sequence_id = any($1)"))
{
    cmd.Parameters.AddWithValue(records.Keys.ToArray());
    await using var r = await cmd.ExecuteReaderAsync(System.Data.CommandBehavior.SequentialAccess);
    while (await r.ReadAsync())
    {
        var record = records[r.GetInt64(0)];
        record.HasGenes = true;
        var payload = r.GetFieldValue<byte[]>(1);
        int count = ZCurveGeneProfile.GeneCountOf(payload);
        if (count > scratch.Length) scratch = new ZCurveGeneProfile.Gene[count];
        ZCurveGeneProfile.ScoreRow(payload, [], scratch, gene =>
        {
            if (dnaAFamilies.Contains(gene.Family)) record.DnaAStarts.Add(gene.Reverse ? gene.End - 1 : gene.Start);
            return false;
        }, (_, _) => { });
    }
}
Console.Error.WriteLine($"gene rows read ({watch.Elapsed.TotalSeconds:F0} s)");

using (var report = File.CreateText(output))
{
    report.WriteLine("SequenceId, Accession, LengthBp, Phylum, Genus, SkewBpPerKbp, OriginBp, DnaACount, DnaAStartBp, OriginToDnaABp, DepositStartToDnaABp");
    foreach (var record in records.Values.OrderBy(r => r.Id))
    {
        var nearest = record.NearestDnaA();
        report.WriteLine(string.Join(", ", record.Id, Csv(record.Accession), record.LengthBp, Csv(record.Phylum), Csv(record.Genus),
            record.SkewBpPerKbp.ToString("F2", CultureInfo.InvariantCulture), record.OriginBp, record.DnaAStarts.Count,
            nearest?.ToString() ?? "", record.OriginToDnaA()?.ToString() ?? "", record.DepositStartToDnaA()?.ToString() ?? ""));
    }
}

var keyed = records.Values.Where(r => r.OriginBp >= 0 && r.DnaAStarts.Count > 0).ToList();
var single = keyed.Where(r => r.DnaAStarts.Count == 1).ToList();
Console.WriteLine($"rotated records {records.Count:N0}; origin placed {records.Values.Count(r => r.OriginBp >= 0):N0}; " +
    $"with a dnaA {keyed.Count:N0} (exactly one {single.Count:N0}); no gene row {records.Values.Count(r => r.OriginBp >= 0 && !r.HasGenes):N0}");
Summarize("one dnaA", single);
Summarize("any dnaA (nearest)", keyed);
Console.WriteLine();
Console.WriteLine("by phylum (one dnaA), within +-1 / +-5 coarse windows, median distance:");
foreach (var phylum in single.GroupBy(r => r.Phylum).OrderByDescending(g => g.Count()).Take(12))
{
    var d = phylum.Select(r => (long)r.OriginToDnaA()!).ToList();
    Console.WriteLine($"  {phylum.Key,-28} {phylum.Count(),6:N0}  {Share(d, 1),6:P1} {Share(d, 5),6:P1}  {Median(d) / 1000.0,8:F1} kbp");
}
Console.WriteLine();
Console.WriteLine("by skew strength (one dnaA), bp of G - C excess per kbp:");
foreach (var (min, max) in new[] { (3.0, 5.0), (5.0, 10.0), (10.0, 20.0), (20.0, double.MaxValue) })
{
    var d = single.Where(r => r.SkewBpPerKbp >= min && r.SkewBpPerKbp < max).Select(r => (long)r.OriginToDnaA()!).ToList();
    Console.WriteLine($"  [{min}, {(max == double.MaxValue ? "inf" : max.ToString())}) {d.Count,6:N0}  {Share(d, 1),6:P1} {Share(d, 5),6:P1}  {Median(d) / 1000.0,8:F1} kbp");
}
Console.Error.WriteLine($"done in {watch.Elapsed.TotalSeconds:F0} s; per-record rows in {output}");

static void Summarize(string label, List<Record> set)
{
    var origin = set.Select(r => (long)r.OriginToDnaA()!).ToList();
    var deposit = set.Select(r => (long)r.DepositStartToDnaA()!).ToList();
    double chance1 = set.Average(r => Math.Min(1, 2.0 * coarseWindow / r.LengthBp));
    double chance5 = set.Average(r => Math.Min(1, 10.0 * coarseWindow / r.LengthBp));
    double ring5 = set.Count(r => (double)r.OriginToDnaA()! / r.LengthBp <= 0.05) / (double)set.Count;
    Console.WriteLine($"{label} (n = {set.Count:N0}):");
    Console.WriteLine($"  skew origin:   within +-1 window {Share(origin, 1):P1}, +-5 {Share(origin, 5):P1}, 5% of ring {ring5:P1}; " +
        $"median {Median(origin) / 1000.0:F1} kbp, 90th percentile {Quantile(origin, 0.9) / 1000.0:F0} kbp; within 5 kbp {origin.Count(d => d <= 5_000) / (double)origin.Count:P1}");
    Console.WriteLine($"  deposit start: within +-1 window {Share(deposit, 1):P1}, +-5 {Share(deposit, 5):P1}; median {Median(deposit) / 1000.0:F1} kbp");
    Console.WriteLine($"  chance:        within +-1 window {chance1:P1}, +-5 {chance5:P1}");
}

static double Share(List<long> distances, int windows) =>
    distances.Count == 0 ? double.NaN : distances.Count(d => d <= (long)windows * coarseWindow) / (double)distances.Count;

static double Median(List<long> distances) => Quantile(distances, 0.5);

static double Quantile(List<long> distances, double q)
{
    if (distances.Count == 0) return double.NaN;
    var sorted = distances.Order().ToList();
    return sorted[(int)Math.Min(sorted.Count - 1, Math.Floor(q * sorted.Count))];
}

static string Csv(string text) => "\"" + text.Replace("\"", "\"\"") + "\"";

sealed class Record(long id, long lengthBp, string phylum, string genus, string accession)
{
    public long Id { get; } = id;
    public string Accession { get; } = accession;
    public long LengthBp { get; } = lengthBp;
    public string Phylum { get; } = phylum;
    public string Genus { get; } = genus;
    public long OriginBp { get; set; } = -1;
    public double SkewBpPerKbp { get; set; } = double.NaN;
    public List<long> DnaAStarts { get; } = [];
    public bool HasGenes { get; set; }

    public long? NearestDnaA() => OriginBp < 0 || DnaAStarts.Count == 0 ? null : DnaAStarts.MinBy(s => Distance(s, OriginBp));

    public long? OriginToDnaA() => NearestDnaA() is long start ? Distance(start, OriginBp) : null;

    // The deposit's own start against its dnaA nearest to position 0: how far the stored cut already is from the key.
    public long? DepositStartToDnaA() => DnaAStarts.Count == 0 ? null : DnaAStarts.Min(s => Distance(s, 0));

    private long Distance(long a, long b)
    {
        long d = Math.Abs(a - b) % LengthBp;
        return Math.Min(d, LengthBp - d);
    }
}
