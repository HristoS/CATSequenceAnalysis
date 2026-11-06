#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Buffers.Binary;
using System.Data;
using System.Globalization;
using Npgsql;
using SequencePlatform.Core;

// Gene layer v2 sizes (2026-11-05), for the short-gene build card: computed from the stored gene headers (no windows needed), for the
// proposed layout. Genes under 50 aa are dropped; 50-499 aa get 5 equal windows over codon positions 1 and 2 (results/ShortProtein_*.csv);
// 500 aa and up keep today's 333-position windows. Windows stay 9 float16 values (18 bytes). Packed header per gene: start as an offset
// from the previous gene's start (genes sorted by start), length in bases, and (PGfam number + 1) * 2 + strand, each a LEB128 varint;
// the window count follows from the length. Short genes are grouped by length into buckets 50 * 1.1^b aa, so a query within +-10% of its
// length reads only the buckets that band touches; reported: the bytes a query of each length reads.
// Args: output path (default GeneLayerV2_Sizes.csv).
const int headerSize = 16, geneSize = 16, fixedWindow = 333, shortWindows = 5, minAa = 50, longAa = 500, windowBytes = 18;
const double bucketRatio = 1.1, band = 0.1;
var output = args.Length > 0 ? args[0] : "GeneLayerV2_Sizes.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

static int Varint(long value)
{
    int bytes = 1;
    while (value >= 128) { value >>= 7; bytes++; }
    return bytes;
}
static int Bucket(int aa) => (int)Math.Floor(Math.Log(aa / (double)minAa, bucketRatio));

string[] binNames = ["under 50 aa", "50-99 aa", "100-166 aa", "167-332 aa", "333-499 aa", "500 aa and up"];
int BinOf(int aa) => aa < 50 ? 0 : aa < 100 ? 1 : aa < 167 ? 2 : aa < 333 ? 3 : aa < 500 ? 4 : 5;
var genesByBin = new long[binNames.Length];
long records = 0, headerBytesNow = 0, packedHeaderBytes = 0, longWindows = 0, shortWindowCount = 0;
var bucketBytes = new Dictionary<int, long>();
long longBytes = 0;
await using (var command = dataSource.CreateCommand("select payload from profile_level where method = $1 and method_version = $2"))
{
    command.Parameters.AddWithValue(GeneProfiles.Method);
    command.Parameters.AddWithValue(GeneProfiles.MethodVersion);
    await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess);
    while (await reader.ReadAsync())
    {
        var bytes = reader.GetFieldValue<byte[]>(0);
        records++;
        int geneCount = BinaryPrimitives.ReadInt32LittleEndian(bytes.AsSpan(8));
        headerBytesNow += headerSize + (long)geneSize * geneCount;
        var genes = new List<(int Start, int End, int Family, bool Reverse)>(geneCount);
        for (int g = 0; g < geneCount; g++)
        {
            var slot = bytes.AsSpan(headerSize + geneSize * g);
            genes.Add((BinaryPrimitives.ReadInt32LittleEndian(slot), BinaryPrimitives.ReadInt32LittleEndian(slot[4..]),
                BinaryPrimitives.ReadInt32LittleEndian(slot[8..]), slot[14] == 1));
        }
        int previous = 0;
        foreach (var gene in genes.OrderBy(g => g.Start))
        {
            // Codons without the stop: BV-BRC's coordinates include it.
            int aa = Math.Max(0, (gene.End - gene.Start) / 3 - 1);
            genesByBin[BinOf(aa)]++;
            if (aa < minAa) continue;
            int header = Varint(gene.Start - previous) + Varint(gene.End - gene.Start) + Varint((gene.Family + 1L) * 2 + (gene.Reverse ? 1 : 0));
            previous = gene.Start;
            packedHeaderBytes += header;
            if (aa >= longAa)
            {
                int windows = 2 * aa / fixedWindow;
                longWindows += windows;
                longBytes += header + (long)windows * windowBytes;
            }
            else
            {
                shortWindowCount += shortWindows;
                int bucket = Bucket(aa);
                bucketBytes[bucket] = bucketBytes.GetValueOrDefault(bucket) + header + (long)shortWindows * windowBytes;
            }
        }
    }
}

string Gb(double value) => (value / 1e9).ToString("F2", CultureInfo.InvariantCulture) + " GB";
long totalGenes = genesByBin.Sum();
var lines = new List<string> { "Measure,Value" };
for (int b = 0; b < binNames.Length; b++)
    lines.Add($"Genes {binNames[b]},{genesByBin[b]} ({(100.0 * genesByBin[b] / totalGenes).ToString("F1", CultureInfo.InvariantCulture)}%)");
lines.Add($"Gene rows,{records}");
lines.Add($"Headers now (16 bytes per gene),{Gb(headerBytesNow)}");
lines.Add($"Packed headers (genes of 50 aa and up),{Gb(packedHeaderBytes)}");
lines.Add($"Long genes (500 aa and up): windows / bytes with headers,{longWindows} / {Gb(longBytes)}");
lines.Add($"Short genes (50-499 aa): windows / bytes with headers,{shortWindowCount} / {Gb(bucketBytes.Values.Sum())}");
lines.Add($"Gene layer v2 total,{Gb(longBytes + bucketBytes.Values.Sum())} (v1 now {Gb(headerBytesNow + 88_243_977L * windowBytes)})");
lines.Add($"Length buckets for short genes,{bucketBytes.Count}");
foreach (var aa in new[] { 50, 75, 100, 150, 200, 300, 400, 499 })
{
    int low = Bucket(Math.Max(minAa, (int)Math.Floor(aa * (1 - band)))), high = Bucket(Math.Min(longAa - 1, (int)Math.Ceiling(aa * (1 + band))));
    long read = Enumerable.Range(low, high - low + 1).Sum(b => bucketBytes.GetValueOrDefault(b));
    lines.Add($"Bytes a {aa} aa query reads (buckets {low}-{high}),{Gb(read)}");
}
lines.Add($"Bytes a 500 aa or longer query reads (long genes only),{Gb(longBytes)}");
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));
