#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Buffers.Binary;
using System.Data;
using System.Globalization;
using Npgsql;
using SequencePlatform.Core;

// Duplicate genes (2026-11-05): the catalogue holds many genomes per species, so many stored gene profiles may be identical. How many,
// and what would storing each distinct profile once save? Reads every gene row (ZCGP v1: 16-byte header, 16 bytes per gene, then 9
// float16 feature blocks of windowCount values) and hashes each gene's window values (64-bit FNV-1a over its 18 bytes per window and
// its window count). Equal bytes mean an equal profile; a synonymous change at codon position 3, which the profile never sees, still
// counts as a duplicate. Read-only. Computed, not measured: the deduplicated size keeps every gene's header and adds a 4-byte
// reference per gene.
// Args: output path (default GeneDuplicates_Check.csv).
const int headerSize = 16, geneSize = 16, featureCount = 9, bytesPerWindow = 2 * featureCount, referenceBytes = 4;
var output = args.Length > 0 ? args[0] : "GeneDuplicates_Check.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var hashes = new List<ulong>(50_000_000);
var windowCounts = new List<ushort>(50_000_000);
var kinds = new List<byte>(50_000_000);
var kindNames = new List<string>();
long genes = 0, windowless = 0, records = 0;
await using (var command = dataSource.CreateCommand(
    "select r.kind, p.payload from profile_level p join sequence_record r on r.id = p.sequence_id where p.method = $1 and p.method_version = $2"))
{
    command.Parameters.AddWithValue(GeneProfiles.Method);
    command.Parameters.AddWithValue(GeneProfiles.MethodVersion);
    await using var reader = await command.ExecuteReaderAsync(CommandBehavior.SequentialAccess);
    while (await reader.ReadAsync())
    {
        var kind = reader.GetString(0);
        if (!kindNames.Contains(kind)) kindNames.Add(kind);
        byte kindIndex = (byte)kindNames.IndexOf(kind);
        var bytes = reader.GetFieldValue<byte[]>(1);
        records++;
        int geneCount = BinaryPrimitives.ReadInt32LittleEndian(bytes.AsSpan(8)), n = BinaryPrimitives.ReadInt32LittleEndian(bytes.AsSpan(12));
        int blocks = headerSize + geneSize * geneCount, first = 0;
        for (int g = 0; g < geneCount; g++)
        {
            int windows = BinaryPrimitives.ReadUInt16LittleEndian(bytes.AsSpan(headerSize + geneSize * g + 12));
            genes++;
            if (windows == 0)
            {
                windowless++;
                continue;
            }
            ulong hash = 14695981039346656037UL;
            void Mix(byte b) => hash = (hash ^ b) * 1099511628211UL;
            Mix((byte)windows);
            Mix((byte)(windows >> 8));
            for (int f = 0; f < featureCount; f++)
                for (int w = 0; w < windows; w++)
                {
                    int at = blocks + 2 * (f * n + first + w);
                    Mix(bytes[at]);
                    Mix(bytes[at + 1]);
                }
            hashes.Add(hash);
            windowCounts.Add((ushort)windows);
            kinds.Add(kindIndex);
            first += windows;
        }
    }
}
Console.WriteLine($"{records:N0} gene rows, {genes:N0} genes, {windowless:N0} without a window.");

var keys = hashes.ToArray();
var order = Enumerable.Range(0, keys.Length).ToArray();
Array.Sort(keys, order);
long distinct = 0, windowsAll = 0, windowsDistinct = 0, inGroupsOf2To10 = 0, inGroupsOver10 = 0, largest = 0;
var duplicatedByKind = new long[kindNames.Count];
var genesByKind = new long[kindNames.Count];
foreach (var k in kinds) genesByKind[k]++;
for (int i = 0; i < keys.Length;)
{
    int j = i;
    while (j + 1 < keys.Length && keys[j + 1] == keys[i]) j++;
    long size = j - i + 1;
    distinct++;
    windowsDistinct += windowCounts[order[i]];
    for (int t = i; t <= j; t++) windowsAll += windowCounts[order[t]];
    if (size > 1)
        for (int t = i; t <= j; t++) duplicatedByKind[kinds[order[t]]]++;
    if (size is >= 2 and <= 10) inGroupsOf2To10 += size;
    else if (size > 10) inGroupsOver10 += size;
    largest = Math.Max(largest, size);
    i = j + 1;
}

string Percent(long part, long whole) => (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
string Gb(double bytes) => (bytes / 1e9).ToString("F2", CultureInfo.InvariantCulture) + " GB";
long withWindows = keys.Length;
var lines = new List<string>
{
    "Measure,Value",
    $"Gene rows,{records}",
    $"Genes,{genes}",
    $"Genes without a window (under 167 aa),{windowless}",
    $"Genes with windows,{withWindows}",
    $"Distinct profiles among them,{distinct} ({Percent(distinct, withWindows)})",
    $"Genes whose profile occurs more than once,{Percent(inGroupsOf2To10 + inGroupsOver10, withWindows)} (2-10 copies {Percent(inGroupsOf2To10, withWindows)}; over 10 {Percent(inGroupsOver10, withWindows)}; largest group {largest})",
    $"Window bytes stored now,{Gb((double)windowsAll * bytesPerWindow)}",
    $"Window bytes with each distinct profile once (computed),{Gb((double)windowsDistinct * bytesPerWindow)} plus {Gb((double)withWindows * referenceBytes)} of references",
};
for (int k = 0; k < kindNames.Count; k++)
    lines.Add($"Genes with a duplicate: {kindNames[k]},{Percent(duplicatedByKind[k], genesByKind[k])} of {genesByKind[k]}");
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));
