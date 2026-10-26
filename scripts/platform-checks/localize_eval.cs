#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using System.Security.Cryptography;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S3.5: the same fragments as localize_check.cs's "min4" run (same seed, sources and draws), through Localization.LocalizeAsync
// with its caches, so its ranks can be compared with the check's (S35_LocalizeCheckMin4.csv). The first few fragments of each
// length are also localized without a cache (streamed from the database, as the CLI does) for timing.
// Pass an output path (default S35_LocalizeEval.csv); a per-length summary goes to the console.
const int seed = 20261022, trialsPerKind = 30, streamedPerLength = 2;
int[] fragmentLengths = { 4_000, 8_000, 10_000, 20_000, 32_000, 50_000, 100_000, 128_000, 200_000 };
var output = args.Length > 0 ? args[0] : "S35_LocalizeEval.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var records = new Dictionary<long, string?>();
var bySha = new Dictionary<string, long>();
await using (var cmd = dataSource.CreateCommand(
    "select r.id, t.tax_genus, r.sequence_sha256 from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        records[r.GetInt64(0)] = r.IsDBNull(1) ? null : r.GetString(1);
        if (!r.IsDBNull(2)) bySha.TryAdd(Convert.ToHexString(r.GetFieldValue<byte[]>(2)), r.GetInt64(0));
    }

string Sha(string dna) => Convert.ToHexString(SHA256.HashData(dna.Where(char.IsAsciiLetter).Select(c => (byte)char.ToUpperInvariant(c)).ToArray()));

var sources = new List<(string File, string Dna, long Id, bool Viral)>();
foreach (var file in Directory.GetFiles("data/genomes", "*.fna.gz").Concat(Directory.GetFiles("data/genomes/viral", "*.fna.gz")).Order())
{
    var dna = FastaReader.ReadPrimaryRecord(file).Sequence;
    if (bySha.TryGetValue(Sha(dna), out var id)) sources.Add((Path.GetFileName(file), dna, id, file.Contains("viral")));
}

var cached = new Localization(dataSource);
var streamed = new Localization(dataSource, useCache: false);
var warm = Stopwatch.StartNew();
foreach (var length in new[] { 4_000, 8_000, 32_000, 128_000 })   // one fragment per level: loads every cached level
    await cached.LocalizeAsync(QueryProfile.FromSequence(sources[0].Dna.Substring(0, length)));
Console.Error.WriteLine($"cache loaded: {cached.WindowBytes / 1e6:F1} MB of windows, {warm.Elapsed.TotalSeconds:F1} s");

var rng = new Random(seed);
var rows = new List<(int Length, bool Viral, int? Rank, bool SameGenus, double Ms)>();
var streamedMs = new Dictionary<int, List<double>>();
using var report = File.CreateText(output);
report.WriteLine("FragmentLength, Trial, Source, SourceKind, Start, Strand, CoarseWindowSize, Candidates, SourceRank, SourceScore, " +
    "Top1SameGenus, OffsetErrorWindows, Log10E, CachedMs, StreamedMs");
foreach (var length in fragmentLengths)
    foreach (var viral in new[] { false, true })
    {
        var eligible = sources.Where(s => s.Viral == viral && s.Dna.Length >= length + 500).ToArray();
        if (eligible.Length == 0) continue;
        for (int trial = 0; trial < trialsPerKind; trial++)
        {
            var source = eligible[rng.Next(eligible.Length)];
            int start = rng.Next(source.Dna.Length - length + 1);
            bool minus = rng.Next(2) == 1;
            var fragment = source.Dna.Substring(start, length);
            if (minus) fragment = QueryProfile.ReverseComplement(fragment);
            var query = QueryProfile.FromSequence(fragment);

            var watch = Stopwatch.StartNew();
            var result = await cached.LocalizeAsync(query, top: Localization.RefineDepth);
            var ms = watch.Elapsed.TotalMilliseconds;
            double? streamMs = null;
            if (trial < streamedPerLength && !viral)
            {
                watch.Restart();
                var fromDb = await streamed.LocalizeAsync(query, top: Localization.RefineDepth);
                streamMs = watch.Elapsed.TotalMilliseconds;
                if (!fromDb.Hits.SequenceEqual(result.Hits)) Console.Error.WriteLine($"streamed result differs: {length} {trial}");
                if (!streamedMs.TryGetValue(length, out var list)) streamedMs[length] = list = new List<double>();
                list.Add(streamMs.Value);
            }

            int index = result.Hits.ToList().FindIndex(h => h.SequenceId == source.Id);
            var hit = index < 0 ? null : result.Hits[index];
            bool sameGenus = records[source.Id] is { } genus && records[result.Hits[0].SequenceId] == genus;
            rows.Add((length, viral, index < 0 ? null : index + 1, sameGenus, ms));
            report.WriteLine($"{length}, {trial}, {source.File}, {(viral ? "virus" : "bacterium")}, {start}, {(minus ? '-' : '+')}, " +
                $"{result.CoarseWindowSize}, {result.Candidates}, {(index < 0 ? "" : index + 1)}, {hit?.Score:F4}, {sameGenus}, " +
                $"{(hit is null ? "" : (hit.OffsetBp - start) / 500.0):F2}, {hit?.Log10E:F2}, {ms:F1}, {streamMs:F1}");
            report.Flush();
        }
    }

double Median(IEnumerable<double> values)
{
    var sorted = values.Order().ToArray();
    return sorted.Length == 0 ? double.NaN : sorted[sorted.Length / 2];
}

Console.WriteLine("FragmentLength, Kind, Trials, InTop100, Top1, Top1SameGenus, MedianCachedMs, MaxCachedMs, MedianStreamedMs");
foreach (var g in rows.GroupBy(r => (r.Length, r.Viral)))
{
    int n = g.Count();
    var streamedText = !g.Key.Viral && streamedMs.TryGetValue(g.Key.Length, out var list) ? $"{Median(list):F0}" : "";
    Console.WriteLine($"{g.Key.Length}, {(g.Key.Viral ? "virus" : "bacterium")}, {n}, {g.Count(r => r.Rank is not null)}/{n}, " +
        $"{g.Count(r => r.Rank == 1)}/{n}, {g.Count(r => r.SameGenus)}/{n}, {Median(g.Select(r => r.Ms)):F0}, {g.Max(r => r.Ms):F0}, {streamedText}");
}
