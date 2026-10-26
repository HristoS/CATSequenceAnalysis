#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using System.Security.Cryptography;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S3.5 pre-check: can a coarse scan of the whole catalogue find where a fragment comes from? Fragments are cut from
// the bundled genomes' DNA at a random start, half of them reverse-complemented, and their source record is found
// in the catalogue by sequence_sha256. Each fragment is scanned at the coarsest stored level its rule allows against
// every record at least as long, both query strands, at query phase 0 only and at all k phases; the top 100 are then
// refined at 500 bp within ±k windows of their coarse offset (as MultiResolutionBenchmark.CoarseToFine).
// Arguments: output path, then the rule. "half" is the first card's k × 500 ≤ L / 2 (S35_LocalizeCheck.csv);
// "min4" (default) keeps at least 4 coarse windows at phase 0 (S35_LocalizeCheckMin4.csv), with lengths at its band edges.
// A per-length summary goes to the console.
const int seed = 20261022, trialsPerKind = 30, refineTop = 100;
int[] factors = { 1, 4, 16, 64 };
var output = args.Length > 0 ? args[0] : "S35_LocalizeCheck.csv";
bool halfRule = args.Length > 1 && args[1] == "half";
int[] fragmentLengths = halfRule
    ? new[] { 2_000, 5_000, 10_000, 20_000, 50_000, 100_000 }
    : new[] { 4_000, 8_000, 10_000, 20_000, 32_000, 50_000, 100_000, 128_000, 200_000 };
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var caches = new Dictionary<int, IReadOnlyDictionary<long, ZCurveWindowedProfile>>();
foreach (var f in factors)
{
    var cache = new ProfileCache(dataSource, 500 * f, lengths: new LengthRange(fragmentLengths.Min(), null));
    await cache.RefreshAsync();
    caches[f] = cache.Profiles;
    Console.Error.WriteLine($"level {500 * f} bp: {cache.Profiles.Count:N0} records, {cache.WindowBytes / 1e6:F1} MB");
}

var records = new Dictionary<long, (long Length, string Kind, string? Genus)>();
var bySha = new Dictionary<string, long>();
await using (var cmd = dataSource.CreateCommand(
    "select r.id, r.length_bp, r.kind, t.tax_genus, r.sequence_sha256 from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        records[r.GetInt64(0)] = (r.GetInt64(1), r.GetString(2), r.IsDBNull(3) ? null : r.GetString(3));
        if (!r.IsDBNull(4)) bySha.TryAdd(Convert.ToHexString(r.GetFieldValue<byte[]>(4)), r.GetInt64(0));
    }

string Sha(string dna)
{
    var letters = dna.Where(char.IsAsciiLetter).Select(c => (byte)char.ToUpperInvariant(c)).ToArray();
    return Convert.ToHexString(SHA256.HashData(letters));
}

string Rc(string s)
{
    var c = new char[s.Length];
    for (int i = 0; i < s.Length; i++)
        c[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c', var o => o };
    return new string(c);
}

var sources = new List<(string File, string Dna, long Id, bool Viral)>();
int unmatched = 0;
foreach (var file in Directory.GetFiles("data/genomes", "*.fna.gz").Concat(Directory.GetFiles("data/genomes/viral", "*.fna.gz")).Order())
{
    var dna = FastaReader.ReadPrimaryRecord(file).Sequence;
    if (bySha.TryGetValue(Sha(dna), out var id)) sources.Add((Path.GetFileName(file), dna, id, file.Contains("viral")));
    else { unmatched++; Console.Error.WriteLine($"not in the catalogue: {file}"); }
}
Console.Error.WriteLine($"{sources.Count} bundled genomes matched by SHA-256, {unmatched} not");

var rng = new Random(seed);
var rows = new List<Row>();
using var report = File.CreateText(output);
report.WriteLine("FragmentLength, Trial, Source, SourceKind, SourceLength, Start, Strand, Factor, Candidates, CoarseRankPhase0, CoarseRankAllPhases, " +
    "SourceCoarseScore, TopCoarseScore, FinalRank, Top1SameGenus, SourceFineScore, OffsetErrorWindows, ScanMsPhase0, ScanMsAllPhases, RefineMs");

foreach (var length in fragmentLengths)
{
    int factor = factors.Where(f => halfRule ? f * 500 <= length / 2 : length / 500 / f >= 4).DefaultIfEmpty(1).Max();
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
            if (minus) fragment = Rc(fragment);
            var fine = new[] { new ZCurveWindowedProfile(fragment), new ZCurveWindowedProfile(Rc(fragment)) };
            var coarse = new ZCurveWindowedProfile[2, factor];
            for (int s = 0; s < 2; s++)
                for (int p = 0; p < factor; p++) coarse[s, p] = fine[s].Downsample(factor, p);

            var level = caches[factor];
            var candidates = level.Keys.Where(id => records[id].Length >= length && level[id].WindowCount >= coarse[0, 0].WindowCount).ToArray();

            var phase0 = new double[candidates.Length];
            var watch = Stopwatch.StartNew();
            Parallel.For(0, candidates.Length, i =>
                phase0[i] = Math.Max(Score(coarse[0, 0], level[candidates[i]]).Score, Score(coarse[1, 0], level[candidates[i]]).Score));
            var phase0Ms = watch.Elapsed.TotalMilliseconds;

            var all = new (double Score, int Strand, int Phase, int Offset)[candidates.Length];
            watch.Restart();
            Parallel.For(0, candidates.Length, i =>
            {
                var best = (Score: double.NegativeInfinity, Strand: 0, Phase: 0, Offset: -1);
                for (int s = 0; s < 2; s++)
                    for (int p = 0; p < factor; p++)
                    {
                        if (coarse[s, p].WindowCount == 0) continue;
                        var (score, offset) = Score(coarse[s, p], level[candidates[i]]);
                        if (score > best.Score) best = (score, s, p, offset);
                    }
                all[i] = best;
            });
            var allMs = watch.Elapsed.TotalMilliseconds;

            int sourceIndex = Array.IndexOf(candidates, source.Id);
            int rank0 = 1 + phase0.Count(x => x > phase0[sourceIndex]);
            int rankAll = 1 + all.Count(x => x.Score > all[sourceIndex].Score);

            var top = Enumerable.Range(0, candidates.Length).OrderByDescending(i => all[i].Score).Take(refineTop).ToArray();
            var refined = new (double Score, int Offset)[top.Length];
            watch.Restart();
            Parallel.For(0, top.Length, j =>
            {
                var (_, strand, phase, coarseOffset) = all[top[j]];
                var reference = caches[1][candidates[top[j]]];
                if (factor == 1) { refined[j] = (all[top[j]].Score, coarseOffset); return; }
                int estimate = coarseOffset * factor - phase;
                int from = Math.Max(0, estimate - factor);
                int to = Math.Min(reference.WindowCount, estimate + factor + fine[strand].WindowCount);
                if (to - from < fine[strand].WindowCount) { refined[j] = (double.NegativeInfinity, -1); return; }
                var neighbourhood = ZCurveWindowedProfile.FromPrecomputedWindows(500, reference.Windows.Skip(from).Take(to - from).ToList());
                var (score, offset) = ZCurveWindowedProfile.BestAlignment(fine[strand], neighbourhood);
                refined[j] = (score, from + offset);
            });
            var refineMs = watch.Elapsed.TotalMilliseconds;

            int inTop = Array.IndexOf(top, sourceIndex);
            int? finalRank = inTop < 0 ? null : 1 + refined.Count(x => x.Score > refined[inTop].Score);
            int winner = Enumerable.Range(0, top.Length).MaxBy(j => refined[j].Score);
            var sourceGenus = records[source.Id].Genus;
            bool sameGenus = sourceGenus is not null && records[candidates[top[winner]]].Genus == sourceGenus;
            double? offsetError = inTop < 0 ? null : refined[inTop].Offset - start / 500.0;
            double? sourceFine = inTop < 0 ? null : refined[inTop].Score;

            var row = new Row(length, viral, factor, rank0, rankAll, finalRank, sameGenus, offsetError, phase0Ms, allMs, refineMs);
            rows.Add(row);
            report.WriteLine($"{length}, {trial}, {source.File}, {(viral ? "virus" : "bacterium")}, {source.Dna.Length}, {start}, {(minus ? '-' : '+')}, " +
                $"{factor}, {candidates.Length}, {rank0}, {rankAll}, {all[sourceIndex].Score:F4}, {all[top[0]].Score:F4}, {finalRank}, {sameGenus}, " +
                $"{sourceFine:F4}, {offsetError:F2}, {phase0Ms:F1}, {allMs:F1}, {refineMs:F1}");
            report.Flush();
        }
    }
}

double Median(IEnumerable<double> values)
{
    var sorted = values.Order().ToArray();
    return sorted.Length == 0 ? double.NaN : sorted[sorted.Length / 2];
}

Console.WriteLine("FragmentLength, Kind, Factor, Trials, Phase0Recall1, Phase0Recall100, AllRecall1, AllRecall10, AllRecall100, FinalTop1, " +
    "FinalTop1AndOffset1, Top1SameGenus, MedianPhase0Ms, MedianAllMs, MedianRefineMs");
foreach (var g in rows.GroupBy(r => (r.Length, r.Viral)))
{
    int n = g.Count();
    string Rate(Func<Row, bool> hit) => $"{g.Count(hit)}/{n}";
    Console.WriteLine($"{g.Key.Length}, {(g.Key.Viral ? "virus" : "bacterium")}, {g.First().Factor}, {n}, {Rate(r => r.Rank0 == 1)}, {Rate(r => r.Rank0 <= 100)}, " +
        $"{Rate(r => r.RankAll == 1)}, {Rate(r => r.RankAll <= 10)}, {Rate(r => r.RankAll <= 100)}, {Rate(r => r.FinalRank == 1)}, " +
        $"{Rate(r => r.FinalRank == 1 && Math.Abs(r.OffsetError!.Value) <= 1)}, {Rate(r => r.SameGenus)}, " +
        $"{Median(g.Select(r => r.Phase0Ms)):F1}, {Median(g.Select(r => r.AllMs)):F1}, {Median(g.Select(r => r.RefineMs)):F1}");
}

// BestAlignment needs the query no longer than the reference; the candidate filter guarantees it for phase 0.
static (double Score, int Offset) Score(ZCurveWindowedProfile query, ZCurveWindowedProfile reference) =>
    query.WindowCount == 0 || query.WindowCount > reference.WindowCount ? (double.NegativeInfinity, -1) : ZCurveWindowedProfile.BestAlignment(query, reference);

record Row(int Length, bool Viral, int Factor, int Rank0, int RankAll, int? FinalRank, bool SameGenus, double? OffsetError,
    double Phase0Ms, double AllMs, double RefineMs);
