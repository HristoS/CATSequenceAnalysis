#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:property PublishAot=false
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// 16-bit genome levels check (2026-11-05), before storing ZCWP levels as float16 on disk (about -1.28 GB, computed): does rounding every
// stored window to half precision change what find-similar ranks or where localization places a fragment?
// Find-similar: a seeded pool of 10,000 public records with a 500 bp level on the local catalogue; the first 300 are queries, each compared
// with every pool record in find-similar's default band [0.5 L, 2 L] by ZCurveProfileLadder.CompareForQuery, once with the stored float
// levels and once with query and candidates rounded to half (the ladder built from the rounded 500 bp level, which differs from rounding each
// stored level by at most one rounding step). Measured: top-1 and top-10 agreement, largest and mean score change.
// Localization: the bundled genomes of at least 100 kbp (data/genomes, run from the repository root), 20 fragments each of 4, 8 and 16 kbp
// from seeded base offsets (not window-aligned, as real fragments are), profiled from the DNA and placed by BestAlignment at 500 bp, float
// against half for both fragment and genome. Measured: offset agreement, and the share within one window of the truth for each.
// Rule fixed before the run: store halves if top-1 agreement and offset agreement are both >= 99%.
// Args: output path (default HalfLevels_Check.csv).
const string method = "zcurve-windowed", version = "1";
const int seed = 20261105, poolSize = 10_000, queryCount = 300, fragmentsPerLength = 20, fine = 500;
int[] fragmentLengths = [4_000, 8_000, 16_000];
var output = args.Length > 0 ? args[0] : "HalfLevels_Check.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

static ZCurveWindowedProfile Rounded(ZCurveWindowedProfile profile) => ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize,
    profile.Windows.Select(w => ((double)(Half)w.X, (double)(Half)w.Y, (double)(Half)w.Z)).ToArray());

var pool = new List<(long Id, long Length)>();
await using (var command = dataSource.CreateCommand(
    "select r.id, r.length_bp from sequence_record r join sequence_source s on s.id = r.source_id where s.owner_id is null and exists (" +
    "select 1 from profile_level p where p.sequence_id = r.id and p.method = $1 and p.method_version = $2 and p.window_size = $3) " +
    "order by md5(r.id::text || $4) limit $5"))
{
    command.Parameters.AddWithValue(method);
    command.Parameters.AddWithValue(version);
    command.Parameters.AddWithValue(fine);
    command.Parameters.AddWithValue(seed.ToString(CultureInfo.InvariantCulture));
    command.Parameters.AddWithValue(poolSize);
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync()) pool.Add((reader.GetInt64(0), reader.GetInt64(1)));
}
var stored = new Dictionary<long, ZCurveWindowedProfile>();
foreach (var chunk in pool.Chunk(1_000))
    foreach (var (id, level) in await ProfileLevels.ReadManyAsync(dataSource, chunk.Select(p => p.Id).ToArray(), method, version, fine))
        stored[id] = level;
var floatLadders = stored.ToDictionary(p => p.Key, p => new ZCurveProfileLadder(p.Value));
var halfLadders = stored.ToDictionary(p => p.Key, p => new ZCurveProfileLadder(Rounded(p.Value)));
Console.WriteLine($"{stored.Count} records loaded.");

var searches = new (bool Ranked, bool Top1, double Top10, double MaxChange, double SumChange, int Pairs)[queryCount];
Parallel.For(0, queryCount, q =>
{
    var (queryId, length) = pool[q];
    var scores = new List<(long Id, double Float, double Half)>();
    foreach (var (id, candidateLength) in pool)
    {
        if (id == queryId || candidateLength < SimilaritySearch.DefaultLengthRatio * length || candidateLength > length / SimilaritySearch.DefaultLengthRatio) continue;
        var exact = ZCurveProfileLadder.CompareForQuery(floatLadders[queryId], floatLadders[id]).Score;
        var rounded = ZCurveProfileLadder.CompareForQuery(halfLadders[queryId], halfLadders[id]).Score;
        if (exact is double e && rounded is double h) scores.Add((id, e, h));
    }
    if (scores.Count < 2) return;
    var byFloat = scores.OrderByDescending(s => s.Float).ThenBy(s => s.Id).Select(s => s.Id).ToList();
    var byHalf = scores.OrderByDescending(s => s.Half).ThenBy(s => s.Id).Select(s => s.Id).ToList();
    int k = Math.Min(10, scores.Count);
    searches[q] = (true, byFloat[0] == byHalf[0], byFloat.Take(k).Intersect(byHalf.Take(k)).Count() / (double)k,
        scores.Max(s => Math.Abs(s.Float - s.Half)), scores.Sum(s => Math.Abs(s.Float - s.Half)), scores.Count);
});
var ranked = searches.Where(s => s.Ranked).ToList();

var placements = new List<(int Length, bool Same, bool FloatRight, bool HalfRight)>();
foreach (var path in Directory.GetFiles(Path.Combine("data", "genomes"), "*.fna.gz").Order())
{
    var genome = FastaReader.ReadPrimaryRecord(path).Sequence;
    if (genome.Length < 100_000) continue;
    var reference = new ZCurveWindowedProfile(genome, fine);
    var roundedReference = Rounded(reference);
    var random = new Random(seed ^ Path.GetFileName(path).Length);
    foreach (var length in fragmentLengths)
        for (int f = 0; f < fragmentsPerLength; f++)
        {
            int start = random.Next(genome.Length - length);
            var fragment = new ZCurveWindowedProfile(genome.Substring(start, length), fine);
            int exact = ZCurveWindowedProfile.BestAlignment(fragment, reference).Offset;
            int rounded = ZCurveWindowedProfile.BestAlignment(Rounded(fragment), roundedReference).Offset;
            double truth = start / (double)fine;
            placements.Add((length, exact == rounded, Math.Abs(exact - truth) <= 1, Math.Abs(rounded - truth) <= 1));
        }
}

string Percent(double share) => (100 * share).ToString("F2", CultureInfo.InvariantCulture) + "%";
var lines = new List<string>
{
    "Part,Measure,Value",
    $"Find-similar,Queries ranked,{ranked.Count}",
    $"Find-similar,Median candidates per query,{ranked.Select(s => s.Pairs).Order().ElementAt(ranked.Count / 2)}",
    $"Find-similar,Top-1 agreement,{Percent(ranked.Average(s => s.Top1 ? 1.0 : 0))}",
    $"Find-similar,Top-10 overlap,{Percent(ranked.Average(s => s.Top10))}",
    string.Create(CultureInfo.InvariantCulture, $"Find-similar,Largest score change,{ranked.Max(s => s.MaxChange):E2}"),
    string.Create(CultureInfo.InvariantCulture, $"Find-similar,Mean score change,{ranked.Sum(s => s.SumChange) / ranked.Sum(s => s.Pairs):E2}"),
};
foreach (var group in placements.GroupBy(p => p.Length))
    lines.Add($"Localization {group.Key / 1000} kbp,Offset agreement / float within 1 window / half within 1 window ({group.Count()} fragments)," +
        $"{Percent(group.Average(p => p.Same ? 1.0 : 0))} / {Percent(group.Average(p => p.FloatRight ? 1.0 : 0))} / {Percent(group.Average(p => p.HalfRight ? 1.0 : 0))}");
double top1 = ranked.Average(s => s.Top1 ? 1.0 : 0), offsets = placements.Average(p => p.Same ? 1.0 : 0);
lines.Add($"Rule,Top-1 agreement >= 99% and offset agreement >= 99%,{Percent(top1)} and {Percent(offsets)}: {(top1 >= 0.99 && offsets >= 0.99 ? "pass: store halves" : "stop")}");
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));
