#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using Npgsql;
using SequencePlatform.Core;

// P9 S9.2 (E2): a phage's host from what the platform stores. The set: viruses whose BV-BRC host_name starts with a stored bacterial
// or archaeal genus (provider lineage, not GTDB, which would merge Shigella into Escherichia and move the answer key). Two predictors:
// - A, composition: rank cellular taxa by the distance from the phage's signature (log rho*, Euclidean) to the nearest genome of each
//   taxon, over every stored genus (secondary: only the host genera). GC alone (|delta GC|) is the baseline; also the most common host
//   and a shuffle of host labels across phages.
// - B, shared sequence (--localize N | all): localize the phage over the cellular records (chromosomes and plasmids) with the windowed
//   profile; the best hit's genus is the prediction, counted as covered when its log10 E <= 0. A prophage in a stored genome is what
//   B can find; composition look-alikes are what the E cut is for.
// Fixed before the run:
// - A passes when top-1 genus accuracy minus GC's has a bootstrap 95 % lower bound above 0, both overall and macro (mean over host
//   genera with >= 5 phages), on the whole set. The set without records under 5 kbp is reported as a check.
// - B is kept for the full run when, on the sample, its top-1 genus precision over the covered phages (at least 10) is above A's
//   accuracy on those same phages.
// Writes <prefix>_PhageHost.csv, and with --localize <prefix>_Localize.csv (per phage).
const int minHostPhages = 5, shuffles = 200, bootstraps = 1000, seed = 20261104;
const long minLengthCheck = 5_000;
const double coveredLog10E = 0;
var prefix = args.FirstOrDefault(a => !a.StartsWith("--", StringComparison.Ordinal) && !char.IsDigit(a[0]) && a != "all") ?? "P9_E2";
int localizeAt = Array.IndexOf(args, "--localize");
int? localizeCount = localizeAt < 0 ? null : localizeAt + 1 < args.Length && int.TryParse(args[localizeAt + 1], out var count) ? count : int.MaxValue;
string[] ranks = ["genus", "family", "order"];

await using var dataSource = ExperimentCorpus.OpenSnapshot();
var corpus = (await ExperimentCorpus.LoadAsync(dataSource)).Where(g => !g.Excluded).ToList();
var cellular = corpus.Where(g => g.Group == SignatureIndex.Cellular && g.Provider.Genus is not null).ToList();
var taxa = ranks.Select((_, r) => new TaxonIndex(cellular.Select(g => Taxon(g.Provider, r)).ToArray())).ToArray();

// A host genus's family and order: the majority over its stored genomes.
var hostLineage = cellular.GroupBy(g => g.Provider.Genus!).ToDictionary(g => g.Key, g => new[]
{
    g.Key,
    g.Select(x => x.Provider.Family).Where(f => f is not null).GroupBy(f => f).OrderByDescending(f => f.Count()).ThenBy(f => f.Key, StringComparer.Ordinal).FirstOrDefault()?.Key,
    g.Select(x => x.Provider.Order).Where(o => o is not null).GroupBy(o => o).OrderByDescending(o => o.Count()).ThenBy(o => o.Key, StringComparer.Ordinal).FirstOrDefault()?.Key,
});
var phages = corpus.Where(g => g.Group == SignatureIndex.Viral && HostGenus(g) is { } h && hostLineage.ContainsKey(h))
    .Select(g => new Phage(g, hostLineage[HostGenus(g)!])).ToList();
var hostGenera = phages.Select(p => taxa[0].Of(p.Host[0])).ToHashSet();
Console.WriteLine($"{phages.Count} phages, {hostGenera.Count} host genera; {cellular.Count} cellular genomes in {taxa[0].Count} genera, " +
    $"{taxa[1].Count} families, {taxa[2].Count} orders.");

var watch = Stopwatch.StartNew();
var signature = Rankings(phages, cellular, taxa, hostGenera, (p, g) => Squared(p.Genome.Signature, g.Signature));
var gc = Rankings(phages, cellular, taxa, hostGenera, (p, g) => Math.Abs(p.Genome.Gc - g.Gc));
Console.WriteLine($"Rankings in {watch.Elapsed.TotalSeconds:F1} s.");

var truth = ranks.Select((_, r) => phages.Select(p => taxa[r].Of(p.Host[r])).ToArray()).ToArray();
var hostOf = phages.Select(p => p.Host[0]!).ToArray();
var macroGenera = hostOf.GroupBy(h => h).Where(h => h.Count() >= minHostPhages).Select(h => h.Key).ToHashSet();
var rows = new List<string>();
bool pass = false;

foreach (var (set, members) in new[]
{
    ("all", Enumerable.Range(0, phages.Count).ToArray()),
    ($"atLeast{minLengthCheck / 1000}kbp", Enumerable.Range(0, phages.Count).Where(i => phages[i].Genome.LengthBp >= minLengthCheck).ToArray()),
})
{
    for (int r = 0; r < ranks.Length; r++)
    {
        var scored = members.Where(i => truth[r][i] >= 0).ToArray();
        Add(set, "all", ranks[r], "signature", scored, i => signature.RankOf[r][i][truth[r][i]]);
        Add(set, "all", ranks[r], "gc", scored, i => gc.RankOf[r][i][truth[r][i]]);
        var common = scored.Select(i => truth[r][i]).GroupBy(t => t).OrderByDescending(t => t.Count()).ThenBy(t => t.Key).Select(t => t.Key).ToArray();
        Add(set, "all", ranks[r], "mostCommon", scored, i => Array.IndexOf(common, truth[r][i]) is var at && at >= 0 ? at + 1 : int.MaxValue);

        // Host labels shuffled across phages: what the signature's ranking scores against a label it has nothing to do with.
        var rng = new Random(seed + r);
        var top1 = new double[shuffles];
        var top5 = new double[shuffles];
        var labels = scored.Select(i => truth[r][i]).ToArray();
        for (int s = 0; s < shuffles; s++)
        {
            rng.Shuffle(labels);
            top1[s] = scored.Select((i, j) => signature.RankOf[r][i][labels[j]] <= 1 ? 1.0 : 0).Average();
            top5[s] = scored.Select((i, j) => signature.RankOf[r][i][labels[j]] <= 5 ? 1.0 : 0).Average();
        }
        Array.Sort(top1);
        Array.Sort(top5);
        rows.Add(Row(set, "all", ranks[r], "shuffledMean", scored.Length, top1.Average(), top5.Average(), null, null, null));
        rows.Add(Row(set, "all", ranks[r], "shuffled97.5", scored.Length, top1[(int)(0.975 * shuffles) - 1], top5[(int)(0.975 * shuffles) - 1], null, null, null));
    }
    Add(set, "hostGenera", "genus", "signature", members, i => signature.PoolRank[i]);
    Add(set, "hostGenera", "genus", "gc", members, i => gc.PoolRank[i]);

    var (overall, overallLow, overallHigh, macro, macroLow, macroHigh) = BootstrapDifference(members,
        members.Select(i => signature.RankOf[0][i][truth[0][i]] == 1).ToArray(), members.Select(i => gc.RankOf[0][i][truth[0][i]] == 1).ToArray(),
        members.Select(i => hostOf[i]).ToArray(), macroGenera);
    bool setPass = overallLow > 0 && macroLow > 0;
    if (set == "all") pass = setPass;
    rows.Add(string.Create(CultureInfo.InvariantCulture,
        $"{set}, all, genus, signatureMinusGc, {members.Length}, {overall:F4}, , {macro:F4}, , {macroGenera.Count}, {overallLow:F4}, {overallHigh:F4}, {macroLow:F4}, {macroHigh:F4}, {setPass}"));
    Console.WriteLine(rows[^1]);

    void Add(string set, string pool, string rank, string predictor, int[] scored, Func<int, int> rankOfTruth)
    {
        var hit1 = scored.Select(i => rankOfTruth(i) <= 1).ToArray();
        var hit5 = scored.Select(i => rankOfTruth(i) <= 5).ToArray();
        var macroScored = scored.Where(i => macroGenera.Contains(hostOf[i])).ToArray();
        double Macro(Func<int, bool> hit) => macroScored.GroupBy(i => hostOf[i]).Average(g => g.Average(i => hit(i) ? 1.0 : 0));
        var position = scored.Select((i, j) => (i, j)).ToDictionary(x => x.i, x => x.j);
        rows.Add(Row(set, pool, rank, predictor, scored.Length, hit1.Average(h => h ? 1.0 : 0), hit5.Average(h => h ? 1.0 : 0),
            Macro(i => hit1[position[i]]), Macro(i => hit5[position[i]]), macroScored.GroupBy(i => hostOf[i]).Count()));
        Console.WriteLine(rows[^1]);
    }
}
Console.WriteLine($"Verdict (A): {(pass ? "pass" : "fail")}.");

using (var csv = File.CreateText($"{prefix}_PhageHost.csv"))
{
    csv.WriteLine("Set, Pool, Rank, Predictor, Phages, Top1, Top5, MacroTop1, MacroTop5, MacroGenera, DiffLow, DiffHigh, MacroDiffLow, MacroDiffHigh, Pass");
    foreach (var row in rows) csv.WriteLine(row);
}

if (localizeCount is int wanted)
    await LocalizeAsync(wanted);

async Task LocalizeAsync(int wanted)
{
    var eligible = Enumerable.Range(0, phages.Count).Where(i => phages[i].Genome.LengthBp >= Localization.MinCatalogueLengthBp).ToArray();
    new Random(seed).Shuffle(eligible);
    var sample = eligible.Take(wanted).OrderBy(i => phages[i].Genome.Id).ToArray();
    var hostRecords = await HostRecordsAsync(dataSource);
    var profiles = await ProfileLevels.ReadManyAsync(dataSource, sample.Select(i => phages[i].Genome.Id).ToArray(),
        RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
    var localization = new Localization(dataSource);
    Console.WriteLine($"B: {sample.Length} of {eligible.Length} phages of at least {Localization.MinCatalogueLengthBp:N0} bp.");

    var results = new List<(int Phage, LocalizationHit? Hit, int Candidates, double Seconds)>();
    var total = Stopwatch.StartNew();
    foreach (var i in sample)
    {
        if (!profiles.TryGetValue(phages[i].Genome.Id, out var fine)) { results.Add((i, null, 0, 0)); continue; }
        var query = QueryProfile.FromStoredRecord(fine, phages[i].Genome.LengthBp);
        var step = Stopwatch.StartNew();
        LocalizationHit? best = null;
        int candidates = 0;
        foreach (var kind in new[] { "bacterium", "archaeon", "plasmid" })
        {
            var found = await localization.LocalizeAsync(query, new CatalogueFilter(Kind: kind), top: 1);
            candidates += found.Candidates;
            if (found.Hits.Count > 0 && (best is null || found.Hits[0].Score > best.Score)) best = found.Hits[0];
        }
        // E over every record scanned in the three calls, not the one call the hit came from.
        if (best?.P is { } p) best = best with { Log10E = p.Log10 + Math.Log10(candidates) };
        results.Add((i, best, candidates, step.Elapsed.TotalSeconds));
        if (results.Count % 25 == 0) Console.WriteLine($"  {results.Count}/{sample.Length}, {total.Elapsed.TotalMinutes:F1} min");
    }

    using var csv = File.CreateText($"{prefix}_Localize.csv");
    csv.WriteLine("PhageId, Name, LengthBp, HostGenus, HostFamily, HitId, HitKind, HitGenus, HitFamily, HitScore, HitLog10E, Covered, BGenusCorrect, BFamilyCorrect, " +
        "AGenus, AGenusCorrect, Candidates, Seconds");
    var covered = new List<(bool B, bool A)>();
    var anyHit = new List<(bool B, bool A)>();
    var combined = new List<bool>();
    foreach (var (i, hit, candidates, seconds) in results)
    {
        var phage = phages[i];
        var record = hit is null ? null : hostRecords.GetValueOrDefault(hit.SequenceId);
        bool isCovered = hit?.Log10E is { } e && e <= coveredLog10E;
        bool bGenus = record?.Genus == phage.Host[0], bFamily = record?.Family is not null && record.Family == phage.Host[1];
        string aGenus = taxa[0].Name(signature.Best[0][i]);
        bool aCorrect = aGenus == phage.Host[0];
        if (record is not null) anyHit.Add((bGenus, aCorrect));
        if (isCovered) covered.Add((bGenus, aCorrect));
        combined.Add(isCovered ? bGenus : aCorrect);
        csv.WriteLine(string.Create(CultureInfo.InvariantCulture,
            $"{phage.Genome.Id}, {Quote(phage.Genome.Name)}, {phage.Genome.LengthBp}, {phage.Host[0]}, {phage.Host[1]}, {hit?.SequenceId}, {record?.Kind}, {record?.Genus}, {record?.Family}, " +
            $"{hit?.Score:F4}, {hit?.Log10E:F2}, {isCovered}, {bGenus}, {bFamily}, {aGenus}, {aCorrect}, {candidates}, {seconds:F2}"));
    }

    double Share(IEnumerable<bool> values) => values.DefaultIfEmpty().Average(v => v ? 1.0 : 0);
    Console.WriteLine($"B: {total.Elapsed.TotalMinutes:F1} min, {results.Average(r => r.Seconds):F2} s per phage.");
    Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
        $"B covered (log10 E <= {coveredLog10E}): {covered.Count}/{results.Count}; B genus precision {Share(covered.Select(c => c.B)):F3}, A on the same {Share(covered.Select(c => c.A)):F3}."));
    Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
        $"B any hit: {anyHit.Count}; B genus {Share(anyHit.Select(c => c.B)):F3}, A on the same {Share(anyHit.Select(c => c.A)):F3}. A-or-B genus top-1 {Share(combined):F3}, A alone {Share(results.Select(r => taxa[0].Name(signature.Best[0][r.Phage]) == phages[r.Phage].Host[0])):F3}."));
    bool keep = covered.Count >= 10 && Share(covered.Select(c => c.B)) > Share(covered.Select(c => c.A));
    Console.WriteLine($"B kept for the full run: {keep}.");
}

static string? HostGenus(CorpusGenome g) => g.Field("host_name")?.Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault();
static string? Taxon(Lineage lineage, int rank) => rank switch { 0 => lineage.Genus, 1 => lineage.Family, _ => lineage.Order };
static string Quote(string text) => "\"" + text.Replace("\"", "\"\"") + "\"";

static double Squared(double[] a, double[] b)
{
    double sum = 0;
    for (int f = 0; f < a.Length; f++) { double diff = a[f] - b[f]; sum += diff * diff; }
    return sum;
}

static string Row(string set, string pool, string rank, string predictor, int phages, double top1, double top5, double? macro1, double? macro5, int? macroGenera) =>
    string.Create(CultureInfo.InvariantCulture,
        $"{set}, {pool}, {rank}, {predictor}, {phages}, {top1:F4}, {top5:F4}, {macro1:F4}, {macro5:F4}, {macroGenera}, , , , , ");

// Per phage and rank: every taxon's 1-based position when taxa are ordered by their nearest genome (Ties in order of taxon index);
// the best taxon; and the true genus's position among the host genera only.
static Ranking Rankings(List<Phage> phages, List<CorpusGenome> cellular, TaxonIndex[] taxa, HashSet<int> hostGenera, Func<Phage, CorpusGenome, double> distance)
{
    var rankOf = taxa.Select(_ => new int[phages.Count][]).ToArray();
    var best = taxa.Select(_ => new int[phages.Count]).ToArray();
    var poolRank = new int[phages.Count];
    Parallel.For(0, phages.Count, i =>
    {
        var nearest = taxa.Select(t => Enumerable.Repeat(double.PositiveInfinity, t.Count).ToArray()).ToArray();
        for (int g = 0; g < cellular.Count; g++)
        {
            double d = distance(phages[i], cellular[g]);
            for (int r = 0; r < taxa.Length; r++)
                if (taxa[r].OfGenome[g] is var t and >= 0 && d < nearest[r][t]) nearest[r][t] = d;
        }
        for (int r = 0; r < taxa.Length; r++)
        {
            var order = Enumerable.Range(0, taxa[r].Count).OrderBy(t => nearest[r][t]).ThenBy(t => t).ToArray();
            var positions = new int[order.Length];
            for (int j = 0; j < order.Length; j++) positions[order[j]] = j + 1;
            rankOf[r][i] = positions;
            best[r][i] = order[0];
        }
        int host = taxa[0].Of(phages[i].Host[0]);
        poolRank[i] = 1 + hostGenera.Count(u => rankOf[0][i][u] < rankOf[0][i][host]);
    });
    return new Ranking(rankOf, best, poolRank);
}

// Signature minus GC in top-1 genus accuracy, overall and macro, with a percentile bootstrap over phages.
static (double Overall, double OverallLow, double OverallHigh, double Macro, double MacroLow, double MacroHigh) BootstrapDifference(int[] members,
    bool[] a, bool[] b, string[] host, HashSet<string> macroGenera)
{
    var draw = new int[members.Length];
    (double Overall, double Macro) Difference(int[] picks)
    {
        double overall = picks.Average(j => (a[j] ? 1.0 : 0) - (b[j] ? 1.0 : 0));
        double macro = picks.Where(j => macroGenera.Contains(host[j])).GroupBy(j => host[j])
            .Average(g => g.Average(j => (a[j] ? 1.0 : 0) - (b[j] ? 1.0 : 0)));
        return (overall, macro);
    }
    var (overallAll, macroAll) = Difference(Enumerable.Range(0, members.Length).ToArray());
    var rng = new Random(20261104);
    var overalls = new double[bootstraps];
    var macros = new double[bootstraps];
    for (int s = 0; s < bootstraps; s++)
    {
        for (int j = 0; j < draw.Length; j++) draw[j] = rng.Next(members.Length);
        (overalls[s], macros[s]) = Difference(draw);
    }
    Array.Sort(overalls);
    Array.Sort(macros);
    int low = (int)(0.025 * bootstraps), high = (int)(0.975 * bootstraps) - 1;
    return (overallAll, overalls[low], overalls[high], macroAll, macros[low], macros[high]);
}

// Lineage of every cellular record and plasmid by id: what a localization hit's genus is.
static async Task<Dictionary<long, HostRecord>> HostRecordsAsync(NpgsqlDataSource dataSource)
{
    var records = new Dictionary<long, HostRecord>();
    await using var command = dataSource.CreateCommand(
        "select id, kind, tax_genus, tax_family from sequence_record where kind in ('bacterium', 'archaeon', 'plasmid')");
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync())
        records[reader.GetInt64(0)] = new HostRecord(reader.GetString(1), reader.IsDBNull(2) ? null : reader.GetString(2), reader.IsDBNull(3) ? null : reader.GetString(3));
    return records;
}

record Phage(CorpusGenome Genome, string?[] Host);

record Ranking(int[][][] RankOf, int[][] Best, int[] PoolRank);

record HostRecord(string Kind, string? Genus, string? Family);

sealed class TaxonIndex
{
    private readonly Dictionary<string, int> index = new(StringComparer.Ordinal);
    private readonly List<string> names = [];

    public int[] OfGenome { get; }
    public int Count => names.Count;

    public TaxonIndex(string?[] genomeTaxa)
    {
        OfGenome = genomeTaxa.Select(t =>
        {
            if (t is null) return -1;
            if (!index.TryGetValue(t, out var i)) { i = names.Count; index[t] = i; names.Add(t); }
            return i;
        }).ToArray();
    }

    public int Of(string? name) => name is not null && index.TryGetValue(name, out var i) ? i : -1;

    public string Name(int i) => names[i];
}
