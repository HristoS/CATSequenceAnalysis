#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using SequenceAnalyses;
using SequencePlatform.Core;

// Check before P9 S9.11 (E2b, phage host), arm C. Round-2 corpus, read only. Fixed before the run (2026-11-12):
// - Phages as E2: viruses whose host_name's first word is a stored cellular genus (provider lineage). Count them.
// - 200 of them, seeded, taken round-robin over host genera so the sample isn't all Escherichia. Their stored genes (gene_profile,
//   genes >= 50 aa): the share with a BV-BRC protein family. Plan's stop: under 20 % ends arm C.
// - Families only: a host genus's score is the most families any of its genomes (the corpus record: largest chromosome) shares with
//   the phage; no shared family is no prediction. Top-1 genus on the 200, against log rho* (nearest genome's genus) on the same 200.
//   Arm C is dropped if families only doesn't beat log rho* there.
// - Arm C timing: in the same pass over every cellular gene row, the first 3 sampled phages with genes have each gene scored against
//   every cellular gene within +-10 % of its length (or both >= 500 aa), as protein search compares; scoring time extrapolated to 200
//   phages, plus one pass of reading. Arm C is dropped if over 2 h.
// Writes <prefix>_Sample.csv (per phage).
const int sampleSize = 200, timedPhages = 3, seed = 20261112;
var prefix = args.FirstOrDefault() ?? "P9_E2b_Check";

await using var dataSource = ExperimentCorpus.OpenLive();
var corpus = (await ExperimentCorpus.LoadAsync(dataSource, ExperimentCorpus.RoundTwoMaxId)).Where(g => !g.Excluded).ToList();
var cellular = corpus.Where(g => g.Group == SignatureIndex.Cellular && g.Provider.Genus is not null).ToList();
var genera = cellular.Select(g => g.Provider.Genus!).ToHashSet(StringComparer.Ordinal);
var phages = corpus.Where(g => g.Group == SignatureIndex.Viral && HostGenus(g) is { } h && genera.Contains(h)).ToList();
Console.WriteLine($"{phages.Count} phages with a stored host genus ({phages.Select(p => HostGenus(p)).Distinct().Count()} host genera; round 1 had 4,419 in 69). " +
    $"{cellular.Count} cellular genomes in {genera.Count} genera.");

var rng = new Random(seed);
var queues = phages.GroupBy(p => HostGenus(p)!).OrderBy(g => g.Key, StringComparer.Ordinal)
    .Select(g => new Queue<CorpusGenome>(g.OrderBy(p => p.Id).OrderBy(_ => rng.Next()))).ToList();
var sample = new List<CorpusGenome>();
while (sample.Count < sampleSize && queues.Any(q => q.Count > 0))
    foreach (var queue in queues.Where(q => q.Count > 0).OrderBy(_ => rng.Next()).ToList())
        if (sample.Count < sampleSize) sample.Add(queue.Dequeue());
Console.WriteLine($"Sample: {sample.Count} phages over {sample.Select(p => HostGenus(p)).Distinct().Count()} host genera.");

// The sample's genes.
var phageGenes = sample.ToDictionary(p => p.Id, _ => new List<(ZCurveGeneProfile Profile, int Index)>());
await using (var command = dataSource.CreateCommand("select sequence_id, payload from gene_profile where sequence_id = any($1)"))
{
    command.Parameters.AddWithValue(sample.Select(p => p.Id).ToArray());
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync())
    {
        var profile = ZCurveGeneProfile.FromBinary(reader.GetFieldValue<byte[]>(1));
        for (int g = 0; g < profile.Genes.Count; g++) phageGenes[reader.GetInt64(0)].Add((profile, g));
    }
}
int genesTotal = phageGenes.Values.Sum(l => l.Count);
int genesWithFamily = phageGenes.Values.Sum(l => l.Count(x => x.Profile.Genes[x.Index].Family != ZCurveGeneProfile.NoFamily));
double familyShare = genesTotal == 0 ? 0 : genesWithFamily / (double)genesTotal;
var phageFamilies = sample.ToDictionary(p => p.Id, p => phageGenes[p.Id].Select(x => x.Profile.Genes[x.Index].Family)
    .Where(f => f != ZCurveGeneProfile.NoFamily).ToHashSet());
var wanted = phageFamilies.Values.SelectMany(f => f).ToHashSet();
Console.WriteLine($"{phageGenes.Count(p => p.Value.Count > 0)} of {sample.Count} phages have stored genes; {genesTotal} genes, {genesWithFamily} with a family " +
    $"({familyShare:P1}; stop under 20 %), {wanted.Count} distinct families.");

// One pass over every cellular gene row: shared families per genome, and arm C's scoring time for a few phages.
var timed = sample.Where(p => phageGenes[p.Id].Count > 0).Take(timedPhages)
    .SelectMany(p => phageGenes[p.Id].Select(x => (Query: x.Profile.WindowsOf(x.Index), x.Profile.Genes[x.Index].Residues))).ToArray();
var cellularIndex = cellular.Select((g, i) => (g.Id, i)).ToDictionary(x => x.Id, x => x.i);
var genomeFamilies = new HashSet<int>[cellular.Count];
for (int i = 0; i < genomeFamilies.Length; i++) genomeFamilies[i] = [];
var total = Stopwatch.StartNew();
var scoring = new Stopwatch();
long bytesRead = 0, comparisons = 0;
var scratch = new ZCurveGeneProfile.Gene[1];
await using (var command = dataSource.CreateCommand("select sequence_id, payload from gene_profile where sequence_id = any($1)"))
{
    command.Parameters.AddWithValue(cellular.Select(g => g.Id).ToArray());
    command.CommandTimeout = 3600;
    await using var reader = await command.ExecuteReaderAsync();
    while (await reader.ReadAsync())
    {
        var bytes = reader.GetFieldValue<byte[]>(1);
        bytesRead += bytes.Length;
        int count = ZCurveGeneProfile.GeneCountOf(bytes);
        if (scratch.Length < count) scratch = new ZCurveGeneProfile.Gene[count];
        var families = genomeFamilies[cellularIndex[reader.GetInt64(0)]];
        ZCurveGeneProfile.ScoreRow(bytes, [], scratch, gene => { if (wanted.Contains(gene.Family)) families.Add(gene.Family); return false; }, (_, _) => { });
        scoring.Start();
        foreach (var (query, residues) in timed)
        {
            bool queryLong = residues >= ZCurveGeneProfile.LongResidues;
            ZCurveGeneProfile.ScoreRow(bytes, query, scratch,
                gene => queryLong ? gene.Residues >= ZCurveGeneProfile.LongResidues : Math.Abs(gene.Residues - residues) <= 0.1 * residues,
                (_, _) => comparisons++);
        }
        scoring.Stop();
    }
}
double passSeconds = total.Elapsed.TotalSeconds - scoring.Elapsed.TotalSeconds;
double armCMinutes = (passSeconds + scoring.Elapsed.TotalSeconds * sampleSize / timedPhages) / 60;
Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
    $"Pass over cellular gene rows: {bytesRead / 1e9:F2} GB in {passSeconds:F0} s. Scoring {timed.Length} genes of {timedPhages} phages: {scoring.Elapsed.TotalSeconds:F1} s, " +
    $"{comparisons:N0} comparisons. Arm C on {sampleSize} phages extrapolated: {armCMinutes:F0} min (limit 120)."));

// Families only against log rho*, top-1 genus on the sample.
var rows = new List<string>();
int familiesHit = 0, familiesPredicted = 0, signatureHit = 0;
foreach (var phage in sample)
{
    var host = HostGenus(phage)!;
    var families = phageFamilies[phage.Id];
    var bestByGenus = new Dictionary<string, int>(StringComparer.Ordinal);
    for (int i = 0; i < cellular.Count; i++)
    {
        int shared = families.Count == 0 ? 0 : genomeFamilies[i].Count(families.Contains);
        if (shared > 0 && shared > bestByGenus.GetValueOrDefault(cellular[i].Provider.Genus!)) bestByGenus[cellular[i].Provider.Genus!] = shared;
    }
    var familyGuess = bestByGenus.OrderByDescending(b => b.Value).ThenBy(b => b.Key, StringComparer.Ordinal).Select(b => b.Key).FirstOrDefault();
    var signatureGuess = cellular.MinBy(g => Squared(g.Signature, phage.Signature))!.Provider.Genus;
    if (familyGuess is not null) familiesPredicted++;
    if (familyGuess == host) familiesHit++;
    if (signatureGuess == host) signatureHit++;
    rows.Add($"{phage.Id}, \"{host}\", {phageGenes[phage.Id].Count}, {families.Count}, \"{familyGuess}\", {(bestByGenus.Count == 0 ? 0 : bestByGenus.Values.Max())}, \"{signatureGuess}\"");
}
Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
    $"Top-1 genus on {sample.Count}: families only {familiesHit / (double)sample.Count:F3} ({familiesPredicted} with any shared family), log rho* {signatureHit / (double)sample.Count:F3}."));
using (var csv = File.CreateText($"{prefix}_Sample.csv"))
{
    csv.WriteLine("PhageId, HostGenus, Genes, Families, FamiliesGuess, SharedFamilies, SignatureGuess");
    foreach (var row in rows) csv.WriteLine(row);
}
Console.WriteLine($"Arm C: families {familyShare >= 0.2}, families beat log rho* {familiesHit > signatureHit}, time {armCMinutes <= 120}.");

static string? HostGenus(CorpusGenome g) => g.Field("host_name")?.Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault();

static double Squared(double[] a, double[] b)
{
    double sum = 0;
    for (int f = 0; f < a.Length; f++) { double d = a[f] - b[f]; sum += d * d; }
    return sum;
}
