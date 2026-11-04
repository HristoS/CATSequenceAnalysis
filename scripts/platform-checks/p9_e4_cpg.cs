#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// P9 S9.4 (E4): known biology in eukaryotes, and a check of the eukaryote counts pipeline. Vertebrate genomes are the textbook case
// of CpG depletion (methylated CpG mutates to TpG), and TpA is depleted broadly. Per genome rho*(CG) and rho*(TA), Karlin's both-strand
// odds ratio (1 = as composition predicts), from the stored counts; per class, the quartiles. Lineage is NCBI's (the provider's if
// missing). Rule fixed before the run: vertebrates pass when their 75th percentile of rho*(CG) is below the 25th percentile of insects
// and of fungi. Bacteria and viruses are summarized as context. Writes <prefix>_Classes.csv and <prefix>_Genomes.csv; with --publish,
// also the run and one finding per eukaryote into the live database (ExperimentStore), for the web UI.
const int cg = 4 * 1 + 2, ta = 4 * 3 + 0;   // index 4 * first + second, A C G T = 0..3
const int minPerClass = 3;
var prefix = args.FirstOrDefault(a => !a.StartsWith("--", StringComparison.Ordinal)) ?? "P9_E4_CpG";
bool publish = args.Contains("--publish");
var summaries = new List<object>();
await using var dataSource = ExperimentCorpus.OpenSnapshot();
var genomes = (await ExperimentCorpus.LoadAsync(dataSource)).Where(g => !g.Excluded).ToList();

var rows = genomes.Select(g =>
{
    var rho = GenomeSignatures.KarlinRho(g.Counts);
    var lineage = g.Ncbi ?? g.Provider;
    return new Row(g, lineage, rho[cg], rho[ta], g.Group == SignatureIndex.Eukaryotic && lineage.Phylum == "Chordata" && !nonVertebrateChordates.Contains(lineage.Class ?? ""));
}).ToList();
var eukaryotes = rows.Where(r => r.Genome.Group == SignatureIndex.Eukaryotic).ToList();

using (var csv = File.CreateText($"{prefix}_Genomes.csv"))
{
    csv.WriteLine("Id, Name, Kingdom, Phylum, Class, Vertebrate, RhoCG, RhoTA, Gc, LengthBp");
    foreach (var r in eukaryotes)
        csv.WriteLine(string.Create(CultureInfo.InvariantCulture,
            $"{r.Genome.Id}, {Quote(r.Genome.Name)}, {r.Lineage.Kingdom}, {r.Lineage.Phylum}, {r.Lineage.Class}, {r.Vertebrate}, {r.RhoCg:F4}, {r.RhoTa:F4}, {r.Genome.Gc:F4}, {r.Genome.LengthBp}"));
}

using var classes = File.CreateText($"{prefix}_Classes.csv");
classes.WriteLine("Set, Kingdom, Phylum, Class, Genomes, RhoCG_Q1, RhoCG_Median, RhoCG_Q3, RhoTA_Q1, RhoTA_Median, RhoTA_Q3");
void Summary(string set, string? kingdom, string? phylum, string? className, IReadOnlyList<Row> members)
{
    var cgs = members.Select(r => r.RhoCg).Order().ToArray();
    var tas = members.Select(r => r.RhoTa).Order().ToArray();
    var line = string.Create(CultureInfo.InvariantCulture,
        $"{set}, {kingdom}, {phylum}, {className}, {members.Count}, {Quantile(cgs, 0.25):F3}, {Quantile(cgs, 0.5):F3}, {Quantile(cgs, 0.75):F3}, " +
        $"{Quantile(tas, 0.25):F3}, {Quantile(tas, 0.5):F3}, {Quantile(tas, 0.75):F3}");
    classes.WriteLine(line);
    Console.WriteLine(line);
    summaries.Add(new { set, kingdom, phylum, @class = className, genomes = members.Count,
        rhoCg = new[] { Quantile(cgs, 0.25), Quantile(cgs, 0.5), Quantile(cgs, 0.75) }, rhoTa = new[] { Quantile(tas, 0.25), Quantile(tas, 0.5), Quantile(tas, 0.75) } });
}
foreach (var c in eukaryotes.Where(r => r.Lineage.Class is not null).GroupBy(r => (r.Lineage.Kingdom, r.Lineage.Phylum, r.Lineage.Class))
    .Where(c => c.Count() >= minPerClass).OrderBy(c => c.Select(r => r.RhoCg).Order().ElementAt(c.Count() / 2)))
    Summary("class", c.Key.Kingdom, c.Key.Phylum, c.Key.Class, c.ToList());

var vertebrates = eukaryotes.Where(r => r.Vertebrate).ToList();
var insects = eukaryotes.Where(r => r.Lineage.Class == "Insecta").ToList();
var fungi = eukaryotes.Where(r => r.Lineage.Kingdom == "Fungi").ToList();
var plants = eukaryotes.Where(r => r.Lineage.Kingdom == "Viridiplantae").ToList();
Summary("vertebrates", "Metazoa", "Chordata", null, vertebrates);
Summary("insects", "Metazoa", "Arthropoda", "Insecta", insects);
Summary("fungi", "Fungi", null, null, fungi);
Summary("plants", "Viridiplantae", null, null, plants);
Summary("bacteria and archaea", null, null, null, rows.Where(r => r.Genome.Group == SignatureIndex.Cellular).ToList());
Summary("viruses", null, null, null, rows.Where(r => r.Genome.Group == SignatureIndex.Viral).ToList());

double vertebrateQ3 = Quantile(vertebrates.Select(r => r.RhoCg).Order().ToArray(), 0.75);
double insectQ1 = Quantile(insects.Select(r => r.RhoCg).Order().ToArray(), 0.25);
double fungalQ1 = Quantile(fungi.Select(r => r.RhoCg).Order().ToArray(), 0.25);
bool pass = vertebrateQ3 < insectQ1 && vertebrateQ3 < fungalQ1;
Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
    $"Rule: vertebrate Q3 {vertebrateQ3:F3} < insect Q1 {insectQ1:F3} and fungal Q1 {fungalQ1:F3}: {(pass ? "pass" : "FAIL")} ({vertebrates.Count} vertebrates, {insects.Count} insects, {fungi.Count} fungi)."));

if (publish)
{
    var json = new System.Text.Json.JsonSerializerOptions { PropertyNamingPolicy = System.Text.Json.JsonNamingPolicy.CamelCase };
    var summary = System.Text.Json.JsonSerializer.Serialize(new
    {
        rule = "Vertebrates pass when their 75th percentile of ρ*(CG) is below the 25th percentile of insects and of fungi.",
        analysed = rows.Count,
        analysedNote = $"{eukaryotes.Count} eukaryotes (every one stored), with {rows.Count - eukaryotes.Count:N0} bacteria, archaea and viruses as context",
        keyNumbers = new object[]
        {
            new { label = "Vertebrates, 75th percentile of ρ*(CG)", value = vertebrateQ3 },
            new { label = "Insects, 25th percentile", value = insectQ1 },
            new { label = "Fungi, 25th percentile", value = fungalQ1 },
        },
        groups = summaries,
    }, json);
    var findings = eukaryotes.Select(r => new ExperimentFinding(r.Genome.Id, "cpg", r.RhoCg, System.Text.Json.JsonSerializer.Serialize(
        new { kingdom = r.Lineage.Kingdom, phylum = r.Lineage.Phylum, @class = r.Lineage.Class, vertebrate = r.Vertebrate, rhoCg = r.RhoCg, rhoTa = r.RhoTa }, json))).ToList();
    await using var live = NpgsqlDataSource.Create(ConnectionString.Resolve());
    var (stored, skipped) = await ExperimentStore.PublishAsync(live, "e4-cpg", "CpG depletion across eukaryotes",
        "Does the stored signature show the known CpG depletion of vertebrates, and does the eukaryote counting pipeline hold up?",
        pass ? ExperimentStore.Pass : ExperimentStore.Fail, ExperimentCorpus.SnapshotDatabase, "scripts/platform-checks/p9_e4_cpg.cs", summary, findings);
    Console.WriteLine($"Published e4-cpg: {stored} findings stored, {skipped} skipped (records gone or not public).");
}

// Linear interpolation between order statistics; the input is sorted.
static double Quantile(double[] sorted, double q)
{
    if (sorted.Length == 0) return double.NaN;
    double position = q * (sorted.Length - 1);
    int low = (int)Math.Floor(position);
    int high = Math.Min(low + 1, sorted.Length - 1);
    return sorted[low] + (position - low) * (sorted[high] - sorted[low]);
}

static string Quote(string text) => "\"" + text.Replace("\"", "\"\"") + "\"";

partial class Program
{
    // Chordates that aren't vertebrates: tunicates and lancelets.
    private static readonly HashSet<string> nonVertebrateChordates = new(StringComparer.Ordinal) { "Ascidiacea", "Appendicularia", "Thaliacea", "Leptocardii" };
}

record Row(CorpusGenome Genome, Lineage Lineage, double RhoCg, double RhoTa, bool Vertebrate);
