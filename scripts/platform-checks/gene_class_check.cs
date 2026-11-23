#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;
using SequenceAnalyses;
using SequencePlatform.Core;

// A "table" for genes (2026-11-22, S9.7 T1 gene variant). Mendeleev's table predicted properties of elements not yet seen from their
// place. The same question for genes: does a gene's place in gene-profile space predict a property of a gene from a family, and a
// phylum, the model never saw? The profile's x + y tracks hydropathy (results/README.md section 2, r -0.84 per window), so the property
// tested is membrane vs soluble protein, labelled from BV-BRC product names (independent of hydropathy). Database only, no downloads.
// Fixed before the run:
// - Corpus: round 2 (sequence_platform_p9b, fingerprint verified), bacteria and archaea, partial left out, one genome per GTDB genus
//   (lowest id), genera ordered by MD5 of the name, the first 1,000 with a gene profile.
// - Labels from protein_family.product, case-insensitive: membrane = permease, transmembrane, integral membrane, symporter, antiporter;
//   soluble = ribosomal protein [LS]<digit>, tRNA synthetase, tRNA ligase, DNA-directed RNA polymerase subunit, elongation factor,
//   chaperone. Names matching both lists, or "hypothetical", are dropped. Gate: fewer than 2,000 genes in either class -> stop unscored.
// - Features per gene: the 9 stored numbers (P0-P2 of x, y, z) averaged over its windows, and log protein length.
// - Arms: full (logistic regression on the 10 features); P0 x + y alone (membrane = lower, direction fixed by the per-window r, no fit);
//   GC + length (logistic on genome GC and log length, the baseline).
// - Splits: families and GTDB phyla each in two halves by MD5. Rule split: train on families A within phyla A, test on families B within
//   phyla B. Context: family split only, phylum split only, and every split with the halves swapped.
// - AUC with a 95% CI from 1,000 bootstrap resamples of test families; the difference to the baseline from the same resamples.
// - Pass: full arm's lower bound >= 0.80 and full minus GC + length >= 0.10 on the rule split. Otherwise stop: "the gene profile
//   organizes genes only by close homology", and no further gene-table work.
// C1, the harder contrast (--contrast enzymes, fixed 2026-11-22 before its run): the same rule with ordinary enzymes as the negative
// class. Enzyme = kinase, dehydrogenase, synthase, synthetase, isomerase, aldolase, mutase, epimerase, transferase, hydrolase,
// reductase; dropped when the name has NADH, quinone, cytochrome, ATP synthase, succinate or formate dehydrogenase, membrane, transport,
// permease, signal peptidase, secretion or hypothetical, or matches either list above (membrane and machinery keep their labels
// exactly as in the first run). --names-only prints each class's 30 commonest products and stops, to read them before scoring.
// Read before scoring (2026-11-22) and also dropped, being membrane-anchored in bacteria: sensor histidine kinase, serine/threonine
// protein kinase, glycosyltransferase / glycosyl transferase, 1-acyl-sn-glycerol-3-phosphate acyltransferase.
// Writes <prefix>_Auc_<contrast>.csv, and for the plots <prefix>_Families.csv (each family's mean over its genes: -(x + y), GC12 and
// log length), <prefix>_FamilyMap.csv (the same for families of at least 100 genes, small enough for the experiments page) and
// <prefix>_Histogram.csv (per-gene -(x + y) by class, then the single cuts), each row with the corpus fingerprint (FamilyMap: one note).
// C2 check (--transfer, fixed 2026-11-22 before its run): can one cut on -(x + y), set on bacteria, mark membrane-like genes elsewhere?
// Bacteria: the 1,000 genera above with a gene profile, bacterial only; archaea: every archaeal GTDB genus with one (first genome);
// viruses: up to 3,000 viral genomes with one, ordered by MD5 of the record id. Soluble = machinery + enzyme. The cut is the best
// balanced accuracy on bacteria. A group gets the mark if its -(x + y) AUC and its balanced accuracy at the bacterial cut are both
// within 0.05 of the bacteria's; a group with fewer than 200 genes in a class isn't testable and gets no mark. The cut must also hold
// within bacteria: balanced accuracy >= 0.80 in every tertile of genome GC, else stop (one fixed cut won't do). Eukaryote genes carry
// no stored names, so they can't be tested here. Writes <prefix>_Transfer.csv.
// C2 on eukaryotes and viruses (--transfer --live, fixed 2026-11-23 before its run): the same bacterial cut, tested on the live
// database's eukaryote and virus genes, labelled from each gene's own stored name (gene_name, migration 0025), not its family. Every
// record of the two kinds with a labelled gene. Names with "hypothetical" are dropped. Two membrane label sets:
// (b) names that state a membrane function: permease, symporter, antiporter, transporter, holin or antiholin, dropped when the name
//     also says ATP-binding, ATPase, substrate-binding, solute-binding, periplasmic or binding protein (soluble parts of transporters);
// (a) (b) plus the generic names of the first run (transmembrane, integral membrane), which annotators may set from a helix predictor,
//     that is from hydropathy itself, so (a) is context only.
// Soluble = machinery + enzyme as above. The rule is judged on (b): a group with fewer than 2,000 genes in a class or fewer than 30
// genomes isn't testable; otherwise it gets the mark if its AUC and its balanced accuracy at the bacterial cut are both within 0.05 of
// the bacteria's. 95% CIs from 1,000 bootstrap resamples of genomes (these genes have no family to resample). Viruses are also shown
// without holins (context). With --names-only, prints each class's 30 commonest names per group and stops, to read before scoring.
// Read before scoring (2026-11-23, --names-only) and also dropped: receptor-type kinases, single-pass membrane proteins common in plant
// annotations (receptor, LRR, leucine-rich repeat ... kinase, resistance kinase, wall-associated), and autotransporters (outer-membrane
// beta-barrels, not helical).
// Writes <prefix>_TransferLive.csv; the live database changes, so each row carries the run date instead of a fingerprint.
// Result: no mark for either (eukaryotes AUC 0.911 but balanced accuracy 0.810 at the bacterial cut). A eukaryote cut of their own
// (--transfer --live --fit-eukaryotes, fixed 2026-11-23 before its run): the same genes and label set (b); eukaryote genomes in two
// halves by MD5 of the accession. On half A, fit the cut of best balanced accuracy on -(x + y), a logistic model on the full profile
// (its own best cut) and the GC12 + length baseline; judge all, and the bacterial cut, on half B, overall and per kingdom (Fungi,
// Viridiplantae, Metazoa; a part with fewer than 2,000 genes in a class isn't scored). Swapped halves shown as context. Pass: the
// eukaryote cut's balanced accuracy on B >= the bacteria's minus 0.05 (0.844) overall and >= 0.80 in every scored kingdom; then the
// mark at that cut for eukaryotes, else none. Also per part (bacteria, eukaryotes, each kingdom): each class's mean and median
// -(x + y), GC12 and median length, with the part's own best cut, to see which side moved. Writes <prefix>_EukaryoteCut.csv.
// Result: 0.843, no mark; eukaryote membrane genes have half the bacterial mean -(x + y) and are 56% longer (dilution, read).
// The most hydrophobic window (--transfer --live --peak, fixed 2026-11-23 before its run): peak = max over a gene's windows of
// -(P0 x + P0 y). Same labels (set (b) for eukaryotes and viruses), same eukaryote halves. Arms in order; the first that passes for a
// group is the one to ship for it:
// 1. peak, one cut for all: the cut of best balanced accuracy on round-2 bacteria; kept for bacteria only if their balanced accuracy is
//    >= 0.874 (today's 0.894 minus 0.02) and >= 0.80 in every GC tertile; archaea, eukaryotes and viruses get it if their AUC and
//    balanced accuracy at that cut are both within 0.05 of the bacteria's (the C2 rule; gate 2,000 genes per class, 30 genomes);
// 2. peak, a eukaryote cut fitted on half A: on half B balanced accuracy >= 0.844 overall and >= 0.80 in every scored kingdom;
// 3. full profile + peak (9 window means, log length, peak; logistic, its best cut on A): the same bar as 2.
// Context: mean -(x + y) with a eukaryote cut, the full profile alone, GC12 + length, swapped halves, viruses without holins.
// Stop: no arm reaches 0.844 on held-out eukaryotes -> no eukaryote mark and the membrane-mark line ends. Also printed, the dilution
// check: each class's mean and median peak per part. Writes <prefix>_Peak.csv.
// Result: arm 1 passes for archaea and viruses, arm 3 for eukaryotes (0.869); the peak didn't close the gap (shipped per kind).
// Length bands (--transfer --live --length-bands, fixed 2026-11-23 before its run): is the window rule (5 equal windows below 500 aa,
// 166-aa windows from 500) behind the bacteria-eukaryote membrane gap? Same genes and labels as --peak (bacteria round 2 by family name;
// eukaryotes and viruses by their own names, set (b)). Bands by residues: 50-149, 150-299, 300-499, 500-699, 700-999, 1000+.
// 1. Per band and part (bacteria, eukaryotes, viruses): genes and mean of the mean -(x + y) and of the peak per class, AUC of each.
//    Confound: the length-matched gap (bacterial minus eukaryote membrane mean -(x + y) per band, weighted by eukaryote membrane genes
//    per band) against the raw 0.082; under half (0.041) -> length explains most of the gap.
// 2. A length-aware single number: per band, the cut of best balanced accuracy on bacteria, for the mean and for the peak; applied to
//    every eukaryote gene unchanged (nothing fitted on eukaryotes). Bar as for the eukaryote model: balanced accuracy >= 0.844 overall
//    and >= 0.80 in every kingdom with 2,000 genes per class -> a card to replace the model with it; else nothing changes. Archaea and
//    viruses at the same band cuts, context.
// 3. The switch at 500 aa: the same means for 450-499 against 500-549 aa, per part. Writes <prefix>_LengthBands.csv.
const int maxGenomes = 1000, minPerClass = 2000, maxViruses = 3000, minTransferPerClass = 200, minLiveGenomes = 30, resamples = 1000, seed = 20261122, bins = 60, namesShown = 30, mapMinGenes = 100;
const double lowerBoundRule = 0.80, marginRule = 0.10, transferTolerance = 0.05, tertileRule = 0.80;
var options = args.Where(a => a.StartsWith("--", StringComparison.Ordinal)).ToList();
var contrast = args.SkipWhile(a => a != "--contrast").Skip(1).FirstOrDefault() ?? "machinery";
if (contrast is not ("machinery" or "enzymes")) throw new ArgumentException("--contrast is machinery or enzymes.");
bool namesOnly = options.Contains("--names-only"), transfer = options.Contains("--transfer"), live = options.Contains("--live"),
    fitEukaryotes = options.Contains("--fit-eukaryotes"), peak = options.Contains("--peak"), lengthBands = options.Contains("--length-bands");
var prefix = args.Where((a, i) => !a.StartsWith("--", StringComparison.Ordinal) && (i == 0 || args[i - 1] != "--contrast")).FirstOrDefault() ?? "GeneClass";
var membrane = new Regex(@"permease|transmembrane|integral membrane|symporter|antiporter", RegexOptions.IgnoreCase);
var soluble = new Regex(@"ribosomal protein [LS]\d|tRNA synthetase|tRNA ligase|DNA-directed RNA polymerase subunit|elongation factor|chaperone",
    RegexOptions.IgnoreCase);
var enzyme = new Regex(@"kinase|dehydrogenase|synthase|synthetase|isomerase|aldolase|mutase|epimerase|transferase|hydrolase|reductase",
    RegexOptions.IgnoreCase);
var notEnzyme = new Regex(@"NADH|quinone|cytochrome|ATP synthase|succinate dehydrogenase|formate dehydrogenase|membrane|transport|permease|" +
    "signal peptidase|secretion|sensor histidine kinase|serine/threonine protein kinase|glycosyl ?transferase|" +
    "1-acyl-sn-glycerol-3-phosphate acyltransferase", RegexOptions.IgnoreCase);
var hypothetical = new Regex("hypothetical", RegexOptions.IgnoreCase);
string[] classNames = ["membrane", "machinery", "enzyme"];

await using var dataSource = ExperimentCorpus.OpenRoundTwo();
var verified = ExperimentCorpus.VerifyRoundTwoAsync(dataSource);
var genomes = (await ExperimentCorpus.LoadAsync(dataSource, ExperimentCorpus.RoundTwoMaxId))
    .Where(g => !g.Excluded && g.Group == SignatureIndex.Cellular && g.Gtdb?.Genus is not null && g.Gtdb.Phylum is not null)
    .GroupBy(g => g.Gtdb!.Genus!).Select(g => g.MinBy(x => x.Id)!)
    .OrderBy(g => Hash(g.Gtdb!.Genus!), StringComparer.Ordinal).ToList();

// Class per family: 0 membrane, 1 machinery (the first run's soluble class), 2 enzyme.
var labels = new Dictionary<int, (int Class, string Product)>();
await using (var command = dataSource.CreateCommand("select id, product from protein_family"))
await using (var reader = await command.ExecuteReaderAsync())
    while (await reader.ReadAsync())
    {
        var product = reader.GetString(1);
        if (hypothetical.IsMatch(product)) continue;
        bool isMembrane = membrane.IsMatch(product), isSoluble = soluble.IsMatch(product);
        if (isMembrane != isSoluble) labels[reader.GetInt32(0)] = (isMembrane ? 0 : 1, product);
        else if (!isMembrane && enzyme.IsMatch(product) && !notEnzyme.IsMatch(product)) labels[reader.GetInt32(0)] = (2, product);
    }
Console.WriteLine(string.Join(", ", classNames.Select((n, c) => $"{labels.Count(l => l.Value.Class == c)} {n}")) + $" families by name; {genomes.Count} GTDB genera.");

// One gene row per labelled gene of a genome: its windows averaged, the group standing in for the phylum.
async Task<List<GeneRow>> LabelledGenesAsync(CorpusGenome genome, string group)
{
    var rows = new List<GeneRow>();
    var profile = await GeneProfiles.ReadAsync(dataSource, genome.Id);
    if (profile is null) return rows;
    for (int g = 0; g < profile.Genes.Count; g++)
    {
        var gene = profile.Genes[g];
        if (gene.Family == ZCurveGeneProfile.NoFamily || !labels.TryGetValue(gene.Family, out var label) || gene.WindowCount == 0) continue;
        var windows = profile.WindowsOf(g);
        var mean = new double[ZCurveGeneProfile.FeatureCount];
        for (int w = 0; w < gene.WindowCount; w++)
            for (int f = 0; f < mean.Length; f++) mean[f] += windows[w * mean.Length + f] / gene.WindowCount;
        rows.Add(new GeneRow(gene.Family, group, label.Class, mean, Math.Log(gene.Residues), genome.Gc, PeakOf(windows, gene.WindowCount),
            genome.Id.ToString(CultureInfo.InvariantCulture)));
    }
    return rows;
}

if (transfer)
{
    var corpus = await ExperimentCorpus.LoadAsync(dataSource, ExperimentCorpus.RoundTwoMaxId);
    bool IsArchaeon(CorpusGenome g) => g.Gtdb?.Superkingdom?.Contains("Archaea", StringComparison.OrdinalIgnoreCase) == true;
    var groups = new List<(string Name, List<GeneRow> Genes, int Genomes)>();
    async Task AddGroupAsync(string name, IEnumerable<CorpusGenome> candidates, int cap)
    {
        var rows = new List<GeneRow>();
        int taken = 0;
        foreach (var genome in candidates)
        {
            if (taken == cap) break;
            var found = await LabelledGenesAsync(genome, name);
            if (found.Count == 0) continue;
            rows.AddRange(found);
            taken++;
        }
        groups.Add((name, rows, taken));
        Console.WriteLine($"{name}: {taken} genomes, {rows.Count(r => r.Membrane)} membrane and {rows.Count(r => !r.Membrane)} soluble genes.");
    }
    await AddGroupAsync("bacteria", genomes.Where(g => !IsArchaeon(g)), maxGenomes);
    await AddGroupAsync("archaea", genomes.Where(IsArchaeon), maxGenomes);
    await AddGroupAsync("viruses", corpus.Where(g => !g.Excluded && g.Group == SignatureIndex.Viral)
        .OrderBy(g => Hash(g.Id.ToString(CultureInfo.InvariantCulture)), StringComparer.Ordinal), maxViruses);
    var transferFingerprint = await verified;

    static double Score(GeneRow g) => -(g.Profile[0] + g.Profile[1]);
    static double BalancedAt(List<GeneRow> rows, double cut) =>
        (rows.Where(r => r.Membrane).Average(r => Score(r) >= cut ? 1.0 : 0) + rows.Where(r => !r.Membrane).Average(r => Score(r) < cut ? 1.0 : 0)) / 2;
    var bacteria = groups[0].Genes;
    double bacteriaAuc = SimpleAuc(bacteria.Select(Score).ToArray(), bacteria.Select(r => r.Membrane).ToArray());
    var (cut, bacteriaBalanced) = BalancedThreshold(bacteria.Select(Score).ToArray(), bacteria.Select(r => r.Membrane).ToArray());
    static string Marked(List<GeneRow> rows, double cut, bool membrane) =>
        F(rows.Where(r => r.Membrane == membrane).Average(r => Score(r) >= cut ? 1.0 : 0));
    var rowsOut = new List<string> { "Group,Genomes,Membrane genes,Soluble genes,AUC of -(x + y),Balanced accuracy at the bacterial cut," +
        "Membrane genes marked,Soluble genes marked,Marked genes that are membrane,Mark,Fingerprint" };
    string Shares(List<GeneRow> rows)
    {
        int tp = rows.Count(r => r.Membrane && Score(r) >= cut), fp = rows.Count(r => !r.Membrane && Score(r) >= cut);
        return string.Join(",", Marked(rows, cut, true), Marked(rows, cut, false), F(tp / (double)Math.Max(1, tp + fp)));
    }
    bool stop = false;
    foreach (var (name, rows, count) in groups)
    {
        int positives = rows.Count(r => r.Membrane), negatives = rows.Count - positives;
        string mark;
        double auc = double.NaN, balanced = double.NaN;
        if (positives < minTransferPerClass || negatives < minTransferPerClass) mark = $"not testable (fewer than {minTransferPerClass} genes in a class)";
        else
        {
            auc = SimpleAuc(rows.Select(Score).ToArray(), rows.Select(r => r.Membrane).ToArray());
            balanced = BalancedAt(rows, cut);
            mark = name == "bacteria" ? "reference" : auc >= bacteriaAuc - transferTolerance && balanced >= bacteriaBalanced - transferTolerance ? "yes" : "no";
        }
        rowsOut.Add(string.Join(",", name, count, positives, negatives, F(auc), F(balanced), double.IsNaN(auc) ? ",," : Shares(rows), mark, transferFingerprint));
    }
    var byGc = bacteria.Select(r => r.Gc).Order().ToArray();
    double[] edges = [byGc[byGc.Length / 3], byGc[2 * byGc.Length / 3]];
    string[] tertiles = ["low GC", "middle GC", "high GC"];
    for (int t = 0; t < 3; t++)
    {
        var part = bacteria.Where(r => (t == 0 ? r.Gc < edges[0] : t == 1 ? r.Gc >= edges[0] && r.Gc < edges[1] : r.Gc >= edges[1])).ToList();
        double balanced = BalancedAt(part, cut), auc = SimpleAuc(part.Select(Score).ToArray(), part.Select(r => r.Membrane).ToArray());
        if (balanced < tertileRule) stop = true;
        rowsOut.Add(string.Join(",", $"bacteria, {tertiles[t]}", "", part.Count(r => r.Membrane), part.Count(r => !r.Membrane), F(auc), F(balanced),
            Shares(part), balanced >= tertileRule ? "cut holds" : "cut fails", transferFingerprint));
    }
    rowsOut.Add("");
    rowsOut.Add($"Bacterial cut,-(x + y) >= {F(cut)} (balanced accuracy {F(bacteriaBalanced)}; GC tertile edges {F(edges[0])} and {F(edges[1])})");
    rowsOut.Add($"Verdict,{(stop ? "STOP: the fixed cut fails in a GC tertile of bacteria" : "the cut holds across bacterial GC; groups marked yes get the mark")}");
    if (!live)
    {
        File.WriteAllLines($"{prefix}_Transfer.csv", rowsOut);
        foreach (var line in rowsOut) Console.WriteLine(line);
        return;
    }
    Console.WriteLine($"Bacterial cut -(x + y) >= {F(cut)}: AUC {F(bacteriaAuc)}, balanced accuracy {F(bacteriaBalanced)}.");

    // Live: each gene's own name. Class codes: 0 membrane (b), 3 membrane (a) only, 1 machinery, 2 enzyme; holins flagged for the
    // context row.
    var membraneFunction = new Regex(@"permease|symporter|antiporter|transporter|(^|[^a-z])(anti)?holin", RegexOptions.IgnoreCase);
    var transporterPart = new Regex(@"ATP-binding|ATPase|substrate-binding|solute-binding|periplasmic|binding protein|autotransporter", RegexOptions.IgnoreCase);
    var receptorKinase = new Regex(@"receptor|LRR|leucine-rich repeat.*kinase|resistance kinase|wall-associated", RegexOptions.IgnoreCase);
    var holin = new Regex(@"(^|[^a-z])(anti)?holin", RegexOptions.IgnoreCase);
    int? NameClass(string name)
    {
        if (hypothetical.IsMatch(name)) return null;
        bool function = membraneFunction.IsMatch(name) && !transporterPart.IsMatch(name);
        bool generic = membrane.IsMatch(name);
        bool anyMembrane = function || generic || membraneFunction.IsMatch(name);
        if (soluble.IsMatch(name)) return anyMembrane ? null : 1;
        if (function) return 0;
        if (generic && !transporterPart.IsMatch(name)) return 3;
        if (anyMembrane) return null;
        return enzyme.IsMatch(name) && !notEnzyme.IsMatch(name) && !receptorKinase.IsMatch(name) ? 2 : null;
    }

    await using var liveSource = ExperimentCorpus.OpenLive();
    var records = new List<(long Id, string Kind, int[] Names)>();
    var accessionOf = new Dictionary<string, string>();
    var kingdomOf = new Dictionary<string, string>();
    await using (var command = liveSource.CreateCommand("select r.id, r.kind, n.names, r.record_id, coalesce(r.tax_kingdom, '-') " +
        "from gene_name n join sequence_record r on r.id = n.sequence_id where r.kind in ('eukaryote', 'virus')"))
    await using (var reader = await command.ExecuteReaderAsync())
        while (await reader.ReadAsync())
        {
            records.Add((reader.GetInt64(0), reader.GetString(1), reader.GetFieldValue<int?[]>(2).Select(n => n ?? 0).ToArray()));
            var key = reader.GetInt64(0).ToString(CultureInfo.InvariantCulture);
            accessionOf[key] = reader.GetString(3);
            kingdomOf[key] = reader.GetString(4);
        }
    var nameIds = records.SelectMany(r => r.Names).Where(n => n != 0).Distinct().ToArray();
    var nameClass = new Dictionary<int, (int Class, bool Holin, string Name)>();
    await using (var command = liveSource.CreateCommand("select id, name from gene_name_text where id = any($1)"))
    {
        command.Parameters.AddWithValue(nameIds);
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync())
        {
            var name = reader.GetString(1);
            if (NameClass(name) is int c) nameClass[reader.GetInt32(0)] = (c, holin.IsMatch(name), name);
        }
    }
    Console.WriteLine($"{records.Count} named eukaryote and virus records, {nameIds.Length} distinct names, {nameClass.Count} labelled.");

    // Per labelled gene: the group, the genome (for the bootstrap), the class code, the holin flag, and the name.
    var liveGenes = new List<(string Kind, GeneRow Row, bool Holin, string Name)>();
    int mismatched = 0;
    foreach (var (id, kind, names) in records)
    {
        if (!names.Any(n => nameClass.ContainsKey(n))) continue;
        var profile = await GeneProfiles.ReadAsync(liveSource, id);
        if (profile is null) continue;
        if (profile.Genes.Count != names.Length) { mismatched++; continue; }
        for (int g = 0; g < names.Length; g++)
        {
            if (!nameClass.TryGetValue(names[g], out var label)) continue;
            var gene = profile.Genes[g];
            if (gene.WindowCount == 0) continue;
            var windows = profile.WindowsOf(g);
            var mean = new double[ZCurveGeneProfile.FeatureCount];
            for (int w = 0; w < gene.WindowCount; w++)
                for (int f = 0; f < mean.Length; f++) mean[f] += windows[w * mean.Length + f] / gene.WindowCount;
            liveGenes.Add((kind, new GeneRow(-1, id.ToString(CultureInfo.InvariantCulture), label.Class, mean, Math.Log(gene.Residues), 0,
                PeakOf(windows, gene.WindowCount)), label.Holin, label.Name));
        }
    }
    Console.WriteLine($"{liveGenes.Count} labelled genes read; {mismatched} records skipped (name count differs from gene count).");

    string[] liveClassNames = ["membrane (b)", "machinery", "enzyme", "membrane (a) only"];
    if (namesOnly)
    {
        foreach (var kind in new[] { "eukaryote", "virus" })
            for (int c = 0; c < 4; c++)
            {
                var ofClass = liveGenes.Where(x => x.Kind == kind && x.Row.Class == c).ToList();
                Console.WriteLine($"{kind}, {liveClassNames[c]}: {ofClass.Count} genes in {ofClass.Select(x => x.Row.Phylum).Distinct().Count()} genomes; commonest names:");
                foreach (var group in ofClass.GroupBy(x => x.Name).OrderByDescending(x => x.Count()).Take(namesShown))
                    Console.WriteLine($"  {group.Count(),6}  {group.Key}");
            }
        return;
    }

    var today = DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture);
    if (fitEukaryotes)
    {
        FitEukaryotes();
        return;
    }
    if (lengthBands)
    {
        LengthBands();
        return;
    }
    if (peak)
    {
        PeakArms();
        return;
    }
    var liveOut =new List<string> { "Group,Label set,Genomes,Membrane genes,Soluble genes,AUC of -(x + y),AUC CI low,AUC CI high," +
        "Balanced accuracy at the bacterial cut,BA CI low,BA CI high,Membrane genes marked,Soluble genes marked,Marked genes that are membrane,Mark,Run date" };
    void ScoreGroup(string group, string set, List<GeneRow> rows, bool rule)
    {
        int positives = rows.Count(r => r.Membrane), negatives = rows.Count - positives;
        var genomeIds = rows.Select(r => r.Phylum).Distinct().ToArray();
        if (positives < minPerClass || negatives < minPerClass || genomeIds.Length < minLiveGenomes)
        {
            liveOut.Add(string.Join(",", group, set, genomeIds.Length, positives, negatives, ",,,,,,,,", rule ? $"not testable (fewer than {minPerClass} genes in a class or {minLiveGenomes} genomes)" : "context", today));
            return;
        }
        var scores = rows.Select(Score).ToArray();
        var truth = rows.Select(r => r.Membrane).ToArray();
        var order = Enumerable.Range(0, scores.Length).OrderBy(i => scores[i]).ToArray();
        var genomeIndex = genomeIds.Select((g, i) => (g, i)).ToDictionary(p => p.g, p => p.i);
        var geneGenome = rows.Select(r => genomeIndex[r.Phylum]).ToArray();
        double WeightedBalanced(double[] weight)
        {
            double tp = 0, p = 0, tn = 0, n = 0;
            for (int i = 0; i < scores.Length; i++)
                if (truth[i]) { p += weight[i]; if (scores[i] >= cut) tp += weight[i]; }
                else { n += weight[i]; if (scores[i] < cut) tn += weight[i]; }
            return (tp / p + tn / n) / 2;
        }
        var ones = Enumerable.Repeat(1.0, scores.Length).ToArray();
        double auc = Auc(scores, order, truth, ones), balanced = WeightedBalanced(ones);
        var bootAuc = new double[resamples];
        var bootBalanced = new double[resamples];
        Parallel.For(0, resamples, r =>
        {
            var rng = new Random(seed + r);
            var multiplicity = new double[genomeIds.Length];
            for (int i = 0; i < genomeIds.Length; i++) multiplicity[rng.Next(genomeIds.Length)]++;
            var weight = geneGenome.Select(g => multiplicity[g]).ToArray();
            bootAuc[r] = Auc(scores, order, truth, weight);
            bootBalanced[r] = WeightedBalanced(weight);
        });
        var (aucLow, aucHigh) = Interval(bootAuc);
        var (baLow, baHigh) = Interval(bootBalanced);
        string mark = !rule ? "context" : auc >= bacteriaAuc - transferTolerance && balanced >= bacteriaBalanced - transferTolerance ? "yes" : "no";
        liveOut.Add(string.Join(",", group, set, genomeIds.Length, positives, negatives, F(auc), F(aucLow), F(aucHigh), F(balanced), F(baLow), F(baHigh),
            Shares(rows), mark, today));
    }
    foreach (var kind in new[] { "eukaryote", "virus" })
    {
        var ofKind = liveGenes.Where(x => x.Kind == kind).ToList();
        ScoreGroup(kind, "(b) named membrane function", ofKind.Where(x => x.Row.Class != 3).Select(x => x.Row).ToList(), rule: true);
        ScoreGroup(kind, "(a) plus generic membrane names", ofKind.Select(x => x.Row.Class == 3 ? x.Row with { Class = 0 } : x.Row).ToList(), rule: false);
        if (kind == "virus")
            ScoreGroup(kind, "(b) without holins", ofKind.Where(x => x.Row.Class != 3 && !x.Holin).Select(x => x.Row).ToList(), rule: false);
    }
    liveOut.Add("");
    liveOut.Add($"Bacterial cut,-(x + y) >= {F(cut)} (round 2 bacteria: AUC {F(bacteriaAuc)}; balanced accuracy {F(bacteriaBalanced)}; fingerprint {transferFingerprint})");
    liveOut.Add($"Rule,(b) only: AUC >= {F(bacteriaAuc - transferTolerance)} and balanced accuracy >= {F(bacteriaBalanced - transferTolerance)}; " +
        $"at least {minPerClass} genes per class and {minLiveGenomes} genomes");
    File.WriteAllLines($"{prefix}_TransferLive.csv", liveOut);
    foreach (var line in liveOut) Console.WriteLine(line);
    return;

    // A cut of their own for eukaryotes, fitted on half the genomes and judged on the other half, per kingdom too; then each class's
    // means per kingdom and for bacteria, to see which side of the cut moved.
    void FitEukaryotes()
    {
        var euk = liveGenes.Where(x => x.Kind == "eukaryote" && x.Row.Class != 3).Select(x => x.Row).ToList();
        int Side(GeneRow r) => Half(accessionOf[r.Phylum]);
        static double Gc12Of(GeneRow g) => (1 + g.Profile[2]) / 2;
        Func<GeneRow, double[]> full = g => [.. g.Profile, g.LogLength], baseline = g => [Gc12Of(g), g.LogLength];
        string[] kingdoms = ["Fungi", "Viridiplantae", "Metazoa"];
        var lines = new List<string> { "Train half,Test part,Arm,Cut,Test genomes,Membrane genes,Soluble genes,AUC,AUC CI low,AUC CI high," +
            "Balanced accuracy,BA CI low,BA CI high,Membrane genes marked,Soluble genes marked,Run date" };
        var verdict = new List<string>();
        foreach (int trainHalf in new[] { 0, 1 })
        {
            var train = euk.Where(r => Side(r) == trainHalf).ToList();
            var test = euk.Where(r => Side(r) != trainHalf).ToList();
            var trainTruth = train.Select(r => r.Membrane).ToArray();
            var fullModel = Logistic.Fit(train, full);
            var baseModel = Logistic.Fit(train, baseline);
            var arms = new (string Name, Func<List<GeneRow>, double[]> Scores, double Cut)[]
            {
                ("-(x + y), bacterial cut", rows => rows.Select(Score).ToArray(), cut),
                ("-(x + y), eukaryote cut", rows => rows.Select(Score).ToArray(), BalancedThreshold(train.Select(Score).ToArray(), trainTruth).Threshold),
                ("full profile", rows => fullModel.Score(rows, full), BalancedThreshold(fullModel.Score(train, full), trainTruth).Threshold),
                ("GC12 + length", rows => baseModel.Score(rows, baseline), BalancedThreshold(baseModel.Score(train, baseline), trainTruth).Threshold),
            };
            foreach (var (part, rows) in new[] { ("all", test) }.Concat(kingdoms.Select(k => (k, test.Where(r => kingdomOf[r.Phylum] == k).ToList()))))
            {
                int positives = rows.Count(r => r.Membrane), negatives = rows.Count - positives;
                bool testable = positives >= minPerClass && negatives >= minPerClass;
                foreach (var (arm, scorer, armCut) in arms)
                {
                    var e = testable ? Evaluate(rows, scorer(rows), armCut) : default;
                    lines.Add(string.Join(",", trainHalf == 0 ? "A" : "B", part, arm, F(armCut), rows.Select(r => r.Phylum).Distinct().Count(), positives, negatives,
                        testable ? string.Join(",", F(e.Auc), F(e.AucLow), F(e.AucHigh), F(e.Balanced), F(e.BalancedLow), F(e.BalancedHigh), F(e.MembraneMarked), F(e.SolubleMarked))
                            : $"not testable (fewer than {minPerClass} genes in a class),,,,,,,", today));
                    if (trainHalf == 0 && arm == "-(x + y), eukaryote cut" && testable)
                    {
                        double needed = part == "all" ? bacteriaBalanced - transferTolerance : tertileRule;
                        verdict.Add($"{part} {F(e.Balanced)} (needs {F(needed)}){(e.Balanced >= needed ? "" : " FAILS")}");
                    }
                }
            }
        }
        bool pass = verdict.Count > 0 && verdict.All(v => !v.EndsWith("FAILS", StringComparison.Ordinal));
        lines.Add("");
        lines.Add($"Verdict (eukaryote cut fitted on A; balanced accuracy on B),{(pass ? "PASS" : "NO MARK")}: {string.Join("; ", verdict)}");

        // Which side moved: each class's -(x + y), GC12 and length, and the part's own best cut, for every eukaryote gene and for bacteria.
        lines.Add("");
        lines.Add("Part,Class,Genes,Mean -(x + y),Median -(x + y),Mean GC12,Median length (aa),Own best cut,Own balanced accuracy,AUC of -(x + y)");
        var parts = new List<(string Name, List<GeneRow> Rows)> { ("bacteria (round 2)", groups[0].Genes), ("eukaryotes", euk) };
        parts.AddRange(kingdoms.Select(k => (k, euk.Where(r => kingdomOf[r.Phylum] == k).ToList())));
        foreach (var (name, rows) in parts)
        {
            var scores = rows.Select(Score).ToArray();
            var truth = rows.Select(r => r.Membrane).ToArray();
            var (ownCut, ownBalanced) = BalancedThreshold(scores, truth);
            double auc = SimpleAuc(scores, truth);
            foreach (bool isMembrane in new[] { true, false })
            {
                var of = rows.Where(r => r.Membrane == isMembrane).ToList();
                var sorted = of.Select(Score).Order().ToArray();
                var lengths = of.Select(r => Math.Exp(r.LogLength)).Order().ToArray();
                lines.Add(string.Join(",", name, isMembrane ? "membrane" : "soluble", of.Count, F(sorted.Average()), F(sorted[sorted.Length / 2]),
                    F(of.Average(Gc12Of)), lengths[lengths.Length / 2].ToString("F0", CultureInfo.InvariantCulture), F(ownCut), F(ownBalanced), F(auc)));
            }
        }
        File.WriteAllLines($"{prefix}_EukaryoteCut.csv", lines);
        foreach (var line in lines) Console.WriteLine(line);
    }

    // The window rule as a confound: classes per length band, per-band cuts from bacteria applied to eukaryotes, and the 500 aa switch.
    void LengthBands()
    {
        static double PeakScore(GeneRow g) => g.Peak;
        static int Residues(GeneRow g) => (int)Math.Round(Math.Exp(g.LogLength));
        int[] edges = [50, 150, 300, 500, 700, 1000, int.MaxValue];
        int BandOf(GeneRow g) { int r = Residues(g), b = 0; while (r >= edges[b + 1]) b++; return b; }
        string BandName(int b) => edges[b + 1] == int.MaxValue ? $"{edges[b]}+" : $"{edges[b]}-{edges[b + 1] - 1}";
        var euk = liveGenes.Where(x => x.Kind == "eukaryote" && x.Row.Class != 3).Select(x => x.Row).ToList();
        var virus = liveGenes.Where(x => x.Kind == "virus" && x.Row.Class != 3).Select(x => x.Row).ToList();
        var parts = new (string Name, List<GeneRow> Rows)[] { ("bacteria", groups[0].Genes), ("eukaryotes", euk), ("viruses", virus) };
        string[] kingdoms = ["Fungi", "Viridiplantae", "Metazoa"];
        double MeanOf(IEnumerable<GeneRow> rows, Func<GeneRow, double> f) => rows.Any() ? rows.Average(f) : double.NaN;
        string FOrBlank(double v) => double.IsNaN(v) ? "" : F(v);
        double AucOf(List<GeneRow> rows, Func<GeneRow, double> f) => rows.Any(r => r.Membrane) && rows.Any(r => !r.Membrane)
            ? SimpleAuc(rows.Select(f).ToArray(), rows.Select(r => r.Membrane).ToArray()) : double.NaN;

        // 1. Per band and part.
        var lines = new List<string> { "Band (aa),Part,Membrane genes,Soluble genes,Membrane mean -(x + y),Soluble mean -(x + y),Membrane peak,Soluble peak," +
            "AUC of the mean,AUC of the peak,Run date" };
        for (int b = 0; b < edges.Length - 1; b++)
            foreach (var (name, all) in parts)
            {
                var rows = all.Where(r => BandOf(r) == b).ToList();
                var m = rows.Where(r => r.Membrane).ToList();
                var s = rows.Where(r => !r.Membrane).ToList();
                lines.Add(string.Join(",", BandName(b), name, m.Count, s.Count, FOrBlank(MeanOf(m, Score)), FOrBlank(MeanOf(s, Score)), FOrBlank(MeanOf(m, PeakScore)),
                    FOrBlank(MeanOf(s, PeakScore)), FOrBlank(AucOf(rows, Score)), FOrBlank(AucOf(rows, PeakScore)), today));
            }
        double rawGap = MeanOf(groups[0].Genes.Where(r => r.Membrane), Score) - MeanOf(euk.Where(r => r.Membrane), Score);
        double matched = 0, weight = 0;
        for (int b = 0; b < edges.Length - 1; b++)
        {
            var bact = groups[0].Genes.Where(r => r.Membrane && BandOf(r) == b).ToList();
            var eukBand = euk.Where(r => r.Membrane && BandOf(r) == b).ToList();
            if (bact.Count == 0 || eukBand.Count == 0) continue;
            matched += eukBand.Count * (MeanOf(bact, Score) - MeanOf(eukBand, Score));
            weight += eukBand.Count;
        }
        matched /= weight;
        lines.Add("");
        lines.Add($"Confound,raw membrane gap {F(rawGap)}; length-matched {F(matched)} ({(matched < rawGap / 2 ? "under half: length explains most of it" : "not under half: length doesn't explain most of it")})");

        // 2. Per-band cuts from bacteria, applied unchanged.
        foreach (var (scoreName, scorer) in new (string, Func<GeneRow, double>)[] { ("mean", Score), ("peak", PeakScore) })
        {
            var cuts = new double[edges.Length - 1];
            for (int b = 0; b < cuts.Length; b++)
            {
                var band = groups[0].Genes.Where(r => BandOf(r) == b).ToList();
                cuts[b] = BalancedThreshold(band.Select(scorer).ToArray(), band.Select(r => r.Membrane).ToArray()).Threshold;
            }
            double BandBalanced(List<GeneRow> rows) =>
                (rows.Where(r => r.Membrane).Average(r => scorer(r) >= cuts[BandOf(r)] ? 1.0 : 0) + rows.Where(r => !r.Membrane).Average(r => scorer(r) < cuts[BandOf(r)] ? 1.0 : 0)) / 2;
            var verdicts = new List<string>();
            var testParts = new List<(string Name, List<GeneRow> Rows, bool Rule)> { ("bacteria", groups[0].Genes, false), ("archaea", groups[1].Genes, false),
                ("eukaryotes", euk, true), ("viruses", virus, false) };
            testParts.AddRange(kingdoms.Select(k => (k, euk.Where(r => kingdomOf[r.Phylum] == k).ToList(), true)));
            foreach (var (name, rows, rule) in testParts)
            {
                int positives = rows.Count(r => r.Membrane), negatives = rows.Count - positives;
                if (positives < minPerClass || negatives < minPerClass) { verdicts.Add($"{name} not scored"); continue; }
                double balanced = BandBalanced(rows);
                double needed = name == "eukaryotes" ? bacteriaBalanced - transferTolerance : tertileRule;
                verdicts.Add(rule ? $"{name} {F(balanced)} ({(balanced >= needed ? "meets" : "FAILS")} {F(needed)})" : $"{name} {F(balanced)}");
            }
            bool pass = verdicts.Where(v => v.Contains("meets") || v.Contains("FAILS")).All(v => v.Contains("meets"));
            lines.Add($"Band cuts on the {scoreName},{string.Join(" ", cuts.Select(F))}");
            lines.Add($"Band cuts on the {scoreName}: {(pass ? "PASS (card to replace the eukaryote model)" : "no")},{string.Join("; ", verdicts)}");
        }

        // 3. The switch at 500 aa.
        lines.Add("");
        lines.Add("Switch,Part,Class,450-499 aa mean -(x + y),500-549 aa mean -(x + y),450-499 aa peak,500-549 aa peak,Genes 450-499,Genes 500-549");
        foreach (var (name, all) in parts)
            foreach (bool isMembrane in new[] { true, false })
            {
                var below = all.Where(r => r.Membrane == isMembrane && Residues(r) is >= 450 and < 500).ToList();
                var above = all.Where(r => r.Membrane == isMembrane && Residues(r) is >= 500 and < 550).ToList();
                lines.Add(string.Join(",", "500 aa", name, isMembrane ? "membrane" : "soluble", FOrBlank(MeanOf(below, Score)), FOrBlank(MeanOf(above, Score)),
                    FOrBlank(MeanOf(below, PeakScore)), FOrBlank(MeanOf(above, PeakScore)), below.Count, above.Count));
            }
        File.WriteAllLines($"{prefix}_LengthBands.csv", lines);
        foreach (var line in lines) Console.WriteLine(line);
    }

    // The most hydrophobic window: one cut from bacteria for every group, then eukaryote-fitted arms on held-out genomes.
    void PeakArms()
    {
        static double PeakScore(GeneRow g) => g.Peak;
        static double Gc12Of(GeneRow g) => (1 + g.Profile[2]) / 2;
        var euk = liveGenes.Where(x => x.Kind == "eukaryote" && x.Row.Class != 3).Select(x => x.Row).ToList();
        var virus = liveGenes.Where(x => x.Kind == "virus" && x.Row.Class != 3).ToList();
        var bacteriaRows = groups[0].Genes;
        string[] kingdoms = ["Fungi", "Viridiplantae", "Metazoa"];
        var lines = new List<string> { "Arm,Group or test part,Train half,Cut,Genomes,Membrane genes,Soluble genes,AUC,AUC CI low,AUC CI high," +
            "Balanced accuracy,BA CI low,BA CI high,Membrane genes marked,Soluble genes marked,Verdict,Run date" };
        // Kingdom parts of a held-out half are gated on genes only, as in --fit-eukaryotes (a half holds 11-20 genomes per kingdom).
        string Row(string arm, string part, string half, double at, List<GeneRow> rows, Func<GeneRow, double> score, Func<double, double, string>? verdictOf)
        {
            int positives = rows.Count(r => r.Membrane), negatives = rows.Count - positives, genomeCount = rows.Select(r => r.Genome ?? r.Phylum).Distinct().Count();
            int genomeGate = kingdoms.Contains(part) ? 0 : minLiveGenomes;
            if (positives < minPerClass || negatives < minPerClass || genomeCount < genomeGate)
            {
                lines.Add(string.Join(",", arm, part, half, F(at), genomeCount, positives, negatives, ",,,,,,,", "not testable", today));
                return "not testable";
            }
            var e = Evaluate(rows, rows.Select(score).ToArray(), at);
            var verdict = verdictOf?.Invoke(e.Auc, e.Balanced) ?? "context";
            lines.Add(string.Join(",", arm, part, half, F(at), genomeCount, positives, negatives, F(e.Auc), F(e.AucLow), F(e.AucHigh), F(e.Balanced),
                F(e.BalancedLow), F(e.BalancedHigh), F(e.MembraneMarked), F(e.SolubleMarked), verdict, today));
            return verdict;
        }

        // Arm 1: one cut, set on bacteria.
        var bacteriaPeaks = bacteriaRows.Select(PeakScore).ToArray();
        var bacteriaTruth = bacteriaRows.Select(r => r.Membrane).ToArray();
        var (peakCut, peakBalanced) = BalancedThreshold(bacteriaPeaks, bacteriaTruth);
        double peakAuc = SimpleAuc(bacteriaPeaks, bacteriaTruth);
        var byGc = bacteriaRows.Select(r => r.Gc).Order().ToArray();
        double[] edges = [byGc[byGc.Length / 3], byGc[2 * byGc.Length / 3]];
        bool tertilesHold = true;
        string Within(double auc, double balanced) => auc >= peakAuc - transferTolerance && balanced >= peakBalanced - transferTolerance ? "yes" : "no";
        Row("1 peak, bacterial cut", "bacteria", "-", peakCut, bacteriaRows, PeakScore, (_, b) => b >= bacteriaBalanced - 0.02 ? "keeps the bar" : "below the bar");
        for (int t = 0; t < 3; t++)
        {
            var part = bacteriaRows.Where(r => t == 0 ? r.Gc < edges[0] : t == 1 ? r.Gc >= edges[0] && r.Gc < edges[1] : r.Gc >= edges[1]).ToList();
            if (Row("1 peak, bacterial cut", $"bacteria GC tertile {t + 1}", "-", peakCut, part, PeakScore, (_, b) => b >= tertileRule ? "holds" : "fails") == "fails")
                tertilesHold = false;
        }
        var arm1 = new Dictionary<string, string>
        {
            ["archaea"] = Row("1 peak, bacterial cut", "archaea", "-", peakCut, groups[1].Genes, PeakScore, Within),
            ["eukaryotes"] = Row("1 peak, bacterial cut", "eukaryotes", "-", peakCut, euk, PeakScore, Within),
            ["viruses"] = Row("1 peak, bacterial cut", "viruses", "-", peakCut, virus.Select(x => x.Row).ToList(), PeakScore, Within),
        };
        Row("1 peak, bacterial cut", "viruses without holins", "-", peakCut, virus.Where(x => !x.Holin).Select(x => x.Row).ToList(), PeakScore, null);
        Row("context: mean, bacterial cut", "viruses", "-", cut, virus.Select(x => x.Row).ToList(), Score, null);

        // Arms 2 and 3, with context arms: fitted on one half of the eukaryote genomes, judged on the other.
        int Side(GeneRow r) => Half(accessionOf[r.Phylum]);
        Func<GeneRow, double[]> full = g => [.. g.Profile, g.LogLength], fullPeak = g => [.. g.Profile, g.LogLength, g.Peak], baseline = g => [Gc12Of(g), g.LogLength];
        var ruleVerdicts = new Dictionary<string, List<string>>();
        foreach (int trainHalf in new[] { 0, 1 })
        {
            var train = euk.Where(r => Side(r) == trainHalf).ToList();
            var test = euk.Where(r => Side(r) != trainHalf).ToList();
            var trainTruth = train.Select(r => r.Membrane).ToArray();
            var fullModel = Logistic.Fit(train, full);
            var fullPeakModel = Logistic.Fit(train, fullPeak);
            var baseModel = Logistic.Fit(train, baseline);
            Func<GeneRow, double> Model(Logistic model, Func<GeneRow, double[]> features) => g => model.Score([g], features)[0];
            var arms = new (string Name, Func<GeneRow, double> Score, bool Rule)[]
            {
                ("2 peak, eukaryote cut", PeakScore, true),
                ("3 full profile + peak", Model(fullPeakModel, fullPeak), true),
                ("context: mean, eukaryote cut", Score, false),
                ("context: full profile", Model(fullModel, full), false),
                ("context: GC12 + length", Model(baseModel, baseline), false),
            };
            foreach (var (name, scorer, rule) in arms)
            {
                double at = BalancedThreshold(train.Select(scorer).ToArray(), trainTruth).Threshold;
                foreach (var (part, rows) in new[] { ("eukaryotes", test) }.Concat(kingdoms.Select(k => (k, test.Where(r => kingdomOf[r.Phylum] == k).ToList()))))
                {
                    double needed = part == "eukaryotes" ? bacteriaBalanced - transferTolerance : tertileRule;
                    bool judged = rule && trainHalf == 0;
                    var verdict = Row(name, part, trainHalf == 0 ? "A" : "B", at, rows, scorer,
                        judged ? (_, b) => b >= needed ? $"meets {F(needed)}" : $"FAILS {F(needed)}" : null);
                    if (judged && verdict != "not testable") (ruleVerdicts.TryGetValue(name, out var list) ? list : ruleVerdicts[name] = []).Add($"{part} {verdict}");
                }
            }
        }

        bool bacteriaKeep = peakBalanced >= bacteriaBalanced - 0.02 && tertilesHold;
        lines.Add("");
        lines.Add($"Arm 1,peak cut {F(peakCut)} on bacteria (AUC {F(peakAuc)}; balanced accuracy {F(peakBalanced)}): bacteria {(bacteriaKeep ? "switch to the peak" : "keep the mean")}; " +
            string.Join("; ", arm1.Select(p => $"{p.Key} {p.Value}")));
        foreach (var (name, verdicts) in ruleVerdicts)
            lines.Add($"Arm {name[0]},{(verdicts.All(v => v.Contains("meets")) ? "PASS" : "no")}: {string.Join("; ", verdicts)}");
        string eukaryoteArm = arm1["eukaryotes"] == "yes" ? "arm 1"
            : ruleVerdicts.TryGetValue("2 peak, eukaryote cut", out var two) && two.All(v => v.Contains("meets")) ? "arm 2"
            : ruleVerdicts.TryGetValue("3 full profile + peak", out var three) && three.All(v => v.Contains("meets")) ? "arm 3" : "none: STOP";
        lines.Add($"Eukaryotes,{eukaryoteArm}");

        // For shipping (2026-11-23, the author's choice per organism): arm 3 refitted on every eukaryote genome, as plain weights on
        // [9 window means, ln residues, peak] with the intercept shifted so the cut is 0 (GeneProfiles.EukaryoteMembraneWeights).
        var allModel = Logistic.Fit(euk, fullPeak);
        var allScores = allModel.Score(euk, fullPeak);
        var (allCut, allBalanced) = BalancedThreshold(allScores, euk.Select(r => r.Membrane).ToArray());
        var weights = allModel.RawWeights();
        weights[0] -= allCut;
        double check = euk.Select((r, i) => Math.Abs(weights[0] + fullPeak(r).Select((v, j) => weights[j + 1] * v).Sum() - (allScores[i] - allCut))).Max();
        lines.Add($"Eukaryote model on all,{euk.Count} genes; balanced accuracy {F(allBalanced)} at 0; largest weight-form difference {check:E1}; weights " +
            string.Join(" ", weights.Select(w => w.ToString("R", CultureInfo.InvariantCulture))));

        // The dilution check: the peak per class and part.
        lines.Add("");
        lines.Add("Part,Class,Genes,Mean peak,Median peak,Mean -(x + y),Own best peak cut,Own balanced accuracy,AUC of the peak");
        var parts = new List<(string Name, List<GeneRow> Rows)> { ("bacteria (round 2)", bacteriaRows), ("archaea (round 2)", groups[1].Genes), ("eukaryotes", euk) };
        parts.AddRange(kingdoms.Select(k => (k, euk.Where(r => kingdomOf[r.Phylum] == k).ToList())));
        parts.Add(("viruses", virus.Select(x => x.Row).ToList()));
        foreach (var (name, rows) in parts)
        {
            var scores = rows.Select(PeakScore).ToArray();
            var truth = rows.Select(r => r.Membrane).ToArray();
            var (ownCut, ownBalanced) = BalancedThreshold(scores, truth);
            double auc = SimpleAuc(scores, truth);
            foreach (bool isMembrane in new[] { true, false })
            {
                var sorted = rows.Where(r => r.Membrane == isMembrane).Select(PeakScore).Order().ToArray();
                lines.Add(string.Join(",", name, isMembrane ? "membrane" : "soluble", sorted.Length, F(sorted.Average()), F(sorted[sorted.Length / 2]),
                    F(rows.Where(r => r.Membrane == isMembrane).Average(Score)), F(ownCut), F(ownBalanced), F(auc)));
            }
        }
        File.WriteAllLines($"{prefix}_Peak.csv", lines);
        foreach (var line in lines) Console.WriteLine(line);
    }

    (double Auc, double AucLow, double AucHigh, double Balanced, double BalancedLow, double BalancedHigh, double MembraneMarked, double SolubleMarked)
        Evaluate(List<GeneRow> rows, double[] scores, double at)
    {
        var truth = rows.Select(r => r.Membrane).ToArray();
        var order = Enumerable.Range(0, scores.Length).OrderBy(i => scores[i]).ToArray();
        var genomeIds = rows.Select(r => r.Genome ?? r.Phylum).Distinct().ToArray();
        var genomeIndex = genomeIds.Select((g, i) => (g, i)).ToDictionary(p => p.g, p => p.i);
        var geneGenome = rows.Select(r => genomeIndex[r.Genome ?? r.Phylum]).ToArray();
        double Balanced(double[] weight)
        {
            double tp = 0, p = 0, tn = 0, n = 0;
            for (int i = 0; i < scores.Length; i++)
                if (truth[i]) { p += weight[i]; if (scores[i] >= at) tp += weight[i]; }
                else { n += weight[i]; if (scores[i] < at) tn += weight[i]; }
            return (tp / p + tn / n) / 2;
        }
        var ones = Enumerable.Repeat(1.0, scores.Length).ToArray();
        var bootAuc = new double[resamples];
        var bootBalanced = new double[resamples];
        Parallel.For(0, resamples, r =>
        {
            var rng = new Random(seed + r);
            var multiplicity = new double[genomeIds.Length];
            for (int i = 0; i < genomeIds.Length; i++) multiplicity[rng.Next(genomeIds.Length)]++;
            var weight = geneGenome.Select(g => multiplicity[g]).ToArray();
            bootAuc[r] = Auc(scores, order, truth, weight);
            bootBalanced[r] = Balanced(weight);
        });
        var (aucLow, aucHigh) = Interval(bootAuc);
        var (balancedLow, balancedHigh) = Interval(bootBalanced);
        double membraneMarked = Enumerable.Range(0, scores.Length).Where(i => truth[i]).Average(i => scores[i] >= at ? 1.0 : 0);
        double solubleMarked = Enumerable.Range(0, scores.Length).Where(i => !truth[i]).Average(i => scores[i] >= at ? 1.0 : 0);
        return (Auc(scores, order, truth, ones), aucLow, aucHigh, Balanced(ones), balancedLow, balancedHigh, membraneMarked, solubleMarked);
    }
}

var genes = new List<GeneRow>();
int used = 0;
foreach (var genome in genomes)
{
    if (used == maxGenomes) break;
    var profile = await GeneProfiles.ReadAsync(dataSource, genome.Id);
    if (profile is null || profile.Genes.Count == 0) continue;
    used++;
    for (int g = 0; g < profile.Genes.Count; g++)
    {
        var gene = profile.Genes[g];
        if (gene.Family == ZCurveGeneProfile.NoFamily || !labels.TryGetValue(gene.Family, out var label) || gene.WindowCount == 0) continue;
        var windows = profile.WindowsOf(g);
        var mean = new double[ZCurveGeneProfile.FeatureCount];
        for (int w = 0; w < gene.WindowCount; w++)
            for (int f = 0; f < mean.Length; f++) mean[f] += windows[w * mean.Length + f] / gene.WindowCount;
        genes.Add(new GeneRow(gene.Family, genome.Gtdb!.Phylum!, label.Class, mean, Math.Log(gene.Residues), genome.Gc));
    }
}
Console.WriteLine($"{used} genomes, {genes.Count} labelled genes: " + string.Join(", ", classNames.Select((n, c) =>
    $"{genes.Count(g => g.Class == c)} {n} ({genes.Where(g => g.Class == c).Select(g => g.Family).Distinct().Count()} families)")));
foreach (int c in new[] { 0, 1, 2 })
{
    Console.WriteLine($"Commonest {classNames[c]} products:");
    foreach (var group in genes.Where(g => g.Class == c).GroupBy(g => labels[g.Family].Product).OrderByDescending(g => g.Count()).Take(namesShown))
        Console.WriteLine($"  {group.Count(),6}  {group.Key}");
}
var fingerprint = await verified;
WritePlotFiles(prefix, genes, labels, classNames, fingerprint, bins, mapMinGenes);
if (namesOnly) return;

int negativeClass = contrast == "enzymes" ? 2 : 1;
genes = genes.Where(g => g.Class == 0 || g.Class == negativeClass).ToList();
int membraneGenes = genes.Count(g => g.Membrane), solubleGenes = genes.Count - membraneGenes;
if (membraneGenes < minPerClass || solubleGenes < minPerClass)
{
    Console.WriteLine($"Stop: fewer than {minPerClass} genes in a class; not scored.");
    return;
}

var arms = new (string Name, Func<GeneRow, double[]>? Features, Func<GeneRow, double>? Fixed)[]
{
    ("full profile", g => [.. g.Profile, g.LogLength], null),
    ("P0 x + y alone", null, g => -(g.Profile[0] + g.Profile[1])),
    ("GC + length", g => [g.Gc, g.LogLength], null),
};
var splits = new (string Name, bool Rule, Func<GeneRow, int> Side)[]
{
    ("families and phyla", true, g => Half(g.Family.ToString(CultureInfo.InvariantCulture)) == Half(g.Phylum) ? Half(g.Phylum) : -1),
    ("families", false, g => Half(g.Family.ToString(CultureInfo.InvariantCulture))),
    ("phyla", false, g => Half(g.Phylum)),
};
var lines = new List<string> { "Split,Train half,Arm,Train genes,Test genes,Test families,AUC,CI low,CI high,Minus GC + length,Diff CI low,Diff CI high,Fingerprint" };
string? verdict = null;
foreach (var (splitName, rule, side) in splits)
    foreach (int trainHalf in new[] { 0, 1 })
    {
        var train = genes.Where(g => side(g) == trainHalf).ToList();
        var test = genes.Where(g => side(g) == 1 - trainHalf).ToList();
        var scores = arms.Select(a => a.Fixed is { } fixedScore ? test.Select(fixedScore).ToArray() : Logistic.Fit(train, a.Features!).Score(test, a.Features!)).ToArray();
        var families = test.Select(g => g.Family).Distinct().ToArray();
        var familyIndex = families.Select((f, i) => (f, i)).ToDictionary(p => p.f, p => p.i);
        var geneFamily = test.Select(g => familyIndex[g.Family]).ToArray();
        var truth = test.Select(g => g.Membrane).ToArray();
        var orders = scores.Select(s => Enumerable.Range(0, s.Length).OrderBy(i => s[i]).ToArray()).ToArray();
        var ones = Enumerable.Repeat(1.0, test.Count).ToArray();
        var point = scores.Select((s, a) => Auc(s, orders[a], truth, ones)).ToArray();
        var boot = new double[arms.Length][];
        for (int a = 0; a < arms.Length; a++) boot[a] = new double[resamples];
        Parallel.For(0, resamples, r =>
        {
            var rng = new Random(seed + r);
            var multiplicity = new double[families.Length];
            for (int i = 0; i < families.Length; i++) multiplicity[rng.Next(families.Length)]++;
            var weight = geneFamily.Select(f => multiplicity[f]).ToArray();
            for (int a = 0; a < arms.Length; a++) boot[a][r] = Auc(scores[a], orders[a], truth, weight);
        });
        int baseline = arms.Length - 1;
        for (int a = 0; a < arms.Length; a++)
        {
            var (low, high) = Interval(boot[a]);
            var (diffLow, diffHigh) = Interval(boot[a].Select((v, r) => v - boot[baseline][r]).ToArray());
            lines.Add(string.Join(",", splitName, trainHalf == 0 ? "A" : "B", arms[a].Name, train.Count, test.Count, families.Length,
                F(point[a]), F(low), F(high), F(point[a] - point[baseline]), F(diffLow), F(diffHigh), fingerprint));
            if (rule && trainHalf == 0 && a == 0)
                verdict = low >= lowerBoundRule && point[a] - point[baseline] >= marginRule
                    ? $"PASS: full profile AUC {F(point[a])} (CI {F(low)}-{F(high)}), {F(point[a] - point[baseline])} above GC + length"
                    : $"STOP: full profile AUC {F(point[a])} (CI {F(low)}-{F(high)}), {F(point[a] - point[baseline])} above GC + length " +
                      $"(rule: lower bound >= {lowerBoundRule}, margin >= {marginRule})";
        }
    }
lines.Add("");
lines.Add($"Verdict (membrane vs {classNames[negativeClass]}; families and phyla, train A),{verdict}");
File.WriteAllLines($"{prefix}_Auc_{contrast}.csv", lines);
foreach (var line in lines) Console.WriteLine(line);

// The plot files: each family's mean over its genes, with the family-level AUC of -(x + y) next to the gene-level one as the check
// that centroids keep the separation; and per-gene -(x + y) binned by class, with the threshold of best balanced accuracy per contrast.
static void WritePlotFiles(string prefix, List<GeneRow> genes, Dictionary<int, (int Class, string Product)> labels, string[] classNames,
    string fingerprint, int bins, int mapMinGenes)
{
    static double Hydropathy(GeneRow g) => -(g.Profile[0] + g.Profile[1]);
    static double Gc12(GeneRow g) => (1 + g.Profile[2]) / 2;
    var families = genes.GroupBy(g => g.Family).Select(f => (Family: f.Key, labels[f.Key].Class, labels[f.Key].Product, Genes: f.Count(),
        Hydropathy: f.Average(Hydropathy), Gc12: f.Average(Gc12), LogLength: f.Average(g => g.LogLength))).OrderBy(f => f.Family).ToList();
    var lines = new List<string> { "Family,Class,Product,Genes,Hydropathy (-(x + y)),GC12,Mean log length,Fingerprint" };
    lines.AddRange(families.Select(f => string.Join(",", f.Family, classNames[f.Class], Quote(f.Product), f.Genes, F(f.Hydropathy), F(f.Gc12),
        F(f.LogLength), fingerprint)));
    File.WriteAllLines($"{prefix}_Families.csv", lines);
    var map = new List<string> { $"# Corpus fingerprint {fingerprint}; families with at least {mapMinGenes} labelled genes",
        "Family,Class,Product,Genes,Hydropathy (-(x + y)),GC12" };
    map.AddRange(families.Where(f => f.Genes >= mapMinGenes).Select(f => string.Join(",", f.Family, classNames[f.Class], Quote(f.Product), f.Genes,
        F(f.Hydropathy), F(f.Gc12))));
    File.WriteAllLines($"{prefix}_FamilyMap.csv", map);

    double lo = genes.Min(Hydropathy), hi = genes.Max(Hydropathy), width = (hi - lo) / bins;
    var histogram = new List<string> { "Bin from,Bin to,Bin centre," + string.Join(",", classNames.Select(n => $"{n} share")) + "," +
        string.Join(",", classNames.Select(n => $"{n} genes")) + ",Fingerprint" };
    var totals = classNames.Select((_, c) => genes.Count(g => g.Class == c)).ToArray();
    // Each gene in exactly one bin, by index: comparing against float edges let a few genes count twice.
    var binCounts = new int[bins, classNames.Length];
    foreach (var g in genes) binCounts[Math.Min(bins - 1, (int)((Hydropathy(g) - lo) / width)), g.Class]++;
    for (int b = 0; b < bins; b++)
        histogram.Add(string.Join(",", F(lo + b * width), F(lo + (b + 1) * width), F(lo + (b + 0.5) * width),
            string.Join(",", classNames.Select((_, c) => (binCounts[b, c] / (double)totals[c]).ToString("F4", CultureInfo.InvariantCulture))),
            string.Join(",", classNames.Select((_, c) => binCounts[b, c])), fingerprint));
    histogram.Add("");
    foreach (int negative in new[] { 1, 2 })
    {
        var pair = genes.Where(g => g.Class == 0 || g.Class == negative).ToList();
        var familyPair = families.Where(f => f.Class == 0 || f.Class == negative).ToList();
        double geneAuc = SimpleAuc(pair.Select(Hydropathy).ToArray(), pair.Select(g => g.Membrane).ToArray());
        double familyAuc = SimpleAuc(familyPair.Select(f => f.Hydropathy).ToArray(), familyPair.Select(f => f.Class == 0).ToArray());
        var (threshold, accuracy) = BalancedThreshold(pair.Select(Hydropathy).ToArray(), pair.Select(g => g.Membrane).ToArray());
        histogram.Add($"Membrane vs {classNames[negative]},-(x + y) AUC over genes {F(geneAuc)}; over family means {F(familyAuc)}; " +
            $"best cut -(x + y) >= {F(threshold)} (balanced accuracy {F(accuracy)})");
        Console.WriteLine($"Membrane vs {classNames[negative]}: -(x + y) AUC over genes {F(geneAuc)}, over family centroids {F(familyAuc)}; " +
            $"balanced accuracy {F(accuracy)} at {F(threshold)}.");
    }
    File.WriteAllLines($"{prefix}_Histogram.csv", histogram);
}

static double PeakOf(double[] windows, int count)
{
    double peak = double.NegativeInfinity;
    for (int w = 0; w < count; w++) peak = Math.Max(peak, -(windows[w * ZCurveGeneProfile.FeatureCount] + windows[w * ZCurveGeneProfile.FeatureCount + 1]));
    return peak;
}

static string Quote(string text) => "\"" + text.Replace("\"", "\"\"") + "\"";

static double SimpleAuc(double[] scores, bool[] truth)
{
    var order = Enumerable.Range(0, scores.Length).OrderBy(i => scores[i]).ToArray();
    return Auc(scores, order, truth, Enumerable.Repeat(1.0, scores.Length).ToArray());
}

// The cut on the score with the best mean of sensitivity and specificity (membrane at or above it).
static (double Threshold, double Accuracy) BalancedThreshold(double[] scores, bool[] truth)
{
    var order = Enumerable.Range(0, scores.Length).OrderBy(i => scores[i]).ToArray();
    double positives = truth.Count(t => t), negatives = truth.Length - positives, negativesBelow = 0, positivesBelow = 0;
    var best = (Threshold: scores[order[0]], Accuracy: 0.5);
    for (int k = 0; k < order.Length; k++)
    {
        if (k > 0 && scores[order[k]] != scores[order[k - 1]])
        {
            double accuracy = (negativesBelow / negatives + (positives - positivesBelow) / positives) / 2;
            if (accuracy > best.Accuracy) best = (scores[order[k]], accuracy);
        }
        if (truth[order[k]]) positivesBelow++; else negativesBelow++;
    }
    return best;
}

static string Hash(string text) => Convert.ToHexStringLower(MD5.HashData(Encoding.UTF8.GetBytes(text)));

static int Half(string text) => MD5.HashData(Encoding.UTF8.GetBytes(text))[0] & 1;

static string F(double value) => value.ToString("F3", CultureInfo.InvariantCulture);

static (double Low, double High) Interval(double[] values)
{
    var sorted = values.Order().ToArray();
    return (sorted[(int)(0.025 * (sorted.Length - 1))], sorted[(int)Math.Ceiling(0.975 * (sorted.Length - 1))]);
}

// Weighted AUC (membrane scoring higher), ties counted half, in one pass over the scores sorted once: a resample only reweights genes.
static double Auc(double[] scores, int[] order, bool[] truth, double[] weight)
{
    double below = 0, sum = 0, positives = 0, negatives = 0;
    for (int start = 0; start < order.Length;)
    {
        int end = start;
        double pos = 0, neg = 0;
        while (end < order.Length && scores[order[end]] == scores[order[start]])
        {
            int i = order[end++];
            if (truth[i]) pos += weight[i]; else neg += weight[i];
        }
        sum += pos * (below + neg / 2);
        below += neg;
        positives += pos;
        negatives += neg;
        start = end;
    }
    return sum / (positives * negatives);
}

// Peak: the most hydrophobic window's -(P0 x + P0 y). Genome: the record, for bootstraps over genomes (live rows keep it in Phylum).
internal sealed record GeneRow(int Family, string Phylum, int Class, double[] Profile, double LogLength, double Gc, double Peak = double.NaN,
    string? Genome = null)
{
    public bool Membrane => Class == 0;
}

// Class-balanced logistic regression on standardized features, Newton steps with a small ridge.
internal sealed class Logistic
{
    private const double ridge = 1e-3;
    private const int steps = 25;

    private readonly double[] mean, sd, beta;

    private Logistic(double[] mean, double[] sd, double[] beta) => (this.mean, this.sd, this.beta) = (mean, sd, beta);

    public static Logistic Fit(List<GeneRow> rows, Func<GeneRow, double[]> features)
    {
        var x = rows.Select(features).ToArray();
        int d = x[0].Length;
        var mean = Enumerable.Range(0, d).Select(j => x.Average(v => v[j])).ToArray();
        var sd = Enumerable.Range(0, d).Select(j => Math.Sqrt(x.Average(v => (v[j] - mean[j]) * (v[j] - mean[j]))) is var s && s > 0 ? s : 1).ToArray();
        int positives = rows.Count(r => r.Membrane);
        double positiveWeight = rows.Count / (2.0 * positives), negativeWeight = rows.Count / (2.0 * (rows.Count - positives));
        var model = new Logistic(mean, sd, new double[d + 1]);
        for (int step = 0; step < steps; step++)
        {
            var gradient = new double[d + 1];
            var hessian = new double[d + 1, d + 1];
            for (int i = 0; i < x.Length; i++)
            {
                var z = model.Design(x[i]);
                double p = 1 / (1 + Math.Exp(-model.Linear(z))), w = rows[i].Membrane ? positiveWeight : negativeWeight;
                double residual = w * ((rows[i].Membrane ? 1 : 0) - p), curvature = w * p * (1 - p);
                for (int a = 0; a <= d; a++)
                {
                    gradient[a] += residual * z[a];
                    for (int b = 0; b <= d; b++) hessian[a, b] += curvature * z[a] * z[b];
                }
            }
            for (int a = 0; a <= d; a++) { gradient[a] -= ridge * model.beta[a]; hessian[a, a] += ridge; }
            var delta = Solve(hessian, gradient);
            for (int a = 0; a <= d; a++) model.beta[a] += delta[a];
            if (delta.Max(Math.Abs) < 1e-9) break;
        }
        return model;
    }

    public double[] Score(List<GeneRow> rows, Func<GeneRow, double[]> features) => rows.Select(r => Linear(Design(features(r)))).ToArray();

    // The same linear score on the raw features: intercept first, the standardization folded in.
    public double[] RawWeights() =>
        [beta[0] - Enumerable.Range(0, mean.Length).Sum(j => beta[j + 1] * mean[j] / sd[j]), .. Enumerable.Range(0, mean.Length).Select(j => beta[j + 1] / sd[j])];

    private double[] Design(double[] values) => [1, .. values.Select((v, j) => (v - mean[j]) / sd[j])];

    private double Linear(double[] z)
    {
        double sum = 0;
        for (int a = 0; a < z.Length; a++) sum += beta[a] * z[a];
        return sum;
    }

    // Gaussian elimination with partial pivoting; the system is small (at most 11 unknowns) and positive definite.
    private static double[] Solve(double[,] matrix, double[] vector)
    {
        int n = vector.Length;
        var a = (double[,])matrix.Clone();
        var b = (double[])vector.Clone();
        for (int col = 0; col < n; col++)
        {
            int pivot = Enumerable.Range(col, n - col).MaxBy(r => Math.Abs(a[r, col]));
            for (int k = 0; k < n; k++) (a[col, k], a[pivot, k]) = (a[pivot, k], a[col, k]);
            (b[col], b[pivot]) = (b[pivot], b[col]);
            for (int r = col + 1; r < n; r++)
            {
                double factor = a[r, col] / a[col, col];
                for (int k = col; k < n; k++) a[r, k] -= factor * a[col, k];
                b[r] -= factor * b[col];
            }
        }
        var solution = new double[n];
        for (int r = n - 1; r >= 0; r--)
        {
            double sum = b[r];
            for (int k = r + 1; k < n; k++) sum -= a[r, k] * solution[k];
            solution[r] = sum / a[r, r];
        }
        return solution;
    }
}
