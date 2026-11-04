#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;

// P9 S9.0: what the shared corpus holds, per group, read from the snapshot database. Every later P9 script starts from the same
// ExperimentCorpus.LoadAsync, so these counts are the denominators behind their results. Writes the given CSV.
var output = args.Length > 0 ? args[0] : "P9_S90_Corpus.csv";
await using var dataSource = ExperimentCorpus.OpenSnapshot();
var genomes = await ExperimentCorpus.LoadAsync(dataSource);

using var csv = File.CreateText(output);
csv.WriteLine("Group, Genomes, Excluded, WithNcbiLineage, WithGtdbLineage, WithProviderGenus, WithHostName, WithCheckm, WithIsolationSource, MedianLengthBp");
foreach (var group in genomes.GroupBy(g => g.Group).OrderBy(g => g.Key, StringComparer.Ordinal))
{
    var lengths = group.Select(g => g.LengthBp).Order().ToArray();
    var line = string.Create(CultureInfo.InvariantCulture,
        $"{group.Key}, {group.Count()}, {group.Count(g => g.Excluded)}, {group.Count(g => g.Ncbi is not null)}, {group.Count(g => g.Gtdb is not null)}, " +
        $"{group.Count(g => g.Provider.Genus is not null)}, {group.Count(g => g.Field("host_name") is not null)}, " +
        $"{group.Count(g => g.Field("checkm_contamination") is not null)}, {group.Count(g => g.Field("isolation_source") is not null)}, {lengths[lengths.Length / 2]}");
    csv.WriteLine(line);
    Console.WriteLine(line);
}
Console.WriteLine($"{genomes.Count} genomes, {genomes.Count(g => g.Excluded)} excluded, from {ExperimentCorpus.SnapshotDatabase}.");
