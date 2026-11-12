#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/ExperimentCorpus.cs
#:property PublishAot=false
using System.Globalization;
using Npgsql;

// P9 corpus: what the shared corpus holds, per group. Every P9 script starts from ExperimentCorpus.LoadAsync, so these counts are the
// denominators behind their results. Default is round 2 (S9.9): the live database cut off at RoundTwoMaxId, its fingerprint verified,
// and the coverage of the inputs round 2 adds (6-mers, amino-acid 4-mers, gene profiles). "--round1" repeats S9.0 on the snapshot.
//   dotnet run p9_corpus.cs [--round1] [output.csv]
//   dotnet run p9_corpus.cs --fingerprint      print the per-table hashes of the live database at the cut-off, check nothing
var roundOne = args.Contains("--round1");
if (args.Contains("--fingerprint"))
{
    await using var live = ExperimentCorpus.OpenLive();
    var tables = await ExperimentCorpus.FingerprintAsync(live, ExperimentCorpus.RoundTwoMaxId);
    foreach (var (table, hash) in tables) Console.WriteLine($"(\"{table}\", \"{hash}\"),");
    Console.WriteLine($"Combined: {ExperimentCorpus.Combine(tables)}");
    return;
}
var output = args.FirstOrDefault(a => !a.StartsWith("--", StringComparison.Ordinal)) ?? (roundOne ? "P9_S90_Corpus.csv" : "P9_S99_Corpus.csv");

await using var dataSource = roundOne ? ExperimentCorpus.OpenSnapshot() : ExperimentCorpus.OpenLive();
var verified = roundOne ? Task.FromResult("") : ExperimentCorpus.VerifyRoundTwoAsync(dataSource);
var genomes = await ExperimentCorpus.LoadAsync(dataSource, roundOne ? long.MaxValue : ExperimentCorpus.RoundTwoMaxId);
(HashSet<long> Kmers6, HashSet<long> AminoKmers, HashSet<long> GeneProfiles)? inputs = roundOne ? null : await RoundTwoInputsAsync(dataSource);
var fingerprint = await verified;

using var csv = File.CreateText(output);
csv.WriteLine("Group, Genomes, Excluded, WithNcbiLineage, WithGtdbLineage, WithProviderGenus, WithHostName, WithCheckm, WithIsolationSource, MedianLengthBp" +
    (roundOne ? "" : ", With6mers, WithAminoKmers, WithGeneProfiles, Fingerprint"));
foreach (var group in genomes.GroupBy(g => g.Group).OrderBy(g => g.Key, StringComparer.Ordinal))
{
    var lengths = group.Select(g => g.LengthBp).Order().ToArray();
    var line = string.Create(CultureInfo.InvariantCulture,
        $"{group.Key}, {group.Count()}, {group.Count(g => g.Excluded)}, {group.Count(g => g.Ncbi is not null)}, {group.Count(g => g.Gtdb is not null)}, " +
        $"{group.Count(g => g.Provider.Genus is not null)}, {group.Count(g => g.Field("host_name") is not null)}, " +
        $"{group.Count(g => g.Field("checkm_contamination") is not null)}, {group.Count(g => g.Field("isolation_source") is not null)}, {lengths[lengths.Length / 2]}");
    if (inputs is not null)
        line += $", {group.Count(g => inputs.Value.Kmers6.Contains(g.Id))}, {group.Count(g => inputs.Value.AminoKmers.Contains(g.Id))}, " +
            $"{group.Count(g => inputs.Value.GeneProfiles.Contains(g.Id))}, {fingerprint}";
    csv.WriteLine(line);
    Console.WriteLine(line);
}
Console.WriteLine(roundOne
    ? $"{genomes.Count} genomes, {genomes.Count(g => g.Excluded)} excluded, from {ExperimentCorpus.SnapshotDatabase}."
    : $"{genomes.Count} genomes, {genomes.Count(g => g.Excluded)} excluded, records <= {ExperimentCorpus.RoundTwoMaxId:N0}, fingerprint {fingerprint}.");

static async Task<(HashSet<long> Kmers6, HashSet<long> AminoKmers, HashSet<long> GeneProfiles)> RoundTwoInputsAsync(NpgsqlDataSource dataSource)
{
    async Task<HashSet<long>> IdsAsync(string sql)
    {
        var ids = new HashSet<long>();
        await using var command = dataSource.CreateCommand(sql);
        command.Parameters.AddWithValue(ExperimentCorpus.RoundTwoMaxId);
        await using var reader = await command.ExecuteReaderAsync();
        while (await reader.ReadAsync()) ids.Add(reader.GetInt64(0));
        return ids;
    }
    return (await IdsAsync("select sequence_id from record_kmer where status = 'done' and sequence_id <= $1"),
        await IdsAsync("select sequence_id from protein_kmer where status = 'done' and sequence_id <= $1"),
        await IdsAsync("select distinct sequence_id from gene_profile where sequence_id <= $1"));
}
