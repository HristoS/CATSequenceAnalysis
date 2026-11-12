#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
#:include shared/EukaryoteFiles.cs
#:property PublishAot=false
using Npgsql;
using SequencePlatform.Core;

// Downloads the eukaryote pairs the cheap checks use (EukaryoteFiles.Pairs) into data/genomes/eukaryotes/, once, each checked against
// the SHA-256 the platform stored. Run from the repo root with SEQUENCE_PLATFORM_DB set; files already there are skipped.
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
foreach (var (id, label) in EukaryoteFiles.Pairs.SelectMany(p => p))
{
    var genome = await EukaryoteFiles.GetAsync(dataSource, id, label, Console.WriteLine);
    long bytes = new[] { genome.Fasta, genome.Cds, genome.Report }.Where(f => f is not null).Sum(f => new FileInfo(f!).Length);
    Console.WriteLine($"{label,-28} {bytes / 1e6,8:F1} MB  {Path.GetFileName(genome.Fasta)}");
}
