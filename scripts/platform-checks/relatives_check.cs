#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S5.5's stop criteria on the local catalogue: E. coli K-12 (bundled, not in the catalogue; its counts are computed from the file
// here and the DNA is not kept) must have Escherichia or Shigella among its 10 relatives; a warm relatives query must take under
// 100 ms; a 500-genome tree under 5 s. Also SARS-CoV-2's viral relatives, index load time and size. Pass a folder to also write
// the 500-genome trees as upgma.nwk and nj.nwk there (to compare implementations).
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var index = new SignatureIndex(dataSource);
var watch = Stopwatch.StartNew();
await index.RefreshAsync();
Console.WriteLine($"index: {index.Entries.Count:N0} genomes ({index.Entries.Values.Count(e => e.Group == SignatureIndex.Cellular):N0} cellular, " +
    $"{index.Entries.Values.Count(e => e.Group == SignatureIndex.Viral):N0} viral), loaded in {watch.Elapsed.TotalMilliseconds:F0} ms");

foreach (var (file, group) in new[] { ("ecoli_k12_mg1655_GCF_000005845.2.fna.gz", SignatureIndex.Cellular), ("viral/sars_cov_2_NC_045512.2.fna.gz", SignatureIndex.Viral) })
{
    var counts = GenomeSignatures.DinucleotideCounts(FastaReader.ReadPrimaryRecord(Path.Combine("data/genomes", file)).Sequence);
    var signature = SignatureIndex.Signature(counts);
    index.Nearest(signature, group, null, 10);   // warm-up
    var times = new List<double>();
    IReadOnlyList<Relative> relatives = [];
    for (int i = 0; i < 20; i++)
    {
        watch.Restart();
        relatives = index.Nearest(signature, group, null, 10);
        times.Add(watch.Elapsed.TotalMilliseconds);
    }
    var records = await Catalogue.ByIdsAsync(dataSource, relatives.Select(r => r.Id).ToList());
    Console.WriteLine($"{file}: median query {times.Order().ElementAt(times.Count / 2):F2} ms, max {times.Max():F2} ms");
    foreach (var r in relatives)
        Console.WriteLine($"  {r.Distance:F4}  {records[r.Id].Lineage.Genus,-22} {records[r.Id].Name}");
}

foreach (var method in new[] { "upgma", "nj" })
{
    var leaves = index.Entries.Values.Where(e => e.Group == SignatureIndex.Cellular && e.OwnerId is null).Take(SignatureIndex.MaxTreeLeaves)
        .Select(e => (e.Id.ToString(), e.Signature)).ToList();
    watch.Restart();
    var tree = SignatureIndex.Tree(leaves, method);
    Console.WriteLine($"{method} tree of {leaves.Count} genomes: {watch.Elapsed.TotalMilliseconds:F0} ms, Newick {tree.ToNewick().Length:N0} chars");
    if (args.Length > 0) File.WriteAllText(Path.Combine(args[0], $"{method}.nwk"), tree.ToNewick());   // to compare implementations
}
