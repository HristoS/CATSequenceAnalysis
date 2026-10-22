#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using Npgsql;
using SequencePlatform.Core;

// Why are all 20 neighbours of the Vibrio fischeri ES114 plasmid (record 7787) coronaviruses? Mean Z-curve axes over
// each record's 500 bp windows: x = purine - pyrimidine, y = amino - keto, z = strong - weak (GC - AT). Groups: the
// records of 15-60 kbp with GC 0.37-0.42, the slice SARS-CoV-2 sits in. For the plasmid both strands are shown,
// because the search takes the better one (the other strand negates x and y).
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());
var groups = new Dictionary<string, List<long>>();
await using (var cmd = dataSource.CreateCommand(
    "select r.id, case when r.id = 7787 then 'the plasmid' when r.id = 19634 then 'SARS-CoV-2' " +
    "when coalesce(t.tax_family, r.tax_family) = 'Coronaviridae' then 'Coronaviridae' " +
    "when r.kind = 'virus' then 'other viruses' else 'other ' || r.kind || 's' end " +
    "from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid " +
    "where (r.length_bp between 14952 and 59806 and r.gc_fraction between 0.37 and 0.42) or r.id in (7787, 19634)"))
await using (var r = await cmd.ExecuteReaderAsync())
    while (await r.ReadAsync())
    {
        var key = r.GetString(1);
        if (!groups.TryGetValue(key, out var list)) groups[key] = list = new List<long>();
        list.Add(r.GetInt64(0));
    }

Console.WriteLine("Group, Records, MeanX, MeanY, MeanZ, SdX, SdY");
foreach (var (group, ids) in groups.OrderBy(g => g.Key))
{
    var means = new List<(double X, double Y, double Z)>();
    foreach (var id in ids)
    {
        var fine = await ProfileLevels.ReadAsync(dataSource, id, RecordStore.Method, RecordStore.MethodVersion, RecordStore.FineWindowSize);
        if (fine is null) continue;
        var w = fine.Windows;
        means.Add((w.Average(v => v.X), w.Average(v => v.Y), w.Average(v => v.Z)));
    }
    double Sd(Func<(double X, double Y, double Z), double> f) =>
        means.Count < 2 ? 0 : Math.Sqrt(means.Sum(m => Math.Pow(f(m) - means.Average(f), 2)) / (means.Count - 1));
    Console.WriteLine($"{group}, {means.Count}, {means.Average(m => m.X):F3}, {means.Average(m => m.Y):F3}, {means.Average(m => m.Z):F3}, {Sd(m => m.X):F3}, {Sd(m => m.Y):F3}");
    if (group == "the plasmid")
        Console.WriteLine($"the plasmid (other strand), 1, {-means[0].X:F3}, {-means[0].Y:F3}, {means[0].Z:F3}, , ");
}
