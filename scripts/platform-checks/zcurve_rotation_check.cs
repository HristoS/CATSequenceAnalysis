#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Globalization;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// P5 (2026-10-30): the map of the complete second-order Z-curve (all 16 Z-curve coefficients of rho*, GenomeSignatures.ToZCurveBasis)
// against the Karlin map, on r4_map.cs's corpora. A rotation of the vectors can only rotate the principal components, so the two maps'
// scores must correlate at |r| = 1, PC by PC (signs are arbitrary). Also the log-scale form, which is not a rotation of rho*.
// Prints the correlations; no output file.
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var cellular = await LoadAsync(
    "select distinct on (r.source_id) r.id, r.name, r.length_bp, r.dinucleotide_counts, g.tax_domain, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species " +
    "from sequence_record r join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null order by r.source_id, r.length_bp desc, r.id");
var viruses = await LoadAsync(
    "select r.id, r.name, r.length_bp, r.dinucleotide_counts, t.tax_realm, t.tax_phylum, t.tax_family, t.tax_genus, t.tax_species " +
    "from sequence_record r left join taxon t on t.taxid = r.ncbi_taxid " +
    "where r.kind = 'virus' and r.dinucleotide_counts is not null and r.length_bp >= 1000 order by r.id");

foreach (var (corpus, genomes) in new[] { ("cellular", cellular), ("viruses", viruses) })
{
    var karlin = Ordination.Pca(genomes.Select(g => GenomeSignatures.KarlinRho(g.Dinucleotides)).ToArray(), 2);
    var full = Ordination.Pca(genomes.Select(g => GenomeSignatures.ToZCurveBasis(GenomeSignatures.KarlinRho(g.Dinucleotides))).ToArray(), 2);
    var log = Ordination.Pca(genomes.Select(g => GenomeSignatures.ToZCurveBasis(GenomeSignatures.KarlinLogRho(g.Dinucleotides))).ToArray(), 2);
    for (int c = 0; c < 2; c++)
        Console.WriteLine(string.Create(CultureInfo.InvariantCulture,
            $"{corpus} PC{c + 1}: complete Z-curve (ratio) vs Karlin r = {Correlation(full.Scores, karlin.Scores, c):F6}; " +
            $"explained {full.ExplainedVariance[c]:P2} vs {karlin.ExplainedVariance[c]:P2}; complete Z-curve (log) vs Karlin r = {Correlation(log.Scores, karlin.Scores, c):F4}"));
}

static double Correlation(double[][] a, double[][] b, int c)
{
    double ma = a.Average(r => r[c]), mb = b.Average(r => r[c]);
    double sab = a.Zip(b, (x, y) => (x[c] - ma) * (y[c] - mb)).Sum();
    return sab / Math.Sqrt(a.Sum(x => (x[c] - ma) * (x[c] - ma)) * b.Sum(y => (y[c] - mb) * (y[c] - mb)));
}

async Task<List<Genome>> LoadAsync(string sql)
{
    var rows = new List<Genome>();
    await using (var cmd = dataSource.CreateCommand(sql))
    await using (var r = await cmd.ExecuteReaderAsync())
        while (await r.ReadAsync())
        {
            string? Text(int column) => r.IsDBNull(column) || r.GetString(column) == "" ? null : r.GetString(column);
            var counts = r.GetFieldValue<long[]>(3);
            if (counts.Sum() == 0) continue;
            double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
            rows.Add(new Genome(r.GetInt64(0), Text(1) ?? "", r.GetInt64(2), gc, counts, Text(4), Text(5), Text(6), Text(7), Text(8)));
        }
    // One genome per species; a genome without a species stands for itself.
    var rng = new Random(20261030);
    return rows.GroupBy(g => g.Species ?? $"id {g.Id}").Select(s => s.OrderBy(_ => rng.Next()).First()).OrderBy(g => g.Id).ToList();
}

record Genome(long Id, string Name, long Length, double Gc, long[] Dinucleotides, string? Top, string? Phylum, string? Family, string? Genus, string? Species);
