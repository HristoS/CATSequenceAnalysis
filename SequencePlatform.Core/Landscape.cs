using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <param name="Domain">Bacteria, Archaea, Viruses or Eukaryota, from the record's kind.</param>
    public sealed record LandscapePoint(long Id, string Name, long LengthBp, double Gc, string Domain, string? Phylum, string? Family, string? Genus, double X, double Y);

    public sealed record LandscapeLoading(string Dinucleotide, double Weight);

    /// <param name="ExplainedVariance">The axis's share of the signatures' total variance.</param>
    /// <param name="TopLoadings">The dinucleotides that weigh most on the axis, heaviest first.</param>
    public sealed record LandscapeAxis(double ExplainedVariance, IReadOnlyList<LandscapeLoading> TopLoadings);

    public sealed record LandscapeMap(string Group, IReadOnlyList<LandscapePoint> Points, IReadOnlyList<LandscapeAxis> Axes);

    /// <summary>
    /// A two-dimensional map of the public catalogue: PCA of the dinucleotide signature (<see cref="SignatureIndex"/>'s log rho*). The
    /// cellular and viral maps hold one genome per species: cellular genomes need a GTDB assignment, whose lineage labels them, and viruses of
    /// at least 1 kbp carry NCBI's; selection follows scripts/platform-checks/r4_map.cs, so the same seed picks the same genomes. The
    /// <see cref="All"/> map holds every public genome <see cref="SignatureIndex"/> would, on one set of axes, labelled by GTDB where
    /// assigned and NCBI otherwise. Two axes keep phyla apart but not genera (R4, results/S53_Map_*.csv): this is a picture of the
    /// catalogue, and nearest relatives come from the full signature.
    ///
    /// Private uploads never enter it, so nobody's upload moves anyone's map. Each group is built on first request and kept until the
    /// record count or the highest record id changes.
    /// </summary>
    public sealed class Landscape
    {
        public const string All = "all";

        private const int minViralLengthBp = 1000;
        private const int speciesSeed = 20261030;
        private const int topLoadings = 4;

        private readonly NpgsqlDataSource dataSource;
        private readonly SemaphoreSlim building = new(1, 1);
        private readonly Dictionary<string, (long Count, long MaxId, LandscapeMap Map)> maps = new();

        public Landscape(NpgsqlDataSource dataSource)
        {
            this.dataSource = dataSource;
        }

        public async Task<LandscapeMap> GetAsync(string group, CancellationToken cancellationToken = default)
        {
            var sql = group switch
            {
                All =>
                    "select r.id, r.name, r.length_bp, r.dinucleotide_counts, coalesce(g.tax_phylum, t.tax_phylum), coalesce(g.tax_family, t.tax_family), " +
                    "coalesce(g.tax_genus, t.tax_genus), coalesce(g.tax_species, t.tax_species), r.kind " +
                    "from (select r.*, row_number() over (partition by r.source_id, r.kind in ('bacterium', 'archaeon') order by r.length_bp desc, r.id) as n " +
                    " from sequence_record r join sequence_source s on s.id = r.source_id " +
                    " where r.dinucleotide_counts is not null and r.kind in ('bacterium', 'archaeon', 'virus', 'eukaryote') and s.owner_id is null) r " +
                    "left join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
                    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
                    "left join taxon t on t.taxid = r.ncbi_taxid " +
                    "where r.kind = 'virus' or r.n = 1 " +
                    "order by r.id",
                SignatureIndex.Cellular =>
                    "select distinct on (r.source_id) r.id, r.name, r.length_bp, r.dinucleotide_counts, g.tax_phylum, g.tax_family, g.tax_genus, g.tax_species, r.kind " +
                    "from sequence_record r join sequence_source s on s.id = r.source_id " +
                    "join gtdb_taxon g on r.annotations->'bvbrc'->>'assembly_accession' ~ '^GC[AF]_' " +
                    "and g.assembly_key = substr(r.annotations->'bvbrc'->>'assembly_accession', 5) " +
                    "where r.kind in ('bacterium', 'archaeon') and r.dinucleotide_counts is not null and s.owner_id is null " +
                    "order by r.source_id, r.length_bp desc, r.id",
                SignatureIndex.Viral =>
                    "select r.id, r.name, r.length_bp, r.dinucleotide_counts, t.tax_phylum, t.tax_family, t.tax_genus, t.tax_species, r.kind " +
                    "from sequence_record r join sequence_source s on s.id = r.source_id " +
                    "left join taxon t on t.taxid = r.ncbi_taxid " +
                    $"where r.kind = 'virus' and r.dinucleotide_counts is not null and r.length_bp >= {minViralLengthBp} and s.owner_id is null " +
                    "order by r.id",
                SignatureIndex.Eukaryotic =>
                    "select r.id, r.name, r.length_bp, r.dinucleotide_counts, coalesce(t.tax_phylum, r.tax_phylum), coalesce(t.tax_family, r.tax_family), " +
                    "coalesce(t.tax_genus, r.tax_genus), coalesce(t.tax_species, r.tax_species), r.kind " +
                    "from sequence_record r join sequence_source s on s.id = r.source_id " +
                    "left join taxon t on t.taxid = r.ncbi_taxid " +
                    "where r.kind = 'eukaryote' and r.dinucleotide_counts is not null and s.owner_id is null " +
                    "order by r.id",
                _ => throw new ArgumentException($"Unknown group '{group}'. Use {All}, {SignatureIndex.Cellular}, {SignatureIndex.Viral} or {SignatureIndex.Eukaryotic}."),
            };

            await building.WaitAsync(cancellationToken);
            try
            {
                var (count, maxId) = await SignatureIndex.StateAsync(dataSource, cancellationToken);
                if (maps.TryGetValue(group, out var known) && known.Count == count && known.MaxId == maxId) return known.Map;

                var map = Build(group, await LoadAsync(sql, onePerSpecies: group != All, cancellationToken));
                maps[group] = (count, maxId, map);
                return map;
            }
            finally
            {
                building.Release();
            }
        }

        private async Task<List<Genome>> LoadAsync(string sql, bool onePerSpecies, CancellationToken cancellationToken)
        {
            var rows = new List<Genome>();
            await using (var command = dataSource.CreateCommand(sql))
            await using (var reader = await command.ExecuteReaderAsync(cancellationToken))
                while (await reader.ReadAsync(cancellationToken))
                {
                    string? Text(int column) => reader.IsDBNull(column) || reader.GetString(column) == "" ? null : reader.GetString(column);
                    var counts = reader.GetFieldValue<long[]>(3);
                    if (counts.Sum() == 0) continue;
                    double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
                    rows.Add(new Genome(reader.GetInt64(0), Text(1) ?? "", reader.GetInt64(2), gc, Domain(reader.GetString(8)), counts, Text(4), Text(5), Text(6), Text(7)));
                }
            if (!onePerSpecies) return rows;

            // One genome per species; a genome without a species stands for itself.
            var rng = new Random(speciesSeed);
            return rows.GroupBy(g => g.Species ?? $"id {g.Id}").Select(s => s.OrderBy(_ => rng.Next()).First()).OrderBy(g => g.Id).ToList();
        }

        private static string Domain(string kind) => kind switch
        {
            "bacterium" => "Bacteria",
            "archaeon" => "Archaea",
            "virus" => "Viruses",
            _ => "Eukaryota",
        };

        private static LandscapeMap Build(string group, List<Genome> genomes)
        {
            if (genomes.Count == 0) return new LandscapeMap(group, [], []);

            var pca = Ordination.Pca(genomes.Select(g => SignatureIndex.Signature(g.Dinucleotides)).ToArray(), 2);
            var points = genomes.Select((g, i) => new LandscapePoint(g.Id, g.Name, g.LengthBp, g.Gc, g.Domain, g.Phylum, g.Family, g.Genus, pca.Scores[i][0], pca.Scores[i][1])).ToList();
            var axes = pca.Loadings.Select((loading, c) => new LandscapeAxis(pca.ExplainedVariance[c], loading
                .Select((weight, i) => new LandscapeLoading($"{"ACGT"[i / 4]}{"ACGT"[i % 4]}", weight))
                .OrderByDescending(l => Math.Abs(l.Weight)).Take(topLoadings).ToList())).ToList();
            return new LandscapeMap(group, points, axes);
        }

        private sealed record Genome(long Id, string Name, long LengthBp, double Gc, string Domain, long[] Dinucleotides, string? Phylum, string? Family, string? Genus, string? Species);
    }
}
