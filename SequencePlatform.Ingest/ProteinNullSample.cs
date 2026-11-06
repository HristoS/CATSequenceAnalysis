using System.Globalization;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// The protein sample behind <see cref="ProteinScoreNull"/> (build-protein-null): real BV-BRC proteins with a protein family, from
    /// stored genomes taken in a seeded order, up to <c>perGenome</c> per length bucket from each genome so no genome dominates, until every
    /// bucket has <c>perBucket</c>. The proteins are held in memory for the build and never stored.
    /// </summary>
    internal static class ProteinNullSample
    {
        private const int perGenome = 2;
        private const int maxGenomes = 600;
        private static readonly string[] mycoplasmaOrders = { "Mycoplasmatales", "Entomoplasmatales", "Mycoplasmoidales" };

        public static async Task<List<ProteinNullQuery>> CollectAsync(NpgsqlDataSource dataSource, BvBrcClient client, int perBucket, int seed)
        {
            var genomes = new List<(string Genome, Dictionary<string, (long Id, string? Order)> Records)>();
            await using (var command = dataSource.CreateCommand(
                "select s.accession, r.record_id, r.id, r.tax_order from sequence_source s join sequence_record r on r.source_id = s.id " +
                "where s.provider = 'BV-BRC' and s.owner_id is null and exists (select 1 from gene_profile g where g.sequence_id = r.id) " +
                "order by md5(s.accession || $1), r.id"))
            {
                command.Parameters.AddWithValue(seed.ToString(CultureInfo.InvariantCulture));
                await using var reader = await command.ExecuteReaderAsync();
                while (await reader.ReadAsync())
                {
                    if (genomes.Count == 0 || genomes[^1].Genome != reader.GetString(0)) genomes.Add((reader.GetString(0), new()));
                    genomes[^1].Records[reader.GetString(1)] = (reader.GetInt64(2), reader.IsDBNull(3) ? null : reader.GetString(3));
                }
            }

            var filled = new int[GeneProfiles.MaxBucket + 1];
            var queries = new List<ProteinNullQuery>();
            var picker = new Random(seed);
            foreach (var (genome, records) in genomes.Take(maxGenomes))
            {
                if (filled.All(f => f >= perBucket)) break;
                var proteins = (await client.ListProteinsAsync(genome))
                    .Where(p => p.LengthAa >= SequenceAnalyses.ZCurveGeneProfile.MinResidues && records.ContainsKey(p.Accession))
                    .GroupBy(p => GeneProfiles.BucketOf(p.LengthAa))
                    .Where(b => filled[b.Key] < perBucket)
                    .SelectMany(b => b.OrderBy(_ => picker.Next()).Take(perGenome)).ToList();
                if (proteins.Count == 0) continue;
                var sequences = await client.ProteinSequencesAsync(proteins.Select(p => p.Md5).ToList());
                foreach (var p in proteins)
                {
                    int bucket = GeneProfiles.BucketOf(p.LengthAa), family = ProteinFamilies.Parse(p.FamilyId);
                    if (filled[bucket] >= perBucket || family == SequenceAnalyses.ZCurveGeneProfile.NoFamily) continue;
                    if (!sequences.TryGetValue(p.Md5, out var protein)) continue;
                    var (id, order) = records[p.Accession];
                    int code = order is not null && mycoplasmaOrders.Contains(order) ? ProteinQuery.MycoplasmaCode : ProteinQuery.StandardCode;
                    ProteinQuery query;
                    try { query = ProteinQuery.FromProtein(protein.TrimEnd('*'), code); }
                    catch (ArgumentException) { continue; }   // an internal stop, letters that aren't amino acids, or too short once trimmed
                    if (GeneProfiles.BucketOf(query.LengthAa) != bucket) continue;
                    queries.Add(new ProteinNullQuery(query, family, id, p.Start, p.End, p.Reverse));
                    filled[bucket]++;
                }
                Console.WriteLine($"  {genome}: {queries.Count:N0} proteins, buckets filled {filled.Count(f => f >= perBucket)} of {filled.Length}");
            }
            return queries;
        }
    }
}
