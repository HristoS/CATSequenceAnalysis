using System.Collections.Concurrent;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>One protein of the null sample: its expected windows, its BV-BRC family, and its own stored gene, which is left out.</summary>
    public sealed record ProteinNullQuery(ProteinQuery Query, int Family, long SequenceId, int GeneStart, int GeneEnd, bool Reverse);

    /// <summary>
    /// p and E for protein hits (P7, migration 0019): per query length bucket (<see cref="GeneProfiles.BucketOf"/>), how often genes of
    /// other BV-BRC families score at least s against real proteins, each compared as <see cref="ProteinSearch"/> compares it. p =
    /// (k + 1) / (N + 1) over those protein-gene pairs, E = p times the genes a search scored. Genes with no family are left out of the
    /// sample. The check before building (results/ProteinNull_Check.csv, 411 proteins): the median protein's best same-family hit is
    /// beyond every other-family score of the rest of its band in every band but 333-499 aa (E 0.025), and the spread between proteins at
    /// a fixed score is 20-41x. Past the sample's top score, p is a bound, as in <see cref="ScoreNullCell"/>.
    /// </summary>
    public sealed class ProteinScoreNull
    {
        public const string MethodVersion = "1";
        public const int Bins = 2000;   // score grid of 1 / Bins from 0 to 1
        private const int minEmpiricalExceedances = 10;

        private readonly Dictionary<int, Cell> cells;

        /// <param name="Exceed">Exceed[i]: pairs scoring at least i / <see cref="Bins"/>.</param>
        public sealed record Cell(int Bucket, int Queries, long Pairs, long[] Exceed);

        private ProteinScoreNull(Dictionary<int, Cell> cells) => this.cells = cells;

        public IReadOnlyCollection<Cell> Cells => cells.Values;

        /// <summary>The stored null; null when none is built.</summary>
        public static async Task<ProteinScoreNull?> LoadAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand("select bucket, queries, pairs, exceed from protein_score_null where method_version = $1");
            command.Parameters.AddWithValue(MethodVersion);
            var cells = new Dictionary<int, Cell>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                cells[reader.GetInt16(0)] = new Cell(reader.GetInt16(0), reader.GetInt32(1), reader.GetInt64(2), reader.GetFieldValue<long[]>(3));
            return cells.Count == 0 ? null : new ProteinScoreNull(cells);
        }

        /// <summary>p and log10 E of a hit scoring <paramref name="score"/> for a protein of this length; null when its bucket has no cell.</summary>
        public (PValue P, double Log10E)? Assess(int lengthAa, double score, long genesScored)
        {
            if (!cells.TryGetValue(GeneProfiles.BucketOf(lengthAa), out var cell)) return null;
            long k = cell.Exceed[Bin(score)];
            double log10P = Math.Log10((k + 1.0) / (cell.Pairs + 1.0));
            var kind = k >= minEmpiricalExceedances ? PValueKind.Empirical : PValueKind.Bound;
            return (new PValue(log10P, kind), log10P + Math.Log10(Math.Max(1, genesScored)));
        }

        // The grid cell a score falls in; a score in a cell counts as reaching that cell's lower edge, which leans conservative.
        private static int Bin(double score) => Math.Clamp((int)Math.Floor(score * Bins), 0, Bins);

        /// <summary>
        /// Scores every query against the genes a search would compare it with, in one pass over gene_profile, and replaces the stored
        /// cells: one per bucket that has queries.
        /// </summary>
        public static async Task<IReadOnlyList<Cell>> BuildAsync(NpgsqlDataSource dataSource, IReadOnlyList<ProteinNullQuery> queries,
            Action<string>? log = null, CancellationToken cancellationToken = default)
        {
            var rules = queries.Select(q => ProteinSearch.CandidateRule(q.Query)).ToArray();
            var buckets = queries.Select(q => GeneProfiles.BucketsFor(q.Query.LengthAa, ProteinSearch.LengthBand).ToHashSet()).ToArray();
            var rows = new BlockingCollection<(long Id, int Bucket, byte[] Payload)>(256);
            var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(() =>
            {
                var counts = queries.Select(_ => new long[Bins + 1]).ToArray();
                foreach (var (id, bucket, payload) in rows.GetConsumingEnumerable(cancellationToken))
                {
                    var profile = ZCurveGeneProfile.FromBinary(payload);
                    for (int q = 0; q < queries.Count; q++)
                    {
                        if (!buckets[q].Contains(bucket)) continue;
                        var query = queries[q];
                        for (int g = 0; g < profile.Genes.Count; g++)
                        {
                            var gene = profile.Genes[g];
                            if (gene.Family == ZCurveGeneProfile.NoFamily || gene.Family == query.Family || !rules[q](gene)) continue;
                            if (id == query.SequenceId && gene.Start == query.GeneStart && gene.End == query.GeneEnd && gene.Reverse == query.Reverse) continue;
                            counts[q][Bin(profile.Score(query.Query.Axes, g))]++;
                        }
                    }
                }
                return counts;
            }, cancellationToken)).ToArray();

            var needed = buckets.SelectMany(b => b).Distinct().Select(b => (short)b).ToArray();
            long read = 0;
            try
            {
                await using var command = dataSource.CreateCommand("select sequence_id, bucket, payload from gene_profile where bucket = any($1) and window_count > 0");
                command.Parameters.AddWithValue(needed);
                await using var reader = await command.ExecuteReaderAsync(cancellationToken);
                while (await reader.ReadAsync(cancellationToken))
                {
                    rows.Add((reader.GetInt64(0), reader.GetInt16(1), reader.GetFieldValue<byte[]>(2)), cancellationToken);
                    if (++read % 100_000 == 0) log?.Invoke($"  {read:N0} gene rows read");
                }
            }
            finally
            {
                rows.CompleteAdding();
            }
            var perWorker = await Task.WhenAll(workers);

            var cells = new List<Cell>();
            foreach (var group in Enumerable.Range(0, queries.Count).GroupBy(q => GeneProfiles.BucketOf(queries[q].Query.LengthAa)).OrderBy(g => g.Key))
            {
                var exceed = new long[Bins + 1];
                for (int b = Bins; b >= 0; b--)
                    exceed[b] = (b < Bins ? exceed[b + 1] : 0) + group.Sum(q => perWorker.Sum(w => w[q][b]));
                if (exceed[0] > 0) cells.Add(new Cell(group.Key, group.Count(), exceed[0], exceed));
            }

            await using var connection = await dataSource.OpenConnectionAsync(cancellationToken);
            await using var transaction = await connection.BeginTransactionAsync(cancellationToken);
            await using (var delete = new NpgsqlCommand("delete from protein_score_null where method_version = $1", connection, transaction))
            {
                delete.Parameters.AddWithValue(MethodVersion);
                await delete.ExecuteNonQueryAsync(cancellationToken);
            }
            foreach (var cell in cells)
            {
                await using var insert = new NpgsqlCommand(
                    "insert into protein_score_null (method_version, bucket, queries, pairs, exceed) values ($1, $2, $3, $4, $5)", connection, transaction);
                insert.Parameters.AddWithValue(MethodVersion);
                insert.Parameters.AddWithValue((short)cell.Bucket);
                insert.Parameters.AddWithValue(cell.Queries);
                insert.Parameters.AddWithValue(cell.Pairs);
                insert.Parameters.AddWithValue(cell.Exceed);
                await insert.ExecuteNonQueryAsync(cancellationToken);
            }
            await transaction.CommitAsync(cancellationToken);
            return cells;
        }
    }
}
