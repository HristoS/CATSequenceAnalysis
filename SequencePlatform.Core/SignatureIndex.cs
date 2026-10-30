using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <param name="Group">Relatives are found within a group: <see cref="SignatureIndex.Cellular"/> or <see cref="SignatureIndex.Viral"/>.</param>
    public sealed record SignatureEntry(long Id, string Group, long? OwnerId, double[] Signature);

    public sealed record Relative(long Id, double Distance);

    /// <summary>One value of a signature view.</summary>
    /// <param name="Ratio">The Karlin view's rho* itself (exp of the log value); null in the Z-curve view.</param>
    public sealed record SignatureValue(string Name, double Value, string Meaning, double? Ratio);

    /// <summary>
    /// The dinucleotide signature of every genome, in memory, for nearest relatives and trees (P5 S5.5). The signature is log rho*
    /// over both strands (<see cref="GenomeSignatures.KarlinLogRho"/>), compared by Euclidean distance: Karlin's genomic signature, which
    /// is also the complete second-order Z-curve, since its Z-curve coefficients (<see cref="GenomeSignatures.ToZCurveBasis"/>) are a
    /// rotation of it; <see cref="Views"/> serves either. Log rather than ratio because it treats depletion and enrichment alike; the
    /// four scale and distance choices measured within about 0.01 of each other (results/S55_*.csv). It comes from the stored
    /// dinucleotide counts (D8), so a profile-only upload has none.
    ///
    /// One record stands for each genome (D9): the largest chromosome of a catalogue bacterium or archaeon, every virus record, and
    /// every record of an upload (an upload may hold many genomes); plasmids are left out. <see cref="RefreshAsync"/> reloads when the
    /// record count or the highest record id changed, which an insert or a deletion always does. Each refresh publishes a new
    /// snapshot, so a query keeps a consistent view while one runs.
    /// </summary>
    public sealed class SignatureIndex
    {
        public const string Cellular = "cellular";
        public const string Viral = "viral";
        public const int MaxTreeLeaves = 500;
        public const string Description = "log rho* over both strands (dinucleotides), Euclidean distance";

        private static readonly string[] axisNames = { "1", "x", "y", "z" };
        private static readonly string[] axisMeanings = { "", "purine/pyrimidine", "amino/keto", "strong/weak" };

        private readonly NpgsqlDataSource dataSource;
        private readonly SemaphoreSlim refreshing = new(1, 1);
        private volatile Snapshot current = Snapshot.Empty;

        public SignatureIndex(NpgsqlDataSource dataSource)
        {
            this.dataSource = dataSource;
        }

        public IReadOnlyDictionary<long, SignatureEntry> Entries => current.Entries;

        /// <returns>Whether the index was reloaded.</returns>
        public async Task<bool> RefreshAsync(CancellationToken cancellationToken = default)
        {
            await refreshing.WaitAsync(cancellationToken);
            try
            {
                var (count, maxId) = await StateAsync(dataSource, cancellationToken);
                if (count == current.RecordCount && maxId == current.MaxRecordId) return false;

                var entries = new Dictionary<long, SignatureEntry>();
                await using var command = dataSource.CreateCommand(
                    "select id, kind, owner_id, dinucleotide_counts from (" +
                    " select r.id, r.kind, s.owner_id, r.dinucleotide_counts," +
                    "  row_number() over (partition by r.source_id, r.kind in ('bacterium', 'archaeon') order by r.length_bp desc, r.id) as n" +
                    " from sequence_record r join sequence_source s on s.id = r.source_id" +
                    " where r.dinucleotide_counts is not null and r.kind in ('bacterium', 'archaeon', 'virus')) x " +
                    "where owner_id is not null or kind = 'virus' or n = 1");
                await using (var reader = await command.ExecuteReaderAsync(cancellationToken))
                    while (await reader.ReadAsync(cancellationToken))
                    {
                        long id = reader.GetInt64(0);
                        entries[id] = new SignatureEntry(id, GroupOf(reader.GetString(1))!, reader.IsDBNull(2) ? null : reader.GetInt64(2),
                            Signature(reader.GetFieldValue<long[]>(3)));
                    }
                current = new Snapshot(entries, count, maxId);
                return true;
            }
            finally
            {
                refreshing.Release();
            }
        }

        internal static async Task<(long Count, long MaxId)> StateAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken)
        {
            await using var command = dataSource.CreateCommand("select count(*), coalesce(max(id), 0) from sequence_record");
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            await reader.ReadAsync(cancellationToken);
            return (reader.GetInt64(0), reader.GetInt64(1));
        }

        /// <summary>The signature of one record's dinucleotide counts: the one definition the index, relatives and trees share.</summary>
        public static double[] Signature(long[] dinucleotideCounts) => GenomeSignatures.KarlinLogRho(dinucleotideCounts);

        /// <returns>The group a record kind belongs to; null for a kind relatives aren't found for (plasmids).</returns>
        public static string? GroupOf(string kind) => kind switch
        {
            "bacterium" or "archaeon" => Cellular,
            "virus" => Viral,
            _ => null,
        };

        /// <summary>
        /// The <paramref name="count"/> genomes of <paramref name="group"/> nearest to <paramref name="signature"/>, nearest first,
        /// among those the viewer may see (public, or the viewer's own uploads); <paramref name="exclude"/> is the query itself.
        /// </summary>
        public IReadOnlyList<Relative> Nearest(double[] signature, string group, long? viewerId, int count, long? exclude = null) =>
            current.Entries.Values
                .Where(e => e.Group == group && e.Id != exclude && (e.OwnerId is null || e.OwnerId == viewerId))
                .Select(e => new Relative(e.Id, Distance(signature, e.Signature)))
                .OrderBy(r => r.Distance).ThenBy(r => r.Id)
                .Take(count)
                .ToList();

        /// <summary>Kind and dinucleotide counts of these records (no visibility check); ids without a row are left out.</summary>
        public static async Task<Dictionary<long, (string Kind, long[]? Counts)>> ReadCountsAsync(NpgsqlDataSource dataSource,
            IReadOnlyCollection<long> ids, CancellationToken cancellationToken = default)
        {
            await using var command = dataSource.CreateCommand("select id, kind, dinucleotide_counts from sequence_record where id = any($1)");
            command.Parameters.AddWithValue(ids.ToArray());
            var rows = new Dictionary<long, (string, long[]?)>();
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
                rows[reader.GetInt64(0)] = (reader.GetString(1), reader.IsDBNull(2) ? null : reader.GetFieldValue<long[]>(2));
            return rows;
        }

        /// <summary>A tree over the leaves' signatures; method is upgma or nj (neighbour-joining, midpoint-rooted for display).</summary>
        public static PhyloTree Tree(IReadOnlyList<(string Label, double[] Signature)> leaves, string method)
        {
            int n = leaves.Count;
            var distances = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    distances[a, b] = distances[b, a] = Distance(leaves[a].Signature, leaves[b].Signature);
            var labels = leaves.Select(l => l.Label).ToArray();
            return method.ToLowerInvariant() switch
            {
                "upgma" => Phylogeny.Upgma(distances, labels),
                "nj" => Phylogeny.NeighbourJoining(distances, labels),
                _ => throw new ArgumentException($"Unknown tree method '{method}'. Use upgma or nj."),
            };
        }

        /// <summary>
        /// The signature as Karlin's per-dinucleotide values (view "karlin", in AA, AC, ... TT order) or as its Z-curve coefficients
        /// (view "zcurve", on the products h_i(n) h_j(n + 1) with h in 1, x, y, z). The same vector either way: the Z-curve view is a
        /// rotation, its Euclidean distances a quarter of the Karlin view's.
        /// </summary>
        public static IReadOnlyList<SignatureValue> Views(long[] dinucleotideCounts, string view)
        {
            var logRho = Signature(dinucleotideCounts);
            switch (view.ToLowerInvariant())
            {
                case "karlin":
                    return Enumerable.Range(0, 16).Select(i => new SignatureValue($"{"ACGT"[i / 4]}{"ACGT"[i % 4]}", logRho[i],
                        "log of observed / expected frequency of this neighbour pair, both strands", Math.Exp(logRho[i]))).ToList();
                case "zcurve":
                    var coefficients = GenomeSignatures.ToZCurveBasis(logRho);
                    return Enumerable.Range(0, 16).Select(k => new SignatureValue(
                        $"{axisNames[k / 4]}(n)·{axisNames[k % 4]}(n+1)", coefficients[k], Meaning(k / 4, k % 4), null)).ToList();
                default:
                    throw new ArgumentException($"Unknown view '{view}'. Use zcurve or karlin.");
            }
        }

        private static string Meaning(int first, int second) => (first, second) switch
        {
            (0, 0) => "mean log ratio over all neighbour pairs",
            (0, _) => $"context of the second base on the {axisMeanings[second]} axis",
            (_, 0) => $"context of the first base on the {axisMeanings[first]} axis",
            _ when first == second => $"{axisMeanings[first]} neighbour run (same class next to same class)",
            _ => $"{axisMeanings[first]} of one base with {axisMeanings[second]} of the next",
        };

        private static double Distance(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; i++) sum += (a[i] - b[i]) * (a[i] - b[i]);
            return Math.Sqrt(sum);
        }

        private sealed record Snapshot(IReadOnlyDictionary<long, SignatureEntry> Entries, long RecordCount, long MaxRecordId)
        {
            public static readonly Snapshot Empty = new(new Dictionary<long, SignatureEntry>(), -1, -1);
        }
    }
}
