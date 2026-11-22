using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <param name="Group">Relatives are found within a group: <see cref="SignatureIndex.Cellular"/> or <see cref="SignatureIndex.Viral"/>.</param>
    /// <param name="CvTree">The genome's <see cref="CvTreeSignature"/> once built; null for a genome without stored amino-acid 4-mers, a
    /// eukaryote, or one still waiting for its vector.</param>
    public sealed record SignatureEntry(long Id, string Group, long? OwnerId, double[] Signature, short[]? CvTree = null);

    public sealed record Relative(long Id, double Distance);

    /// <param name="Signature">What was compared: <see cref="SignatureIndex.Description"/> or <see cref="CvTreeSignature.Description"/>.</param>
    public sealed record Relatives(string Signature, IReadOnlyList<Relative> Items);

    /// <summary>One value of a signature view.</summary>
    /// <param name="Ratio">The Karlin view's rho* itself (exp of the log value); null in the Z-curve view.</param>
    public sealed record SignatureValue(string Name, double Value, string Meaning, double? Ratio);

    /// <summary>
    /// Two signatures of every genome, in memory, for nearest relatives and trees (P5 S5.5): amino-acid CVTree
    /// (<see cref="CvTreeSignature"/>) for bacteria, archaea and viruses with stored 4-mers, and the dinucleotide signature for every
    /// genome. CVTree finds a genome's genus about 13 to 16 points more often among its 10 nearest (results/README.md §4, 2026-11-22),
    /// so a query that has a vector is ranked by it among the genomes that have one; any other query (a eukaryote, an upload, a record
    /// without a gene layer) by the dinucleotide signature among all. Vectors take minutes to build, so they are built in the
    /// background, on <c>vectorParallelism</c> threads, after the refresh that finds them missing; until every vector is in, every query
    /// uses the dinucleotide signature, and each answer says which signature it used.
    ///
    /// The dinucleotide signature The signature is log rho*
    /// over both strands (<see cref="GenomeSignatures.KarlinLogRho"/>), compared by Euclidean distance: Karlin's genomic signature, which
    /// is also the complete second-order Z-curve, since its Z-curve coefficients (<see cref="GenomeSignatures.ToZCurveBasis"/>) are a
    /// rotation of it; <see cref="Views"/> serves either. Log rather than ratio because it treats depletion and enrichment alike; the
    /// four scale and distance choices measured within about 0.01 of each other (results/S55_*.csv). It comes from the stored
    /// dinucleotide counts (D8), so a profile-only upload has none.
    ///
    /// One record stands for each genome (D9): the largest chromosome of a catalogue bacterium or archaeon, every virus record, a
    /// eukaryote's one record (its nuclear counts summed), and every record of an upload (an upload may hold many genomes); plasmids
    /// are left out. <see cref="RefreshAsync"/> reloads when the
    /// record count, the highest record id or the number of stored 4-mer rows changed, which an insert, a deletion or a gene backfill
    /// always does; vectors already built are kept (a rewritten 4-mer row of the same record keeps its old vector until a restart). Each
    /// refresh publishes a new snapshot, so a query keeps a consistent view while one runs.
    /// </summary>
    public sealed class SignatureIndex
    {
        public const string Cellular = "cellular";
        public const string Viral = "viral";
        public const string Eukaryotic = "eukaryotic";
        public const int MaxTreeLeaves = 500;
        public const string Description = "log rho* over both strands (dinucleotides), Euclidean distance";
        // Every request refreshes; counting the 4-mer rows costs about 30 ms, and only a gene backfill changes them.
        private static readonly TimeSpan kmerCheckInterval = TimeSpan.FromMinutes(1);

        private static readonly string[] axisNames = { "1", "x", "y", "z" };
        private static readonly string[] axisMeanings = { "", "purine/pyrimidine", "amino/keto", "strong/weak" };

        private readonly NpgsqlDataSource dataSource;
        private readonly int vectorParallelism;
        private readonly SemaphoreSlim refreshing = new(1, 1);
        private volatile Snapshot current = Snapshot.Empty;
        private readonly object buildLock = new();
        private Task building = Task.CompletedTask;
        private DateTime kmersCheckedAt = DateTime.MinValue;
        // Records whose stored 4-mers hold no word at all, so they get no vector and aren't waited for again.
        private readonly HashSet<long> withoutVector = new();

        /// <param name="vectorParallelism">Threads that build CVTree vectors; 0: half the cores, at least one, so searches keep the rest.</param>
        public SignatureIndex(NpgsqlDataSource dataSource, int vectorParallelism = 0)
        {
            this.dataSource = dataSource;
            this.vectorParallelism = vectorParallelism > 0 ? vectorParallelism : Math.Max(1, Environment.ProcessorCount / 2);
        }

        public IReadOnlyDictionary<long, SignatureEntry> Entries => current.Entries;

        /// <summary>Whether every genome that can have a CVTree vector has it, so relatives and trees use CVTree.</summary>
        public bool VectorsReady => current.Pending.Count == 0;

        /// <summary>Completes when the vectors the refreshes so far found missing are built, or their build failed.</summary>
        public Task VectorsBuiltAsync()
        {
            lock (buildLock) return building;
        }

        /// <returns>Whether the index was reloaded.</returns>
        public async Task<bool> RefreshAsync(CancellationToken cancellationToken = default)
        {
            await refreshing.WaitAsync(cancellationToken);
            try
            {
                var (count, maxId) = await StateAsync(dataSource, cancellationToken);
                var previous = current;
                bool recordsChanged = count != previous.RecordCount || maxId != previous.MaxRecordId;
                long kmerRows = previous.KmerRows;
                if (recordsChanged || DateTime.UtcNow - kmersCheckedAt > kmerCheckInterval)
                {
                    kmerRows = await KmerRowsAsync(cancellationToken);
                    kmersCheckedAt = DateTime.UtcNow;
                }
                bool changed = recordsChanged || kmerRows != previous.KmerRows;
                if (changed)
                    current = await LoadAsync(previous, count, maxId, kmerRows, cancellationToken);
                // Also after a failed build: the next refresh tries again.
                lock (buildLock)
                    if (current.Pending.Count > 0 && building.IsCompleted) building = Task.Run(BuildVectorsAsync, CancellationToken.None);
                return changed;
            }
            finally
            {
                refreshing.Release();
            }
        }

        private async Task<Snapshot> LoadAsync(Snapshot previous, long count, long maxId, long kmerRows, CancellationToken cancellationToken)
        {
            var entries = new Dictionary<long, SignatureEntry>();
            var pending = new HashSet<long>();
            await using var command = dataSource.CreateCommand(
                "select id, kind, owner_id, dinucleotide_counts, has_kmers from (" +
                " select r.id, r.kind, s.owner_id, r.dinucleotide_counts," +
                "  exists (select 1 from protein_kmer p where p.sequence_id = r.id and p.k = $1 and p.status = $2) as has_kmers," +
                "  row_number() over (partition by r.source_id, r.kind in ('bacterium', 'archaeon') order by r.length_bp desc, r.id) as n" +
                " from sequence_record r join sequence_source s on s.id = r.source_id" +
                " where r.dinucleotide_counts is not null and r.kind in ('bacterium', 'archaeon', 'virus', 'eukaryote')) x " +
                "where owner_id is not null or kind = 'virus' or n = 1");
            command.Parameters.AddWithValue((short)CvTreeSignature.K);
            command.Parameters.AddWithValue(ProteinKmers.Done);
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
            {
                long id = reader.GetInt64(0);
                var group = GroupOf(reader.GetString(1))!;
                bool wanted = reader.GetBoolean(4) && HasCvTree(group);
                var vector = wanted && previous.Entries.TryGetValue(id, out var old) ? old.CvTree : null;
                if (wanted && vector is null)
                    lock (withoutVector)
                        if (!withoutVector.Contains(id)) pending.Add(id);
                entries[id] = new SignatureEntry(id, group, reader.IsDBNull(2) ? null : reader.GetInt64(2),
                    Signature(reader.GetFieldValue<long[]>(3)), vector);
            }
            return new Snapshot(entries, count, maxId, kmerRows, pending);
        }

        // Builds the vectors the current snapshot is missing and publishes them into whatever snapshot is current by then (a refresh may
        // have run meanwhile); loops while that one is missing more.
        private async Task BuildVectorsAsync()
        {
            while (current.Pending.Count > 0)
            {
                var requested = current.Pending.ToHashSet();
                var built = await CvTreeSignature.ReadAsync(dataSource, requested, vectorParallelism);
                await refreshing.WaitAsync();
                try
                {
                    var snapshot = current;
                    var entries = new Dictionary<long, SignatureEntry>(snapshot.Entries.Count);
                    foreach (var (id, entry) in snapshot.Entries)
                        entries[id] = entry.CvTree is null && built.TryGetValue(id, out var vector) ? entry with { CvTree = vector } : entry;
                    lock (withoutVector)
                        withoutVector.UnionWith(requested.Where(id => !built.ContainsKey(id)));
                    current = new Snapshot(entries, snapshot.RecordCount, snapshot.MaxRecordId, snapshot.KmerRows,
                        snapshot.Pending.Where(id => !requested.Contains(id)).ToHashSet());
                }
                finally
                {
                    refreshing.Release();
                }
            }
        }

        private async Task<long> KmerRowsAsync(CancellationToken cancellationToken)
        {
            await using var command = dataSource.CreateCommand("select count(*) from protein_kmer where k = $1 and status = $2");
            command.Parameters.AddWithValue((short)CvTreeSignature.K);
            command.Parameters.AddWithValue(ProteinKmers.Done);
            return (long)(await command.ExecuteScalarAsync(cancellationToken))!;
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
            "eukaryote" => Eukaryotic,
            _ => null,
        };

        /// <summary>
        /// The <paramref name="count"/> genomes of <paramref name="group"/> nearest to <paramref name="signature"/>, nearest first,
        /// among those the viewer may see (public, or the viewer's own uploads); <paramref name="exclude"/> is the query itself.
        /// </summary>
        public IReadOnlyList<Relative> Nearest(double[] signature, string group, long? viewerId, int count, long? exclude = null) =>
            Nearest(current, signature, group, viewerId, count, exclude);

        private static IReadOnlyList<Relative> Nearest(Snapshot snapshot, double[] signature, string group, long? viewerId, int count, long? exclude) =>
            snapshot.Entries.Values
                .Where(e => e.Group == group && e.Id != exclude && (e.OwnerId is null || e.OwnerId == viewerId))
                .Select(e => new Relative(e.Id, Distance(signature, e.Signature)))
                .OrderBy(r => r.Distance).ThenBy(r => r.Id)
                .Take(count)
                .ToList();

        /// <summary>
        /// The <paramref name="count"/> relatives of record <paramref name="id"/> of <paramref name="group"/> the viewer may see, nearest
        /// first: by CVTree among the group's genomes with a vector when every vector is built and the record has one (its index
        /// entry's, or built from its stored 4-mers when the index doesn't hold it, as a second chromosome); otherwise by the
        /// dinucleotide signature of <paramref name="dinucleotideCounts"/> among all.
        /// </summary>
        public async Task<Relatives> RelativesAsync(long id, long[] dinucleotideCounts, string group, long? viewerId, int count,
            CancellationToken cancellationToken = default)
        {
            var snapshot = current;
            var vector = await VectorOfAsync(snapshot, id, group, cancellationToken);
            if (vector is null) return new Relatives(Description, Nearest(snapshot, Signature(dinucleotideCounts), group, viewerId, count, id));
            return new Relatives(CvTreeSignature.Description, NearestByCvTree(snapshot, vector, group, viewerId, count, id));
        }

        // One pass keeping the best count in order (ties by id), over the group's genomes with a vector in id order.
        private static List<Relative> NearestByCvTree(Snapshot snapshot, short[] vector, string group, long? viewerId, int count, long exclude)
        {
            var best = new List<Relative>(count + 1);
            if (!snapshot.WithVector.TryGetValue(group, out var candidates)) return best;
            foreach (var e in candidates)
            {
                if (e.Id == exclude || (e.OwnerId is not null && e.OwnerId != viewerId)) continue;
                double distance = CvTreeSignature.Distance(vector, e.CvTree!);
                if (best.Count == count && distance >= best[^1].Distance) continue;
                int slot = best.Count;
                while (slot > 0 && best[slot - 1].Distance > distance) slot--;
                best.Insert(slot, new Relative(e.Id, distance));
                if (best.Count > count) best.RemoveAt(count);
            }
            return best;
        }

        private async Task<short[]?> VectorOfAsync(Snapshot snapshot, long id, string group, CancellationToken cancellationToken)
        {
            if (!HasCvTree(group) || snapshot.Pending.Count > 0) return null;
            if (snapshot.Entries.TryGetValue(id, out var entry)) return entry.CvTree;
            return (await CvTreeSignature.ReadAsync(dataSource, [id], 1, cancellationToken)).GetValueOrDefault(id);
        }

        // Eukaryotes keep the dinucleotide signature: 131 of 587 had stored 4-mers on 2026-11-22.
        private static bool HasCvTree(string group) => group is Cellular or Viral;

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

        /// <summary>
        /// A tree over records of one group, labelled by id: by CVTree when every vector is built and every leaf's genome has one in the
        /// index, otherwise by the dinucleotide signatures of the leaves' counts.
        /// </summary>
        public (PhyloTree Tree, string Signature) Tree(IReadOnlyList<(long Id, long[] Counts)> leaves, string group, string method)
        {
            var snapshot = current;
            var labels = leaves.Select(l => l.Id.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray();
            var vectors = leaves.Select(l => snapshot.Entries.TryGetValue(l.Id, out var e) ? e.CvTree : null).ToArray();
            if (HasCvTree(group) && snapshot.Pending.Count == 0 && vectors.All(v => v is not null))
                return (Tree(labels, (a, b) => CvTreeSignature.Distance(vectors[a]!, vectors[b]!), method), CvTreeSignature.Description);
            var signatures = leaves.Select(l => Signature(l.Counts)).ToArray();
            return (Tree(labels, (a, b) => Distance(signatures[a], signatures[b]), method), Description);
        }

        /// <summary>A tree over the leaves' dinucleotide signatures; method is upgma or nj (neighbour-joining, midpoint-rooted for display).</summary>
        public static PhyloTree Tree(IReadOnlyList<(string Label, double[] Signature)> leaves, string method) =>
            Tree(leaves.Select(l => l.Label).ToArray(), (a, b) => Distance(leaves[a].Signature, leaves[b].Signature), method);

        private static PhyloTree Tree(string[] labels, Func<int, int, double> distance, string method)
        {
            int n = labels.Length;
            var distances = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    distances[a, b] = distances[b, a] = distance(a, b);
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

        /// <param name="pending">Genomes whose CVTree vector is wanted and not built yet.</param>
        private sealed class Snapshot(IReadOnlyDictionary<long, SignatureEntry> entries, long recordCount, long maxRecordId, long kmerRows,
            IReadOnlySet<long> pending)
        {
            public static readonly Snapshot Empty = new(new Dictionary<long, SignatureEntry>(), -1, -1, -1, new HashSet<long>());

            public IReadOnlyDictionary<long, SignatureEntry> Entries { get; } = entries;
            public long RecordCount { get; } = recordCount;
            public long MaxRecordId { get; } = maxRecordId;
            public long KmerRows { get; } = kmerRows;
            public IReadOnlySet<long> Pending { get; } = pending;

            /// <summary>Each group's genomes with a CVTree vector, in id order: what a CVTree query scans.</summary>
            public IReadOnlyDictionary<string, SignatureEntry[]> WithVector { get; } = entries.Values.Where(e => e.CvTree is not null)
                .GroupBy(e => e.Group).ToDictionary(g => g.Key, g => g.OrderBy(e => e.Id).ToArray());
        }
    }
}
