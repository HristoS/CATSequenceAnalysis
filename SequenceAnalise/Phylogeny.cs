using System.Globalization;
using System.Text;

namespace SequenceAnalyses
{
    /// <summary>Per taxon with at least two leaves: whether its leaves form a split of the tree, and how close they come.</summary>
    /// <param name="Consistency">Taxon size over the size of the smallest split side holding all its leaves: 1 when monophyletic.</param>
    public readonly record struct TaxonRecovery(string Taxon, int Leaves, bool Monophyletic, double Consistency);

    public sealed class PhyloNode
    {
        public PhyloNode(int leaf, IReadOnlyList<PhyloNode> children)
        {
            Leaf = leaf;
            Children = children;
        }

        /// <summary>The leaf's index into the distance matrix; -1 for an internal node.</summary>
        public int Leaf { get; }
        public IReadOnlyList<PhyloNode> Children { get; }
        public double BranchLength { get; internal set; }
    }

    public sealed class PhyloTree
    {
        private List<ulong[]>? splits;

        public PhyloTree(PhyloNode root, IReadOnlyList<string> labels)
        {
            Root = root;
            Labels = labels;
        }

        public PhyloNode Root { get; }
        public IReadOnlyList<string> Labels { get; }

        public string ToNewick()
        {
            var builder = new StringBuilder();
            Write(Root);
            return builder.Append(';').ToString();

            void Write(PhyloNode node)
            {
                if (node.Leaf >= 0)
                    builder.Append(NewickLabel(Labels[node.Leaf]));
                else
                {
                    builder.Append('(');
                    for (int i = 0; i < node.Children.Count; i++)
                    {
                        if (i > 0) builder.Append(',');
                        Write(node.Children[i]);
                    }
                    builder.Append(')');
                }
                if (node != Root) builder.Append(':').Append(node.BranchLength.ToString("G6", CultureInfo.InvariantCulture));
            }
        }

        private static string NewickLabel(string label) =>
            new string(label.Select(c => "()[]':;, \t".Contains(c) ? '_' : c).ToArray());

        /// <summary>Leaf indexes left to right: the order the tree places the organisms in.</summary>
        public List<int> LeafOrder()
        {
            var order = new List<int>();
            Visit(Root);
            return order;

            void Visit(PhyloNode node)
            {
                if (node.Leaf >= 0) order.Add(node.Leaf);
                foreach (var child in node.Children) Visit(child);
            }
        }

        /// <summary>Path length between every pair of leaves.</summary>
        public double[,] PatristicDistances()
        {
            int n = Labels.Count;
            var result = new double[n, n];
            var depth = new double[n];
            var ancestors = new List<PhyloNode>[n];
            Walk(Root, 0, new List<PhyloNode>());
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                {
                    int shared = 0;
                    while (shared < ancestors[a].Count && shared < ancestors[b].Count && ancestors[a][shared] == ancestors[b][shared]) shared++;
                    double meet = ancestors[a].Take(shared).Skip(1).Sum(node => node.BranchLength);
                    result[a, b] = result[b, a] = depth[a] + depth[b] - 2 * meet;
                }
            return result;

            void Walk(PhyloNode node, double length, List<PhyloNode> path)
            {
                path.Add(node);
                if (node.Leaf >= 0)
                {
                    depth[node.Leaf] = length;
                    ancestors[node.Leaf] = new List<PhyloNode>(path);
                }
                foreach (var child in node.Children) Walk(child, length + child.BranchLength, path);
                path.RemoveAt(path.Count - 1);
            }
        }

        /// <summary>
        /// Recovery of every taxon with at least two leaves (singletons are trivially monophyletic, so they're left out).
        /// Judged on the unrooted tree: a taxon is monophyletic when some edge separates exactly its leaves, so the result
        /// doesn't depend on where the tree was rooted.
        /// </summary>
        public List<TaxonRecovery> Recovery(IReadOnlyList<string> taxa)
        {
            var bipartitions = Splits();
            int n = Labels.Count, words = (n + 63) / 64;
            var result = new List<TaxonRecovery>();
            foreach (var group in Enumerable.Range(0, n).GroupBy(i => taxa[i]).Where(g => g.Count() >= 2))
            {
                var members = new ulong[words];
                foreach (int leaf in group) members[leaf / 64] |= 1UL << (leaf % 64);
                int size = group.Count(), smallest = n;
                foreach (var side in bipartitions)
                {
                    int sideSize = PopCount(side);
                    if (sideSize >= size && sideSize < smallest && Contains(side, members)) smallest = sideSize;
                    var other = Complement(side, n);
                    if (n - sideSize >= size && n - sideSize < smallest && Contains(other, members)) smallest = n - sideSize;
                }
                result.Add(new TaxonRecovery(group.Key, size, smallest == size, (double)size / smallest));
            }
            return result;
        }

        /// <summary>
        /// Normalized Robinson-Foulds distance to a tree over the same leaves: the share of non-trivial splits found in only
        /// one of the two (0 identical, 1 nothing shared). Unrooted, like <see cref="Recovery"/>.
        /// </summary>
        public double NormalizedRobinsonFoulds(PhyloTree other)
        {
            int n = Labels.Count;
            if (other.Labels.Count != n) throw new ArgumentException("The trees have different leaf counts.", nameof(other));
            var mine = Splits().Select(s => Canonical(s, n)).Distinct().ToHashSet();
            var theirs = other.Splits().Select(s => Canonical(s, n)).Distinct().ToHashSet();
            int total = mine.Count + theirs.Count;
            return total == 0 ? 0 : (total - 2.0 * mine.Count(theirs.Contains)) / total;
        }

        // The side without leaf 0, as a string key, so a split and its complement compare equal.
        private static string Canonical(ulong[] side, int n) =>
            string.Join(",", ((side[0] & 1) == 0 ? side : Complement(side, n)).Select(w => w.ToString("x")));

        // One side of every internal edge, as a leaf bit set.
        private List<ulong[]> Splits()
        {
            if (splits != null) return splits;
            int n = Labels.Count, words = (n + 63) / 64;
            splits = new List<ulong[]>();
            Collect(Root);
            return splits;

            ulong[] Collect(PhyloNode node)
            {
                var set = new ulong[words];
                if (node.Leaf >= 0) set[node.Leaf / 64] |= 1UL << (node.Leaf % 64);
                foreach (var child in node.Children)
                {
                    var below = Collect(child);
                    for (int w = 0; w < words; w++) set[w] |= below[w];
                }
                int size = PopCount(set);
                if (node != Root && size >= 2 && size <= n - 2) splits.Add(set);
                return set;
            }
        }

        private static int PopCount(ulong[] set) => set.Sum(w => System.Numerics.BitOperations.PopCount(w));

        private static bool Contains(ulong[] set, ulong[] members)
        {
            for (int w = 0; w < set.Length; w++)
                if ((members[w] & ~set[w]) != 0) return false;
            return true;
        }

        private static ulong[] Complement(ulong[] set, int n)
        {
            var result = new ulong[set.Length];
            for (int w = 0; w < set.Length; w++) result[w] = ~set[w];
            if (n % 64 != 0) result[^1] &= (1UL << (n % 64)) - 1;
            return result;
        }
    }

    public static class Phylogeny
    {
        // Saitou & Nei 1987; the formulas as in https://en.wikipedia.org/wiki/Neighbor_joining. A negative branch length is set to
        // 0 and the difference moved to its sister, so the pair's distance is kept. The unrooted result is rooted at the midpoint
        // of its longest leaf-to-leaf path, for display only: recovery is judged unrooted.
        public static PhyloTree NeighbourJoining(double[,] distances, IReadOnlyList<string> labels)
        {
            int n = labels.Count;
            if (n < 2) throw new ArgumentException("A tree needs at least two leaves.", nameof(labels));
            var adjacency = new List<List<(int Node, double Length)>>();
            for (int i = 0; i < n; i++) adjacency.Add(new List<(int, double)>());

            // A plain matrix, read from the upper triangle with a zero diagonal. A joined pair's new node takes over the first one's
            // row, so the matrix never grows; node[row] is the tree node a row stands for. Sums run over active rows in list order and
            // the new node goes last, as a growing node list would put it, so ties break as they always have.
            var d = new double[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    d[a, b] = d[b, a] = distances[a, b];
            var node = Enumerable.Range(0, n).ToArray();
            var sums = new double[n];

            void Connect(int a, int b, double length)
            {
                adjacency[a].Add((b, length));
                adjacency[b].Add((a, length));
            }

            var active = Enumerable.Range(0, n).ToList();
            while (active.Count > 2)
            {
                int r = active.Count;
                foreach (int i in active)
                {
                    double sum = 0;
                    foreach (int k in active) sum += d[i, k];
                    sums[i] = sum;
                }
                (int I, int J) best = (-1, -1);
                double bestQ = double.PositiveInfinity;
                for (int x = 0; x < r; x++)
                    for (int y = x + 1; y < r; y++)
                    {
                        int i = active[x], j = active[y];
                        double q = (r - 2) * d[i, j] - sums[i] - sums[j];
                        if (q < bestQ) (bestQ, best) = (q, (i, j));
                    }

                var (bi, bj) = best;
                double dij = d[bi, bj];
                double li = dij / 2 + (sums[bi] - sums[bj]) / (2.0 * (r - 2));
                if (li < 0) li = 0;
                if (li > dij) li = dij;
                double lj = dij - li;

                int u = adjacency.Count;
                adjacency.Add(new List<(int, double)>());
                Connect(u, node[bi], li);
                Connect(u, node[bj], lj);
                foreach (int k in active)
                    if (k != bi && k != bj) d[bi, k] = d[k, bi] = (d[bi, k] + d[bj, k] - dij) / 2;
                node[bi] = u;
                active.Remove(bi);
                active.Remove(bj);
                active.Add(bi);
            }
            Connect(node[active[0]], node[active[1]], Math.Max(0, d[active[0], active[1]]));

            return new PhyloTree(MidpointRoot(adjacency, n), labels);
        }

        private static PhyloNode MidpointRoot(List<List<(int Node, double Length)>> adjacency, int leaves)
        {
            (int Far, double Distance, int[] Parent) Farthest(int from)
            {
                var distance = new double[adjacency.Count];
                var parent = new int[adjacency.Count];
                Array.Fill(distance, -1);
                distance[from] = 0;
                parent[from] = -1;
                var stack = new Stack<int>();
                stack.Push(from);
                while (stack.Count > 0)
                {
                    int node = stack.Pop();
                    foreach (var (next, length) in adjacency[node])
                        if (distance[next] < 0)
                        {
                            distance[next] = distance[node] + length;
                            parent[next] = node;
                            stack.Push(next);
                        }
                }
                int far = Enumerable.Range(0, leaves).MaxBy(i => distance[i]);
                return (far, distance[far], parent);
            }

            var (a, _, _) = Farthest(0);
            var (b, total, parentToA) = Farthest(a);
            // Walk from b towards a until the midpoint falls on an edge; root there.
            double walked = 0;
            int lower = b;
            while (true)
            {
                int upper = parentToA[lower];
                double length = adjacency[lower].First(e => e.Node == upper).Length;
                if (walked + length >= total / 2 || parentToA[upper] < 0)
                {
                    double belowRoot = Math.Clamp(total / 2 - walked, 0, length);
                    var left = Build(lower, upper);
                    left.BranchLength = belowRoot;
                    var right = Build(upper, lower);
                    right.BranchLength = length - belowRoot;
                    return new PhyloNode(-1, new[] { left, right });
                }
                walked += length;
                lower = upper;
            }

            PhyloNode Build(int node, int from)
            {
                var children = new List<PhyloNode>();
                foreach (var (next, length) in adjacency[node])
                {
                    if (next == from) continue;
                    var child = Build(next, node);
                    child.BranchLength = length;
                    children.Add(child);
                }
                return new PhyloNode(node < leaves ? node : -1, children);
            }
        }

        // Average linkage (Sokal & Michener 1958): clusters merge at half their mean pairwise distance, so the tree is ultrametric.
        public static PhyloTree Upgma(double[,] distances, IReadOnlyList<string> labels)
        {
            int n = labels.Count;
            if (n < 2) throw new ArgumentException("A tree needs at least two leaves.", nameof(labels));
            var clusters = Enumerable.Range(0, n).Select(i => (Node: new PhyloNode(i, Array.Empty<PhyloNode>()), Height: 0.0, Size: 1)).ToList();
            var d = new double[n, n];
            Array.Copy(distances, d, distances.Length);
            var ids = Enumerable.Range(0, n).ToList();   // row of d for each live cluster

            while (clusters.Count > 1)
            {
                int bestX = 0, bestY = 1;
                for (int x = 0; x < clusters.Count; x++)
                    for (int y = x + 1; y < clusters.Count; y++)
                        if (d[ids[x], ids[y]] < d[ids[bestX], ids[bestY]]) (bestX, bestY) = (x, y);

                var (left, right) = (clusters[bestX], clusters[bestY]);
                double height = d[ids[bestX], ids[bestY]] / 2;
                left.Node.BranchLength = Math.Max(0, height - left.Height);
                right.Node.BranchLength = Math.Max(0, height - right.Height);
                int size = left.Size + right.Size;
                int row = ids[bestX];
                for (int k = 0; k < clusters.Count; k++)
                {
                    if (k == bestX || k == bestY) continue;
                    double merged = (d[ids[bestX], ids[k]] * left.Size + d[ids[bestY], ids[k]] * right.Size) / size;
                    d[row, ids[k]] = d[ids[k], row] = merged;
                }
                clusters[bestX] = (new PhyloNode(-1, new[] { left.Node, right.Node }), height, size);
                clusters.RemoveAt(bestY);
                ids.RemoveAt(bestY);
            }
            return new PhyloTree(clusters[0].Node, labels);
        }
    }
}
