using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class PhylogenyTests
    {
        private static readonly string[] fiveTaxa = { "a", "b", "c", "d", "e" };

        private static double[,] Matrix(double[][] rows)
        {
            var result = new double[rows.Length, rows.Length];
            for (int i = 0; i < rows.Length; i++)
                for (int j = 0; j < rows.Length; j++)
                    result[i, j] = rows[i][j];
            return result;
        }

        private static void AssertDistances(double[,] actual, double[,] expected)
        {
            for (int i = 0; i < expected.GetLength(0); i++)
                for (int j = 0; j < expected.GetLength(1); j++)
                    Assert.That(actual[i, j], Is.EqualTo(expected[i, j]).Within(1e-9), $"leaves {i}, {j}");
        }

        // https://en.wikipedia.org/wiki/Neighbor_joining#Example: the matrix is additive (a:2, b:3 on u; u:3, c:4 on v; v:1, d:2,
        // e:1 on w), so the tree NJ builds must reproduce every distance as a path length, which fixes topology and branch lengths.
        [Test]
        public void NeighbourJoiningRecoversTheTextbookTree()
        {
            var distances = Matrix(new[]
            {
                new double[] { 0, 5, 9, 9, 8 },
                new double[] { 5, 0, 10, 10, 9 },
                new double[] { 9, 10, 0, 8, 7 },
                new double[] { 9, 10, 8, 0, 3 },
                new double[] { 8, 9, 7, 3, 0 },
            });

            var tree = Phylogeny.NeighbourJoining(distances, fiveTaxa);

            AssertDistances(tree.PatristicDistances(), distances);
            var recovery = tree.Recovery(new[] { "ab", "ab", "cde", "cde", "cde" });
            Assert.That(recovery.All(r => r.Monophyletic && r.Consistency == 1));
        }

        // https://en.wikipedia.org/wiki/UPGMA#Working_example (5S rRNA): (a,b) at 8.5, then e at 11, (c,d) at 14, root at 16.5.
        [Test]
        public void UpgmaReproducesTheTextbookHeights()
        {
            var distances = Matrix(new[]
            {
                new double[] { 0, 17, 21, 31, 23 },
                new double[] { 17, 0, 30, 34, 21 },
                new double[] { 21, 30, 0, 28, 39 },
                new double[] { 31, 34, 28, 0, 43 },
                new double[] { 23, 21, 39, 43, 0 },
            });
            var ultrametric = Matrix(new[]
            {
                new double[] { 0, 17, 33, 33, 22 },
                new double[] { 17, 0, 33, 33, 22 },
                new double[] { 33, 33, 0, 28, 33 },
                new double[] { 33, 33, 28, 0, 33 },
                new double[] { 22, 22, 33, 33, 0 },
            });

            var tree = Phylogeny.Upgma(distances, fiveTaxa);

            AssertDistances(tree.PatristicDistances(), ultrametric);
            Assert.That(tree.ToNewick(), Is.EqualTo("(((a:8.5,b:8.5):2.5,e:11):5.5,(c:14,d:14):2.5);"));
        }

        [Test]
        public void NeighbourJoiningIsExactOnARandomAdditiveTree()
        {
            var rng = new Random(7);
            const int leaves = 40;
            var reference = Phylogeny.Upgma(RandomMatrix(leaves, rng), Enumerable.Range(0, leaves).Select(i => $"t{i}").ToArray());
            // Make it non-ultrametric (NJ must not rely on a clock) by lengthening every leaf edge differently.
            var distances = reference.PatristicDistances();
            var extra = Enumerable.Range(0, leaves).Select(_ => rng.NextDouble()).ToArray();
            for (int i = 0; i < leaves; i++)
                for (int j = 0; j < leaves; j++)
                    if (i != j) distances[i, j] += extra[i] + extra[j];

            var tree = Phylogeny.NeighbourJoining(distances, reference.Labels);

            AssertDistances(tree.PatristicDistances(), distances);
        }

        [Test]
        public void RecoveryIsJudgedOnTheUnrootedTreeAndScoresNearMisses()
        {
            // ((a,b),(c,(d,e))) unrooted: {c,d,e} is a split even though the root sits beside it; {a,c} is not, and the smallest
            // side holding both is {a,b,c} (the other side of {d,e}), so its consistency is 2/3. The singleton taxon is left out.
            var distances = Matrix(new[]
            {
                new double[] { 0, 2, 4, 4, 4 },
                new double[] { 2, 0, 4, 4, 4 },
                new double[] { 4, 4, 0, 3, 3 },
                new double[] { 4, 4, 3, 0, 2 },
                new double[] { 4, 4, 3, 2, 0 },
            });
            var tree = Phylogeny.Upgma(distances, fiveTaxa);

            var byTaxon = tree.Recovery(new[] { "x", "y", "x", "z", "z" }).ToDictionary(r => r.Taxon);
            Assert.That(byTaxon.Keys, Is.EquivalentTo(new[] { "x", "z" }));
            Assert.That(byTaxon["z"].Monophyletic, Is.True);
            Assert.That(byTaxon["x"].Monophyletic, Is.False);
            Assert.That(byTaxon["x"].Consistency, Is.EqualTo(2.0 / 3));

            var cde = tree.Recovery(new[] { "p", "p", "q", "q", "q" });
            Assert.That(cde.All(r => r.Monophyletic));
            Assert.That(tree.LeafOrder(), Is.EquivalentTo(Enumerable.Range(0, 5)));

            // NJ on the same additive matrix finds the same two splits, rooted elsewhere; ((a,c),(b,(d,e))) shares only {d,e}.
            Assert.That(Phylogeny.NeighbourJoining(distances, fiveTaxa).NormalizedRobinsonFoulds(tree), Is.EqualTo(0));
            var swapped = Matrix(new[]
            {
                new double[] { 0, 4, 2, 4, 4 },
                new double[] { 4, 0, 4, 3, 3 },
                new double[] { 2, 4, 0, 4, 4 },
                new double[] { 4, 3, 4, 0, 2 },
                new double[] { 4, 3, 4, 2, 0 },
            });
            Assert.That(Phylogeny.Upgma(swapped, fiveTaxa).NormalizedRobinsonFoulds(tree), Is.EqualTo(0.5));
        }

        private static double[,] RandomMatrix(int n, Random rng)
        {
            var result = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    result[i, j] = result[j, i] = 1 + rng.NextDouble() * 10;
            return result;
        }
    }
}
