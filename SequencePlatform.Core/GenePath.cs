using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// A record's Z-curve path built from its gene layer alone: its stored genes in chromosome order (from the replication origin when
    /// asked; a eukaryote's sequences one after another, or one sequence alone), each window adding its window positions times its means of x, y and z over codon positions 1 and 2, a reverse-strand gene
    /// flipped onto the forward strand (its windows in reverse order, x and y negated). Intergenic DNA and third codon positions are left
    /// out. On K-12 against two relatives and B. subtilis it lines relatives up a little less well than the genome-window path and
    /// separates the unrelated genome a little more (results/GenePath_Check.csv), so it's a view of the coding strand, not a better path.
    /// </summary>
    public static class GenePath
    {
        /// <param name="FirstPoint">The point the gene's stretch starts at; it runs <see cref="ZCurveGeneProfile.Gene.WindowCount"/> points on.</param>
        public sealed record Span(int Gene, int FirstPoint);

        /// <param name="Points">Running sums, starting at the origin of the axes; one point per gene window, after a leading zero.</param>
        /// <param name="Genes">Each gene's stretch, in path order.</param>
        /// <param name="StartBp">Where the path starts on the record: its replication origin, or 0.</param>
        public sealed record Result(IReadOnlyList<(double X, double Y, double Z)> Points, IReadOnlyList<Span> Genes, long StartBp);

        /// <param name="sequence">Only the genes on this sequence of the record (<see cref="ZCurveGeneProfile.Gene.Chromosome"/>), the way a
        /// eukaryote's chromosomes are drawn one at a time; null for all of them, sequence by sequence.</param>
        public static Result Build(ZCurveGeneProfile profile, long lengthBp, long startBp = 0, int? sequence = null)
        {
            long Offset(int g) => ((profile.Genes[g].Start - startBp) % lengthBp + lengthBp) % lengthBp;
            var order = Enumerable.Range(0, profile.Genes.Count).Where(g => sequence is null || profile.Genes[g].Chromosome == sequence)
                .OrderBy(g => profile.Genes[g].Chromosome).ThenBy(Offset).ThenBy(g => g).ToList();
            var points = new List<(double X, double Y, double Z)>(order.Sum(g => profile.Genes[g].WindowCount) + 1) { (0, 0, 0) };
            var spans = new List<Span>(order.Count);
            foreach (int g in order)
            {
                var gene = profile.Genes[g];
                int positions = ZCurveGeneProfile.WindowOf(gene.Residues);
                var windows = profile.WindowsOf(g);
                spans.Add(new Span(g, points.Count - 1));
                for (int i = 0; i < gene.WindowCount; i++)
                {
                    int w = gene.Reverse ? gene.WindowCount - 1 - i : i, f = w * ZCurveGeneProfile.FeatureCount;
                    double sign = gene.Reverse ? -1 : 1;
                    var last = points[^1];
                    points.Add((last.X + sign * positions * windows[f], last.Y + sign * positions * windows[f + 1], last.Z + positions * windows[f + 2]));
                }
            }
            return new Result(points, spans, startBp);
        }
    }
}
