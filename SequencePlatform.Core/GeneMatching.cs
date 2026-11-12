using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// A gene's best matches among another record's stored genes, by the Z-curve gene score alone (<see cref="ZCurveGeneProfile.Score"/>,
    /// the query being the gene's own stored windows), each gene compared only with genes of its own kind: from 500 aa the other long
    /// genes (windows of 333 positions from the start codon), below that the genes within ±10% of its length (5 windows each). A match is
    /// reciprocal when it is the best and its own best match back in the first record is the gene: the usual way orthologs are paired.
    /// Between two Podospora and two Nelumbo genomes, 99.8% and 92.9% of reciprocal pairs with named proteins on both sides carry the same NCBI protein name, against 0.5% and 0.2% by chance, as good as amino-acid 4-mers
    /// (gene_match_check.cs, results/GeneMatch_Check.csv); between distant species it hasn't been measured. Used where no family is
    /// annotated, as on NCBI eukaryotes. Within one record a gene isn't its own match, so the matches are its closest paralogs.
    /// </summary>
    public static class GeneMatching
    {
        /// <param name="Gene">Its index in the other record's genes.</param>
        /// <param name="Reciprocal">A reciprocal best hit: the best match, and its own best match back in the first record is the queried gene.</param>
        public sealed record Match(int Gene, double Score, bool Reciprocal);

        /// <returns>The <paramref name="top"/> best matches, best first, and how many genes of <paramref name="other"/> were compared.</returns>
        public static (IReadOnlyList<Match> Matches, int Compared) Find(ZCurveGeneProfile profile, int gene, ZCurveGeneProfile other, int top,
            bool sameRecord = false)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(gene);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(gene, profile.Genes.Count);
            var ranked = Rank(profile, gene, other, sameRecord ? gene : -1, out int compared);
            var matches = ranked.Take(top).Select((m, rank) =>
                new Match(m.Gene, m.Score, rank == 0 && Best(other, m.Gene, profile, sameRecord ? m.Gene : -1) == gene)).ToList();
            return (matches, compared);
        }

        // Every candidate of `other` with its score, best first (ties by index); `skip` is a gene of `other` left out.
        private static List<(int Gene, double Score)> Rank(ZCurveGeneProfile profile, int gene, ZCurveGeneProfile other, int skip, out int compared)
        {
            var own = profile.Genes[gene];
            var query = profile.WindowsOf(gene);
            var candidate = Comparable(own.Residues);
            var scored = new List<(int Gene, double Score)>();
            for (int g = 0; g < other.Genes.Count; g++)
                if (g != skip && candidate(other.Genes[g]) && other.Score(query, g) is var score && !double.IsNaN(score)) scored.Add((g, score));
            compared = scored.Count;
            scored.Sort((a, b) => b.Score != a.Score ? b.Score.CompareTo(a.Score) : a.Gene.CompareTo(b.Gene));
            return scored;
        }

        private static Func<ZCurveGeneProfile.Gene, bool> Comparable(int residues) => residues >= ZCurveGeneProfile.LongResidues
            ? gene => gene.Residues >= ZCurveGeneProfile.LongResidues
            : gene => gene.Residues < ZCurveGeneProfile.LongResidues && Math.Abs(gene.Residues - residues) <= ProteinSearch.LengthBand * residues;

        private static int Best(ZCurveGeneProfile profile, int gene, ZCurveGeneProfile other, int skip)
        {
            var own = profile.Genes[gene];
            var query = profile.WindowsOf(gene);
            var candidate = Comparable(own.Residues);
            int best = -1;
            double top = double.NegativeInfinity;
            for (int g = 0; g < other.Genes.Count; g++)
                if (g != skip && candidate(other.Genes[g]) && other.Score(query, g) is var score && score > top) (top, best) = (score, g);
            return best;
        }
    }
}
