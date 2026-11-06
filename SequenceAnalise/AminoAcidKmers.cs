using System;
using System.Collections.Generic;

namespace SequenceAnalyses
{
    /// <summary>
    /// Amino-acid k-mer counts of a genome's proteins, the input of Hao's composition vector (CVTree). On catalogue genomes the
    /// Markov-corrected 4-mers track GTDB's tree far better than dinucleotides: within-family Spearman 0.842 against Karlin's 0.586
    /// (results/T2_CvTreeGtdb.csv). Words never cross from one protein to the next, and a letter outside <see cref="Alphabet"/>
    /// (a stop '*', an unknown 'X') breaks a word.
    /// </summary>
    public static class AminoAcidKmers
    {
        public const string Alphabet = "ACDEFGHIKLMNPQRSTVWY";

        /// <summary>20^k counts, index = the word read as a base-20 number over <see cref="Alphabet"/>, first residue most significant.</summary>
        public static int[] Count(IEnumerable<string> proteins, int k)
        {
            if (k < 1 || k > 6) throw new ArgumentOutOfRangeException(nameof(k));
            int space = 1;
            for (int i = 0; i < k; i++)
                space *= Alphabet.Length;
            var lookup = new int[128];
            Array.Fill(lookup, -1);
            for (int i = 0; i < Alphabet.Length; i++)
                lookup[Alphabet[i]] = i;

            var counts = new int[space];
            foreach (var protein in proteins)
            {
                int word = 0, valid = 0;
                foreach (var residue in protein)
                {
                    int letter = residue < 128 ? lookup[char.ToUpperInvariant(residue)] : -1;
                    if (letter < 0)
                    {
                        valid = 0;
                        continue;
                    }
                    word = word * Alphabet.Length % space + letter;
                    if (++valid >= k) counts[word]++;
                }
            }
            return counts;
        }
    }
}
