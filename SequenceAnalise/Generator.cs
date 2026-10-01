using System;

namespace SequenceAnalyses
{
    public static class Generator
    {
        // Order is fixed (A, T, C, G) - changing it would change every seeded sequence.
        private static readonly char[] bases = { 'A', 'T', 'C', 'G' };

        public static string GetRandomDNA(int length) => GetRandomDNA(length, Random.Shared);

        public static string GetRandomDNA(int length, int seed) => GetRandomDNA(length, new Random(seed));

        public static string GetRandomDNA(int length, Random random)
        {
            var dna = new char[length];
            for (int i = 0; i < length; i++)
            {
                dna[i] = bases[random.Next(bases.Length)];
            }

            return new string(dna);
        }
    }
}
