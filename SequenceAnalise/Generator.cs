using System;

namespace SequenceAnalyses
{
    public static class Generator
    {
        // Index matches Base's declared order (A, T, C, G), so seeded output is stable.
        private static readonly char[] Bases = { 'A', 'T', 'C', 'G' };

        public static string GetRandomDNA(int length) => GetRandomDNA(length, Random.Shared);

        public static string GetRandomDNA(int length, int seed) => GetRandomDNA(length, new Random(seed));

        public static string GetRandomDNA(int length, Random random)
        {
            var dna = new char[length];
            for (int i = 0; i < length; i++)
            {
                dna[i] = Bases[random.Next(Bases.Length)];
            }

            return new string(dna);
        }
    }
}
