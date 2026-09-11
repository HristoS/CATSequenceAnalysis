namespace SequenceAnaliseTests
{
    internal static class AccuracyBenchmarkShared
    {
        public const int SeedCount = 5;

        public static readonly double[] MutationRates = { 0.00, 0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.50 };

        // Mutate a fraction of bases deterministically, to simulate a diverged homolog. Same
        // algorithm already duplicated in ZCurveProfileTests.cs / ZCurveWindowedProfileTests.cs;
        // shared here so the small-scale (NW-anchored) and large-scale (known-by-construction)
        // accuracy tables are provably driven by the identical mutation procedure - important
        // since both tables are presented side by side as one reliability argument.
        public static string Mutate(string s, double rate, int seed)
        {
            var rng = new Random(seed);
            char[] bases = { 'A', 'C', 'G', 'T' };
            var chars = s.ToCharArray();
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < rate) chars[i] = bases[rng.Next(bases.Length)];
            return new string(chars);
        }
    }
}
