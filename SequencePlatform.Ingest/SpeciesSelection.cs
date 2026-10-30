namespace SequencePlatform.Ingest
{
    /// <summary>
    /// What --per-species adds: one genome for each species the catalogue doesn't hold. BV-BRC lists thousands of near-identical strains of
    /// a few species, so taking everything would grow the catalogue without widening it. The order is shuffled with a fixed seed, so a
    /// run with --limit is a fair sample rather than the start of the alphabet.
    /// </summary>
    internal static class SpeciesSelection
    {
        private const int seed = 20261030;
        private static readonly string[] referenceClasses = { "Reference", "Representative" };

        public static List<string> Missing(IEnumerable<string> species, IReadOnlySet<string> stored, int limit)
        {
            var rng = new Random(seed);
            return species.Where(s => !stored.Contains(s)).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal)
                .OrderBy(_ => rng.Next()).Take(limit).ToList();
        }

        /// <summary>Per missing species a Reference or Representative genome if there is one, else the lowest genome id.</summary>
        public static List<BvBrcGenome> PickNew(IEnumerable<BvBrcGenome> genomes, IReadOnlySet<string> stored, int limit)
        {
            var bySpecies = genomes.Where(g => !string.IsNullOrEmpty(g.Field("species")))
                .GroupBy(g => g.Field("species")!, StringComparer.Ordinal)
                .ToDictionary(group => group.Key, group => group
                    .OrderBy(g => referenceClasses.Contains(g.Field("reference_genome")) ? 0 : 1)
                    .ThenBy(g => g.GenomeId, StringComparer.Ordinal)
                    .First(), StringComparer.Ordinal);
            return Missing(bySpecies.Keys, stored, limit).Select(species => bySpecies[species]).ToList();
        }
    }
}
