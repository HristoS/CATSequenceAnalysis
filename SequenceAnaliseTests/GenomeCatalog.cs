namespace SequenceAnaliseTests
{
    /// <summary>
    /// Locates the repository's data/genomes/ directory (30 real RefSeq bacterial genomes - see
    /// data/genomes/README.md) relative to wherever the test assembly happens to be running from,
    /// so tests don't hardcode a path tied to Debug/Release or a specific target framework moniker.
    /// </summary>
    internal static class GenomeCatalog
    {
        public static string DirectoryPath { get; } = FindGenomesDirectory();

        private static string FindGenomesDirectory()
        {
            for (var dir = new DirectoryInfo(AppContext.BaseDirectory); dir is not null; dir = dir.Parent)
            {
                var candidate = Path.Combine(dir.FullName, "data", "genomes");
                if (Directory.Exists(candidate) && Directory.EnumerateFiles(candidate, "*.fna.gz").Any())
                    return candidate;
            }

            throw new DirectoryNotFoundException(
                $"Could not find a populated 'data/genomes' directory above '{AppContext.BaseDirectory}'. " +
                "Expected to find it by walking up to the repository root - see data/genomes/README.md.");
        }

        /// <summary>Full path to one genome file, by the filename listed in data/genomes/README.md.</summary>
        public static string PathTo(string fileName) => Path.Combine(DirectoryPath, fileName);

        /// <summary>Every downloaded genome file's path, for tests that want to sweep all of them.</summary>
        public static IEnumerable<string> AllFiles =>
            Directory.EnumerateFiles(DirectoryPath, "*.fna.gz").OrderBy(f => f);
    }
}
