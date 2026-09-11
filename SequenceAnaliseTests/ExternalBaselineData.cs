using System.Globalization;
using System.Runtime.InteropServices;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Parsers for the committed Mash/FastANI raw output under results/raw/, shared by
    /// ExternalBaselineComparisonBenchmark.cs (baseline comparison, §8.7) and
    /// RetrievalClassificationBenchmark.cs (retrieval framing, §8.8) so both read the same source
    /// files through the same code - the same reason GenomeTaxonomy.cs exists rather than each
    /// benchmark carrying its own copy of the labels.
    ///
    /// Mash and FastANI were run once, externally, via Docker (see results/raw/README.md for exact
    /// versions/commands); their raw TSV output is committed as source data, so everything built on it
    /// is regenerable via `dotnet test` without requiring Docker.
    /// </summary>
    internal static class ExternalBaselineData
    {
        private static string FindRepoRoot()
        {
            // GenomeCatalog.DirectoryPath is <repoRoot>/data/genomes - walk up two levels from there.
            var dataDir = Directory.GetParent(GenomeCatalog.DirectoryPath) ?? throw new DirectoryNotFoundException();
            var repoRoot = dataDir.Parent ?? throw new DirectoryNotFoundException();
            return repoRoot.FullName;
        }

        public static string RawResultsPath(string fileName) =>
            Path.Combine(FindRepoRoot(), "results", "raw", fileName);

        /// <summary>Order-independent key for an unordered genome pair.</summary>
        public static string PairKey(string fileA, string fileB) =>
            string.CompareOrdinal(fileA, fileB) <= 0 ? $"{fileA}|{fileB}" : $"{fileB}|{fileA}";

        /// <summary>Reference, Query, Distance, PValue, SharedHashes/SketchSize - mash dist output.</summary>
        public static Dictionary<string, double> ParseMashDistances(string path)
        {
            var result = new Dictionary<string, double>();
            foreach (var line in File.ReadLines(path))
            {
                var cols = line.Split('\t');
                var (fileA, fileB) = (cols[0], cols[1]);
                var distance = double.Parse(cols[2], CultureInfo.InvariantCulture);
                var key = PairKey(fileA, fileB);
                // Mash's Jaccard-derived distance is symmetric by construction; keep the first value
                // seen per unordered pair rather than assuming and skip verifying the reverse direction.
                result.TryAdd(key, distance);
            }
            return result;
        }

        /// <summary>Query, Reference, ANI%, bidirectional fragment mappings, total query fragments.
        /// Only pairs above FastANI's effective identity threshold appear at all - averages both
        /// directions when both are present (ANI is not perfectly symmetric; one direction otherwise).</summary>
        public static Dictionary<string, List<double>> ParseFastAni(string path)
        {
            var result = new Dictionary<string, List<double>>();
            foreach (var line in File.ReadLines(path))
            {
                var cols = line.Split('\t');
                var fileA = Path.GetFileName(cols[0]);
                var fileB = Path.GetFileName(cols[1]);
                if (fileA == fileB) continue; // skip self-pairs (always ANI=100, not a comparison)
                var ani = double.Parse(cols[2], CultureInfo.InvariantCulture);
                var key = PairKey(fileA, fileB);

                // CollectionsMarshal.GetValueRefOrAddDefault: one dictionary probe that returns a
                // mutable ref to the slot, instead of the usual TryGetValue-then-indexer-set pattern
                // (two probes). exists is false only the first time a given pair's key is seen.
                ref var list = ref CollectionsMarshal.GetValueRefOrAddDefault(result, key, out bool exists);
                if (!exists) list = new List<double>();
                list!.Add(ani);
            }
            return result;
        }
    }
}
