using System.Globalization;
using System.Text.Json;
using SequenceAnalyses;

namespace SequenceAnalise.Benchmarks
{
    public static class WindowSizeRobustnessBenchmark
    {
        public static void Run(
            string genomesDirectory,
            string manifestPath,
            string outputPath,
            int[] windowSizes,
            int[] queryLengths)
        {
            if (windowSizes is null || windowSizes.Length == 0)
                throw new ArgumentException("At least one window size is required.", nameof(windowSizes));
            if (queryLengths is null || queryLengths.Length == 0)
                throw new ArgumentException("At least one query length is required.", nameof(queryLengths));

            var manifest = JsonSerializer.Deserialize<LocalizationDatasetManifest>(File.ReadAllText(manifestPath))
                ?? throw new InvalidDataException("Manifest is empty or invalid JSON.");

            var fullOutputPath = Path.GetFullPath(outputPath);
            var outputDirectory = Path.GetDirectoryName(fullOutputPath);
            if (!string.IsNullOrEmpty(outputDirectory)) Directory.CreateDirectory(outputDirectory);

            using var writer = new StreamWriter(fullOutputPath, false, System.Text.Encoding.UTF8);
            writer.WriteLine(
                "ReferenceFileName,Accession,LocusId,QueryId,Method,WindowSizeBp,QueryLengthBp,TrueStartBp,BestOffsetBp,AbsoluteErrorBp,ExactHit,TrueHitRank,TrueScore,BestFalseScore,ScoreMargin,BestScore");

            foreach (var reference in manifest.References)
            {
                var referencePath = Path.Combine(genomesDirectory, reference.FileName);
                var record = FastaReader.ReadPrimaryRecord(referencePath);

                foreach (var windowSize in windowSizes)
                foreach (var method in new[] { "ZCurve", "WindowedKmer" })
                foreach (var locus in reference.Loci)
                foreach (var queryRecipe in locus.Queries.Where(q => queryLengths.Contains(q.LengthBp)))
                {
                    var query = record.Sequence.Substring(locus.StartBp, queryRecipe.LengthBp);
                    var result = WindowSizeLocalizationBenchmark.EvaluateLocalization(record.Sequence, query, locus.StartBp, windowSize, method);

                    writer.WriteLine(string.Join(",",
                        CsvEscape(reference.FileName),
                        CsvEscape(reference.Accession),
                        CsvEscape(locus.LocusId),
                        CsvEscape(queryRecipe.QueryId),
                        CsvEscape(result.Method),
                        result.WindowSizeBp,
                        result.QueryLengthBp,
                        result.TrueStartBp,
                        result.BestOffsetBp,
                        result.AbsoluteErrorBp,
                        result.ExactHit.ToString().ToLowerInvariant(),
                        result.TrueHitRank,
                        result.TrueScore.ToString("G17", CultureInfo.InvariantCulture),
                        result.BestFalseScore.ToString("G17", CultureInfo.InvariantCulture),
                        result.ScoreMargin.ToString("G17", CultureInfo.InvariantCulture),
                        result.BestScore.ToString("G17", CultureInfo.InvariantCulture)));
                }
            }
        }

        private static string CsvEscape(string value)
        {
            if (value is null) return string.Empty;
            if (value.Contains(',') || value.Contains('"') || value.Contains('\n') || value.Contains('\r'))
                return '"' + value.Replace("\"", "\"\"") + '"';
            return value;
        }
    }
}
