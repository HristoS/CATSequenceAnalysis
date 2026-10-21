using System.Globalization;
using System.Text.Json;
using SequenceAnalyses;

namespace SequenceAnalise.Benchmarks
{
    public static class DivergenceLocalizationBenchmark
    {
        public static void Run(
            string genomesDirectory,
            string manifestPath,
            string outputPath,
            int windowSizeBp,
            int[] queryLengths,
            double[] divergenceRates)
        {
            if (windowSizeBp <= 0) throw new ArgumentOutOfRangeException(nameof(windowSizeBp));
            if (queryLengths is null || queryLengths.Length == 0)
                throw new ArgumentException("At least one query length is required.", nameof(queryLengths));
            if (divergenceRates is null || divergenceRates.Length == 0)
                throw new ArgumentException("At least one divergence rate is required.", nameof(divergenceRates));

            var manifest = JsonSerializer.Deserialize<LocalizationDatasetManifest>(File.ReadAllText(manifestPath))
                ?? throw new InvalidDataException("Manifest is empty or invalid JSON.");

            var fullOutputPath = Path.GetFullPath(outputPath);
            var outputDirectory = Path.GetDirectoryName(fullOutputPath);
            if (!string.IsNullOrEmpty(outputDirectory)) Directory.CreateDirectory(outputDirectory);

            using var writer = new StreamWriter(fullOutputPath, false, System.Text.Encoding.UTF8);
            writer.WriteLine(
                "ReferenceFileName,Accession,LocusId,QueryId,Method,DivergenceRate,WindowSizeBp,QueryLengthBp,TrueStartBp,BestOffsetBp,AbsoluteErrorBp,ExactHit,TrueHitRank,TrueScore,BestFalseScore,ScoreMargin,BestScore");

            foreach (var reference in manifest.References)
            {
                var referencePath = Path.Combine(genomesDirectory, reference.FileName);
                var record = FastaReader.ReadPrimaryRecord(referencePath);

                foreach (var locus in reference.Loci)
                foreach (var queryRecipe in locus.Queries.Where(q => queryLengths.Contains(q.LengthBp)))
                {
                    var originalQuery = record.Sequence.Substring(locus.StartBp, queryRecipe.LengthBp);

                    foreach (var divergenceRate in divergenceRates)
                    foreach (var method in new[] { "ZCurve", "WindowedKmer" })
                    {
                        var mutatedQuery = MutateSequence(originalQuery, divergenceRate, DeriveSeed(reference.FileName, locus.LocusId, queryRecipe.QueryId, divergenceRate));
                        var result = WindowSizeLocalizationBenchmark.EvaluateLocalization(record.Sequence, mutatedQuery, locus.StartBp, windowSizeBp, method);
                        writer.WriteLine(string.Join(",",
                            CsvEscape(reference.FileName),
                            CsvEscape(reference.Accession),
                            CsvEscape(locus.LocusId),
                            CsvEscape(queryRecipe.QueryId),
                            CsvEscape(result.Method),
                            divergenceRate.ToString("G17", CultureInfo.InvariantCulture),
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
        }

        private static string MutateSequence(string query, double divergenceRate, int seed)
        {
            if (divergenceRate < 0.0 || divergenceRate > 1.0)
                throw new ArgumentOutOfRangeException(nameof(divergenceRate), divergenceRate, "Divergence rate must be between 0.0 and 1.0.");

            if (divergenceRate == 0.0)
                return query;

            var buffer = query.ToCharArray();
            var random = new Random(seed);
            foreach (var index in Enumerable.Range(0, buffer.Length))
            {
                if (random.NextDouble() < divergenceRate)
                {
                    var choices = "ACGT".Where(c => c != buffer[index]).ToArray();
                    buffer[index] = choices[random.Next(choices.Length)];
                }
            }

            return new string(buffer);
        }

        private static int DeriveSeed(string fileName, string locusId, string queryId, double divergenceRate)
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 31 + fileName.GetHashCode(StringComparison.Ordinal);
                hash = hash * 31 + locusId.GetHashCode(StringComparison.Ordinal);
                hash = hash * 31 + queryId.GetHashCode(StringComparison.Ordinal);
                hash = hash * 31 + BitConverter.SingleToInt32Bits((float)divergenceRate);
                return hash;
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
