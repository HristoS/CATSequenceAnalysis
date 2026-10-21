using System.Globalization;

namespace SequenceAnalise.Benchmarks
{
    public static class BenchmarkSummary
    {
        public static void WriteDivergenceSummary(string inputCsvPath, string outputCsvPath)
        {
            var rows = ReadRows(inputCsvPath);
            var grouped = rows
                .GroupBy(r => new { r.Method, r.DivergenceRate, r.WindowSizeBp, r.QueryLengthBp })
                .OrderBy(g => g.Key.Method, StringComparer.Ordinal)
                .ThenBy(g => g.Key.DivergenceRate)
                .ThenBy(g => g.Key.WindowSizeBp)
                .ThenBy(g => g.Key.QueryLengthBp)
                .ToList();

            var lines = new List<string>
            {
                "Method,DivergenceRate,WindowSizeBp,QueryLengthBp,Trials,Top1Accuracy,Top5Accuracy,MedianAbsoluteErrorBp,MeanAbsoluteErrorBp,MeanScoreMargin,MeanBestScore"
            };

            foreach (var group in grouped)
            {
                var exactCount = group.Count(r => r.ExactHit);
                var top5Count = group.Count(r => r.TrueHitRank <= 5);
                var absoluteErrors = group.Select(r => r.AbsoluteErrorBp).OrderBy(v => v).ToList();
                var medianAbsoluteError = Median(absoluteErrors);
                var meanAbsoluteError = absoluteErrors.Count == 0 ? 0.0 : absoluteErrors.Average();
                var meanScoreMargin = group.Select(r => r.ScoreMargin).DefaultIfEmpty(0.0).Average();
                var meanBestScore = group.Select(r => r.BestScore).DefaultIfEmpty(0.0).Average();

                lines.Add(string.Join(",",
                    group.Key.Method,
                    group.Key.DivergenceRate.ToString("G17", CultureInfo.InvariantCulture),
                    group.Key.WindowSizeBp,
                    group.Key.QueryLengthBp,
                    group.Count(),
                    (exactCount / (double)group.Count()).ToString("0.######", CultureInfo.InvariantCulture),
                    (top5Count / (double)group.Count()).ToString("0.######", CultureInfo.InvariantCulture),
                    medianAbsoluteError.ToString("0.######", CultureInfo.InvariantCulture),
                    meanAbsoluteError.ToString("0.######", CultureInfo.InvariantCulture),
                    meanScoreMargin.ToString("0.######", CultureInfo.InvariantCulture),
                    meanBestScore.ToString("0.######", CultureInfo.InvariantCulture)));
            }

            var outputDirectory = Path.GetDirectoryName(Path.GetFullPath(outputCsvPath));
            if (!string.IsNullOrEmpty(outputDirectory)) Directory.CreateDirectory(outputDirectory);
            File.WriteAllLines(outputCsvPath, lines);
        }

        public static void WriteWindowSizeSummary(string inputCsvPath, string outputCsvPath)
        {
            var rows = ReadRows(inputCsvPath);
            var grouped = rows
                .GroupBy(r => new { r.Method, r.WindowSizeBp, r.QueryLengthBp })
                .OrderBy(g => g.Key.Method, StringComparer.Ordinal)
                .ThenBy(g => g.Key.WindowSizeBp)
                .ThenBy(g => g.Key.QueryLengthBp)
                .ToList();

            var lines = new List<string>
            {
                "Method,WindowSizeBp,QueryLengthBp,Trials,Top1Accuracy,Top5Accuracy,MedianAbsoluteErrorBp,MeanAbsoluteErrorBp,MeanScoreMargin,MeanBestScore"
            };

            foreach (var group in grouped)
            {
                var exactCount = group.Count(r => r.ExactHit);
                var top5Count = group.Count(r => r.TrueHitRank <= 5);
                var absoluteErrors = group.Select(r => r.AbsoluteErrorBp).OrderBy(v => v).ToList();
                var medianAbsoluteError = Median(absoluteErrors);
                var meanAbsoluteError = absoluteErrors.Count == 0 ? 0.0 : absoluteErrors.Average();
                var meanScoreMargin = group.Select(r => r.ScoreMargin).DefaultIfEmpty(0.0).Average();
                var meanBestScore = group.Select(r => r.BestScore).DefaultIfEmpty(0.0).Average();

                lines.Add(string.Join(",",
                    group.Key.Method,
                    group.Key.WindowSizeBp,
                    group.Key.QueryLengthBp,
                    group.Count(),
                    (exactCount / (double)group.Count()).ToString("0.######", CultureInfo.InvariantCulture),
                    (top5Count / (double)group.Count()).ToString("0.######", CultureInfo.InvariantCulture),
                    medianAbsoluteError.ToString("0.######", CultureInfo.InvariantCulture),
                    meanAbsoluteError.ToString("0.######", CultureInfo.InvariantCulture),
                    meanScoreMargin.ToString("0.######", CultureInfo.InvariantCulture),
                    meanBestScore.ToString("0.######", CultureInfo.InvariantCulture)));
            }

            var outputDirectory = Path.GetDirectoryName(Path.GetFullPath(outputCsvPath));
            if (!string.IsNullOrEmpty(outputDirectory)) Directory.CreateDirectory(outputDirectory);
            File.WriteAllLines(outputCsvPath, lines);
        }

        private static List<BenchmarkRow> ReadRows(string csvPath)
        {
            if (!File.Exists(csvPath)) throw new FileNotFoundException($"Benchmark CSV not found: {csvPath}", csvPath);

            var lines = File.ReadAllLines(csvPath);
            if (lines.Length < 2) return new List<BenchmarkRow>();

            var header = ParseCsvLine(lines[0]);
            var index = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
            for (int i = 0; i < header.Length; i++) index[header[i]] = i;

            var rows = new List<BenchmarkRow>(lines.Length - 1);
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var values = ParseCsvLine(lines[i]);
                rows.Add(new BenchmarkRow(
                    Get(values, index, "Method"),
                    ParseDouble(Get(values, index, "DivergenceRate"), 0.0),
                    ParseInt(Get(values, index, "WindowSizeBp"), 0),
                    ParseInt(Get(values, index, "QueryLengthBp"), 0),
                    ParseInt(Get(values, index, "TrueHitRank"), 0),
                    ParseBool(Get(values, index, "ExactHit")),
                    ParseDouble(Get(values, index, "ScoreMargin"), 0.0),
                    ParseDouble(Get(values, index, "BestScore"), 0.0),
                    ParseInt(Get(values, index, "AbsoluteErrorBp"), 0)));
            }

            return rows;
        }

        private static string Get(string[] values, Dictionary<string, int> index, string fieldName)
        {
            if (!index.TryGetValue(fieldName, out var position)) return string.Empty;
            return position < values.Length ? values[position] : string.Empty;
        }

        private static string[] ParseCsvLine(string line)
        {
            var values = new List<string>();
            var current = new System.Text.StringBuilder();
            var inQuotes = false;
            for (int i = 0; i < line.Length; i++)
            {
                var ch = line[i];
                if (ch == '"')
                {
                    if (inQuotes && i + 1 < line.Length && line[i + 1] == '"')
                    {
                        current.Append('"');
                        i++;
                    }
                    else
                    {
                        inQuotes = !inQuotes;
                    }
                }
                else if (ch == ',' && !inQuotes)
                {
                    values.Add(current.ToString());
                    current.Clear();
                }
                else
                {
                    current.Append(ch);
                }
            }

            values.Add(current.ToString());
            return values.ToArray();
        }

        private static double Median(IReadOnlyList<int> values)
        {
            if (values.Count == 0) return 0.0;
            var ordered = values.OrderBy(v => v).ToList();
            int mid = ordered.Count / 2;
            if (ordered.Count % 2 == 0)
                return (ordered[mid - 1] + ordered[mid]) / 2.0;
            return ordered[mid];
        }

        private static bool ParseBool(string value) => bool.TryParse(value, out var b) ? b : value.Equals("true", StringComparison.OrdinalIgnoreCase);

        private static double ParseDouble(string value, double defaultValue)
            => double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out var parsed) ? parsed : defaultValue;

        private static int ParseInt(string value, int defaultValue)
            => int.TryParse(value, NumberStyles.Integer, CultureInfo.InvariantCulture, out var parsed) ? parsed : defaultValue;

        private sealed record BenchmarkRow(
            string Method,
            double DivergenceRate,
            int WindowSizeBp,
            int QueryLengthBp,
            int TrueHitRank,
            bool ExactHit,
            double ScoreMargin,
            double BestScore,
            int AbsoluteErrorBp);
    }
}
