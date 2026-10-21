using System.Globalization;
using System.Text.Json;
using SequenceAnalyses;

namespace SequenceAnalise.Benchmarks
{
    public sealed record LocalizationResult(
        string Method,
        int WindowSizeBp,
        int QueryLengthBp,
        string ReferenceFileName,
        string Accession,
        string LocusId,
        string QueryId,
        int TrueStartBp,
        int BestOffsetBp,
        int AbsoluteErrorBp,
        bool ExactHit,
        int TrueHitRank,
        double TrueScore,
        double BestFalseScore,
        double ScoreMargin,
        double BestScore);

    public static class WindowSizeLocalizationBenchmark
    {
        public static LocalizationResult EvaluateLocalization(
            string reference,
            string query,
            int trueStartBp,
            int windowSizeBp,
            string method)
        {
            if (reference is null) throw new ArgumentNullException(nameof(reference));
            if (query is null) throw new ArgumentNullException(nameof(query));
            if (windowSizeBp <= 0) throw new ArgumentOutOfRangeException(nameof(windowSizeBp));
            if (trueStartBp < 0) throw new ArgumentOutOfRangeException(nameof(trueStartBp));
            if (trueStartBp + query.Length > reference.Length)
                throw new ArgumentException("The query extends beyond the provided reference sequence.");

            var trueOffset = trueStartBp / windowSizeBp;
            var ranked = ScoreAllOffsets(reference, query, windowSizeBp, method)
                .OrderByDescending(r => r.Score)
                .ThenBy(r => r.Offset)
                .ToList();

            if (ranked.Count == 0)
                throw new InvalidOperationException("No valid window offsets were evaluated.");

            var best = ranked[0];
            var bestOffset = best.Offset * windowSizeBp;
            var exactHit = bestOffset == trueStartBp;
            var trueHitRank = ranked.FindIndex(r => r.Offset == trueOffset);
            var trueScore = ranked[trueHitRank >= 0 ? trueHitRank : 0].Score;
            var bestFalseScore = ranked
                .Where(r => r.Offset != trueOffset)
                .DefaultIfEmpty(new RankedOffset(best.Offset, double.NegativeInfinity))
                .MaxBy(r => r.Score)
                .Score;
            if (double.IsNegativeInfinity(bestFalseScore))
                bestFalseScore = trueScore;

            return new LocalizationResult(
                Method: method,
                WindowSizeBp: windowSizeBp,
                QueryLengthBp: query.Length,
                ReferenceFileName: string.Empty,
                Accession: string.Empty,
                LocusId: string.Empty,
                QueryId: string.Empty,
                TrueStartBp: trueStartBp,
                BestOffsetBp: bestOffset,
                AbsoluteErrorBp: Math.Abs(bestOffset - trueStartBp),
                ExactHit: exactHit,
                TrueHitRank: trueHitRank >= 0 ? trueHitRank + 1 : 0,
                TrueScore: trueScore,
                BestFalseScore: bestFalseScore,
                ScoreMargin: trueScore - bestFalseScore,
                BestScore: best.Score);
        }

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

            var outputDirectory = Path.GetDirectoryName(Path.GetFullPath(outputPath));
            if (!string.IsNullOrEmpty(outputDirectory)) Directory.CreateDirectory(outputDirectory);

            using var writer = new StreamWriter(outputPath, false, System.Text.Encoding.UTF8);
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
                    var result = EvaluateLocalization(record.Sequence, query, locus.StartBp, windowSize, method);
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

        private static List<RankedOffset> ScoreAllOffsets(string reference, string query, int windowSizeBp, string method)
        {
            var queryLength = query.Length;
            var windowCount = queryLength / windowSizeBp;
            if (windowCount == 0)
                return new List<RankedOffset>();

            var referenceWindowCount = reference.Length / windowSizeBp;
            var candidateOffsets = referenceWindowCount - windowCount + 1;
            var results = new List<RankedOffset>(candidateOffsets);

            for (int offset = 0; offset < candidateOffsets; offset++)
            {
                var sliceStart = offset * windowSizeBp;
                var slice = reference.Substring(sliceStart, queryLength);
                var score = ScoreForSlice(query, slice, windowSizeBp, method);
                results.Add(new RankedOffset(offset, score));
            }

            return results;
        }

        private static double ScoreForSlice(string query, string referenceSlice, int windowSizeBp, string method)
        {
            return method switch
            {
                "ZCurve" => CompareZCurve(query, referenceSlice, windowSizeBp),
                "WindowedKmer" => CompareWindowedKMer(query, referenceSlice, windowSizeBp),
                _ => throw new ArgumentOutOfRangeException(nameof(method), method, "Unsupported method name.")
            };
        }

        private static double CompareWindowedKMer(string query, string referenceSlice, int windowSizeBp)
        {
            var q = new WindowedProfile(query, windowSizeBp, 4);
            var r = new WindowedProfile(referenceSlice, windowSizeBp, 4);
            return WindowedProfile.Compare(q, r);
        }

        private static double CompareZCurve(string query, string referenceSlice, int windowSizeBp)
        {
            var q = new ZCurveWindowedProfile(query, windowSizeBp);
            var r = new ZCurveWindowedProfile(referenceSlice, windowSizeBp);
            return ZCurveWindowedProfile.Compare(q, r);
        }

        private static string CsvEscape(string value)
        {
            if (value is null) return string.Empty;
            if (value.Contains(',') || value.Contains('"') || value.Contains('\n') || value.Contains('\r'))
                return '"' + value.Replace("\"", "\"\"") + '"';
            return value;
        }

        private readonly record struct RankedOffset(int Offset, double Score);
    }
}
