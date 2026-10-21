using System.Diagnostics;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Every profile here is built once at the fine 500 bp window and every coarser resolution is
    /// derived with <see cref="ZCurveWindowedProfile.Downsample"/> - exact up to float rounding
    /// (ZCurveWindowedProfileTests.DownsampleMatchesProfileBuiltAtTheCoarserWindowSize) - which is
    /// the point being tested: a platform could store only the fine profile and derive the rest.
    ///
    /// Six questions, one test each, on the 30-genome bacterial corpus of
    /// RealGenomeDistanceMatrixBenchmark/RetrievalClassificationBenchmark and
    /// the 31-genome viral corpus of ViralTaxonomy (WindowCountSweep and LadderRuleAgainstFineOnly):
    ///  - ResolutionSweep: how far can window size grow before relationship ranking (retrieval AUC)
    ///    degrades? WindowSizeSensitivityBenchmark only went to 2,000 bp.
    ///  - CoarseToFineLocalization: does searching a coarse level first, then refining at 500 bp
    ///    near the coarse hit, find the same offset as the full fine search, and how much faster?
    ///  - WindowCountSweep: ranking depends on how many windows the *shorter* sequence of a pair
    ///    keeps, not on window size itself - how few can it get before retrieval AUC breaks? Each
    ///    corpus is pushed to coarser and coarser levels so the same pairs are measured throughout
    ///    (bucketing mixed pairs by window count would confound it with which pairs land where).
    ///  - LadderRuleAgainstFineOnly: does ZCurveProfileLadder.Compare (coarsest level keeping the
    ///    shorter sequence at >= 8 windows) keep the fine-only AUC on both corpora, and how much faster?
    ///  - PerQueryLevelAgainstPerPair: does one level per query (what a search must use, so a ranking
    ///    never mixes levels) rank as well as LadderRuleAgainstFineOnly's one level per pair?
    ///  - GcAxisWeightSweep: does weighting the GC (z) axis differently from the other two - the
    ///    only thing a cylindrical frame around z could change about the metric - help ranking?
    ///
    /// An int16-sum storage prototype (6 bytes per window instead of 12, AVX2 integer distance) was
    /// also measured here and removed: scores matched to 4.5e-7 but comparison was 0.94x - slower
    /// than the float path - so it bought storage only. See results/README.md, 2026-09-26.
    /// </summary>
    public class MultiResolutionBenchmark
    {
        private const int fineWindowSize = 500;
        private const int localizationTrials = 10;
        private const int timingRepetitions = 3;

        private static readonly int[] downsampleFactors = { 1, 2, 4, 8, 16, 32, 64, 128 };
        private static readonly int[] fragmentLengths = { 2_000, 10_000, 50_000 };
        private static readonly int[] coarseFactors = { 4, 8, 16 };
        private static readonly double[] gcAxisWeights = { 0, 0.5, 1, 2, 4, 16 };
        private static readonly int[] bacterialCountSweepFactors = { 1, 4, 16, 64, 256, 1024, 2048, 4096 };
        private static readonly int[] viralCountSweepFactors = { 1, 2, 4, 8, 16, 32 };

        private static readonly RelationshipLevel[] levelsInOrder =
        {
            RelationshipLevel.SameSpecies, RelationshipLevel.SameGenus,
            RelationshipLevel.SamePhylum, RelationshipLevel.DifferentPhylum,
        };

        // Read and profiled once per test run; each genome's primary record only.
        private static List<(GenomeMeta Meta, string Sequence, ZCurveWindowedProfile Profile)>? corpus;

        private static List<(GenomeMeta Meta, string Sequence, ZCurveWindowedProfile Profile)> Corpus =>
            corpus ??= GenomeTaxonomy.Genomes.Select(g =>
            {
                var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence;
                return (g, sequence, new ZCurveWindowedProfile(sequence, fineWindowSize));
            }).ToList();

        private static List<(ViralGenomeMeta Meta, ZCurveWindowedProfile Profile)>? viralCorpus;

        private static List<(ViralGenomeMeta Meta, ZCurveWindowedProfile Profile)> ViralCorpus =>
            viralCorpus ??= ViralTaxonomy.Genomes.Select(g => (g, new ZCurveWindowedProfile(
                FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(g.FileName)).Sequence, fineWindowSize))).ToList();

        private static List<(RelationshipLevel Level, double Score)> ScoreAllPairs(
            IReadOnlyList<(GenomeMeta Meta, ZCurveWindowedProfile Profile)> profiles)
        {
            var rows = new List<(RelationshipLevel Level, double Score)>();
            for (int a = 0; a < profiles.Count; a++)
                for (int b = a + 1; b < profiles.Count; b++)
                    rows.Add((GenomeTaxonomy.Classify(profiles[a].Meta, profiles[b].Meta),
                        ZCurveWindowedProfile.Compare(profiles[a].Profile, profiles[b].Profile)));
            return rows;
        }

        // Same two thresholds as RetrievalClassificationBenchmark, so AUCs are directly comparable.
        private static (double SameGenusOrCloser, double SamePhylumOrCloser) Aucs(List<(RelationshipLevel Level, double Score)> rows)
        {
            double Auc(Func<RelationshipLevel, bool> isPositive) => RetrievalClassificationBenchmark.ComputeAuc(
                rows.Where(r => isPositive(r.Level)).Select(r => r.Score).ToList(),
                rows.Where(r => !isPositive(r.Level)).Select(r => r.Score).ToList());

            return (
                Auc(l => l is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus),
                Auc(l => l is not RelationshipLevel.DifferentPhylum));
        }

        private static (double Mean, double StdDev) MeanAndStdDev(IEnumerable<double> values)
        {
            var v = values.ToList();
            var mean = v.Average();
            var stdDev = v.Count < 2 ? 0 : Math.Sqrt(v.Sum(x => (x - mean) * (x - mean)) / (v.Count - 1));
            return (mean, stdDev);
        }

        [Test]
        public void ResolutionSweep()
        {
            var fine = Corpus.Select(c => (c.Meta, c.Profile)).ToList();

            using var report = File.CreateText("MultiResolution_Sweep.csv");
            report.WriteLine("WindowSize, Factor, " +
                string.Join(", ", levelsInOrder.Select(l => $"{l} Mean, {l} StdDev")) +
                ", AucSameGenusOrCloser, AucSamePhylumOrCloser, AllPairsCompareMs");

            foreach (var factor in downsampleFactors)
            {
                var profiles = fine.Select(f => (f.Meta, Profile: f.Profile.Downsample(factor))).ToList();

                var stopwatch = Stopwatch.StartNew();
                var rows = ScoreAllPairs(profiles);
                stopwatch.Stop();

                var levelStats = levelsInOrder.Select(l => MeanAndStdDev(rows.Where(r => r.Level == l).Select(r => r.Score))).ToList();
                var (aucGenus, aucPhylum) = Aucs(rows);

                report.WriteLine($"{fineWindowSize * factor}, {factor}, " +
                    string.Join(", ", levelStats.Select(s => $"{s.Mean}, {s.StdDev}")) +
                    $", {aucGenus}, {aucPhylum}, {stopwatch.Elapsed.TotalMilliseconds}");
                report.Flush();

                TestContext.Out.WriteLine(
                    $"w={fineWindowSize * factor,6}  AUC genus={aucGenus:F4} phylum={aucPhylum:F4}  " +
                    $"compare={stopwatch.Elapsed.TotalMilliseconds,9:F0} ms");
            }

            Assert.Pass();
        }

        [Test]
        public void CoarseToFineLocalization()
        {
            using var report = File.CreateText("MultiResolution_CoarseToFine.csv");
            report.WriteLine("FragmentLength, Trial, Genome, FragmentStart, ExpectedWindow, Factor, " +
                "FineOffset, FineScore, CoarseToFineOffset, CoarseToFineScore, SameOffset, FineMs, CoarseToFineMs");

            foreach (var fragmentLength in fragmentLengths)
                for (int trial = 0; trial < localizationTrials; trial++)
                {
                    // Random base offset, deliberately not window-aligned.
                    var rng = new Random(BenchmarkTiming.Seed("coarse-to-fine", fragmentLength, trial));
                    var entry = Corpus[rng.Next(Corpus.Count)];
                    int start = rng.Next(0, entry.Sequence.Length - fragmentLength);

                    var query = new ZCurveWindowedProfile(entry.Sequence.Substring(start, fragmentLength), fineWindowSize);
                    var reference = entry.Profile;
                    var referenceWindows = reference.Windows;

                    var (fineScore, fineOffset) = ZCurveWindowedProfile.BestAlignment(query, reference);
                    var fineMs = BenchmarkTiming.Measure(timingRepetitions, () => ZCurveWindowedProfile.BestAlignment(query, reference)).Average();

                    foreach (var factor in coarseFactors)
                    {
                        // Coarse reference is what a platform would store alongside the fine one,
                        // so building it is not part of the timed search.
                        var coarseReference = reference.Downsample(factor);
                        var result = CoarseToFine(query, coarseReference, referenceWindows, factor);
                        var prefix = $"{fragmentLength}, {trial}, {entry.Meta.FileName}, {start}, {(double)start / fineWindowSize}, {factor}, {fineOffset}, {fineScore}";

                        if (result is null)
                        {
                            // Fragment too short to have a single coarse window at this factor.
                            report.WriteLine($"{prefix}, , , NotApplicable, {fineMs}, ");
                            continue;
                        }

                        var coarseMs = BenchmarkTiming.Measure(timingRepetitions, () => CoarseToFine(query, coarseReference, referenceWindows, factor)).Average();
                        var (score, offset) = result.Value;
                        report.WriteLine($"{prefix}, {offset}, {score}, {offset == fineOffset}, {fineMs}, {coarseMs}");
                    }
                    report.Flush();
                }

            Assert.Pass();
        }

        // Coarse pass: one BestAlignment per query phase (so a fragment starting mid-group still
        // lines up with some coarse grid), each at ~1/factor^2 of the fine cost. Refine pass: fine
        // BestAlignment over only a +-factor-window neighbourhood of each coarse hit. Heuristic - the
        // coarse best can differ from the fine best, which is exactly what the report measures.
        private static (double Score, int Offset)? CoarseToFine(
            ZCurveWindowedProfile query,
            ZCurveWindowedProfile coarseReference,
            IReadOnlyList<(double X, double Y, double Z)> referenceWindows,
            int factor)
        {
            (double Score, int Offset)? best = null;
            for (int phase = 0; phase < factor; phase++)
            {
                var coarseQuery = query.Downsample(factor, phase);
                if (coarseQuery.WindowCount == 0) continue;

                var (_, coarseOffset) = ZCurveWindowedProfile.BestAlignment(coarseQuery, coarseReference);
                if (coarseOffset < 0) continue;

                // Query fine window `phase` sits at reference fine window coarseOffset * factor.
                int estimate = coarseOffset * factor - phase;
                int from = Math.Max(0, estimate - factor);
                int to = Math.Min(referenceWindows.Count, estimate + factor + query.WindowCount);
                if (to - from < query.WindowCount) continue;

                var neighbourhood = ZCurveWindowedProfile.FromPrecomputedWindows(
                    query.WindowSize, referenceWindows.Skip(from).Take(to - from).ToList());
                var (score, offset) = ZCurveWindowedProfile.BestAlignment(query, neighbourhood);
                if (best is null || score > best.Value.Score) best = (score, from + offset);
            }
            return best;
        }

        [Test]
        public void WindowCountSweep()
        {
            // Two positive classes per corpus, at the closest and the next rank each taxonomy supports.
            var bacterialPairs = PairsOf(Corpus.Count, (a, b) =>
            {
                var level = GenomeTaxonomy.Classify(Corpus[a].Meta, Corpus[b].Meta);
                return (level is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus, level is not RelationshipLevel.DifferentPhylum);
            });
            var viralPairs = PairsOf(ViralCorpus.Count, (a, b) =>
            {
                var level = ViralTaxonomy.Classify(ViralCorpus[a].Meta, ViralCorpus[b].Meta);
                return (level is ViralRelationshipLevel.SameGenus, level is not ViralRelationshipLevel.DifferentFamily);
            });

            using var report = File.CreateText("MultiResolution_WindowCount.csv");
            report.WriteLine("Corpus, WindowSize, Factor, PairCount, IncomparablePairs, ShorterWindowsMin, ShorterWindowsMedian, " +
                "AucClose_AllPairs, AucRelated_AllPairs, AucClose_ComparableOnly, AucRelated_ComparableOnly, CompareMs");

            SweepCorpus("Bacterial (close = same genus or closer, related = same phylum or closer)",
                Corpus.Select(c => c.Profile).ToList(), bacterialPairs, bacterialCountSweepFactors, report);
            SweepCorpus("Viral (close = same genus, related = same family or closer)",
                ViralCorpus.Select(c => c.Profile).ToList(), viralPairs, viralCountSweepFactors, report);

            Assert.Pass();
        }

        private static List<(int A, int B, bool Close, bool Related)> PairsOf(int count, Func<int, int, (bool Close, bool Related)> label)
        {
            var pairs = new List<(int A, int B, bool Close, bool Related)>();
            for (int a = 0; a < count; a++)
                for (int b = a + 1; b < count; b++)
                {
                    var (close, related) = label(a, b);
                    pairs.Add((a, b, close, related));
                }
            return pairs;
        }

        private static void SweepCorpus(
            string corpusName,
            List<ZCurveWindowedProfile> fine,
            List<(int A, int B, bool Close, bool Related)> pairs,
            int[] factors,
            StreamWriter report)
        {
            foreach (var factor in factors)
            {
                var profiles = fine.Select(p => p.Downsample(factor)).ToList();

                var stopwatch = Stopwatch.StartNew();
                // Compare returns 0 when either profile has no windows - production behaviour, kept
                // here so the "AllPairs" AUC shows what an unguarded caller would actually see.
                var scores = pairs.Select(p => ZCurveWindowedProfile.Compare(profiles[p.A], profiles[p.B])).ToList();
                stopwatch.Stop();

                var shorterWindows = pairs.Select(p => Math.Min(profiles[p.A].WindowCount, profiles[p.B].WindowCount)).ToList();
                var comparable = Enumerable.Range(0, pairs.Count).Where(i => shorterWindows[i] > 0).ToList();
                var sortedShorter = shorterWindows.OrderBy(w => w).ToList();

                double Auc(IEnumerable<int> indices, Func<int, bool> isPositive)
                {
                    var list = indices.ToList();
                    var positives = list.Where(isPositive).Select(i => scores[i]).ToList();
                    var negatives = list.Where(i => !isPositive(i)).Select(i => scores[i]).ToList();
                    return positives.Count == 0 || negatives.Count == 0
                        ? double.NaN
                        : RetrievalClassificationBenchmark.ComputeAuc(positives, negatives);
                }

                var all = Enumerable.Range(0, pairs.Count);
                var row = (
                    CloseAll: Auc(all, i => pairs[i].Close),
                    RelatedAll: Auc(all, i => pairs[i].Related),
                    CloseComparable: Auc(comparable, i => pairs[i].Close),
                    RelatedComparable: Auc(comparable, i => pairs[i].Related));

                report.WriteLine($"\"{corpusName}\", {fineWindowSize * factor}, {factor}, {pairs.Count}, {pairs.Count - comparable.Count}, " +
                    $"{sortedShorter[0]}, {sortedShorter[sortedShorter.Count / 2]}, " +
                    $"{row.CloseAll}, {row.RelatedAll}, {row.CloseComparable}, {row.RelatedComparable}, {stopwatch.Elapsed.TotalMilliseconds}");
                report.Flush();

                TestContext.Out.WriteLine(
                    $"{corpusName.Split(' ')[0],-9} w={fineWindowSize * factor,8}  shorter windows min={sortedShorter[0],5} median={sortedShorter[sortedShorter.Count / 2],5}  " +
                    $"incomparable={pairs.Count - comparable.Count,3}  AUC close={row.CloseAll:F3}/{row.CloseComparable:F3}  related={row.RelatedAll:F3}/{row.RelatedComparable:F3}");
            }
        }

        [Test]
        public void LadderRuleAgainstFineOnly()
        {
            using var report = File.CreateText("MultiResolution_LadderRule.csv");
            report.WriteLine("Corpus, Method, PairCount, " +
                string.Join(", ", ZCurveProfileLadder.DefaultFactors.Select(f => $"PairsAt{fineWindowSize * f}bp")) +
                ", NullPairs, AucClose, AucRelated, CompareMs");

            var bacterial = Corpus.Select(c => c.Profile).ToList();
            var viral = ViralCorpus.Select(c => c.Profile).ToList();
            LadderVersusFine("Bacterial", bacterial, PairsOf(bacterial.Count, (a, b) =>
            {
                var level = GenomeTaxonomy.Classify(Corpus[a].Meta, Corpus[b].Meta);
                return (level is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus, level is not RelationshipLevel.DifferentPhylum);
            }), report);
            LadderVersusFine("Viral", viral, PairsOf(viral.Count, (a, b) =>
            {
                var level = ViralTaxonomy.Classify(ViralCorpus[a].Meta, ViralCorpus[b].Meta);
                return (level is ViralRelationshipLevel.SameGenus, level is not ViralRelationshipLevel.DifferentFamily);
            }), report);

            Assert.Pass();
        }

        private static void LadderVersusFine(
            string corpusName,
            List<ZCurveWindowedProfile> fine,
            List<(int A, int B, bool Close, bool Related)> pairs,
            StreamWriter report)
        {
            // Ladders are what the platform would store, so building them isn't timed.
            var ladders = fine.Select(p => new ZCurveProfileLadder(p)).ToList();

            var fineWatch = Stopwatch.StartNew();
            var fineScores = pairs.Select(p => (double?)ZCurveWindowedProfile.Compare(fine[p.A], fine[p.B])).ToList();
            fineWatch.Stop();

            var ladderWatch = Stopwatch.StartNew();
            var ladderScores = pairs.Select(p => ZCurveProfileLadder.Compare(ladders[p.A], ladders[p.B])).ToList();
            ladderWatch.Stop();

            var chosenFactors = pairs.Select(p => ZCurveProfileLadder.SelectFactor(ladders[p.A], ladders[p.B])).ToList();

            void Write(string method, List<double?> scores, double ms, bool countLevels)
            {
                var comparable = Enumerable.Range(0, pairs.Count).Where(i => scores[i].HasValue).ToList();
                double Auc(Func<int, bool> isPositive) => RetrievalClassificationBenchmark.ComputeAuc(
                    comparable.Where(isPositive).Select(i => scores[i]!.Value).ToList(),
                    comparable.Where(i => !isPositive(i)).Select(i => scores[i]!.Value).ToList());

                var levelCounts = ZCurveProfileLadder.DefaultFactors.Select(f =>
                    countLevels ? chosenFactors.Count(c => c == f) : (f == 1 ? pairs.Count : 0));
                var aucClose = Auc(i => pairs[i].Close);
                var aucRelated = Auc(i => pairs[i].Related);

                report.WriteLine($"{corpusName}, {method}, {pairs.Count}, {string.Join(", ", levelCounts)}, " +
                    $"{pairs.Count - comparable.Count}, {aucClose}, {aucRelated}, {ms}");
                TestContext.Out.WriteLine($"{corpusName,-9} {method,-6} AUC close={aucClose:F4} related={aucRelated:F4}  " +
                    $"levels=[{string.Join(", ", levelCounts)}]  compare={ms,9:F1} ms");
            }

            Write("Fine", fineScores, fineWatch.Elapsed.TotalMilliseconds, countLevels: false);
            Write("Ladder", ladderScores, ladderWatch.Elapsed.TotalMilliseconds, countLevels: true);
            report.Flush();
        }

        // C4 (platform P3): a search scores every candidate at one level chosen by the query alone
        // (ZCurveProfileLadder.CompareForQuery) and moves candidates too short for it to a separate
        // short list. Pairs become ordered (query, candidate). Ranked scores can't differ from the
        // per-pair ladder: a candidate with >= 8 windows at the query's level gets that level per pair
        // too, since the query has < 8 at every coarser one (ScoreChangedPairs checks it). What changes
        // is that short-listed pairs leave the ranking, so both rules are also compared on the same
        // ranked subset, pooled and as the mean of each query's own AUC (what one search returns).
        [Test]
        public void PerQueryLevelAgainstPerPair()
        {
            using var report = File.CreateText("MultiResolution_PerQueryLevel.csv");
            report.WriteLine("Corpus, Method, Subset, OrderedPairs, RankedPairs, ShortListPairs, NotComparablePairs, ScoreChangedPairs, " +
                "AucClose, AucRelated, QueriesClose, MeanQueryAucClose, QueriesRelated, MeanQueryAucRelated");

            PerQueryVersusPerPair("Bacterial", Corpus.Select(c => new ZCurveProfileLadder(c.Profile)).ToList(), (a, b) =>
            {
                var level = GenomeTaxonomy.Classify(Corpus[a].Meta, Corpus[b].Meta);
                return (level is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus, level is not RelationshipLevel.DifferentPhylum);
            }, report);
            PerQueryVersusPerPair("Viral", ViralCorpus.Select(c => new ZCurveProfileLadder(c.Profile)).ToList(), (a, b) =>
            {
                var level = ViralTaxonomy.Classify(ViralCorpus[a].Meta, ViralCorpus[b].Meta);
                return (level is ViralRelationshipLevel.SameGenus, level is not ViralRelationshipLevel.DifferentFamily);
            }, report);

            Assert.Pass();
        }

        private static void PerQueryVersusPerPair(
            string corpusName,
            List<ZCurveProfileLadder> ladders,
            Func<int, int, (bool Close, bool Related)> label,
            StreamWriter report)
        {
            var pairs = new List<(int Query, int Candidate, bool Close, bool Related, double? PerPair, SearchMatch PerQuery)>();
            for (int q = 0; q < ladders.Count; q++)
                for (int c = 0; c < ladders.Count; c++)
                {
                    if (q == c) continue;
                    var (close, related) = label(q, c);
                    pairs.Add((q, c, close, related,
                        ZCurveProfileLadder.Compare(ladders[q], ladders[c]),
                        ZCurveProfileLadder.CompareForQuery(ladders[q], ladders[c])));
                }

            var ranked = pairs.Where(p => p.PerQuery.List == SearchList.Ranked).ToList();
            int shortList = pairs.Count(p => p.PerQuery.List == SearchList.ShortList);
            int notComparable = pairs.Count(p => p.PerQuery.List == SearchList.NotComparable);
            int changed = ranked.Count(p => p.PerQuery.Score != p.PerPair);

            static double Auc(IEnumerable<(bool Positive, double Score)> rows)
            {
                var list = rows.ToList();
                var positives = list.Where(r => r.Positive).Select(r => r.Score).ToList();
                var negatives = list.Where(r => !r.Positive).Select(r => r.Score).ToList();
                return positives.Count == 0 || negatives.Count == 0
                    ? double.NaN
                    : RetrievalClassificationBenchmark.ComputeAuc(positives, negatives);
            }

            void Write(string method, string subset, List<(int Query, bool Close, bool Related, double Score)> rows)
            {
                (int Queries, double Mean) PerQuery(Func<(int Query, bool Close, bool Related, double Score), bool> isPositive)
                {
                    var aucs = rows.GroupBy(r => r.Query)
                        .Select(g => Auc(g.Select(r => (isPositive(r), r.Score))))
                        .Where(a => !double.IsNaN(a))
                        .ToList();
                    return (aucs.Count, aucs.Count == 0 ? double.NaN : aucs.Average());
                }

                var aucClose = Auc(rows.Select(r => (r.Close, r.Score)));
                var aucRelated = Auc(rows.Select(r => (r.Related, r.Score)));
                var queryClose = PerQuery(r => r.Close);
                var queryRelated = PerQuery(r => r.Related);

                report.WriteLine($"{corpusName}, {method}, {subset}, {pairs.Count}, {ranked.Count}, {shortList}, {notComparable}, {changed}, " +
                    $"{aucClose}, {aucRelated}, {queryClose.Queries}, {queryClose.Mean}, {queryRelated.Queries}, {queryRelated.Mean}");
                TestContext.Out.WriteLine($"{corpusName,-9} {method,-8} {subset,-10} pooled AUC close={aucClose:F4} related={aucRelated:F4}  " +
                    $"per-query close={queryClose.Mean:F4} (n={queryClose.Queries}) related={queryRelated.Mean:F4} (n={queryRelated.Queries})");
            }

            Write("PerPair", "All", pairs.Where(p => p.PerPair.HasValue)
                .Select(p => (p.Query, p.Close, p.Related, p.PerPair!.Value)).ToList());
            Write("PerPair", "RankedOnly", ranked.Select(p => (p.Query, p.Close, p.Related, p.PerPair!.Value)).ToList());
            Write("PerQuery", "RankedOnly", ranked.Select(p => (p.Query, p.Close, p.Related, p.PerQuery.Score!.Value)).ToList());
            TestContext.Out.WriteLine($"{corpusName,-9} ordered pairs={pairs.Count} ranked={ranked.Count} short={shortList} " +
                $"not comparable={notComparable} score changed={changed}");
            report.Flush();
        }

        [Test]
        public void GcAxisWeightSweep()
        {
            using var report = File.CreateText("MultiResolution_GcAxisWeight.csv");
            report.WriteLine("GcAxisWeight, AucSameGenusOrCloser, AucSamePhylumOrCloser, SameGenus Mean, DifferentPhylum Mean");

            foreach (var weight in gcAxisWeights)
            {
                // sqrt(dx^2 + dy^2 + weight*dz^2) is the plain Euclidean distance after scaling z by
                // sqrt(weight), so the production SIMD Compare scores the weighted metric unchanged.
                var scale = Math.Sqrt(weight);
                var profiles = Corpus.Select(c => (c.Meta, Profile: ZCurveWindowedProfile.FromPrecomputedWindows(
                    c.Profile.WindowSize, c.Profile.Windows.Select(w => (w.X, w.Y, w.Z * scale)).ToList()))).ToList();

                var rows = ScoreAllPairs(profiles);
                var (aucGenus, aucPhylum) = Aucs(rows);
                var sameGenus = rows.Where(r => r.Level == RelationshipLevel.SameGenus).Average(r => r.Score);
                var differentPhylum = rows.Where(r => r.Level == RelationshipLevel.DifferentPhylum).Average(r => r.Score);

                report.WriteLine($"{weight}, {aucGenus}, {aucPhylum}, {sameGenus}, {differentPhylum}");
                report.Flush();
                TestContext.Out.WriteLine($"weight={weight,5}  AUC genus={aucGenus:F4} phylum={aucPhylum:F4}");
            }

            Assert.Pass();
        }
    }
}
