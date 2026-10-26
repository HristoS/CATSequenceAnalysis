using System.Globalization;
using Npgsql;
using SequenceAnaliseTests;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    /// <summary>
    /// P3 step S3.7: the engine, run through the database on the bundled corpora, against the library results the manuscript
    /// reports (criteria restated 2026-10-26, your choice, after the first run; results/S37_*.csv).
    ///
    /// Find-similar: the library's length-band computation (MultiResolutionBenchmark.BandCorpus, r = 0.5) redone from the DNA must
    /// reproduce the committed mean per-query AUCs (results/MultiResolution_LengthBand.csv), and, with the platform's two designed
    /// differences added (the better of both strands, S3.2; the band on length_bp rather than window counts, S3.4), must equal what
    /// the platform returns through the database. Against the committed values that costs at most 0.0102 (viral close relatives:
    /// half strand, half band edge).
    ///
    /// Localization: the fragments of results/MultiResolution_CoarseToFine.csv, forward and reverse-complemented, are judged against
    /// their true start (±1 window), not the committed full fine search, which is itself wrong for one 10 kbp Vibrio fragment.
    /// Every 10 and 50 kbp one, localized over the whole bundled corpus, must come first in its own record on the right strand there.
    /// The 2 kbp ones (under Localization.MinCatalogueLengthBp, placed inside their record) are reported, not required: 4 windows
    /// are a weak signal. Writes S37_*.csv to the test directory.
    /// </summary>
    [Explicit("Profiles and stores all 61 bundled genomes (about 100 Mbp); run it by name.")]
    public class PlatformValidationTests
    {
        private static readonly double[] committedBand = { 0.8529028075467704, 0.7446777108560474, 0.9115482556659028, 0.7349647266313935 };

        private static Task<IReadOnlyDictionary<string, long>>? loading;

        // Every bundled genome loaded once per run; the primary record's id per file.
        private static Task<IReadOnlyDictionary<string, long>> PrimaryRecordsAsync() => loading ??= LoadAsync();

        private static async Task<IReadOnlyDictionary<string, long>> LoadAsync()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            foreach (var sql in new[]
            {
                "delete from profile_level where sequence_id in (select r.id from sequence_record r join sequence_source s on s.id = r.source_id where s.provider = 'NCBI RefSeq')",
                "delete from sequence_record where source_id in (select id from sequence_source where provider = 'NCBI RefSeq')",
                "delete from sequence_source where provider = 'NCBI RefSeq'",
            })
            {
                await using var command = dataSource.CreateCommand(sql);
                await command.ExecuteNonQueryAsync();
            }

            var primaries = new Dictionary<string, long>();
            foreach (var genome in BundledCorpus.Genomes)
            {
                var sourceId = await BundledCorpus.LoadAsync(dataSource, genome);
                var primary = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genome.FileName)).Id;
                await using var command = dataSource.CreateCommand("select id from sequence_record where source_id = $1 and record_id = $2");
                command.Parameters.AddWithValue(sourceId);
                command.Parameters.AddWithValue(primary);
                primaries[genome.FileName] = (long)(await command.ExecuteScalarAsync())!;
            }
            return primaries;
        }

        [Test]
        public async Task FindSimilarReproducesTheLibraryLengthBand()
        {
            var primaries = await PrimaryRecordsAsync();
            var search = new SimilaritySearch(ThrowawayDatabase.DataSource);

            var bacteria = GenomeTaxonomy.Genomes.Select(g => (primaries[g.FileName], (object)g)).ToList();
            var viruses = ViralTaxonomy.Genomes.Select(g => (primaries[g.FileName], (object)g)).ToList();
            var measured = new List<(string Corpus, int Queries, double Close, double Related)>
            {
                await MeanAucsAsync("Bacterial", bacteria, (a, b) =>
                {
                    var level = GenomeTaxonomy.Classify((GenomeMeta)a, (GenomeMeta)b);
                    return (level is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus, level is not RelationshipLevel.DifferentPhylum);
                }, search),
                await MeanAucsAsync("Viral", viruses, (a, b) =>
                {
                    var level = ViralTaxonomy.Classify((ViralGenomeMeta)a, (ViralGenomeMeta)b);
                    return (level is ViralRelationshipLevel.SameGenus, level is not ViralRelationshipLevel.DifferentFamily);
                }, search),
            };

            // The library computation redone from the DNA, one difference at a time: as committed (forward strand, band on window
            // counts), with the better of both strands (the reverse a profile flip, as for a stored record), and with the band on
            // length_bp as the platform has it. The last should equal the platform.
            var recomputed = new List<(double Close, double Related)[]>
            {
                LibraryVariants(GenomeTaxonomy.Genomes.Select(g => (object)g).ToList(), g => ((GenomeMeta)g).FileName, (a, b) =>
                {
                    var level = GenomeTaxonomy.Classify((GenomeMeta)a, (GenomeMeta)b);
                    return (level is RelationshipLevel.SameSpecies or RelationshipLevel.SameGenus, level is not RelationshipLevel.DifferentPhylum);
                }),
                LibraryVariants(ViralTaxonomy.Genomes.Select(g => (object)g).ToList(), g => ((ViralGenomeMeta)g).FileName, (a, b) =>
                {
                    var level = ViralTaxonomy.Classify((ViralGenomeMeta)a, (ViralGenomeMeta)b);
                    return (level is ViralRelationshipLevel.SameGenus, level is not ViralRelationshipLevel.DifferentFamily);
                }),
            };

            using var report = File.CreateText("S37_FindSimilarValidation.csv");
            report.WriteLine("Corpus, Relationship, PlatformMeanQueryAuc, LibraryMeanQueryAuc, Difference, " +
                "RecomputedForwardWindowBand, RecomputedBothStrandsWindowBand, RecomputedBothStrandsLengthBand");
            var rows = new List<(string Corpus, string Relationship, double Platform, double Library)>();
            for (int i = 0; i < measured.Count; i++)
            {
                rows.Add((measured[i].Corpus, "close", measured[i].Close, committedBand[2 * i]));
                rows.Add((measured[i].Corpus, "related", measured[i].Related, committedBand[2 * i + 1]));
            }
            for (int i = 0; i < rows.Count; i++)
            {
                var (corpus, relationship, platform, committed) = rows[i];
                var variants = recomputed[i / 2].Select(v => i % 2 == 0 ? v.Close : v.Related);
                report.WriteLine(string.Create(CultureInfo.InvariantCulture,
                    $"{corpus}, {relationship}, {platform}, {committed}, {platform - committed}, {string.Join(", ", variants)}"));
            }
            report.Flush();

            for (int i = 0; i < rows.Count; i++)
            {
                var (corpus, relationship, platform, committed) = rows[i];
                var variants = recomputed[i / 2].Select(v => i % 2 == 0 ? v.Close : v.Related).ToArray();
                Assert.That(variants[0], Is.EqualTo(committed).Within(1e-12), $"{corpus} {relationship}: the recomputation reproduces the committed value");
                Assert.That(platform, Is.EqualTo(variants[2]).Within(1e-12), $"{corpus} {relationship}: the platform equals the library under its rules");
                Assert.That(platform, Is.EqualTo(committed).Within(0.011), $"{corpus} {relationship}: the designed differences stay small");
            }
        }

        // As MultiResolutionBenchmark.BandCorpus: each primary record searched (r = 0.5, by score), AUC over the corpus's other
        // primary records in its band, averaged over the queries that have both positives and negatives.
        private static async Task<(string Corpus, int Queries, double Close, double Related)> MeanAucsAsync(string corpus,
            List<(long Id, object Meta)> members, Func<object, object, (bool Close, bool Related)> label, SimilaritySearch search)
        {
            var metaById = members.ToDictionary(m => m.Id, m => m.Meta);
            var close = new List<double>();
            var related = new List<double>();
            foreach (var (id, meta) in members)
            {
                var result = await search.SearchRecordAsync(id, top: 100_000, useNeighbours: false);
                var scored = result.Ranked.Where(h => metaById.ContainsKey(h.SequenceId))
                    .Select(h => (Label: label(meta, metaById[h.SequenceId]), h.Score)).ToList();
                if (Auc(scored.Select(s => (s.Label.Close, s.Score))) is double c) close.Add(c);
                if (Auc(scored.Select(s => (s.Label.Related, s.Score))) is double r) related.Add(r);
            }
            return (corpus, members.Count, close.Average(), related.Average());
        }

        // MultiResolutionBenchmark.BandCorpus's per-query AUC from the bundled DNA, in three variants (see the caller).
        private static (double Close, double Related)[] LibraryVariants(List<object> members, Func<object, string> fileName,
            Func<object, object, (bool Close, bool Related)> label)
        {
            var sequences = members.Select(m => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName(m))).Sequence).ToList();
            var forward = sequences.Select(s => new ZCurveProfileLadder(new ZCurveWindowedProfile(s))).ToList();
            var reverse = forward.Select(l => new ZCurveProfileLadder(ZCurveWindowedProfile.FromPrecomputedWindows(l.Fine.WindowSize,
                l.Fine.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList()))).ToList();

            (double, double) Variant(bool bothStrands, bool lengthBand)
            {
                var close = new List<double>();
                var related = new List<double>();
                for (int q = 0; q < members.Count; q++)
                {
                    if (ZCurveProfileLadder.SelectBandFactor(forward[q], 0.5) is not int factor) continue;
                    var band = SimilaritySearch.Band(sequences[q].Length, 0.5);
                    var scored = new List<((bool Close, bool Related) Label, double Score)>();
                    for (int c = 0; c < members.Count; c++)
                    {
                        if (c == q) continue;
                        bool kept = lengthBand
                            ? band.Contains(sequences[c].Length)
                            : Math.Min(forward[q].Fine.WindowCount, forward[c].Fine.WindowCount) >= 0.5 * Math.Max(forward[q].Fine.WindowCount, forward[c].Fine.WindowCount);
                        if (!kept || forward[c].Level(factor).WindowCount == 0) continue;
                        double score = ZCurveWindowedProfile.Compare(forward[q].Level(factor), forward[c].Level(factor));
                        if (bothStrands) score = Math.Max(score, ZCurveWindowedProfile.Compare(reverse[q].Level(factor), forward[c].Level(factor)));
                        scored.Add((label(members[q], members[c]), score));
                    }
                    if (Auc(scored.Select(s => (s.Label.Close, s.Score))) is double a) close.Add(a);
                    if (Auc(scored.Select(s => (s.Label.Related, s.Score))) is double b) related.Add(b);
                }
                return (close.Average(), related.Average());
            }

            return new[] { Variant(false, false), Variant(true, false), Variant(true, true) };
        }

        // Mann-Whitney AUC with ties counted half, as RetrievalClassificationBenchmark.ComputeAuc; null without both classes.
        private static double? Auc(IEnumerable<(bool Positive, double Score)> rows)
        {
            var list = rows.ToList();
            var positives = list.Where(r => r.Positive).Select(r => r.Score).ToList();
            var negatives = list.Where(r => !r.Positive).Select(r => r.Score).ToList();
            if (positives.Count == 0 || negatives.Count == 0) return null;
            double wins = 0;
            foreach (var p in positives)
                foreach (var n in negatives)
                    wins += p > n ? 1 : p == n ? 0.5 : 0;
            return wins / (positives.Count * (double)negatives.Count);
        }

        [Test]
        public async Task LocalizationReproducesTheCoarseToFineTrials()
        {
            var primaries = await PrimaryRecordsAsync();
            var localization = new Localization(ThrowawayDatabase.DataSource);
            var trialsPath = Path.GetFullPath(Path.Combine(GenomeCatalog.DirectoryPath, "..", "..", "results", "MultiResolution_CoarseToFine.csv"));
            var trials = File.ReadLines(trialsPath).Skip(1)
                .Select(line => line.Split(',').Select(f => f.Trim()).ToArray())
                .Select(f => (Length: int.Parse(f[0]), Trial: int.Parse(f[1]), Genome: f[2], Start: int.Parse(f[3]), FineOffset: int.Parse(f[6])))
                .Distinct()
                .ToList();

            using var report = File.CreateText("S37_LocalizationValidation.csv");
            report.WriteLine("FragmentLength, Trial, Genome, Start, Strand, Mode, CommittedFineOffset, TopIsSource, TopStrand, OffsetWindows, Score, " +
                "AtCommittedFineOffset, AtTrueStart");
            var failures = new List<string>();
            int shortPlaced = 0, shortTotal = 0;
            foreach (var (length, trial, genome, start, fineOffset) in trials)
            {
                var dna = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genome)).Sequence;
                var piece = dna.Substring(start, length);
                long source = primaries[genome];
                foreach (var minus in new[] { false, true })
                {
                    var query = QueryProfile.FromSequence(minus ? QueryProfile.ReverseComplement(piece) : piece);
                    bool inRecord = length < Localization.MinCatalogueLengthBp;
                    var result = inRecord
                        ? await localization.LocalizeInRecordAsync(query, source)
                        : await localization.LocalizeAsync(query, top: 5);
                    var top = result.Hits[0];
                    long offsetWindows = top.OffsetBp / RecordStore.FineWindowSize;
                    var expectedStrand = minus ? Strand.Reverse : Strand.Forward;
                    bool placed = top.SequenceId == source && top.Strand == expectedStrand;
                    bool atFine = placed && (minus ? Math.Abs(offsetWindows - fineOffset) <= 1 : offsetWindows == fineOffset);
                    bool atTruth = placed && Math.Abs(offsetWindows - start / (double)RecordStore.FineWindowSize) <= 1;
                    report.WriteLine($"{length}, {trial}, {genome}, {start}, {(minus ? '-' : '+')}, {(inRecord ? "named record" : "catalogue")}, " +
                        $"{fineOffset}, {top.SequenceId == source}, {top.Strand}, {offsetWindows}, {top.Score.ToString("F4", CultureInfo.InvariantCulture)}, " +
                        $"{atFine}, {atTruth}");
                    if (inRecord)
                    {
                        shortTotal++;
                        if (atTruth) shortPlaced++;
                    }
                    else if (!atTruth)
                        failures.Add($"{length} bp trial {trial} ({genome}, {(minus ? "minus" : "plus")}): top {top.SequenceId} {top.Strand} " +
                            $"at window {offsetWindows}, expected {source} at {start / (double)RecordStore.FineWindowSize:F1}");
                }
            }
            report.Flush();
            TestContext.Out.WriteLine($"2 kbp inside their own record: {shortPlaced} of {shortTotal} at the true start (reported, not required)");

            Assert.That(trials, Has.Count.EqualTo(30));
            Assert.That(failures, Is.Empty);
        }
    }
}
