using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Platform question R6 (P7): can a protein query use Z-curve profiles at all? A protein has no
    /// nucleotides, so it has to be back-translated first, and the codons chosen are a guess. This
    /// takes in-frame regions of the bundled genomes, translates them, back-translates them several
    /// ways, and localizes each result in its own source genome with BestAlignment - the easiest case
    /// a protein search could face. The Dna variant (the region itself) is the control, so the gap
    /// between it and the back-translated variants is the synonymous-codon effect alone.
    ///
    /// Query starts sit on window boundaries and only forward-strand regions are used, which keeps
    /// window phase and strand out of the measurement; both would cost a real protein query more.
    /// </summary>
    public class ProteinBackTranslationFeasibility
    {
        private const int windowSize = ZCurveWindowedProfile.DefaultWindowSize;
        private const int regionsPerGenome = 20;
        private const int codonUsageMinLength = 900;
        // Two thirds of the default 500 bp: a positions-1-and-2 window spans the same ~500 bp of gene.
        private const int pos12WindowSize = 333;
        private const string codonBases = "TCAG";
        // NCBI translation tables in TCAG order (index = 16*first + 4*second + third). Table 11
        // (bacteria) translates internal codons exactly like table 1; table 4 (Mycoplasma) reads TGA
        // as Trp instead of a stop.
        private const string standardCode = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        private const string mycoplasmaCode = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        private const string genericUsageFile = "ecoli_k12_mg1655_GCF_000005845.2.fna.gz";
        // R6d: a target gene is the query's ortholog if it holds at least this share of the query's amino-acid 5-mers.
        private const double orthologMinShare = 0.3;
        private const int homologQueryLength = 1500;
        private const int homologQueryCount = 200;
        private const int orthologSearchMinLength = 300;
        private static readonly int[] queryLengths = { 1500, 3000 };
        private static readonly string[] variants =
            { "Dna", "MostFrequentOwn", "WeightedRandomOwn", "MostFrequentGeneric", "UniformRandom", "RandomDna" };
        private static readonly string[] homologTargets =
        {
            "ecoli_o157h7_sakai_GCF_000008865.2.fna.gz",
            "shigella_flexneri_2a_301_GCF_000006925.2.fna.gz",
            "salmonella_typhimurium_lt2_GCF_000006945.2.fna.gz",
            "klebsiella_pneumoniae_mgh78578_GCF_000016305.1.fna.gz",
        };

        [Test]
        public void BackTranslatedRegionsLocalizeInTheirSourceGenome()
        {
            var genericSequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genericUsageFile)).Sequence;
            var genericUsage = CodonUsage(genericSequence, OpenStretches(genericSequence, standardCode, codonUsageMinLength));
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();

            var rows = new List<QueryResult>();
            using (var raw = File.CreateText("R6_ProteinBackTranslation_Raw.csv"))
            {
                raw.WriteLine("Group,FileName,GeneticCode,ReferenceWindows,QueryLengthBp,RegionIndex,TrueStartBp,Variant," +
                              "DifferingFraction,BestOffsetBp,ExactHit,WithinOneWindow,TrueRank,TrueScore,BestFalseScore,ScoreMargin");

                foreach (var (group, fileName) in Corpus())
                {
                    var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence;
                    var code = Path.GetFileName(fileName).StartsWith("mycoplasma_", StringComparison.Ordinal) ? mycoplasmaCode : standardCode;
                    var synonyms = Synonyms(code);
                    var ownUsage = CodonUsage(sequence, OpenStretches(sequence, code, codonUsageMinLength));
                    var reference = new ZCurveWindowedProfile(sequence);
                    var referenceWindows = reference.Windows.ToArray();

                    foreach (var length in queryLengths)
                    {
                        var seed = StableSeed(fileName) ^ length;
                        var picker = new Random(seed);
                        var picked = AlignedStarts(OpenStretches(sequence, code, length), length)
                            .OrderBy(_ => picker.Next()).Take(regionsPerGenome).OrderBy(s => s).ToList();

                        for (int regionIndex = 0; regionIndex < picked.Count; regionIndex++)
                        {
                            var start = picked[regionIndex];
                            var original = sequence.Substring(start, length);
                            var protein = Translate(original, code);
                            var random = new Random(seed * 31 + regionIndex);

                            foreach (var variant in variants)
                            {
                                var query = BuildVariant(variant, original, protein, synonyms, ownUsage, genericUsage, uniformUsage, random);
                                var result = Localize(query, reference, referenceWindows, start / windowSize);
                                var differing = DifferingFraction(original, query);
                                rows.Add(new QueryResult(group, length, variant, differing, result));
                                raw.WriteLine(string.Join(",",
                                    group, fileName, code == mycoplasmaCode ? 4 : 11, referenceWindows.Length, length, regionIndex, start, variant,
                                    Format(differing), result.BestOffset * windowSize,
                                    result.ExactHit.ToString().ToLowerInvariant(), result.WithinOneWindow.ToString().ToLowerInvariant(),
                                    result.TrueRank, Format(result.TrueScore), Format(result.BestFalseScore),
                                    Format(result.TrueScore - result.BestFalseScore)));
                            }
                        }
                    }
                }
            }

            using (var summary = File.CreateText("R6_ProteinBackTranslation_Summary.csv"))
            {
                summary.WriteLine("Group,QueryLengthBp,Variant,Queries,Top1,Top1WithinOneWindow,MedianTrueRank,MedianScoreMargin,MeanDifferingFraction");
                foreach (var cell in rows
                    .GroupBy(r => (r.Group, r.Length, r.Variant))
                    .OrderBy(g => g.Key.Group).ThenBy(g => g.Key.Length).ThenBy(g => Array.IndexOf(variants, g.Key.Variant)))
                {
                    var line = string.Join(",",
                        cell.Key.Group, cell.Key.Length, cell.Key.Variant, cell.Count(),
                        Format(cell.Average(r => r.Result.ExactHit ? 1.0 : 0.0)),
                        Format(cell.Average(r => r.Result.WithinOneWindow ? 1.0 : 0.0)),
                        Format(Median(cell.Select(r => (double)r.Result.TrueRank))),
                        Format(Median(cell.Select(r => r.Result.TrueScore - r.Result.BestFalseScore))),
                        Format(cell.Average(r => r.DifferingFraction)));
                    summary.WriteLine(line);
                    TestContext.Out.WriteLine(line);
                }
            }

            // Validity only: if the unmodified region doesn't find itself, the harness is wrong and the
            // back-translated rates mean nothing.
            foreach (var group in rows.Where(r => r.Variant == "Dna").GroupBy(r => r.Group))
                Assert.That(group.Average(r => r.Result.ExactHit ? 1.0 : 0.0), Is.GreaterThanOrEqualTo(0.9), group.Key);
        }

        /// <summary>
        /// R6 follow-up: the two things R6 left out (reverse-strand genes, query starts off the window
        /// grid) and gene-aligned windows as a possible rescue. FixedWindows is what a protein query
        /// faces without annotations: the back-translated gene start lands at any phase against the
        /// stored 500 bp windows, on either strand. GeneAligned is the best case: one profile per gene,
        /// windows starting at the gene start in coding orientation, so phase and strand drop out and the
        /// query only competes against the genome's other genes. A gene start here is the first codon
        /// after a stop, a stand-in for annotations the bundled data doesn't carry.
        /// </summary>
        [Test]
        public void GeneAlignedWindowsRescueBackTranslation()
        {
            var genericSequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genericUsageFile)).Sequence;
            var genericUsage = CodonUsage(genericSequence, OpenStretches(genericSequence, standardCode, codonUsageMinLength));
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();

            var rows = new List<GeneQueryResult>();
            using (var raw = File.CreateText("R6b_GeneAligned_Raw.csv"))
            {
                raw.WriteLine("Group,FileName,Mode,QueryLengthBp,GeneIndex,Strand,ForwardStartBp,Variant,DifferingFraction," +
                              "Candidates,TrueRank,Top1,TrueScore,BestFalseScore,MatchedReverseComplement");

                foreach (var (group, fileName) in Corpus())
                {
                    var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence;
                    var reverse = ReverseComplement(forward);
                    var code = Path.GetFileName(fileName).StartsWith("mycoplasma_", StringComparison.Ordinal) ? mycoplasmaCode : standardCode;
                    var synonyms = Synonyms(code);
                    var ownUsage = CodonUsage(forward, OpenStretches(forward, code, codonUsageMinLength));
                    var genomeWindows = new ZCurveWindowedProfile(forward).Windows.ToArray();

                    foreach (var length in queryLengths)
                    {
                        // Every gene long enough for an L bp prefix, on both strands; Start is in that strand's coordinates.
                        var genes = OpenStretches(forward, code, length).Select(s => (Reverse: false, s.Start))
                            .Concat(OpenStretches(reverse, code, length).Select(s => (Reverse: true, s.Start)))
                            .ToList();
                        var index = genes
                            .Select(g => new ZCurveWindowedProfile((g.Reverse ? reverse : forward).Substring(g.Start, length)).Windows.ToArray())
                            .ToList();

                        var seed = StableSeed(fileName) ^ length ^ 0x5bd1e995;
                        var picker = new Random(seed);
                        var picked = Enumerable.Range(0, genes.Count).OrderBy(_ => picker.Next()).Take(regionsPerGenome).OrderBy(i => i).ToList();

                        foreach (var geneIndex in picked)
                        {
                            var gene = genes[geneIndex];
                            var original = (gene.Reverse ? reverse : forward).Substring(gene.Start, length);
                            var forwardStart = gene.Reverse ? forward.Length - gene.Start - length : gene.Start;
                            var protein = Translate(original, code);
                            var random = new Random(seed * 31 + geneIndex);

                            foreach (var variant in variants)
                            {
                                var query = BuildVariant(variant, original, protein, synonyms, ownUsage, genericUsage, uniformUsage, random);
                                var queryWindows = new ZCurveWindowedProfile(query).Windows.ToArray();
                                var differing = DifferingFraction(original, query);

                                var aligned = RankAmongGenes(queryWindows, index, geneIndex);
                                var fixedWindows = RankOnFixedWindows(query, genomeWindows, forwardStart, gene.Reverse);
                                foreach (var (mode, result) in new[] { ("GeneAligned", aligned), ("FixedWindows", fixedWindows) })
                                {
                                    rows.Add(new GeneQueryResult(group, mode, length, variant, differing, result));
                                    raw.WriteLine(string.Join(",",
                                        group, fileName, mode, length, geneIndex, gene.Reverse ? "-" : "+", forwardStart, variant,
                                        Format(differing), result.Candidates, result.TrueRank,
                                        (result.TrueRank == 1).ToString().ToLowerInvariant(),
                                        Format(result.TrueScore), Format(result.BestFalseScore),
                                        result.MatchedReverseComplement.ToString().ToLowerInvariant()));
                                }
                            }
                        }
                    }
                }
            }

            using (var summary = File.CreateText("R6b_GeneAligned_Summary.csv"))
            {
                summary.WriteLine("Group,Mode,QueryLengthBp,Variant,Queries,Top1,MedianTrueRank,MedianCandidates,MeanDifferingFraction");
                foreach (var cell in rows
                    .GroupBy(r => (r.Group, r.Mode, r.Length, r.Variant))
                    .OrderBy(g => g.Key.Group).ThenBy(g => g.Key.Mode).ThenBy(g => g.Key.Length)
                    .ThenBy(g => Array.IndexOf(variants, g.Key.Variant)))
                {
                    var line = string.Join(",",
                        cell.Key.Group, cell.Key.Mode, cell.Key.Length, cell.Key.Variant, cell.Count(),
                        Format(cell.Average(r => r.Result.TrueRank == 1 ? 1.0 : 0.0)),
                        Format(Median(cell.Select(r => (double)r.Result.TrueRank))),
                        Format(Median(cell.Select(r => (double)r.Result.Candidates))),
                        Format(cell.Average(r => r.DifferingFraction)));
                    summary.WriteLine(line);
                    TestContext.Out.WriteLine(line);
                }
            }

            foreach (var group in rows.Where(r => r.Mode == "GeneAligned" && r.Variant == "Dna").GroupBy(r => r.Group))
                Assert.That(group.Average(r => r.Result.TrueRank == 1 ? 1.0 : 0.0), Is.GreaterThanOrEqualTo(0.9), group.Key);
        }

        /// <summary>
        /// R6 second follow-up: windows that follow the gene instead of a fixed size. GeneWhole is one
        /// window spanning the whole gene (exons and introns; the bundled genomes have no detectable
        /// introns, so a gene is one stop-free stretch). GenePos12 keeps GeneAligned's ~500 bp windows
        /// but profiles only codon positions 1 and 2, where synonymous changes almost never land (Leu,
        /// Arg and Ser are the exceptions). GeneWholePos12 combines the two. GeneAligned is repeated on
        /// the same genes as the baseline. All four rank the query among every gene of its genome that
        /// is at least MinGeneLengthBp long; the back-translation is done once per gene and shared.
        /// </summary>
        [Test]
        public void PerGeneWindowsRescueBackTranslation()
        {
            var genericSequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genericUsageFile)).Sequence;
            var genericUsage = CodonUsage(genericSequence, OpenStretches(genericSequence, standardCode, codonUsageMinLength));
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();
            var modes = new (string Name, Func<string, int, (double X, double Y, double Z)[]> Windows)[]
            {
                ("GeneAligned", (gene, length) => new ZCurveWindowedProfile(gene.Substring(0, length)).Windows.ToArray()),
                ("GeneWhole", (gene, _) => new ZCurveWindowedProfile(gene, gene.Length).Windows.ToArray()),
                ("GenePos12", (gene, length) => new ZCurveWindowedProfile(CodonPositions12(gene.Substring(0, length)), pos12WindowSize).Windows.ToArray()),
                ("GeneWholePos12", (gene, _) => new ZCurveWindowedProfile(CodonPositions12(gene), gene.Length / 3 * 2).Windows.ToArray()),
            };

            var rows = new List<GeneQueryResult>();
            using (var raw = File.CreateText("R6c_PerGene_Raw.csv"))
            {
                raw.WriteLine("Group,FileName,Mode,MinGeneLengthBp,GeneIndex,Strand,GeneLengthBp,Variant,DifferingFraction," +
                              "Candidates,TrueRank,Top1,TrueScore,BestFalseScore");

                foreach (var (group, fileName) in Corpus())
                {
                    var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence;
                    var reverse = ReverseComplement(forward);
                    var code = Path.GetFileName(fileName).StartsWith("mycoplasma_", StringComparison.Ordinal) ? mycoplasmaCode : standardCode;
                    var synonyms = Synonyms(code);
                    var ownUsage = CodonUsage(forward, OpenStretches(forward, code, codonUsageMinLength));

                    foreach (var length in queryLengths)
                    {
                        var genes = OpenStretches(forward, code, length).Select(s => (Reverse: false, Sequence: forward.Substring(s.Start, s.End - s.Start)))
                            .Concat(OpenStretches(reverse, code, length).Select(s => (Reverse: true, Sequence: reverse.Substring(s.Start, s.End - s.Start))))
                            .ToList();
                        var indexes = modes.Select(m => genes.Select(g => m.Windows(g.Sequence, length)).ToList()).ToArray();

                        // Same seed and picker as GeneAlignedWindowsRescueBackTranslation, so the same genes are queried.
                        var seed = StableSeed(fileName) ^ length ^ 0x5bd1e995;
                        var picker = new Random(seed);
                        var picked = Enumerable.Range(0, genes.Count).OrderBy(_ => picker.Next()).Take(regionsPerGenome).OrderBy(i => i).ToList();

                        foreach (var geneIndex in picked)
                        {
                            var gene = genes[geneIndex];
                            var protein = Translate(gene.Sequence, code);
                            var random = new Random(seed * 31 + geneIndex);

                            foreach (var variant in variants)
                            {
                                var query = BuildVariant(variant, gene.Sequence, protein, synonyms, ownUsage, genericUsage, uniformUsage, random);
                                var differing = DifferingFraction(gene.Sequence, query);

                                for (int m = 0; m < modes.Length; m++)
                                {
                                    var result = RankAmongGenes(modes[m].Windows(query, length), indexes[m], geneIndex);
                                    rows.Add(new GeneQueryResult(group, modes[m].Name, length, variant, differing, result));
                                    raw.WriteLine(string.Join(",",
                                        group, fileName, modes[m].Name, length, geneIndex, gene.Reverse ? "-" : "+", gene.Sequence.Length, variant,
                                        Format(differing), result.Candidates, result.TrueRank,
                                        (result.TrueRank == 1).ToString().ToLowerInvariant(),
                                        Format(result.TrueScore), Format(result.BestFalseScore)));
                                }
                            }
                        }
                    }
                }
            }

            using (var summary = File.CreateText("R6c_PerGene_Summary.csv"))
            {
                summary.WriteLine("Group,Mode,MinGeneLengthBp,Variant,Queries,Top1,MedianTrueRank,MedianCandidates,MeanDifferingFraction");
                foreach (var cell in rows
                    .GroupBy(r => (r.Group, r.Mode, r.Length, r.Variant))
                    .OrderBy(g => g.Key.Group).ThenBy(g => Array.FindIndex(modes, m => m.Name == g.Key.Mode)).ThenBy(g => g.Key.Length)
                    .ThenBy(g => Array.IndexOf(variants, g.Key.Variant)))
                {
                    var line = string.Join(",",
                        cell.Key.Group, cell.Key.Mode, cell.Key.Length, cell.Key.Variant, cell.Count(),
                        Format(cell.Average(r => r.Result.TrueRank == 1 ? 1.0 : 0.0)),
                        Format(Median(cell.Select(r => (double)r.Result.TrueRank))),
                        Format(Median(cell.Select(r => (double)r.Result.Candidates))),
                        Format(cell.Average(r => r.DifferingFraction)));
                    summary.WriteLine(line);
                    TestContext.Out.WriteLine(line);
                }
            }

            foreach (var cell in rows.Where(r => r.Variant == "Dna").GroupBy(r => (r.Group, r.Mode)))
                Assert.That(cell.Average(r => r.Result.TrueRank == 1 ? 1.0 : 0.0), Is.GreaterThanOrEqualTo(0.9), cell.Key.ToString());
        }

        /// <summary>
        /// R6d: does R6c's GenePos12 signal survive the real job? Part 1 pools every gene of every genome
        /// in a group, so a query competes with about 30x more genes than in R6c. HomologTop1 also counts
        /// the top hit as right when it is another genome's copy of the protein (≥ 0.3 of the query's
        /// amino-acid 5-mers), since "organisms containing a protein" wants exactly those. Part 2 searches
        /// E. coli K-12 proteins against four related genomes and asks whether the ortholog ranks first among
        /// that genome's genes; orthologs are picked alignment-free by shared 5-mers. Genes start at their
        /// first in-frame ATG here (the longest-ORF convention), so a protein starting at Met lines up with
        /// the stored gene's first window.
        /// </summary>
        [Test]
        public void PooledAndHomologRetrievalWithPos12()
        {
            var genericSequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genericUsageFile)).Sequence;
            var genericUsage = CodonUsage(genericSequence, OpenStretches(genericSequence, standardCode, codonUsageMinLength));
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();
            var pooledLengths = new[] { 999, 1500, 3000 };

            var pooledRows = new List<(string Group, int Length, string Variant, int PoolSize, int TrueRank, bool HomologTop1, double Differing)>();
            using (var raw = File.CreateText("R6d_Pooled_Raw.csv"))
            {
                raw.WriteLine("Group,FileName,MinGeneLengthBp,GeneIndex,Strand,GeneLengthBp,Variant,DifferingFraction,PoolSize,TrueRank,Top1,HomologTop1,TopHitFileName,TopHitKmerShare");

                foreach (var group in new[] { "Bacteria", "Virus" })
                {
                    var fileNames = Corpus().Where(c => c.Group == group).Select(c => c.FileName).ToList();
                    var pools = pooledLengths.ToDictionary(l => l, _ => new List<PooledGene>());
                    var queries = pooledLengths.ToDictionary(l => l, _ => new List<PooledQuery>());

                    for (int genome = 0; genome < fileNames.Count; genome++)
                    {
                        var fileName = fileNames[genome];
                        var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence;
                        var reverse = ReverseComplement(forward);
                        var code = Path.GetFileName(fileName).StartsWith("mycoplasma_", StringComparison.Ordinal) ? mycoplasmaCode : standardCode;
                        var synonyms = Synonyms(code);
                        var ownUsage = CodonUsage(forward, OpenStretches(forward, code, codonUsageMinLength));

                        foreach (var length in pooledLengths)
                        {
                            var genes = AtgGenes(forward, code, length).Select(g => (Reverse: false, Sequence: forward.Substring(g.Start, g.End - g.Start)))
                                .Concat(AtgGenes(reverse, code, length).Select(g => (Reverse: true, Sequence: reverse.Substring(g.Start, g.End - g.Start))))
                                .ToList();
                            int firstPoolIndex = pools[length].Count;
                            foreach (var gene in genes)
                                pools[length].Add(new PooledGene(genome, Translate(gene.Sequence, code), Pos12Prefix(gene.Sequence, length)));

                            var seed = StableSeed(fileName) ^ length ^ 0x2545f491;
                            var picker = new Random(seed);
                            foreach (var geneIndex in Enumerable.Range(0, genes.Count).OrderBy(_ => picker.Next()).Take(regionsPerGenome).OrderBy(i => i))
                            {
                                var gene = genes[geneIndex];
                                var protein = pools[length][firstPoolIndex + geneIndex].Protein;
                                var random = new Random(seed * 31 + geneIndex);
                                foreach (var variant in variants)
                                {
                                    var query = BuildVariant(variant, gene.Sequence, protein, synonyms, ownUsage, genericUsage, uniformUsage, random);
                                    queries[length].Add(new PooledQuery(fileName, geneIndex, gene.Reverse, gene.Sequence.Length, variant,
                                        DifferingFraction(gene.Sequence, query), firstPoolIndex + geneIndex, Pos12Prefix(query, length)));
                                }
                            }
                        }
                    }

                    foreach (var length in pooledLengths)
                    {
                        var pool = pools[length];
                        foreach (var query in queries[length])
                        {
                            var trueScore = EqualLengthScore(query.Windows, pool[query.PoolIndex].Windows);
                            int higher = 0, top = query.PoolIndex;
                            var topScore = trueScore;
                            for (int i = 0; i < pool.Count; i++)
                            {
                                if (i == query.PoolIndex) continue;
                                var score = EqualLengthScore(query.Windows, pool[i].Windows);
                                if (score > trueScore) higher++;
                                if (score > topScore)
                                {
                                    topScore = score;
                                    top = i;
                                }
                            }

                            var share = top == query.PoolIndex ? 1.0 : KmerShare(pool[query.PoolIndex].Protein, pool[top].Protein);
                            var homologTop1 = share >= orthologMinShare;
                            pooledRows.Add((group, length, query.Variant, pool.Count, higher + 1, homologTop1, query.DifferingFraction));
                            raw.WriteLine(string.Join(",",
                                group, query.FileName, length, query.GeneIndex, query.Reverse ? "-" : "+", query.GeneLength, query.Variant,
                                Format(query.DifferingFraction), pool.Count, higher + 1, (higher == 0).ToString().ToLowerInvariant(),
                                homologTop1.ToString().ToLowerInvariant(), fileNames[pool[top].Genome], Format(share)));
                        }
                    }
                }
            }

            using (var summary = File.CreateText("R6d_Pooled_Summary.csv"))
            {
                summary.WriteLine("Group,MinGeneLengthBp,Variant,Queries,PoolSize,Top1,Top10,HomologTop1,MedianTrueRank,MeanDifferingFraction");
                foreach (var cell in pooledRows
                    .GroupBy(r => (r.Group, r.Length, r.Variant))
                    .OrderBy(g => g.Key.Group).ThenBy(g => g.Key.Length).ThenBy(g => Array.IndexOf(variants, g.Key.Variant)))
                {
                    var line = string.Join(",",
                        cell.Key.Group, cell.Key.Length, cell.Key.Variant, cell.Count(), cell.First().PoolSize,
                        Format(cell.Average(r => r.TrueRank == 1 ? 1.0 : 0.0)),
                        Format(cell.Average(r => r.TrueRank <= 10 ? 1.0 : 0.0)),
                        Format(cell.Average(r => r.HomologTop1 ? 1.0 : 0.0)),
                        Format(Median(cell.Select(r => (double)r.TrueRank))),
                        Format(cell.Average(r => r.Differing)));
                    summary.WriteLine(line);
                    TestContext.Out.WriteLine(line);
                }
            }

            HomologRetrieval(genericUsage, uniformUsage);

            foreach (var cell in pooledRows.Where(r => r.Variant == "Dna").GroupBy(r => (r.Group, r.Length)))
                Assert.That(cell.Average(r => r.TrueRank == 1 ? 1.0 : 0.0), Is.GreaterThanOrEqualTo(0.9), cell.Key.ToString());
        }

        private static void HomologRetrieval(long[] genericUsage, long[] uniformUsage)
        {
            var homologVariants = new[] { "Dna", "MostFrequentGeneric", "UniformRandom" };
            var querySequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genericUsageFile)).Sequence;
            var queryReverse = ReverseComplement(querySequence);
            var synonyms = Synonyms(standardCode);
            var picker = new Random(StableSeed(genericUsageFile) ^ 0x68e31da4);
            var queryGenes = AtgGenes(querySequence, standardCode, homologQueryLength).Select(g => querySequence.Substring(g.Start, g.End - g.Start))
                .Concat(AtgGenes(queryReverse, standardCode, homologQueryLength).Select(g => queryReverse.Substring(g.Start, g.End - g.Start)))
                .OrderBy(_ => picker.Next())
                .Take(homologQueryCount)
                .ToList();

            var rows = new List<(string Target, string Variant, string Bin, int TrueRank, int Candidates)>();
            var skipped = new Dictionary<string, (int NoOrtholog, int TooShort)>();
            using var raw = File.CreateText("R6d_Homolog_Raw.csv");
            raw.WriteLine("Target,QueryIndex,QueryGeneLengthBp,OrthologGeneLengthBp,KmerShare,Variant,Candidates,TrueRank,Top1");

            foreach (var target in homologTargets)
            {
                var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(target)).Sequence;
                var reverse = ReverseComplement(forward);
                var genes = AtgGenes(forward, standardCode, orthologSearchMinLength).Select(g => forward.Substring(g.Start, g.End - g.Start))
                    .Concat(AtgGenes(reverse, standardCode, orthologSearchMinLength).Select(g => reverse.Substring(g.Start, g.End - g.Start)))
                    .ToList();
                var proteins = genes.Select(g => Translate(g, standardCode)).ToList();
                var kmerIndex = new Dictionary<int, List<int>>();
                for (int gene = 0; gene < proteins.Count; gene++)
                    foreach (var kmer in PeptideKmers(proteins[gene]))
                    {
                        if (!kmerIndex.TryGetValue(kmer, out var list)) kmerIndex[kmer] = list = new List<int>();
                        list.Add(gene);
                    }

                // Only genes long enough for a homologQueryLength prefix are candidates; slotOf maps a gene to its slot.
                var candidates = Enumerable.Range(0, genes.Count).Where(g => genes[g].Length >= homologQueryLength).ToList();
                var slotOf = candidates.Select((gene, slot) => (gene, slot)).ToDictionary(p => p.gene, p => p.slot);
                var candidateWindows = candidates.Select(g => Pos12Prefix(genes[g], homologQueryLength)).ToList();
                int noOrtholog = 0, tooShort = 0;

                for (int queryIndex = 0; queryIndex < queryGenes.Count; queryIndex++)
                {
                    var queryGene = queryGenes[queryIndex];
                    var protein = Translate(queryGene, standardCode);
                    var queryKmers = PeptideKmers(protein);
                    var counts = new Dictionary<int, int>();
                    foreach (var kmer in queryKmers)
                        if (kmerIndex.TryGetValue(kmer, out var list))
                            foreach (var gene in list)
                                counts[gene] = counts.GetValueOrDefault(gene) + 1;

                    if (counts.Count == 0)
                    {
                        noOrtholog++;
                        continue;
                    }
                    var best = counts.MaxBy(c => c.Value);
                    var ortholog = best.Key;
                    var share = best.Value / (double)queryKmers.Count;
                    if (share < orthologMinShare)
                    {
                        noOrtholog++;
                        continue;
                    }
                    if (!slotOf.TryGetValue(ortholog, out var trueSlot))
                    {
                        tooShort++;
                        continue;
                    }

                    var bin = share >= 0.8 ? "0.8-1.0" : share >= 0.5 ? "0.5-0.8" : "0.3-0.5";
                    var random = new Random(StableSeed(target) * 31 + queryIndex);
                    foreach (var variant in homologVariants)
                    {
                        var query = BuildVariant(variant, queryGene, protein, synonyms, genericUsage, genericUsage, uniformUsage, random);
                        var queryWindows = Pos12Prefix(query, homologQueryLength);
                        var trueScore = EqualLengthScore(queryWindows, candidateWindows[trueSlot]);
                        int higher = 0;
                        for (int slot = 0; slot < candidateWindows.Count; slot++)
                            if (slot != trueSlot && EqualLengthScore(queryWindows, candidateWindows[slot]) > trueScore) higher++;

                        rows.Add((target, variant, bin, higher + 1, candidateWindows.Count));
                        raw.WriteLine(string.Join(",",
                            target, queryIndex, queryGene.Length, genes[ortholog].Length, Format(share), variant,
                            candidateWindows.Count, higher + 1, (higher == 0).ToString().ToLowerInvariant()));
                    }
                }
                skipped[target] = (noOrtholog, tooShort);
            }

            using var summary = File.CreateText("R6d_Homolog_Summary.csv");
            summary.WriteLine("Target,Variant,KmerShareBin,Queries,Candidates,Top1,Top10,MedianTrueRank,SkippedNoOrtholog,SkippedOrthologTooShort");
            foreach (var cell in rows
                .Concat(rows.Select(r => r with { Bin = "all" }))
                .GroupBy(r => (r.Target, r.Variant, r.Bin))
                .OrderBy(g => Array.IndexOf(homologTargets, g.Key.Target)).ThenBy(g => Array.IndexOf(homologVariants, g.Key.Variant)).ThenBy(g => g.Key.Bin))
            {
                var line = string.Join(",",
                    cell.Key.Target, cell.Key.Variant, cell.Key.Bin, cell.Count(), cell.First().Candidates,
                    Format(cell.Average(r => r.TrueRank == 1 ? 1.0 : 0.0)),
                    Format(cell.Average(r => r.TrueRank <= 10 ? 1.0 : 0.0)),
                    Format(Median(cell.Select(r => (double)r.TrueRank))),
                    skipped[cell.Key.Target].NoOrtholog, skipped[cell.Key.Target].TooShort);
                summary.WriteLine(line);
                TestContext.Out.WriteLine(line);
            }
        }

        // A stop-free stretch trimmed to its first in-frame ATG, kept if at least minLength remains.
        private static IEnumerable<(int Start, int End)> AtgGenes(string strand, string code, int minLength)
        {
            foreach (var (start, end) in OpenStretches(strand, code, minLength))
                for (int i = start; i + minLength <= end; i += 3)
                    if (strand[i] == 'A' && strand[i + 1] == 'T' && strand[i + 2] == 'G')
                    {
                        yield return (i, end);
                        break;
                    }
        }

        private static (double X, double Y, double Z)[] Pos12Prefix(string gene, int length) =>
            new ZCurveWindowedProfile(CodonPositions12(gene.Substring(0, length)), pos12WindowSize).Windows.ToArray();

        // Allocation-free twin of ScoreAllOffsets for two profiles with the same window count.
        private static double EqualLengthScore((double X, double Y, double Z)[] a, (double X, double Y, double Z)[] b)
        {
            double sum = 0;
            for (int j = 0; j < a.Length; j++)
            {
                var dx = a[j].X - b[j].X;
                var dy = a[j].Y - b[j].Y;
                var dz = a[j].Z - b[j].Z;
                sum += Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }
            return 1 - sum / a.Length;
        }

        // Amino-acid 5-mers packed 5 bits per residue; letters are A-Z, so each fits.
        private static HashSet<int> PeptideKmers(string protein)
        {
            var kmers = new HashSet<int>();
            for (int i = 0; i + 5 <= protein.Length; i++)
            {
                int kmer = 0;
                for (int k = 0; k < 5; k++)
                    kmer = (kmer << 5) | (protein[i + k] - 'A');
                kmers.Add(kmer);
            }
            return kmers;
        }

        private static double KmerShare(string query, string other)
        {
            var queryKmers = PeptideKmers(query);
            if (queryKmers.Count == 0) return 0;
            var otherKmers = PeptideKmers(other);
            return queryKmers.Count(otherKmers.Contains) / (double)queryKmers.Count;
        }

        // Every codon's first two bases, in order: the positions a synonymous change almost never touches.
        private static string CodonPositions12(string gene)
        {
            var positions = new char[gene.Length / 3 * 2];
            for (int codon = 0; codon < gene.Length / 3; codon++)
            {
                positions[2 * codon] = gene[3 * codon];
                positions[2 * codon + 1] = gene[3 * codon + 1];
            }
            return new string(positions);
        }

        // Equal window counts, so each gene gets one "1 - mean distance" score and no sliding. Ties with
        // the true gene (duplicated genes) don't count against it.
        private static GeneRank RankAmongGenes((double X, double Y, double Z)[] query, List<(double X, double Y, double Z)[]> index, int trueGene)
        {
            var trueScore = ScoreAllOffsets(query, index[trueGene])[0];
            var bestFalseScore = double.NegativeInfinity;
            int higher = 0;
            for (int gene = 0; gene < index.Count; gene++)
            {
                if (gene == trueGene) continue;
                var score = ScoreAllOffsets(query, index[gene])[0];
                if (score > bestFalseScore) bestFalseScore = score;
                if (score > trueScore) higher++;
            }
            return new GeneRank(index.Count, higher + 1, trueScore, bestFalseScore, false);
        }

        // The query is scanned as given and reverse-complemented, since a protein carries no strand. The
        // true hit is the best offset within one window of the gene's forward-coordinate start, in the
        // orientation that matches the gene's strand; every other offset in either orientation competes.
        private static GeneRank RankOnFixedWindows(string query, (double X, double Y, double Z)[] genomeWindows, int forwardStart, bool reverseGene)
        {
            var asGiven = ScoreAllOffsets(new ZCurveWindowedProfile(query).Windows.ToArray(), genomeWindows);
            var complemented = ScoreAllOffsets(new ZCurveWindowedProfile(ReverseComplement(query)).Windows.ToArray(), genomeWindows);
            var trueScores = reverseGene ? complemented : asGiven;

            var trueScore = double.NegativeInfinity;
            for (int offset = 0; offset < trueScores.Length; offset++)
                if (IsNearTrueStart(offset, forwardStart) && trueScores[offset] > trueScore) trueScore = trueScores[offset];

            var bestFalseScore = double.NegativeInfinity;
            int higher = 0;
            foreach (var (scores, isTrueOrientation) in new[] { (asGiven, !reverseGene), (complemented, reverseGene) })
                for (int offset = 0; offset < scores.Length; offset++)
                {
                    if (isTrueOrientation && IsNearTrueStart(offset, forwardStart)) continue;
                    if (scores[offset] > bestFalseScore) bestFalseScore = scores[offset];
                    if (scores[offset] > trueScore) higher++;
                }

            var bestAsGiven = asGiven.Max();
            var bestComplemented = complemented.Max();
            return new GeneRank(asGiven.Length + complemented.Length, higher + 1, trueScore, bestFalseScore, bestComplemented > bestAsGiven);
        }

        private static bool IsNearTrueStart(int offset, int forwardStart) => Math.Abs(offset * windowSize - forwardStart) <= windowSize;

        private static string ReverseComplement(string sequence)
        {
            var complement = new char[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
            {
                complement[sequence.Length - 1 - i] = sequence[i] switch
                {
                    'A' => 'T',
                    'T' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    _ => 'N',
                };
            }
            return new string(complement);
        }

        [Test]
        public void AllOffsetScoresAgreeWithBestAlignment()
        {
            var reference = new ZCurveWindowedProfile(Generator.GetRandomDNA(200_000, 61));
            var referenceWindows = reference.Windows.ToArray();
            foreach (var length in queryLengths)
            for (int seed = 0; seed < 5; seed++)
            {
                var query = new ZCurveWindowedProfile(Generator.GetRandomDNA(length, 1000 + seed));
                var scores = ScoreAllOffsets(query.Windows.ToArray(), referenceWindows);
                var bestByScan = Array.IndexOf(scores, scores.Max());
                var (bestScore, bestOffset) = ZCurveWindowedProfile.BestAlignment(query, reference);

                Assert.That(bestByScan, Is.EqualTo(bestOffset));
                Assert.That(scores[bestByScan], Is.EqualTo(bestScore).Within(1e-5));
            }
        }

        private static IEnumerable<(string Group, string FileName)> Corpus()
        {
            var viral = ViralTaxonomy.Genomes.Select(g => g.FileName).ToList();
            foreach (var path in GenomeCatalog.AllFiles)
            {
                var fileName = Path.GetFileName(path);
                if (!viral.Contains(fileName)) yield return ("Bacteria", fileName);
            }
            foreach (var fileName in viral)
                yield return ("Virus", fileName);
        }

        // string.GetHashCode is randomized per process, so seeds come from FNV-1a instead.
        private static int StableSeed(string text)
        {
            unchecked
            {
                uint hash = 2166136261;
                foreach (var c in text)
                    hash = (hash ^ c) * 16777619;
                return (int)hash;
            }
        }

        private static int CodonIndex(string sequence, int i)
        {
            int first = codonBases.IndexOf(sequence[i]);
            int second = codonBases.IndexOf(sequence[i + 1]);
            int third = codonBases.IndexOf(sequence[i + 2]);
            return first < 0 || second < 0 || third < 0 ? -1 : 16 * first + 4 * second + third;
        }

        private static string CodonAt(int index) =>
            new string(new[] { codonBases[index / 16], codonBases[index / 4 % 4], codonBases[index % 4] });

        // Forward-strand stop-free stretches, per frame, as [Start, End) in bases. A non-ACGT codon ends a
        // stretch like a stop does, so every stretch is ACGT-only. At these lengths a stop-free stretch is
        // essentially always a gene, since random sequence hits a stop about every 21 codons.
        private static List<(int Start, int End)> OpenStretches(string sequence, string code, int minLength)
        {
            var stretches = new List<(int Start, int End)>();
            for (int frame = 0; frame < 3; frame++)
            {
                int start = frame;
                int i = frame;
                for (; i + 3 <= sequence.Length; i += 3)
                {
                    int codon = CodonIndex(sequence, i);
                    if (codon >= 0 && code[codon] != '*') continue;
                    if (i - start >= minLength) stretches.Add((start, i));
                    start = i + 3;
                }
                if (i - start >= minLength) stretches.Add((start, i));
            }
            return stretches;
        }

        // The first position in each stretch that is both a window boundary and in frame; one of any
        // three consecutive window boundaries is in frame, because 500 mod 3 = 2.
        private static IEnumerable<int> AlignedStarts(List<(int Start, int End)> stretches, int length)
        {
            foreach (var (start, end) in stretches)
            {
                int aligned = (start + windowSize - 1) / windowSize * windowSize;
                while ((aligned - start) % 3 != 0) aligned += windowSize;
                if (aligned + length <= end) yield return aligned;
            }
        }

        private static long[] CodonUsage(string sequence, List<(int Start, int End)> stretches)
        {
            var counts = new long[64];
            foreach (var (start, end) in stretches)
                for (int i = start; i < end; i += 3)
                    counts[CodonIndex(sequence, i)]++;
            return counts;
        }

        private static Dictionary<char, int[]> Synonyms(string code) =>
            Enumerable.Range(0, 64)
                .Where(i => code[i] != '*')
                .GroupBy(i => code[i])
                .ToDictionary(g => g.Key, g => g.ToArray());

        private static string Translate(string region, string code)
        {
            var protein = new char[region.Length / 3];
            for (int i = 0; i < protein.Length; i++)
                protein[i] = code[CodonIndex(region, 3 * i)];
            return new string(protein);
        }

        private static string BuildVariant(
            string variant, string original, string protein, Dictionary<char, int[]> synonyms,
            long[] ownUsage, long[] genericUsage, long[] uniformUsage, Random random) => variant switch
        {
            "Dna" => original,
            "MostFrequentOwn" => BackTranslateMostFrequent(protein, synonyms, ownUsage),
            "WeightedRandomOwn" => BackTranslateWeighted(protein, synonyms, ownUsage, random),
            "MostFrequentGeneric" => BackTranslateMostFrequent(protein, synonyms, genericUsage),
            "UniformRandom" => BackTranslateWeighted(protein, synonyms, uniformUsage, random),
            "RandomDna" => Generator.GetRandomDNA(original.Length, random),
            _ => throw new InvalidOperationException(variant),
        };

        private static string BackTranslateMostFrequent(string protein, Dictionary<char, int[]> synonyms, long[] usage)
        {
            var dna = new System.Text.StringBuilder(protein.Length * 3);
            foreach (var aminoAcid in protein)
            {
                var codons = synonyms[aminoAcid];
                int best = codons[0];
                foreach (var codon in codons)
                    if (usage[codon] > usage[best]) best = codon;
                dna.Append(CodonAt(best));
            }
            return dna.ToString();
        }

        // Draws each codon in proportion to usage among the amino acid's synonyms; falls back to uniform
        // when the table has never seen any of them.
        private static string BackTranslateWeighted(string protein, Dictionary<char, int[]> synonyms, long[] usage, Random random)
        {
            var dna = new System.Text.StringBuilder(protein.Length * 3);
            foreach (var aminoAcid in protein)
            {
                var codons = synonyms[aminoAcid];
                long total = codons.Sum(c => usage[c]);
                int chosen = codons[random.Next(codons.Length)];
                if (total > 0)
                {
                    long draw = random.NextInt64(total);
                    foreach (var codon in codons)
                    {
                        if (draw < usage[codon])
                        {
                            chosen = codon;
                            break;
                        }
                        draw -= usage[codon];
                    }
                }
                dna.Append(CodonAt(chosen));
            }
            return dna.ToString();
        }

        private static double DifferingFraction(string original, string variant)
        {
            int differing = 0;
            for (int i = 0; i < original.Length; i++)
                if (original[i] != variant[i]) differing++;
            return differing / (double)original.Length;
        }

        private static Localization Localize(
            string query, ZCurveWindowedProfile reference, (double X, double Y, double Z)[] referenceWindows, int trueOffset)
        {
            var queryProfile = new ZCurveWindowedProfile(query);
            var (_, bestOffset) = ZCurveWindowedProfile.BestAlignment(queryProfile, reference);
            var scores = ScoreAllOffsets(queryProfile.Windows.ToArray(), referenceWindows);

            var trueScore = scores[trueOffset];
            var bestFalseScore = double.NegativeInfinity;
            int higher = 0;
            for (int offset = 0; offset < scores.Length; offset++)
            {
                if (offset == trueOffset) continue;
                if (scores[offset] > bestFalseScore) bestFalseScore = scores[offset];
                if (scores[offset] > trueScore) higher++;
            }

            return new Localization(bestOffset, bestOffset == trueOffset, Math.Abs(bestOffset - trueOffset) <= 1,
                higher + 1, trueScore, bestFalseScore);
        }

        // Same "1 - mean Euclidean distance" as BestAlignment, but keeps every offset's score so the true
        // offset's rank and margin can be read off. AllOffsetScoresAgreeWithBestAlignment checks the two agree.
        private static double[] ScoreAllOffsets((double X, double Y, double Z)[] query, (double X, double Y, double Z)[] reference)
        {
            var scores = new double[reference.Length - query.Length + 1];
            for (int offset = 0; offset < scores.Length; offset++)
            {
                double sum = 0;
                for (int j = 0; j < query.Length; j++)
                {
                    var dx = query[j].X - reference[offset + j].X;
                    var dy = query[j].Y - reference[offset + j].Y;
                    var dz = query[j].Z - reference[offset + j].Z;
                    sum += Math.Sqrt(dx * dx + dy * dy + dz * dz);
                }
                scores[offset] = 1 - sum / query.Length;
            }
            return scores;
        }

        private static double Median(IEnumerable<double> values)
        {
            var sorted = values.OrderBy(v => v).ToArray();
            int middle = sorted.Length / 2;
            return sorted.Length % 2 == 1 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
        }

        private static string Format(double value) => value.ToString("G6", CultureInfo.InvariantCulture);

        private readonly record struct Localization(
            int BestOffset, bool ExactHit, bool WithinOneWindow, int TrueRank, double TrueScore, double BestFalseScore);

        private sealed record QueryResult(string Group, int Length, string Variant, double DifferingFraction, Localization Result);

        private readonly record struct GeneRank(
            int Candidates, int TrueRank, double TrueScore, double BestFalseScore, bool MatchedReverseComplement);

        private sealed record PooledGene(int Genome, string Protein, (double X, double Y, double Z)[] Windows);

        private sealed record PooledQuery(
            string FileName, int GeneIndex, bool Reverse, int GeneLength, string Variant, double DifferingFraction,
            int PoolIndex, (double X, double Y, double Z)[] Windows);

        private sealed record GeneQueryResult(string Group, string Mode, int Length, string Variant, double DifferingFraction, GeneRank Result);
    }
}
