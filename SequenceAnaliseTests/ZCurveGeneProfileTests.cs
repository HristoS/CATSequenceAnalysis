using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class ZCurveGeneProfileTests
    {
        private const string codonBases = "TCAG";

        // A start codon, `codons` random sense codons, and a stop: a gene in coding orientation.
        private static string RandomGene(int codons, int seed)
        {
            var rng = new Random(seed);
            var sense = Enumerable.Range(0, 64).Where(c => ZCurveGeneProfile.standardCode[c] != '*').ToArray();
            var gene = new System.Text.StringBuilder("ATG");
            for (int i = 0; i < codons; i++)
                gene.Append(Codon(sense[rng.Next(sense.Length)]));
            return gene.Append("TAA").ToString();
        }

        private static string Codon(int index) =>
            new string(new[] { codonBases[index / 16], codonBases[index / 4 % 4], codonBases[index % 4] });

        private static string Translate(string gene) =>
            new string(Enumerable.Range(0, gene.Length / 3)
                .Select(c => ZCurveGeneProfile.standardCode[16 * codonBases.IndexOf(gene[3 * c]) + 4 * codonBases.IndexOf(gene[3 * c + 1]) + codonBases.IndexOf(gene[3 * c + 2])])
                .Where(a => a != '*').ToArray());

        private static string ReverseComplement(string s) =>
            new string(s.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' }).ToArray());

        [Test]
        public void MeansEqualTheWindowedProfileOfCodonPositionsOneAndTwo()
        {
            var gene = RandomGene(600, 1);
            var positions12 = new string(Enumerable.Range(0, gene.Length / 3 - 1).SelectMany(c => new[] { gene[3 * c], gene[3 * c + 1] }).ToArray());
            var windows = new ZCurveWindowedProfile(positions12, ZCurveGeneProfile.WindowBases).Windows;
            var query = ZCurveGeneProfile.QueryFromGene(gene);

            Assert.That(query.Length / ZCurveGeneProfile.FeatureCount, Is.EqualTo(windows.Count));
            for (int w = 0; w < windows.Count; w++)
            {
                Assert.That(query[w * 9], Is.EqualTo(windows[w].X).Within(1e-6));
                Assert.That(query[w * 9 + 1], Is.EqualTo(windows[w].Y).Within(1e-6));
                Assert.That(query[w * 9 + 2], Is.EqualTo(windows[w].Z).Within(1e-6));
            }
        }

        [Test]
        public void TrendSeesOrderThatTheMeanCannot()
        {
            // 100 aa, so 5 windows of 40 positions-1-and-2 bases (20 codons): in each, purines (GA) in the first half then pyrimidines
            // (CT), or the reverse. The means are equal, while the trends are about -sqrt(3) / 2 and +sqrt(3) / 2.
            var purinesFirst = string.Concat(Enumerable.Repeat(string.Concat(Enumerable.Repeat("GAA", 10)) + string.Concat(Enumerable.Repeat("CTA", 10)), 5)) + "TAA";
            var pyrimidinesFirst = string.Concat(Enumerable.Repeat(string.Concat(Enumerable.Repeat("CTA", 10)) + string.Concat(Enumerable.Repeat("GAA", 10)), 5)) + "TAA";
            var a = ZCurveGeneProfile.QueryFromGene(purinesFirst);
            var b = ZCurveGeneProfile.QueryFromGene(pyrimidinesFirst);

            Assert.That(a.Length / ZCurveGeneProfile.FeatureCount, Is.EqualTo(ZCurveGeneProfile.SplitWindows));
            Assert.That(a[0], Is.EqualTo(b[0]).Within(1e-12), "the same mean x");
            Assert.That(a[3], Is.LessThan(-0.8), "purines (+x) early: falling x");
            Assert.That(b[3], Is.GreaterThan(0.8), "purines late: rising x");
        }

        [Test]
        public void ShortGeneMeansEqualFiveEqualWindowsOfItsCodonPositionsOneAndTwo()
        {
            var gene = RandomGene(149, 11);
            var positions12 = new string(Enumerable.Range(0, gene.Length / 3 - 1).SelectMany(c => new[] { gene[3 * c], gene[3 * c + 1] }).ToArray());
            var windows = new ZCurveWindowedProfile(positions12, positions12.Length / ZCurveGeneProfile.SplitWindows).Windows;
            var query = ZCurveGeneProfile.QueryFromGene(gene);

            Assert.That((windows.Count, query.Length / ZCurveGeneProfile.FeatureCount), Is.EqualTo((5, 5)), "150 aa: 5 windows of 60");
            for (int w = 0; w < windows.Count; w++)
                Assert.That((query[w * 9], query[w * 9 + 1], query[w * 9 + 2]), Is.EqualTo(windows[w]).Using<(double, double, double)>(
                    (p, q) => Math.Abs(p.Item1 - q.Item1) < 1e-6 && Math.Abs(p.Item2 - q.Item2) < 1e-6 && Math.Abs(p.Item3 - q.Item3) < 1e-6));
        }

        [Test]
        public void WindowsOfGiveEachGenesStoredWindowsAndMeanWindowAveragesThemAll()
        {
            var (short1, long1) = (RandomGene(120, 21), RandomGene(700, 22));
            var record = short1 + long1;
            var profile = ZCurveGeneProfile.FromRecord(record, new[] { (0, short1.Length, false, 1), (short1.Length, record.Length, false, 2) });

            var stored = new[] { short1, long1 }.Select(g => ZCurveGeneProfile.QueryFromGene(g).Select(v => (double)(Half)v).ToArray()).ToArray();
            Assert.That(profile.WindowsOf(0), Is.EqualTo(stored[0]));
            Assert.That(profile.WindowsOf(1), Is.EqualTo(stored[1]));
            Assert.That((ZCurveGeneProfile.WindowOf(120), ZCurveGeneProfile.WindowOf(700), ZCurveGeneProfile.WindowOf(49)), Is.EqualTo((48, 333, 0)));

            var all = stored.SelectMany(s => s).ToArray();
            var expected = Enumerable.Range(0, ZCurveGeneProfile.FeatureCount)
                .Select(f => Enumerable.Range(0, all.Length / ZCurveGeneProfile.FeatureCount).Average(w => all[w * ZCurveGeneProfile.FeatureCount + f]));
            Assert.That(profile.MeanWindow(), Is.EqualTo(expected).Within(1e-12));
        }

        [Test]
        public void FinalStopCodonIsLeftOut()
        {
            var gene = RandomGene(700, 2);
            Assert.That(ZCurveGeneProfile.QueryFromGene(gene), Is.EqualTo(ZCurveGeneProfile.QueryFromGene(gene[..^3])));
        }

        [Test]
        public void ReverseStrandGeneIsReadInCodingOrientation()
        {
            var gene = RandomGene(800, 3);
            var flank = Generator.GetRandomDNA(1000, 4);
            var record = flank + ReverseComplement(gene) + flank;
            var profile = ZCurveGeneProfile.FromRecord(record, new[] { (1000, 1000 + gene.Length, true, 42) });

            Assert.That(profile.Genes.Single(), Is.EqualTo(new ZCurveGeneProfile.Gene(1000, 1000 + gene.Length, true, 42, 0, 4, 801)));
            Assert.That(profile.Score(ZCurveGeneProfile.QueryFromGene(gene), 0), Is.GreaterThan(0.998), "only half-precision rounding remains");
        }

        [Test]
        public void ProteinQueryWithOneCodonPerAminoAcidEqualsItsBackTranslation()
        {
            var gene = RandomGene(900, 5);
            var protein = Translate(gene);
            var usage = new long[64];
            var chosen = Enumerable.Range(0, 64).Where(c => ZCurveGeneProfile.standardCode[c] != '*')
                .GroupBy(c => ZCurveGeneProfile.standardCode[c]).ToDictionary(g => g.Key, g => g.Last());
            foreach (var codon in chosen.Values) usage[codon] = 1;
            var backTranslated = string.Concat(protein.Select(a => Codon(chosen[a])));

            var fromProtein = ZCurveGeneProfile.QueryFromProtein(protein, usage);
            var fromDna = ZCurveGeneProfile.QueryFromGene(backTranslated);
            Assert.That(fromProtein.Length, Is.EqualTo(fromDna.Length));
            for (int i = 0; i < fromDna.Length; i++)
                Assert.That(fromProtein[i], Is.EqualTo(fromDna[i]).Within(1e-12));
        }

        [Test]
        public void ProteinQueryFindsItsOwnGeneAmongRandomGenes()
        {
            var genes = Enumerable.Range(0, 200).Select(i => RandomGene(500, 100 + i)).ToArray();
            var record = string.Concat(genes);
            int offset = 0;
            var coordinates = genes.Select(g => { var c = (offset, offset + g.Length, false, ZCurveGeneProfile.NoFamily); offset += g.Length; return c; }).ToArray();
            var profile = ZCurveGeneProfile.FromRecord(record, coordinates);

            for (int q = 0; q < 20; q++)
            {
                var query = ZCurveGeneProfile.QueryFromProtein(Translate(genes[q]));
                var best = Enumerable.Range(0, genes.Length).MaxBy(g => profile.Score(query, g));
                Assert.That(best, Is.EqualTo(q), $"query {q}");
            }
        }

        [Test]
        public void GenesUnder50AaAreDroppedAndShortOnesGetFiveWindows()
        {
            var tiny = RandomGene(40, 6);
            var shortGene = RandomGene(150, 12);
            var record = shortGene + tiny;
            var profile = ZCurveGeneProfile.FromRecord(record, new[] { (shortGene.Length, record.Length, false, 7), (0, shortGene.Length, false, 8) });

            Assert.That(profile.Genes.Single(), Is.EqualTo(new ZCurveGeneProfile.Gene(0, shortGene.Length, false, 8, 0, 5, 151)));
            var query = ZCurveGeneProfile.QueryFromProtein(Translate(shortGene));
            var others = Enumerable.Range(0, 30).Select(i => RandomGene(150, 200 + i)).ToArray();
            var rivals = ZCurveGeneProfile.FromRecord(string.Concat(others),
                others.Select((g, i) => (i * g.Length, (i + 1) * g.Length, false, ZCurveGeneProfile.NoFamily)));
            Assert.That(profile.Score(query, 0), Is.GreaterThan(Enumerable.Range(0, 30).Max(g => rivals.Score(query, g))),
                "its own gene beats random genes of the same length");
            Assert.That(ZCurveGeneProfile.QueryFromGene(tiny), Is.Empty, "41 aa: no window");
        }

        [Test]
        public void SelectAndCombineSplitAndRejoinARecordsGenes()
        {
            var genes = new[] { RandomGene(700, 13), RandomGene(100, 14), RandomGene(300, 15) };
            var record = string.Concat(genes);
            int a = genes[0].Length, b = a + genes[1].Length;
            var profile = ZCurveGeneProfile.FromRecord(record, new[] { (0, a, false, 1), (a, b, true, 2), (b, record.Length, false, 3) });

            var parts = new[] { profile.Select(new[] { 1 }), profile.Select(new[] { 2, 0 }) };
            var joined = ZCurveGeneProfile.Combine(parts.Select(p => ZCurveGeneProfile.FromBinary(p.ToBinary())));

            Assert.That(parts.Select(p => p.Genes.Select(g => g.Family)), Is.EqualTo(new[] { new[] { 2 }, new[] { 1, 3 } }));
            Assert.That(joined.Genes, Is.EqualTo(profile.Genes));
            var query = ZCurveGeneProfile.QueryFromGene(genes[2]);
            Assert.That(joined.Score(query, 2), Is.EqualTo(profile.Score(query, 2)));
        }

        [Test]
        public void BinaryRoundTripKeepsGenesAndScores()
        {
            var genes = new[] { RandomGene(700, 7), RandomGene(100, 8), RandomGene(1200, 9) };
            var record = genes[0] + ReverseComplement(genes[1]) + genes[2];
            int a = genes[0].Length, b = a + genes[1].Length;
            var profile = ZCurveGeneProfile.FromRecord(record, new[] { (0, a, false, 11), (a, b, true, ZCurveGeneProfile.NoFamily), (b, record.Length, false, 13) });
            var copy = ZCurveGeneProfile.FromBinary(profile.ToBinary());

            Assert.That(copy.Genes, Is.EqualTo(profile.Genes));
            Assert.That(copy.WindowCount, Is.EqualTo(profile.WindowCount));
            var query = ZCurveGeneProfile.QueryFromGene(genes[2]);
            Assert.That(copy.Score(query, 2), Is.EqualTo(profile.Score(query, 2)));
            Assert.That(copy.Score(query, 0), Is.EqualTo(profile.Score(query, 0)));
        }

        // ScoreRow is what protein search runs on stored rows; it must give FromBinary + Score's numbers exactly, for both formats, every
        // query shape (long windows, split short windows) and with the candidate rule applied.
        [Test]
        public void ScoreRowEqualsDecodingThenScoringBitForBit()
        {
            var genes = new[] { RandomGene(700, 20), RandomGene(120, 21), RandomGene(1200, 22), RandomGene(60, 23), RandomGene(130, 24) };
            var unspliced = ZCurveGeneProfile.FromCodingSequences(genes.Select((g, i) => (0, 10_000 * i, 10_000 * i + g.Length, i % 2 == 1, i, g)));
            var spliced = ZCurveGeneProfile.FromCodingSequences(genes.Select((g, i) => (i % 3, 7_000 * i, 7_000 * i + g.Length + 100 * i, i % 2 == 0, i, g)));
            var queries = new[]
            {
                ZCurveGeneProfile.QueryFromGene(genes[2]), ZCurveGeneProfile.QueryFromProtein(Translate(genes[0])),
                ZCurveGeneProfile.QueryFromProtein(Translate(genes[1])), ZCurveGeneProfile.QueryFromGene(genes[3]),
            };
            foreach (var profile in new[] { unspliced, spliced })
            {
                var bytes = profile.ToBinary();
                var decoded = ZCurveGeneProfile.FromBinary(bytes);
                foreach (var query in queries)
                {
                    var expected = Enumerable.Range(0, decoded.Genes.Count).Where(g => decoded.Genes[g].Residues >= 100)
                        .Select(g => (decoded.Genes[g], decoded.Score(query, g))).ToList();
                    var actual = new List<(ZCurveGeneProfile.Gene, double)>();
                    ZCurveGeneProfile.ScoreRow(bytes, query, new ZCurveGeneProfile.Gene[ZCurveGeneProfile.GeneCountOf(bytes)],
                        g => g.Residues >= 100, (g, score) => actual.Add((g, score)));

                    Assert.That(actual.Select(a => a.Item1), Is.EqualTo(expected.Select(e => e.Item1)));
                    Assert.That(actual.Select(a => BitConverter.DoubleToInt64Bits(a.Item2)), Is.EqualTo(expected.Select(e => BitConverter.DoubleToInt64Bits(e.Item2))));
                }
            }
            Assert.That(BitConverter.ToUInt16(spliced.ToBinary(), 4), Is.EqualTo(3));
        }

        [Test]
        public void GenesOfTheRecordItselfStayFormatVersion2()
        {
            var profile = ZCurveGeneProfile.FromRecord(RandomGene(700, 16), new[] { (0, 2106, false, 1) });
            Assert.That(profile.BinaryFormatVersion, Is.EqualTo(2));
            Assert.That(BitConverter.ToUInt16(profile.ToBinary(), 4), Is.EqualTo(2));
        }

        [Test]
        public void SplicedGenesOnSeveralChromosomesRoundTripAsVersion3()
        {
            var genes = new[] { RandomGene(700, 17), RandomGene(120, 18), RandomGene(900, 19) };
            // Spans as an annotation gives them: wider than the coding sequence by the introns, two genes on chromosome 2.
            var spliced = new[]
            {
                (Chromosome: 2, Start: 5_000, End: 5_000 + genes[0].Length + 640, Reverse: true, Family: 7, Gene: genes[0]),
                (Chromosome: 0, Start: 90, End: 90 + genes[1].Length, Reverse: false, Family: ZCurveGeneProfile.NoFamily, Gene: genes[1]),
                (Chromosome: 2, Start: 1_200, End: 1_200 + genes[2].Length + 3_333, Reverse: false, Family: ZCurveGeneProfile.NoFamily, Gene: genes[2]),
            };
            var profile = ZCurveGeneProfile.FromCodingSequences(spliced);
            var bytes = profile.ToBinary();
            var copy = ZCurveGeneProfile.FromBinary(bytes);

            Assert.That(BitConverter.ToUInt16(bytes, 4), Is.EqualTo(3));
            Assert.That(profile.Genes.Select(g => (g.Chromosome, g.Start, g.Introns)), Is.EqualTo(new[] { (0, 90, 0), (2, 1_200, 3_333), (2, 5_000, 640) }));
            Assert.That(profile.Genes.Select(g => g.Residues), Is.EqualTo(new[] { 121, 901, 701 }), "residues from the coding length, not the span");
            Assert.That(copy.Genes, Is.EqualTo(profile.Genes));
            var alone = ZCurveGeneProfile.FromRecord(genes[2], new[] { (0, genes[2].Length, false, 1) });
            Assert.That(copy.WindowsOf(1), Is.EqualTo(alone.WindowsOf(0)), "introns change where a gene is, not its windows");
        }

        [Test]
        public void SplicedGeneShorterThanItsSpanIsRequired()
        {
            var gene = RandomGene(100, 20);
            Assert.Throws<ArgumentException>(() => ZCurveGeneProfile.FromCodingSequences(new[] { (0, 0, gene.Length - 3, false, 1, gene) }));
        }

        [Test]
        public void AlternativeYeastCodeReadsCtgAsSerine()
        {
            Assert.That(ZCurveGeneProfile.Translate("ATGCTGTTGTAA", ZCurveGeneProfile.alternativeYeastCode), Is.EqualTo("MSL"));
            Assert.That(ZCurveGeneProfile.Translate("ATGCTGTTGTAA"), Is.EqualTo("MLL"));
        }

        [Test]
        public void RejectsForeignOrTruncatedPayloads()
        {
            var bytes = ZCurveGeneProfile.FromRecord(RandomGene(700, 10), new[] { (0, 2106, false, 1) }).ToBinary();
            Assert.Throws<FormatException>(() => ZCurveGeneProfile.FromBinary(bytes.AsSpan(0, bytes.Length - 1)));
            bytes[0] ^= 0xFF;
            Assert.Throws<FormatException>(() => ZCurveGeneProfile.FromBinary(bytes));
        }
    }
}
