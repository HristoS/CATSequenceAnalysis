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
