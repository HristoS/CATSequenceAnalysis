using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Could the platform tell a tumour genome from its normal one (2026-11-03, before looking for cancer data)? Somatic point
    /// mutations run at about 1-10 per Mbp in most tumours and up to about 100 in hypermutated ones; two unrelated people differ at
    /// about 1 SNV per kbp. E. coli K-12's chromosome stands in for a genome (no human DNA): substitutions are planted at those rates,
    /// either uniformly or as C>T at CpG (SBS1, the commonest somatic signature, and the one that moves the CG dinucleotide the
    /// signature is most sensitive to). Each variant is scored against the original the way the platform scores it: 1 - Compare at
    /// 500 bp, 1 - Compare through the ladder, and the Euclidean distance between log rho* signatures (SignatureIndex). Real strain
    /// pairs give the scale. Stop criterion agreed beforehand: if 100 per Mbp moves the score less than the closest strain pair,
    /// cancer genomes aren't worth loading for this method.
    /// </summary>
    public class SomaticMutationRateBenchmark
    {
        private const string ecoliK12 = "ecoli_k12_mg1655_GCF_000005845.2.fna.gz";
        private const int seedCount = 5;

        private static readonly double[] ratesPerMbp = { 1, 10, 100, 1000 };

        private static readonly (string Name, string A, string B)[] strainPairs =
        {
            ("E. coli K-12 vs O157:H7", ecoliK12, "ecoli_o157h7_sakai_GCF_000008865.2.fna.gz"),
            ("E. coli K-12 vs Shigella flexneri 2a", ecoliK12, "shigella_flexneri_2a_301_GCF_000006925.2.fna.gz"),
            ("Y. pestis CO92 vs Y. pseudotuberculosis", "yersinia_pestis_co92_GCF_000009065.1.fna.gz", "yersinia_pseudotuberculosis_ip32953_GCF_000047365.1.fna.gz"),
            ("B. anthracis Ames vs B. cereus ATCC 14579", "bacillus_anthracis_ames_GCF_000007845.1.fna.gz", "bacillus_cereus_atcc14579_GCF_006094295.1.fna.gz"),
        };

        [Test]
        public void SomaticRates()
        {
            var original = Primary(ecoliK12);
            var reference = Measure(original);

            using var report = File.CreateText("SomaticMutationRate.csv");
            report.WriteLine("Comparison, Spectrum, RatePerMbp, Seed, Substitutions, Distance500bp, LadderFactor, DistanceLadder, SignatureDistance");

            foreach (var (name, a, b) in strainPairs)
            {
                var first = Measure(Primary(a));
                var second = Measure(Primary(b));
                var (fine, factor, ladder, signature) = Distances(first, second);
                Write(report, $"\"{name}\", strain pair, , , , {Format(fine)}, {factor}, {Format(ladder)}, {Format(signature)}");
            }

            foreach (var spectrum in new[] { Spectrum.Uniform, Spectrum.CpgTransition })
                foreach (var rate in ratesPerMbp)
                    for (int seed = 0; seed < seedCount; seed++)
                    {
                        var count = (int)Math.Round(rate * original.Length / 1e6);
                        var mutated = Plant(original, count, spectrum, BenchmarkTiming.Seed(nameof(SomaticRates), spectrum, rate, seed));
                        var (fine, factor, ladder, signature) = Distances(reference, Measure(mutated));
                        Write(report, $"\"K-12 vs K-12 mutated\", {spectrum}, {rate}, {seed}, {count}, {Format(fine)}, {factor}, {Format(ladder)}, {Format(signature)}");
                    }

            Assert.Pass();
        }

        private enum Spectrum { Uniform, CpgTransition }

        private sealed record Measured(ZCurveProfileLadder Ladder, double[] Signature);

        private static Measured Measure(string sequence) =>
            new(new ZCurveProfileLadder(new ZCurveWindowedProfile(sequence)), GenomeSignatures.KarlinLogRho(GenomeSignatures.DinucleotideCounts(sequence)));

        private static (double Fine, int? Factor, double? Ladder, double Signature) Distances(Measured x, Measured y)
        {
            var fine = 1 - ZCurveWindowedProfile.Compare(x.Ladder.Fine, y.Ladder.Fine);
            var ladder = 1 - ZCurveProfileLadder.Compare(x.Ladder, y.Ladder);
            var signature = Math.Sqrt(x.Signature.Zip(y.Signature, (a, b) => (a - b) * (a - b)).Sum());
            return (fine, ZCurveProfileLadder.SelectFactor(x.Ladder, y.Ladder), ladder, signature);
        }

        // Exactly count distinct positions, each changed to a different base. CpG transitions take the C of a CG (C>T) or, for the
        // other strand's C, the G (G>A), half each.
        private static string Plant(string sequence, int count, Spectrum spectrum, int seed)
        {
            var rng = new Random(seed);
            var chars = sequence.ToCharArray();
            var candidates = spectrum == Spectrum.Uniform
                ? null
                : Enumerable.Range(0, chars.Length - 1).Where(i => chars[i] == 'C' && chars[i + 1] == 'G').ToArray();
            var chosen = new HashSet<int>();
            while (chosen.Count < count)
            {
                if (candidates is null)
                {
                    var position = rng.Next(chars.Length);
                    if (!chosen.Add(position)) continue;
                    var bases = "ACGT".Where(c => c != chars[position]).ToArray();
                    chars[position] = bases[rng.Next(bases.Length)];
                }
                else
                {
                    var site = candidates[rng.Next(candidates.Length)];
                    var position = rng.Next(2) == 0 ? site : site + 1;
                    if (chosen.Contains(site) || chosen.Contains(site + 1) || !chosen.Add(position)) continue;
                    chars[position] = chars[position] == 'C' ? 'T' : 'A';
                }
            }
            return new string(chars);
        }

        private static void Write(StreamWriter report, string line)
        {
            report.WriteLine(line);
            TestContext.Out.WriteLine(line);
        }

        private static string Format(double? value) => value?.ToString("G6", CultureInfo.InvariantCulture) ?? "";

        private static string Primary(string fileName) => FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence.ToUpperInvariant();
    }
}
