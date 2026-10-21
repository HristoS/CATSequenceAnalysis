using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// C3 (platform P3): do strand and rotation break whole-genome comparison? Neither the library nor
    /// the Ingest handles reverse complement or circular topology, so a genome deposited on the other
    /// strand or with a different origin is compared as-is. E. coli K-12 is profiled forward,
    /// reverse-complemented, rotated by half its length, and both, then compared with the original at
    /// 500 bp and through <see cref="ZCurveProfileLadder"/>. O157:H7 (same species) and B. subtilis
    /// (different phylum) give the scale the variants should be read against. Near 1: S3.2 can be
    /// skipped. Low: S3.2 is required before P5's trees.
    /// </summary>
    public class StrandAndRotationBenchmark
    {
        private const string ecoliK12 = "ecoli_k12_mg1655_GCF_000005845.2.fna.gz";
        private const string ecoliO157 = "ecoli_o157h7_sakai_GCF_000008865.2.fna.gz";
        private const string bsubtilis = "bsubtilis_168_GCF_000009045.1.fna.gz";

        [Test]
        public void StrandAndRotation()
        {
            var original = Primary(ecoliK12);
            var o157 = Primary(ecoliO157);

            var rows = new (string Name, string Sequence)[]
            {
                ("K-12 forward (self)", original),
                ("K-12 reverse complement", ReverseComplement(original)),
                ("K-12 rotated by half", RotateBy(original, original.Length / 2)),
                ("K-12 reverse complement, rotated by half", RotateBy(ReverseComplement(original), original.Length / 2)),
                ("O157:H7 forward", o157),
                ("O157:H7 reverse complement", ReverseComplement(o157)),
                ("B. subtilis forward", Primary(bsubtilis)),
            };

            var reference = new ZCurveProfileLadder(new ZCurveWindowedProfile(original));

            using var report = File.CreateText("C3_StrandRotation.csv");
            report.WriteLine("Variant, LengthBp, Compare500bp, LadderFactor, CompareLadder");
            foreach (var (name, sequence) in rows)
            {
                var ladder = new ZCurveProfileLadder(new ZCurveWindowedProfile(sequence));
                var fine = ZCurveWindowedProfile.Compare(reference.Fine, ladder.Fine);
                var factor = ZCurveProfileLadder.SelectFactor(reference, ladder);
                var coarse = ZCurveProfileLadder.Compare(reference, ladder);

                report.WriteLine($"\"{name}\", {sequence.Length}, {fine}, {factor}, {coarse}");
                TestContext.Out.WriteLine($"{name,-42} 500 bp={fine:F4}  ladder (x{factor})={coarse:F4}");
            }

            Assert.Pass();
        }

        // S3.2's gain rests on how often real deposits differ in orientation. Every pair of the
        // 30-genome corpus is scored as stored and with the second genome reverse-complemented; a pair
        // stored on opposite strands should gain about as much as K-12 lost in StrandAndRotation (~0.19
        // at 500 bp). Unrelated pairs gain or lose by noise, so the gains are read per relationship level.
        [Test]
        public void OrientationInTheBundledCorpus()
        {
            var genomes = GenomeTaxonomy.Genomes;
            var sequences = genomes.Select(g => Primary(g.FileName)).ToList();
            var forward = sequences.Select(s => new ZCurveProfileLadder(new ZCurveWindowedProfile(s))).ToList();
            var reversed = sequences.Select(s => new ZCurveProfileLadder(new ZCurveWindowedProfile(ReverseComplement(s)))).ToList();

            using var report = File.CreateText("C3_CorpusOrientation.csv");
            report.WriteLine("GenomeA, GenomeB, Level, Forward500bp, Reversed500bp, Gain500bp, ForwardLadder, ReversedLadder, GainLadder");
            var rows = new List<(RelationshipLevel Level, double Gain, double GainLadder)>();
            for (int a = 0; a < genomes.Length; a++)
                for (int b = a + 1; b < genomes.Length; b++)
                {
                    var level = GenomeTaxonomy.Classify(genomes[a], genomes[b]);
                    var fineForward = ZCurveWindowedProfile.Compare(forward[a].Fine, forward[b].Fine);
                    var fineReversed = ZCurveWindowedProfile.Compare(forward[a].Fine, reversed[b].Fine);
                    var ladderForward = ZCurveProfileLadder.Compare(forward[a], forward[b])!.Value;
                    var ladderReversed = ZCurveProfileLadder.Compare(forward[a], reversed[b])!.Value;

                    rows.Add((level, fineReversed - fineForward, ladderReversed - ladderForward));
                    report.WriteLine($"{genomes[a].FileName}, {genomes[b].FileName}, {level}, {fineForward}, {fineReversed}, " +
                        $"{fineReversed - fineForward}, {ladderForward}, {ladderReversed}, {ladderReversed - ladderForward}");
                }

            foreach (var group in rows.GroupBy(r => r.Level).OrderBy(g => g.Key))
                TestContext.Out.WriteLine($"{group.Key,-16} n={group.Count(),3}  500 bp: reversed higher {group.Count(r => r.Gain > 0),3}, " +
                    $"by > 0.02 {group.Count(r => r.Gain > 0.02),3}, max {group.Max(r => r.Gain):F4}  " +
                    $"32 kbp: higher {group.Count(r => r.GainLadder > 0),3}, by > 0.02 {group.Count(r => r.GainLadder > 0.02),3}, max {group.Max(r => r.GainLadder):F4}");

            Assert.Pass();
        }

        // Can the reverse-complement profile be derived from the stored one, now that the DNA is gone?
        // Flipping a profile (windows reversed, x and y negated, z kept) is exact only when the length
        // is a multiple of the window: otherwise the reverse-complement windows start length mod 500
        // bases off the stored grid, and the trailing partial window was never stored. The trimmed
        // column checks the flip itself; the other two measure the grid error on real lengths.
        [Test]
        public void ProfileFlipAgainstReverseComplementDna()
        {
            using var report = File.CreateText("C3_ProfileFlip.csv");
            report.WriteLine("Genome, LengthBp, LengthMod500, Flip500bp, FlipLadder, LadderFactor, FlipTrimmedToMultiple500bp");
            foreach (var genome in GenomeTaxonomy.Genomes)
            {
                var sequence = Primary(genome.FileName);
                var fromDna = new ZCurveProfileLadder(new ZCurveWindowedProfile(ReverseComplement(sequence)));
                var flipped = new ZCurveProfileLadder(Flip(new ZCurveWindowedProfile(sequence)));
                var trimmed = sequence[..(sequence.Length - sequence.Length % 500)];
                var trimmedScore = ZCurveWindowedProfile.Compare(
                    new ZCurveWindowedProfile(ReverseComplement(trimmed)), Flip(new ZCurveWindowedProfile(trimmed)));

                report.WriteLine($"{genome.FileName}, {sequence.Length}, {sequence.Length % 500}, " +
                    $"{ZCurveWindowedProfile.Compare(fromDna.Fine, flipped.Fine)}, {ZCurveProfileLadder.Compare(fromDna, flipped)}, " +
                    $"{ZCurveProfileLadder.SelectFactor(fromDna, flipped)}, {trimmedScore}");
            }

            Assert.Pass();
        }

        private static ZCurveWindowedProfile Flip(ZCurveWindowedProfile profile) =>
            ZCurveWindowedProfile.FromPrecomputedWindows(profile.WindowSize,
                profile.Windows.Reverse().Select(w => (-w.X, -w.Y, w.Z)).ToList());

        private static string Primary(string fileName) =>
            FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence;

        // IUPAC codes other than the four bases are kept unchanged; both profile types score them as zero anyway.
        private static string ReverseComplement(string sequence)
        {
            var chars = new char[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
                chars[sequence.Length - 1 - i] = sequence[i] switch
                {
                    'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
                    'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c',
                    var other => other,
                };
            return new string(chars);
        }

        private static string RotateBy(string sequence, int shift) => sequence[shift..] + sequence[..shift];
    }
}
