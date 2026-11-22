using System.Globalization;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// What the gene profile's axes mean for the protein (2026-11-22). A protein query is built from each amino acid's expected Z-curve
    /// axes at codon positions 1 and 2 (ZCurveGeneProfile.ExpectedAxes, uniform over synonymous codons); this correlates those 6 numbers,
    /// over the 20 amino acids, with standard residue properties: Kyte and Doolittle's hydropathy (1982), Zamyatnin's residue volume
    /// (1972, cubic angstroms) and side-chain charge at pH 7 (His counted as 0). Also position 2's x + y, which is +2 for A, -2 for T and
    /// 0 for C and G. Pearson and Spearman r, with a two-sided permutation p (10,000 shuffles, seeded; n = 20). Writes
    /// AminoAcidAxisProperties.csv. A description, not a rule: it says what a profile difference means biologically.
    /// </summary>
    public class AminoAcidAxisProperties
    {
        private const string aminoAcids = "ARNDCQEGHILKMFPSTWYV";
        private const int shuffles = 10_000, seed = 20261122;
        private static readonly double[] hydropathy = { 1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2 };
        private static readonly double[] volume = { 88.6, 173.4, 114.1, 111.1, 108.5, 143.8, 138.4, 60.1, 153.2, 166.7, 166.7, 168.6, 162.9, 189.9, 112.7, 89.0, 116.1, 227.8, 193.6, 140.0 };
        private static readonly double[] charge = { 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 };

        [Test]
        public void ExpectedAxesAgainstResidueProperties()
        {
            var table = ZCurveGeneProfile.ExpectedAxes();
            (double X, double Y, double Z) At(char aminoAcid, int position) => table[2 * (aminoAcid - 'A') + position];
            // Spot checks: Met is ATG only (A at position 1: x, y, z = 1, 1, -1); Trp is TGG only (G at position 2: 1, -1, 1).
            Assert.That(At('M', 0), Is.EqualTo((1.0, 1.0, -1.0)));
            Assert.That(At('W', 1), Is.EqualTo((1.0, -1.0, 1.0)));

            var features = new List<(string Name, double[] Values)>();
            foreach (int position in new[] { 0, 1 })
            {
                features.Add(($"position {position + 1} x (purine)", aminoAcids.Select(a => At(a, position).X).ToArray()));
                features.Add(($"position {position + 1} y (amino)", aminoAcids.Select(a => At(a, position).Y).ToArray()));
                features.Add(($"position {position + 1} z (strong)", aminoAcids.Select(a => At(a, position).Z).ToArray()));
            }
            features.Add(("position 2 x + y (A minus T)", aminoAcids.Select(a => At(a, 1).X + At(a, 1).Y).ToArray()));
            var properties = new (string Name, double[] Values)[] { ("hydropathy", hydropathy), ("volume", volume), ("charge", charge) };

            var lines = new List<string> { "Feature,Property,Pearson r,Spearman r,Permutation p" };
            var rng = new Random(seed);
            foreach (var (featureName, feature) in features)
                foreach (var (propertyName, property) in properties)
                {
                    double r = Pearson(feature, property);
                    var shuffled = property.ToArray();
                    int extreme = 0;
                    for (int s = 0; s < shuffles; s++)
                    {
                        rng.Shuffle(shuffled);
                        if (Math.Abs(Pearson(feature, shuffled)) >= Math.Abs(r) - 1e-12) extreme++;
                    }
                    lines.Add(string.Join(",", featureName, propertyName, Format(r), Format(Pearson(Ranks(feature), Ranks(property))),
                        Format((extreme + 1) / (double)(shuffles + 1))));
                }
            lines.Add("");
            lines.Add("Amino acid," + string.Join(",", features.Select(f => f.Name)) + ",hydropathy,volume,charge");
            for (int i = 0; i < aminoAcids.Length; i++)
                lines.Add(string.Join(",", new[] { aminoAcids[i].ToString() }.Concat(features.Select(f => Format(f.Values[i])))
                    .Concat(new[] { hydropathy[i], volume[i], charge[i] }.Select(Format))));

            File.WriteAllLines("AminoAcidAxisProperties.csv", lines);
            foreach (var line in lines) TestContext.Out.WriteLine(line);
        }

        /// <summary>
        /// The same on real windows (2026-11-22): genes called on the first 6 bundled bacteria (GeneLevelSignatureBenchmark.CallGenes,
        /// from 50 aa), windows from the DNA as stored (QueryFromGene). Each window's Kyte-Doolittle series, one value per codon position
        /// 1 and 2 base (its residue's hydropathy), gets the same three moments as the axes, so window and hydropathy moments are taken
        /// over exactly the same bases. Pearson r across windows, separately for genes from 500 aa (333-base windows) and shorter ones (5
        /// windows), per moment, for x, y, z and x + y. Stored windows mix positions 1 and 2, and position 1 carries no hydropathy, so this
        /// is diluted against the amino-acid table. Reading fixed before the run: P0 x + y against mean hydropathy |r| >= 0.5 -> "the
        /// profile is a DNA-level hydropathy profile" holds for real windows; < 0.3 -> the claim stays at the amino-acid level.
        /// Writes AminoAcidAxisProperties_Windows.csv.
        /// </summary>
        [Test]
        public void WindowAxesTrackWindowHydropathy()
        {
            var files = ProteinBackTranslationFeasibility.Corpus().Where(c => c.Group == "Bacteria").Select(c => c.FileName).Take(6).ToList();
            const int F = ZCurveGeneProfile.FeatureCount;
            var pairs = new List<(bool Long, double[] Axes, double[] Hydropathy)>();
            int genes = 0;
            foreach (var file in files)
            {
                var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(file)).Sequence;
                var code = Path.GetFileName(file).StartsWith("mycoplasma_", StringComparison.Ordinal)
                    ? ProteinBackTranslationFeasibility.mycoplasmaCode : ProteinBackTranslationFeasibility.standardCode;
                foreach (var dna in GeneLevelSignatureBenchmark.CallGenes(sequence, code, 3 * ZCurveGeneProfile.MinResidues))
                {
                    var protein = ProteinBackTranslationFeasibility.Translate(dna, code).TrimEnd('*');
                    if (protein.Length < ZCurveGeneProfile.MinResidues || protein.Contains('*')) continue;
                    var windows = ZCurveGeneProfile.QueryFromGene(dna, code);
                    int window = ZCurveGeneProfile.WindowOf(protein.Length), count = windows.Length / F;
                    Assert.That(count, Is.EqualTo(2 * protein.Length / window), "window count");
                    genes++;
                    for (int w = 0; w < count; w++)
                    {
                        var moments = new double[3];
                        for (int i = 0; i < window; i++)
                        {
                            int index = aminoAcids.IndexOf(protein[(w * window + i) / 2]);
                            double value = index < 0 ? 0 : hydropathy[index], t = -1 + 2 * (i + 0.5) / window;
                            moments[0] += value;
                            moments[1] += value * t;
                            moments[2] += value * (3 * t * t - 1) / 2;
                        }
                        for (int k = 0; k < 3; k++) moments[k] *= Math.Sqrt(2 * k + 1) / window;
                        pairs.Add((protein.Length >= ZCurveGeneProfile.LongResidues, windows[(w * F)..((w + 1) * F)], moments));
                    }
                }
            }

            string[] momentNames = { "mean (P0)", "trend (P1)", "curvature (P2)" };
            var lines = new List<string> { $"Genes,Moment,Axis,Windows,Pearson r ({genes:N0} genes on {files.Count} genomes)" };
            foreach (bool isLong in new[] { true, false })
            {
                var cell = pairs.Where(p => p.Long == isLong).ToList();
                for (int k = 0; k < 3; k++)
                {
                    var target = cell.Select(p => p.Hydropathy[k]).ToArray();
                    foreach (var (axisName, axis) in new (string, Func<double[], double>)[]
                        { ("x", a => a[3 * k]), ("y", a => a[3 * k + 1]), ("z", a => a[3 * k + 2]), ("x + y", a => a[3 * k] + a[3 * k + 1]) })
                        lines.Add(string.Join(",", isLong ? "500+ aa (333-base windows)" : "50-499 aa (5 windows)", momentNames[k], axisName, cell.Count,
                            Format(Pearson(cell.Select(p => axis(p.Axes)).ToArray(), target))));
                }
            }
            // Reference: the amino-acid table with positions 1 and 2 averaged, as a stored window mixes them.
            var table = ZCurveGeneProfile.ExpectedAxes();
            var mixed = aminoAcids.Select(a => (table[2 * (a - 'A')].X + table[2 * (a - 'A')].Y + table[2 * (a - 'A') + 1].X + table[2 * (a - 'A') + 1].Y) / 2).ToArray();
            lines.Add($"Amino-acid table: positions 1 and 2 averaged x + y against hydropathy,,,20,{Format(Pearson(mixed, hydropathy))}");
            double reading = Math.Abs(Pearson(pairs.Where(p => p.Long).Select(p => p.Axes[0] + p.Axes[1]).ToArray(), pairs.Where(p => p.Long).Select(p => p.Hydropathy[0]).ToArray()));
            lines.Add($"Reading (P0 x + y against mean hydropathy, 500+ aa: |r| >= 0.5 holds, < 0.3 amino-acid level only),,,,{Format(reading)}: " +
                (reading >= 0.5 ? "holds for real windows" : reading < 0.3 ? "amino-acid level only" : "in between"));

            File.WriteAllLines("AminoAcidAxisProperties_Windows.csv", lines);
            foreach (var line in lines) TestContext.Out.WriteLine(line);
        }

        private static double Pearson(double[] a, double[] b)
        {
            double ma = a.Average(), mb = b.Average(), ab = 0, aa = 0, bb = 0;
            for (int i = 0; i < a.Length; i++) { ab += (a[i] - ma) * (b[i] - mb); aa += (a[i] - ma) * (a[i] - ma); bb += (b[i] - mb) * (b[i] - mb); }
            return aa == 0 || bb == 0 ? double.NaN : ab / Math.Sqrt(aa * bb);
        }

        // Average ranks, ties sharing their mean rank.
        private static double[] Ranks(double[] values)
        {
            var order = Enumerable.Range(0, values.Length).OrderBy(i => values[i]).ToArray();
            var ranks = new double[values.Length];
            for (int start = 0; start < order.Length;)
            {
                int end = start;
                while (end + 1 < order.Length && values[order[end + 1]] == values[order[start]]) end++;
                for (int k = start; k <= end; k++) ranks[order[k]] = (start + end) / 2.0 + 1;
                start = end + 1;
            }
            return ranks;
        }

        private static string Format(double value) => value.ToString("F3", CultureInfo.InvariantCulture);
    }
}
