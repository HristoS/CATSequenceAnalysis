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
