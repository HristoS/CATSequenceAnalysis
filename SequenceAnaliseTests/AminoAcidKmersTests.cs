using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    public class AminoAcidKmersTests
    {
        [Test]
        public void CountsEveryWordInsideEachProtein()
        {
            var proteins = new[] { "MKVLAAGIVG", "ACDE", "WY" };
            var counts = AminoAcidKmers.Count(proteins, 4);

            Assert.That(counts.Length, Is.EqualTo(160_000));
            Assert.That(counts.Sum(), Is.EqualTo(7 + 1 + 0), "words don't cross from one protein to the next");
            Assert.That(counts[1 * 400 + 2 * 20 + 3], Is.EqualTo(1), "ACDE is 0123 in base 20");
        }

        [Test]
        public void StopsAndUnknownResiduesBreakWords()
        {
            Assert.That(AminoAcidKmers.Count(new[] { "MKV*LAA" }, 3).Sum(), Is.EqualTo(1 + 1));
            Assert.That(AminoAcidKmers.Count(new[] { "MKVXLAA" }, 3).Sum(), Is.EqualTo(1 + 1));
            Assert.That(AminoAcidKmers.Count(new[] { "mkvlaa" }, 3).Sum(), Is.EqualTo(4), "lower case counts as upper case");
        }

        [Test]
        public void TranslationLeavesOutTheFinalStopOnlyAndMarksUnknownCodons()
        {
            Assert.That(ZCurveGeneProfile.Translate("ATGAAATAA"), Is.EqualTo("MK"));
            Assert.That(ZCurveGeneProfile.Translate("ATGTAAAAATAG"), Is.EqualTo("M*K"));
            Assert.That(ZCurveGeneProfile.Translate("ATGANAAAA"), Is.EqualTo("MXK"));
            Assert.That(ZCurveGeneProfile.Translate("ATGTGAAAA", ZCurveGeneProfile.mycoplasmaCode), Is.EqualTo("MWK"));
        }
    }
}
