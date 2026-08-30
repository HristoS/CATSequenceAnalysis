using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    /// <summary>
    /// Checks CATProfile/ZCurveWindowedProfile similarity against genuine phylogenetic distance
    /// using the real RefSeq genomes in data/genomes/, instead of only synthetic random sequences
    /// and independent per-base substitution mutants (the gap recorded in
    /// docs/zcurve-windowed-profile-factcheck.md, "Real biological data"). Each genome's primary
    /// record (the chromosome, picked by FastaReader.ReadPrimaryRecord - see data/genomes/README.md
    /// for why several assemblies also carry plasmids that are deliberately excluded here) is
    /// profiled once in [OneTimeSetUp] and reused across every test in this class.
    ///
    /// The distance ladder tested - same species -> same family -> different phylum - mirrors
    /// data/genomes/README.md's "Built-in close pairs" section directly.
    /// </summary>
    public class RealGenomeProfileTests
    {
        private ZCurveWindowedProfile _ecoliK12 = null!;
        private ZCurveWindowedProfile _ecoliO157H7 = null!;
        private ZCurveWindowedProfile _salmonellaTyphimurium = null!;
        private ZCurveWindowedProfile _yersiniaPestis = null!;
        private ZCurveWindowedProfile _yersiniaPseudotuberculosis = null!;
        private ZCurveWindowedProfile _bacillusSubtilis = null!;
        private ZCurveWindowedProfile _deinococcusRadiodurans = null!;

        private static ZCurveWindowedProfile LoadProfile(string fileName)
        {
            var primary = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName));
            return new ZCurveWindowedProfile(primary.Sequence);
        }

        [OneTimeSetUp]
        public void LoadRealGenomes()
        {
            _ecoliK12 = LoadProfile("ecoli_k12_mg1655_GCF_000005845.2.fna.gz");
            _ecoliO157H7 = LoadProfile("ecoli_o157h7_sakai_GCF_000008865.2.fna.gz");
            _salmonellaTyphimurium = LoadProfile("salmonella_typhimurium_lt2_GCF_000006945.2.fna.gz");
            _yersiniaPestis = LoadProfile("yersinia_pestis_co92_GCF_000009065.1.fna.gz");
            _yersiniaPseudotuberculosis = LoadProfile("yersinia_pseudotuberculosis_ip32953_GCF_000047365.1.fna.gz");
            _bacillusSubtilis = LoadProfile("bsubtilis_168_GCF_000009045.1.fna.gz");
            _deinococcusRadiodurans = LoadProfile("deinococcus_radiodurans_r1_GCF_020546685.1.fna.gz");
        }

        [Test]
        public void SameSpeciesStrainsScoreHigherThanCrossFamilyPair()
        {
            // E. coli K-12 vs its own O157:H7 strain (same species) vs K-12 against Salmonella
            // (same family, different genus - Enterobacteriaceae).
            var sameSpecies = ZCurveWindowedProfile.Compare(_ecoliK12, _ecoliO157H7);
            var crossFamily = ZCurveWindowedProfile.Compare(_ecoliK12, _salmonellaTyphimurium);

            TestContext.Out.WriteLine($"E. coli K-12 vs O157:H7 (same species): {sameSpecies}");
            TestContext.Out.WriteLine($"E. coli K-12 vs Salmonella Typhimurium (same family): {crossFamily}");

            Assert.That(sameSpecies, Is.GreaterThan(crossFamily));
        }

        [Test]
        public void CrossFamilyPairScoresHigherThanCrossPhylumPair()
        {
            // Salmonella (Enterobacteriaceae, same family as E. coli) vs Bacillus subtilis
            // (Bacillota - a different phylum entirely from E. coli's Pseudomonadota).
            var crossFamily = ZCurveWindowedProfile.Compare(_ecoliK12, _salmonellaTyphimurium);
            var crossPhylum = ZCurveWindowedProfile.Compare(_ecoliK12, _bacillusSubtilis);

            TestContext.Out.WriteLine($"E. coli K-12 vs Salmonella Typhimurium (same family): {crossFamily}");
            TestContext.Out.WriteLine($"E. coli K-12 vs B. subtilis (different phylum): {crossPhylum}");

            Assert.That(crossFamily, Is.GreaterThan(crossPhylum));
        }

        [Test]
        public void RecentlyDivergedGenusPairScoresHigherThanCrossFamilyPair()
        {
            // Y. pestis is a clonal derivative of Y. pseudotuberculosis on an evolutionary
            // timescale of roughly 1,500-20,000 years (see data/genomes/README.md) - about as
            // close as two named bacterial species get. That should score higher than either
            // Yersinia against E. coli (same family, but nowhere near as recently diverged).
            var sameGenus = ZCurveWindowedProfile.Compare(_yersiniaPestis, _yersiniaPseudotuberculosis);
            var crossGenus = ZCurveWindowedProfile.Compare(_yersiniaPestis, _ecoliK12);

            TestContext.Out.WriteLine($"Y. pestis vs Y. pseudotuberculosis (recently diverged): {sameGenus}");
            TestContext.Out.WriteLine($"Y. pestis vs E. coli K-12 (same family, older split): {crossGenus}");

            Assert.That(sameGenus, Is.GreaterThan(crossGenus));
        }

        [Test]
        public void FullDistanceLadderFromEColiK12IsMonotonic()
        {
            // The comprehensive version of the two tests above, as a single ordered ladder:
            // same species > same family > different phylum entirely, all measured against the
            // same reference genome (E. coli K-12) so the comparison is apples-to-apples.
            var sameSpecies = ZCurveWindowedProfile.Compare(_ecoliK12, _ecoliO157H7);
            var sameFamily = ZCurveWindowedProfile.Compare(_ecoliK12, _salmonellaTyphimurium);
            var differentPhylum = ZCurveWindowedProfile.Compare(_ecoliK12, _bacillusSubtilis);
            var mostDistantOutgroup = ZCurveWindowedProfile.Compare(_ecoliK12, _deinococcusRadiodurans);

            TestContext.Out.WriteLine($"same species (E. coli K-12 vs O157:H7):      {sameSpecies}");
            TestContext.Out.WriteLine($"same family  (E. coli K-12 vs Salmonella):   {sameFamily}");
            TestContext.Out.WriteLine($"diff. phylum (E. coli K-12 vs B. subtilis):  {differentPhylum}");
            TestContext.Out.WriteLine($"outgroup     (E. coli K-12 vs Deinococcus):  {mostDistantOutgroup}");

            Assert.That(sameSpecies, Is.GreaterThan(sameFamily));
            Assert.That(sameFamily, Is.GreaterThan(differentPhylum));
            Assert.That(sameFamily, Is.GreaterThan(mostDistantOutgroup));
        }

        [Test]
        public void IdenticalGenomeComparesAsOne()
        {
            // Sanity check that loading the same file twice through the full FASTA -> profile
            // pipeline reproduces ZCurveWindowedProfileTests.IdenticalSequencesCompareAsOne's
            // guarantee, not just for Generator's synthetic sequences.
            var reloaded = LoadProfile("ecoli_k12_mg1655_GCF_000005845.2.fna.gz");

            Assert.That(ZCurveWindowedProfile.Compare(_ecoliK12, reloaded), Is.EqualTo(1.0).Within(1e-9));
        }
    }
}
