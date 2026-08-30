namespace SequenceAnaliseTests
{
    internal enum RelationshipLevel
    {
        SameSpecies,       // different strain of the same named species
        SameGenus,         // different species, same genus
        SamePhylum,        // different genus, same phylum
        DifferentPhylum,
    }

    internal sealed record GenomeMeta(string FileName, string OrganismName, string Genus, string Phylum);

    /// <summary>
    /// The 30-genome real bacterial corpus's taxonomy labels, shared by
    /// RealGenomeDistanceMatrixBenchmark.cs (ZCurveWindowedProfile vs. taxonomy) and
    /// ExternalBaselineComparisonBenchmark.cs (Mash/FastANI vs. taxonomy) so both benchmarks compare
    /// against the exact same relationship labels rather than two independently-typed copies that
    /// could silently drift apart. Genus/phylum are the only taxonomic claims made anywhere - see
    /// RealGenomeDistanceMatrixBenchmark's remarks for why (family/order placements are unsettled for
    /// several of these organisms).
    /// </summary>
    internal static class GenomeTaxonomy
    {
        public static readonly GenomeMeta[] Genomes =
        {
            new("ecoli_k12_mg1655_GCF_000005845.2.fna.gz", "Escherichia coli K-12 MG1655", "Escherichia", "Pseudomonadota"),
            new("ecoli_o157h7_sakai_GCF_000008865.2.fna.gz", "Escherichia coli O157:H7 Sakai", "Escherichia", "Pseudomonadota"),
            new("shigella_flexneri_2a_301_GCF_000006925.2.fna.gz", "Shigella flexneri 2a 301", "Shigella", "Pseudomonadota"),
            new("klebsiella_pneumoniae_mgh78578_GCF_000016305.1.fna.gz", "Klebsiella pneumoniae MGH 78578", "Klebsiella", "Pseudomonadota"),
            new("salmonella_typhimurium_lt2_GCF_000006945.2.fna.gz", "Salmonella enterica Typhimurium LT2", "Salmonella", "Pseudomonadota"),
            new("yersinia_pestis_co92_GCF_000009065.1.fna.gz", "Yersinia pestis CO92", "Yersinia", "Pseudomonadota"),
            new("yersinia_pseudotuberculosis_ip32953_GCF_000047365.1.fna.gz", "Yersinia pseudotuberculosis IP 32953", "Yersinia", "Pseudomonadota"),
            new("pseudomonas_aeruginosa_pao1_GCF_000006765.1.fna.gz", "Pseudomonas aeruginosa PAO1", "Pseudomonas", "Pseudomonadota"),
            new("vibrio_cholerae_n16961_GCF_900205735.1.fna.gz", "Vibrio cholerae N16961", "Vibrio", "Pseudomonadota"),
            new("haemophilus_influenzae_rd_kw20_GCF_000027305.1.fna.gz", "Haemophilus influenzae Rd KW20", "Haemophilus", "Pseudomonadota"),
            new("acinetobacter_baumannii_atcc17978_GCF_013372085.1.fna.gz", "Acinetobacter baumannii ATCC 17978", "Acinetobacter", "Pseudomonadota"),
            new("bsubtilis_168_GCF_000009045.1.fna.gz", "Bacillus subtilis 168", "Bacillus", "Bacillota"),
            new("bacillus_anthracis_ames_GCF_000007845.1.fna.gz", "Bacillus anthracis Ames", "Bacillus", "Bacillota"),
            new("bacillus_cereus_atcc14579_GCF_006094295.1.fna.gz", "Bacillus cereus ATCC 14579", "Bacillus", "Bacillota"),
            new("staphylococcus_aureus_nctc8325_GCF_000013425.1.fna.gz", "Staphylococcus aureus NCTC 8325", "Staphylococcus", "Bacillota"),
            new("staphylococcus_epidermidis_atcc12228_GCF_000007645.1.fna.gz", "Staphylococcus epidermidis ATCC 12228", "Staphylococcus", "Bacillota"),
            new("listeria_monocytogenes_egde_GCF_000196035.1.fna.gz", "Listeria monocytogenes EGD-e", "Listeria", "Bacillota"),
            new("streptococcus_pneumoniae_tigr4_GCF_000273445.1.fna.gz", "Streptococcus pneumoniae TIGR4", "Streptococcus", "Bacillota"),
            new("streptococcus_pyogenes_m1gas_GCF_043231225.1.fna.gz", "Streptococcus pyogenes M1 GAS", "Streptococcus", "Bacillota"),
            new("enterococcus_faecalis_v583_GCF_000007785.1.fna.gz", "Enterococcus faecalis V583", "Enterococcus", "Bacillota"),
            new("clostridioides_difficile_630_GCF_000009205.2.fna.gz", "Clostridioides difficile 630", "Clostridioides", "Bacillota"),
            new("mycoplasma_genitalium_g37_GCF_000027325.1.fna.gz", "Mycoplasmoides genitalium G37", "Mycoplasmoides", "Mycoplasmatota"),
            new("mycoplasma_pneumoniae_m129_GCF_910574535.1.fna.gz", "Mycoplasmoides pneumoniae M129", "Mycoplasmoides", "Mycoplasmatota"),
            new("mycobacterium_tuberculosis_h37rv_GCF_000195955.2.fna.gz", "Mycobacterium tuberculosis H37Rv", "Mycobacterium", "Actinomycetota"),
            new("mycobacterium_smegmatis_mc2155_GCF_000767605.1.fna.gz", "Mycolicibacterium smegmatis MC2 155", "Mycolicibacterium", "Actinomycetota"),
            new("streptomyces_coelicolor_a3_GCF_008931305.1.fna.gz", "Streptomyces coelicolor A3(2)", "Streptomyces", "Actinomycetota"),
            new("helicobacter_pylori_26695_GCF_000307795.1.fna.gz", "Helicobacter pylori 26695", "Helicobacter", "Campylobacterota"),
            new("campylobacter_jejuni_nctc11168_GCF_000009085.1.fna.gz", "Campylobacter jejuni NCTC 11168", "Campylobacter", "Campylobacterota"),
            new("borrelia_burgdorferi_b31_GCF_000008685.2.fna.gz", "Borreliella burgdorferi B31", "Borreliella", "Spirochaetota"),
            new("deinococcus_radiodurans_r1_GCF_020546685.1.fna.gz", "Deinococcus radiodurans R1", "Deinococcus", "Deinococcota"),
        };

        public static RelationshipLevel Classify(GenomeMeta a, GenomeMeta b)
        {
            if (a.Genus == b.Genus)
                return a.OrganismName.Split(' ')[0] == b.OrganismName.Split(' ')[0]
                    && a.OrganismName.Split(' ')[1] == b.OrganismName.Split(' ')[1]
                    ? RelationshipLevel.SameSpecies
                    : RelationshipLevel.SameGenus;
            return a.Phylum == b.Phylum ? RelationshipLevel.SamePhylum : RelationshipLevel.DifferentPhylum;
        }
    }
}
