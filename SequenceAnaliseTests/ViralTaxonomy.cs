namespace SequenceAnaliseTests
{
    internal enum ViralRelationshipLevel
    {
        SameGenus,
        SameFamily,       // different genus, same family
        DifferentFamily,
    }

    internal sealed record ViralGenomeMeta(string FileName, string Name, string Genus, string Family);

    /// <summary>
    /// The 31-genome viral corpus (data/genomes/viral/, plus the phiX174 and lambda files already in
    /// data/genomes/) with the genus and family NCBI Taxonomy returned for each record's taxid on
    /// 2026-09-26 - see data/genomes/viral/README.md. FileName is relative to data/genomes/, so
    /// GenomeCatalog.PathTo resolves both locations.
    ///
    /// Family is the coarsest rank used: viral taxonomy above family (order/class/phylum/realm) is
    /// sparse and was reorganized repeatedly, whereas family is populated for every genome here.
    /// </summary>
    internal static class ViralTaxonomy
    {
        public static readonly ViralGenomeMeta[] Genomes =
        {
            new("viral/hsv1_NC_001806.2.fna.gz", "Human herpesvirus 1 strain 17", "Simplexvirus", "Orthoherpesviridae"),
            new("viral/hsv2_NC_001798.2.fna.gz", "Human herpesvirus 2 strain HG52", "Simplexvirus", "Orthoherpesviridae"),
            new("viral/vzv_NC_001348.1.fna.gz", "Human herpesvirus 3", "Varicellovirus", "Orthoherpesviridae"),
            new("viral/ebv_NC_007605.1.fna.gz", "Human gammaherpesvirus 4", "Lymphocryptovirus", "Orthoherpesviridae"),
            new("viral/hcmv_NC_006273.2.fna.gz", "Human herpesvirus 5 strain Merlin", "Cytomegalovirus", "Orthoherpesviridae"),
            new("viral/vaccinia_NC_006998.1.fna.gz", "Vaccinia virus", "Orthopoxvirus", "Poxviridae"),
            new("viral/cowpox_NC_003663.2.fna.gz", "Cowpox virus", "Orthopoxvirus", "Poxviridae"),
            new("viral/molluscum_contagiosum_NC_001731.1.fna.gz", "Molluscum contagiosum virus subtype 1", "Molluscipoxvirus", "Poxviridae"),
            new("viral/sars_cov_2_NC_045512.2.fna.gz", "Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1", "Betacoronavirus", "Coronaviridae"),
            new("viral/mers_cov_NC_019843.3.fna.gz", "Middle East respiratory syndrome-related coronavirus isolate HCoV-EMC/2012", "Betacoronavirus", "Coronaviridae"),
            new("viral/hcov_oc43_NC_006213.1.fna.gz", "Human coronavirus OC43 strain ATCC VR-759", "Betacoronavirus", "Coronaviridae"),
            new("viral/hcov_229e_NC_002645.1.fna.gz", "Human coronavirus 229E", "Alphacoronavirus", "Coronaviridae"),
            new("viral/dengue2_NC_001474.2.fna.gz", "Dengue virus 2", "Orthoflavivirus", "Flaviviridae"),
            new("viral/zika_NC_012532.1.fna.gz", "Zika virus strain MR 766", "Orthoflavivirus", "Flaviviridae"),
            new("viral/hcv1_NC_004102.1.fna.gz", "Hepatitis C virus genotype 1", "Orthohepacivirus", "Hepaciviridae"),
            new("viral/hiv1_NC_001802.1.fna.gz", "Human immunodeficiency virus 1", "Lentivirus", "Retroviridae"),
            new("viral/hiv2_NC_001722.1.fna.gz", "Human immunodeficiency virus 2", "Lentivirus", "Retroviridae"),
            new("viral/htlv1_NC_001436.1.fna.gz", "Human T-lymphotropic virus 1", "Deltaretrovirus", "Retroviridae"),
            new("viral/hadv5_AC_000008.1.fna.gz", "Human adenovirus 5", "Mastadenovirus", "Adenoviridae"),
            new("viral/hadv12_NC_001460.1.fna.gz", "Human adenovirus A", "Mastadenovirus", "Adenoviridae"),
            new("viral/hpv16_NC_001526.4.fna.gz", "Human papillomavirus type 16", "Alphapapillomavirus", "Papillomaviridae"),
            new("viral/hpv18_NC_001357.1.fna.gz", "Human papillomavirus - 18", "Alphapapillomavirus", "Papillomaviridae"),
            new("viral/hpv1_NC_001356.1.fna.gz", "Human papillomavirus - 1", "Mupapillomavirus", "Papillomaviridae"),
            new("viral/poliovirus1_NC_002058.3.fna.gz", "Poliovirus", "Enterovirus", "Picornaviridae"),
            new("viral/coxsackievirus_a16_NC_001612.1.fna.gz", "Human enterovirus A", "Enterovirus", "Picornaviridae"),
            new("viral/hepatitis_a_NC_001489.1.fna.gz", "Hepatitis A virus", "Hepatovirus", "Picornaviridae"),
            new("viral/phage_t4_NC_000866.4.fna.gz", "Enterobacteria phage T4", "Tequatrovirus", "Straboviridae"),
            new("viral/phage_t7_NC_001604.1.fna.gz", "Enterobacteria phage T7", "Teseptimavirus", "Autotranscriptaviridae"),
            new("viral/phage_t3_NC_003298.1.fna.gz", "Enterobacteria phage T3", "Teetrevirus", "Autotranscriptaviridae"),
            new("phix174_GCF_000819615.1.fna.gz", "Escherichia phage phiX174", "Sinsheimervirus", "Eubullaviridae"),
            new("lambda_phage_GCF_000840245.1.fna.gz", "Enterobacteria phage lambda", "Lambdavirus", "Zimmerviridae"),
        };

        public static ViralRelationshipLevel Classify(ViralGenomeMeta a, ViralGenomeMeta b)
        {
            if (a.Genus == b.Genus) return ViralRelationshipLevel.SameGenus;
            return a.Family == b.Family ? ViralRelationshipLevel.SameFamily : ViralRelationshipLevel.DifferentFamily;
        }
    }
}
