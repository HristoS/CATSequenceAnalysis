# Viral genomes for window-count validation

31 complete RefSeq viral genomes: 29 downloaded here on 2026-09-26, plus the phiX174 and lambda files already in
`data/genomes/`. They exist to test `ZCurveWindowedProfile` and `ZCurveProfileLadder` on sequences short enough
(5.4–236 kbp) that coarse windows leave them with few or zero windows. The 30-genome bacterial corpus never gets that
short. They're used by `MultiResolutionBenchmark.WindowCountSweep` and `.LadderRuleAgainstFineOnly`, with taxonomy in
`SequenceAnaliseTests/ViralTaxonomy.cs`.

| File | Organism (RefSeq title) | Accession | Length (bp) | Genus | Family | Realm |
|---|---|---|---|---|---|---|
| `hsv1_NC_001806.2.fna.gz` | Human herpesvirus 1 strain 17 | NC_001806.2 | 152,222 | *Simplexvirus* | Orthoherpesviridae | Duplodnaviria |
| `hsv2_NC_001798.2.fna.gz` | Human herpesvirus 2 strain HG52 | NC_001798.2 | 154,675 | *Simplexvirus* | Orthoherpesviridae | Duplodnaviria |
| `vzv_NC_001348.1.fna.gz` | Human herpesvirus 3 | NC_001348.1 | 124,884 | *Varicellovirus* | Orthoherpesviridae | Duplodnaviria |
| `ebv_NC_007605.1.fna.gz` | Human gammaherpesvirus 4 | NC_007605.1 | 171,823 | *Lymphocryptovirus* | Orthoherpesviridae | Duplodnaviria |
| `hcmv_NC_006273.2.fna.gz` | Human herpesvirus 5 strain Merlin | NC_006273.2 | 235,646 | *Cytomegalovirus* | Orthoherpesviridae | Duplodnaviria |
| `vaccinia_NC_006998.1.fna.gz` | Vaccinia virus | NC_006998.1 | 194,711 | *Orthopoxvirus* | Poxviridae | Varidnaviria |
| `cowpox_NC_003663.2.fna.gz` | Cowpox virus | NC_003663.2 | 224,499 | *Orthopoxvirus* | Poxviridae | Varidnaviria |
| `molluscum_contagiosum_NC_001731.1.fna.gz` | Molluscum contagiosum virus subtype 1 | NC_001731.1 | 190,289 | *Molluscipoxvirus* | Poxviridae | Varidnaviria |
| `sars_cov_2_NC_045512.2.fna.gz` | Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1 | NC_045512.2 | 29,903 | *Betacoronavirus* | Coronaviridae | Riboviria |
| `mers_cov_NC_019843.3.fna.gz` | Middle East respiratory syndrome-related coronavirus isolate HCoV-EMC/2012 | NC_019843.3 | 30,119 | *Betacoronavirus* | Coronaviridae | Riboviria |
| `hcov_oc43_NC_006213.1.fna.gz` | Human coronavirus OC43 strain ATCC VR-759 | NC_006213.1 | 30,741 | *Betacoronavirus* | Coronaviridae | Riboviria |
| `hcov_229e_NC_002645.1.fna.gz` | Human coronavirus 229E | NC_002645.1 | 27,317 | *Alphacoronavirus* | Coronaviridae | Riboviria |
| `dengue2_NC_001474.2.fna.gz` | Dengue virus 2 | NC_001474.2 | 10,723 | *Orthoflavivirus* | Flaviviridae | Riboviria |
| `zika_NC_012532.1.fna.gz` | Zika virus strain MR 766 | NC_012532.1 | 10,794 | *Orthoflavivirus* | Flaviviridae | Riboviria |
| `hcv1_NC_004102.1.fna.gz` | Hepatitis C virus genotype 1 | NC_004102.1 | 9,646 | *Orthohepacivirus* | Hepaciviridae | Riboviria |
| `hiv1_NC_001802.1.fna.gz` | Human immunodeficiency virus 1 | NC_001802.1 | 9,181 | *Lentivirus* | Retroviridae | Riboviria |
| `hiv2_NC_001722.1.fna.gz` | Human immunodeficiency virus 2 | NC_001722.1 | 10,359 | *Lentivirus* | Retroviridae | Riboviria |
| `htlv1_NC_001436.1.fna.gz` | Human T-lymphotropic virus 1 | NC_001436.1 | 8,507 | *Deltaretrovirus* | Retroviridae | Riboviria |
| `hadv5_AC_000008.1.fna.gz` | Human adenovirus 5 | AC_000008.1 | 35,938 | *Mastadenovirus* | Adenoviridae | Varidnaviria |
| `hadv12_NC_001460.1.fna.gz` | Human adenovirus A | NC_001460.1 | 34,125 | *Mastadenovirus* | Adenoviridae | Varidnaviria |
| `hpv16_NC_001526.4.fna.gz` | Human papillomavirus type 16 | NC_001526.4 | 7,906 | *Alphapapillomavirus* | Papillomaviridae | Floreoviria |
| `hpv18_NC_001357.1.fna.gz` | Human papillomavirus - 18 | NC_001357.1 | 7,857 | *Alphapapillomavirus* | Papillomaviridae | Floreoviria |
| `hpv1_NC_001356.1.fna.gz` | Human papillomavirus - 1 | NC_001356.1 | 7,815 | *Mupapillomavirus* | Papillomaviridae | Floreoviria |
| `poliovirus1_NC_002058.3.fna.gz` | Poliovirus | NC_002058.3 | 7,440 | *Enterovirus* | Picornaviridae | Riboviria |
| `coxsackievirus_a16_NC_001612.1.fna.gz` | Human enterovirus A | NC_001612.1 | 7,413 | *Enterovirus* | Picornaviridae | Riboviria |
| `hepatitis_a_NC_001489.1.fna.gz` | Hepatitis A virus | NC_001489.1 | 7,478 | *Hepatovirus* | Picornaviridae | Riboviria |
| `phage_t4_NC_000866.4.fna.gz` | Enterobacteria phage T4 | NC_000866.4 | 168,903 | *Tequatrovirus* | Straboviridae | Duplodnaviria |
| `phage_t7_NC_001604.1.fna.gz` | Enterobacteria phage T7 | NC_001604.1 | 39,937 | *Teseptimavirus* | Autotranscriptaviridae | Duplodnaviria |
| `phage_t3_NC_003298.1.fna.gz` | Enterobacteria phage T3 | NC_003298.1 | 38,208 | *Teetrevirus* | Autotranscriptaviridae | Duplodnaviria |
| `../phix174_GCF_000819615.1.fna.gz` | Escherichia phage phiX174 | NC_001422.1 | 5,386 | *Sinsheimervirus* | Eubullaviridae | Floreoviria |
| `../lambda_phage_GCF_000840245.1.fna.gz` | Enterobacteria phage lambda | NC_001416.1 | 48,502 | *Lambdavirus* | Zimmerviridae | Duplodnaviria |

Relationship ladder: **10 same-genus pairs** (for example HSV-1/HSV-2, vaccinia/cowpox, SARS-CoV-2/MERS-CoV/HCoV-OC43,
dengue/Zika, HIV-1/HIV-2, HPV16/HPV18), **21 same-family, different-genus pairs** (herpesviruses, poxviruses,
coronaviruses, retroviruses, papillomaviruses, picornaviruses, T7/T3), and 434 different-family pairs. Genus and family are
whatever NCBI Taxonomy returned on the download date. That includes recent ICTV renames: phage T7 and T3 are now separate
genera (*Teseptimavirus*/*Teetrevirus*) in *Autotranscriptaviridae*, and HCV sits in *Hepaciviridae*, not *Flaviviridae*.
Select-agent genomes (variola, monkeypox, SARS-CoV-1) were deliberately left out.

## Provenance / re-downloading

For each accession:
1. `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=<accession>&retmode=json` gives the title,
   sequence length and taxid.
2. `https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/<taxid>/dataset_report` gives genus, family and realm.
3. `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<accession>&rettype=fasta&retmode=text`
   gives the sequence, written gzip-compressed as-is.

Every file was verified after download: exactly one FASTA record, a base count equal to NCBI's reported length, and a
byte-identical gzip round trip. RNA viruses are stored as RefSeq provides them (DNA alphabet, T for U).
