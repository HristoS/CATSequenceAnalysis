# Real bacterial genomes for validation

30 complete RefSeq reference genomes, gzip-compressed FASTA as provided by NCBI, downloaded 2026-08-21 from
the NCBI Genomes FTP site (`https://ftp.ncbi.nlm.nih.gov/genomes/all/`). Chosen to give a real-world
phylogenetic distance ladder — several same-species/same-genus pairs, several same-family relatives, and
several different-phylum outgroups — so `ZCurveWindowedProfile`/`CATProfile` similarity scores can be
checked against genuine evolutionary relationships instead of only synthetic random-substitution mutants
(the gap flagged in `docs/zcurve-windowed-profile-factcheck.md`, "Real biological data").

Every accession below was resolved programmatically against the
[NCBI Datasets v2 API](https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/) (not guessed), and every file
was verified after download: gzip integrity, sequence count, and total base count all recomputed directly
from the downloaded bytes rather than taken from the source page.

## Genomes, grouped by taxonomy

### Enterobacteriaceae (Gammaproteobacteria) — the dense cluster

| File | Organism / strain | Assembly | Seqs | Length (bp) |
|---|---|---|---|---|
| `ecoli_k12_mg1655_GCF_000005845.2.fna.gz` | *Escherichia coli* K-12 MG1655 | GCF_000005845.2 | 1 | 4,641,652 |
| `ecoli_o157h7_sakai_GCF_000008865.2.fna.gz` | *Escherichia coli* O157:H7 Sakai | GCF_000008865.2 | 3 | 5,594,605 |
| `shigella_flexneri_2a_301_GCF_000006925.2.fna.gz` | *Shigella flexneri* 2a str. 301 | GCF_000006925.2 | 2 | 4,828,820 |
| `klebsiella_pneumoniae_mgh78578_GCF_000016305.1.fna.gz` | *Klebsiella pneumoniae* MGH 78578 | GCF_000016305.1 | 6 | 5,694,894 |
| `salmonella_typhimurium_lt2_GCF_000006945.2.fna.gz` | *Salmonella enterica* Typhimurium LT2 | GCF_000006945.2 | 2 | 4,951,383 |
| `yersinia_pestis_co92_GCF_000009065.1.fna.gz` | *Yersinia pestis* CO92 | GCF_000009065.1 | 4 | 4,829,855 |
| `yersinia_pseudotuberculosis_ip32953_GCF_000047365.1.fna.gz` | *Yersinia pseudotuberculosis* IP 32953 | GCF_000047365.1 | 3 | 4,840,898 |

### Other Gammaproteobacteria

| File | Organism / strain | Assembly | Seqs | Length (bp) |
|---|---|---|---|---|
| `pseudomonas_aeruginosa_pao1_GCF_000006765.1.fna.gz` | *Pseudomonas aeruginosa* PAO1 | GCF_000006765.1 | 1 | 6,264,404 |
| `vibrio_cholerae_n16961_GCF_900205735.1.fna.gz` | *Vibrio cholerae* O1 El Tor N16961 | GCF_900205735.1 | 2 | 4,033,501 |
| `haemophilus_influenzae_rd_kw20_GCF_000027305.1.fna.gz` | *Haemophilus influenzae* Rd KW20 | GCF_000027305.1 | 1 | 1,830,138 |
| `acinetobacter_baumannii_atcc17978_GCF_013372085.1.fna.gz` | *Acinetobacter baumannii* ATCC 17978 | GCF_013372085.1 | 3 | 4,030,054 |

### Bacillota (Firmicutes)

| File | Organism / strain | Assembly | Seqs | Length (bp) |
|---|---|---|---|---|
| `bsubtilis_168_GCF_000009045.1.fna.gz` | *Bacillus subtilis* 168 | GCF_000009045.1 | 1 | 4,215,606 |
| `bacillus_anthracis_ames_GCF_000007845.1.fna.gz` | *Bacillus anthracis* Ames | GCF_000007845.1 | 1 | 5,227,293 |
| `bacillus_cereus_atcc14579_GCF_006094295.1.fna.gz` | *Bacillus cereus* ATCC 14579 | GCF_006094295.1 | 2 | 5,431,377 |
| `staphylococcus_aureus_nctc8325_GCF_000013425.1.fna.gz` | *Staphylococcus aureus* NCTC 8325 | GCF_000013425.1 | 1 | 2,821,361 |
| `staphylococcus_epidermidis_atcc12228_GCF_000007645.1.fna.gz` | *Staphylococcus epidermidis* ATCC 12228 | GCF_000007645.1 | 7 | 2,564,615 |
| `listeria_monocytogenes_egde_GCF_000196035.1.fna.gz` | *Listeria monocytogenes* EGD-e | GCF_000196035.1 | 1 | 2,944,528 |
| `streptococcus_pneumoniae_tigr4_GCF_000273445.1.fna.gz` | *Streptococcus pneumoniae* TIGR4 | GCF_000273445.1 | 1 | 2,163,340 |
| `streptococcus_pyogenes_m1gas_GCF_043231225.1.fna.gz` | *Streptococcus pyogenes* M1 GAS | GCF_043231225.1 | 1 | 1,853,025 |
| `enterococcus_faecalis_v583_GCF_000007785.1.fna.gz` | *Enterococcus faecalis* V583 | GCF_000007785.1 | 4 | 3,359,974 |
| `clostridioides_difficile_630_GCF_000009205.2.fna.gz` | *Clostridioides difficile* 630 | GCF_000009205.2 | 2 | 4,298,133 |
| `mycoplasma_genitalium_g37_GCF_000027325.1.fna.gz` | *Mycoplasmoides genitalium* G37 | GCF_000027325.1 | 1 | 580,076 |
| `mycoplasma_pneumoniae_m129_GCF_910574535.1.fna.gz` | *Mycoplasmoides pneumoniae* M129 | GCF_910574535.1 | 1 | 816,550 |

### Actinomycetota (Actinobacteria)

| File | Organism / strain | Assembly | Seqs | Length (bp) |
|---|---|---|---|---|
| `mycobacterium_tuberculosis_h37rv_GCF_000195955.2.fna.gz` | *Mycobacterium tuberculosis* H37Rv | GCF_000195955.2 | 1 | 4,411,532 |
| `mycobacterium_smegmatis_mc2155_GCF_000767605.1.fna.gz` | *Mycolicibacterium smegmatis* MC2 155 | GCF_000767605.1 | 1 | 6,988,269 |
| `streptomyces_coelicolor_a3_GCF_008931305.1.fna.gz` | *Streptomyces coelicolor* A3(2) | GCF_008931305.1 | 1 | 8,667,664 |

### Other phyla (distant outgroups)

| File | Organism / strain | Assembly | Seqs | Length (bp) |
|---|---|---|---|---|
| `helicobacter_pylori_26695_GCF_000307795.1.fna.gz` | *Helicobacter pylori* 26695 | GCF_000307795.1 | 1 | 1,667,892 |
| `campylobacter_jejuni_nctc11168_GCF_000009085.1.fna.gz` | *Campylobacter jejuni* NCTC 11168 | GCF_000009085.1 | 1 | 1,641,481 |
| `borrelia_burgdorferi_b31_GCF_000008685.2.fna.gz` | *Borreliella burgdorferi* B31 | GCF_000008685.2 | 22 | 1,521,208 |
| `deinococcus_radiodurans_r1_GCF_020546685.1.fna.gz` | *Deinococcus radiodurans* R1 | GCF_020546685.1 | 4 | 3,279,598 |

### The manuscript's motivating pair (bacteriophages, not bacteria)

| File | Organism | Assembly | Seqs | Length (bp) |
|---|---|---|---|---|
| `phix174_GCF_000819615.1.fna.gz` | *Escherichia phage phiX174* (NC_001422.1) | GCF_000819615.1 | 1 | 5,386 |
| `lambda_phage_GCF_000840245.1.fna.gz` | *Escherichia phage Lambda* (NC_001416.1) | GCF_000840245.1 | 1 | 48,502 |

These two are the real genomes behind §2.3's motivating claim — two unrelated bacteriophages of
similar composition (45% GC vs 50% GC) that `CATProfile.Compare` scores indistinguishably from a
near-duplicate — computed by `PhageMotivatingExampleTests.cs` (see
`docs/zcurve-windowed-profile-factcheck.md`, discrepancy 3, for why this needed fixing: the claim
previously had no backing genome data or test anywhere in the repository). Not part of the
30-genome bacterial corpus below or its taxonomic labeling in `RealGenomeDistanceMatrixBenchmark.cs`.

Two organisms surfaced under names newer than the ones most papers still use — NCBI's current taxonomy has
reclassified *Mycobacterium smegmatis* → *Mycolicibacterium smegmatis*, *Mycoplasma* → *Mycoplasmoides*
(genitalium/pneumoniae), and *Borrelia burgdorferi* → *Borreliella burgdorferi*. Filenames and the
common/legacy names are kept above for readability; `organism_name` in each API lookup (see Provenance)
carries the current name NCBI actually returned.

Total: **30 genomes, 33 MB gzip-compressed**, 580 kbp (*M. genitalium*, one of the smallest known free-living
genomes) to 8.67 Mbp (*S. coelicolor*) — comfortably inside the 20,000–10,240,000 bp range already exercised
by `TrilaterationPerformanceBenchmarks.cs`/`TrilaterationReliabilityAtScaleBenchmarks.cs`, so plain
`string`-based construction is sufficient; none of these need the `FromStream` streaming path.

## Built-in close pairs, for a graded distance ladder

Deliberately included multiple genomes per genus so real phylogenetic distance can be checked at several
depths in one dataset, not just "related" vs. "unrelated":

- **Near-identical (same species, different strain):** *E. coli* K-12 vs. O157:H7.
- **Very close (different genus, but famously near-indistinguishable at genome level):** *E. coli* /
  *Shigella flexneri* — textbook case where 16S/whole-genome similarity blurs the species boundary despite
  the different genus name.
- **Close (same genus, recent divergence):** *Yersinia pestis* vs. *Y. pseudotuberculosis* (*Y. pestis* is a
  clonal derivative of *Y. pseudotuberculosis* on an evolutionary timescale of only ~1,500–20,000 years —
  about as close as two named species get); *Bacillus anthracis* vs. *B. cereus* (same relationship: *B.
  anthracis* is essentially a specialized *B. cereus* lineage).
- **Same genus, moderate divergence:** *Staphylococcus aureus* vs. *S. epidermidis*; *Streptococcus
  pneumoniae* vs. *S. pyogenes*; *Mycoplasmoides genitalium* vs. *M. pneumoniae*; *Mycobacterium
  tuberculosis* vs. *Mycolicibacterium smegmatis* (same family, but *smegmatis* is a fast-growing
  non-pathogen — one of the more diverged same-family pairs here).
- **Same family, different genus:** the full Enterobacteriaceae cluster above (*E. coli*/*Shigella*,
  *Klebsiella*, *Salmonella*, *Yersinia*) against each other.
- **Same phylum, different class/genus:** *Bacillus*/*Staphylococcus*/*Streptococcus*/*Listeria*/
  *Enterococcus*/*Clostridioides* all sit in Bacillota but span a wide internal spread.
- **Different phylum (distant outgroups):** *Mycobacterium*/*Streptomyces* (Actinomycetota),
  *Helicobacter*/*Campylobacter* (Campylobacterota), *Borreliella* (Spirochaetota), *Deinococcus*
  (Deinococcota) against everything in Gammaproteobacteria/Bacillota — the real-genome analogue of the
  phiX174-vs-lambda-phage "unrelated" check already in the manuscript's introduction.

## Provenance / re-downloading

Each `GCF_*` accession and its RefSeq assembly name were resolved via:

```
https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/<organism query>/dataset_report
```

preferring entries marked `"refseq_category": "reference genome"` (or `"representative genome"`) at
`"assembly_level": "Complete Genome"`. The FASTA itself is the `_genomic.fna.gz` product from that assembly's
directory, e.g. for `GCF_000005845.2_ASM584v2`:

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

SHA-256 checksums for every file are not reproduced in full here to keep this table readable; they're
recorded in the download session and can be regenerated with `sha256sum *.fna.gz` — re-verify against a fresh
NCBI download if a checksum ever needs auditing, since RefSeq occasionally revises a genomic FASTA in place
under the same accession.

## Not yet done

These are the raw multi-FASTA downloads only — nothing in this repo parses them yet. Loading (strip headers,
concatenate or select the primary chromosome, decide how to handle the extra plasmid/contig sequences several
of these assemblies carry) and wiring them into an actual `ZCurveWindowedProfile`/`CATProfile` comparison
experiment is separate follow-up work.
