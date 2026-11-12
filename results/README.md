# Results

What the windowed Z-curve profile, the gene profile and the platform built on them have shown, with the numbers. The CSVs kept
here back a manuscript table, a claim in `docs/article/plan.md`, a limit in `docs/article/method-scope.md`, or a statement in
`CLAUDE.md`.

**Removed files.** On 2026-11-12, 114 CSVs, the Newick trees under `trees/`, and the long per-file README were deleted. Their
numbers are summarized below. They are still in git at commit `f5cb655`:

```
git show f5cb655:results/README.md            # the full per-file descriptions
git show f5cb655:results/S34_BandedSearch.csv # any deleted file
git ls-tree --name-only f5cb655 results/      # what was there
```

Code comments that cite a removed CSV refer to that commit.

**Reproducing.** Library results come from `SequenceAnaliseTests` (output in `SequenceAnaliseTests/bin/Debug/net10.0/`, copied
here by hand). Platform results come from `scripts/platform-checks/` or `SequencePlatform.Ingest` on the local database, so they
depend on what the database held that day, and many can no longer be reproduced exactly (the downloads were deleted after ingest).
`raw/` holds the Mash and FastANI outputs that `ExternalBaselineComparisonBenchmark` reads (see `raw/README.md`).

Machine: i7-10510U, 4 cores / 8 threads. Corpora: 30 bundled bacteria (`data/genomes/`), 31 bundled viruses
(`data/genomes/viral/`), and the local catalogue (BV-BRC, about 61,600 records by 2026-11).

## 1. The method on genomes (Paper A, genome part)

Files behind the TCBB manuscript's tables (section numbers are the manuscript's):

| File | § | Content |
|---|---|---|
| `ReportTrilaterationBuildTime.csv` | 8.1 | Profile build time, 20 kbp to 10.24 Mbp. Construction runs at 39.6–49.7 M bases/s |
| `ReportTrilaterationCompareTime_WholeSequence.csv` | 8.2 | Whole-vs-whole compare: 0.064 ms at 10.24 Mbp |
| `ReportTrilaterationLocalizationTime.csv` | 8.2 | `BestAlignment` fragment search: 2.38 ms at 10.24 Mbp |
| `ReportReliabilityAtScale_*.csv` | 8.3 | Score against known mutation rate, 5 seeds per length |
| `ReportTrilaterationAccuracy_*.csv`, `ReportNWCalibration.csv` | 8.4 | Against Needleman–Wunsch, 500 bp to 20 kbp |
| `GenomeScale*_3100000000.csv` | 8.2, 8.3, 8.5 | Synthetic 3.1 Gbp streaming build, planted-fragment search 141.7 ms, profile 74.4 MB (83:1 against the sequence) |
| `RealGenomeDistanceMatrix_*.csv` | 8.6 | 435 pairs of the 30 bacteria; means 0.8187 / 0.7914 / 0.7573 / 0.6063 for same species / genus / phylum / different phylum |
| `ExternalBaseline_*.csv` | 8.7 | Mash and FastANI on the same pairs. Mash saturates for distant pairs while the Z-curve stays graded; r = 0.39 on the 22 pairs above Mash's floor. FastANI returns nothing for most pairs |
| `ReportRetrievalClassification.csv`, `ReportRetrievalOperatingPoint.csv` | 8.8 | **Mash ranks better**: AUC 0.969 vs 0.835 (same genus) and 0.802 vs 0.789 (same phylum); best F1 0.480 vs 0.229 |
| `WindowSizeSensitivity.csv` | 8.9 | 250 / 500 / 1,000 / 2,000 bp windows |
| `ReportPhageMotivatingExample.csv` | 2.3 | CAT scores phiX174 against lambda at about 0.99: the single-point ceiling (CAT removed; commit `e152aa4`) |
| `ReportWindowingStrategy.csv`, `ReportRegionalStructureDetection.csv` | 5.1, 5.3 | Fixed window size, not fixed count, keeps unrelated pairs from drifting together with length |

Caveats carried into the article plan: the timing tables were taken in a Debug build (A4: rerun in Release), and the 3.1 Gbp
compare and reliability rows predate the SIMD rewrite (upper bounds).

Measured after the manuscript:
- **Window count, not window size, decides quality** (`MultiResolution_WindowCount.csv`). Bacteria lose no AUC down to a median of
  about 4 windows in the shorter genome (500 bp → 512 kbp windows: 0.835 → 0.851). Viral same-genus AUC falls gently, 0.881 (18
  windows) → 0.850 (2). Viral families are weak at every level (about 0.60). With 0 windows `Compare` returns 0.
- **The ladder** (`MultiResolution_LadderRule.csv`): 500 bp / 2 / 8 / 32 kbp derived by `Downsample`, compared at the coarsest level
  keeping ≥ 8 windows. Bacteria: same AUC (0.853 / 0.783), all pairs in 8.6 ms instead of 23.5 s. Scores from different levels
  aren't on one scale, so a search ranks at one level per query (0 of 1,485 ranked scores changed).
- Coarse-to-fine localization works when the fragment spans at least 2 coarse windows. Weighting GC by 0.5 nudged AUC up (0.837 /
  0.806); dropping GC collapses it (0.606 / 0.709). int16 window sums were 0.94× (slower) and were removed.
- **GC alone gets close** (`GcBaseline_Retrieval.csv`): GC-only AUC 0.852 vs 0.835 at genus (only 8 positive pairs, CI about ±0.17),
  0.769 vs 0.789 at phylum. On 437 catalogue genomes (`S51_GcStratumCheck.csv`), 1 − `Compare` has Spearman 0.993 with |ΔGC|;
  on same-GC pairs it ties GC at genus (0.614 vs 0.610) and beats it at family (+0.074) and phylum (+0.10 to +0.14).
- **Strand and rotation** (`C3_StrandRotation.csv`): E. coli K-12 against its own reverse complement, half rotation, or both scores
  0.815 / 0.814 / 0.819 at 500 bp, the same as against O157:H7 (0.819). For diverged pairs orientation moves the score by ≤ 0.016.
  The profile-level flip is exact only at lengths that are multiples of 500.
- **Replication origin from the profile** (`C3_OriginRotation.csv`): the minimum of the cumulative GC skew (x − y) finds K-12's oriC.
  Rotating there restores K-12 against itself rotated by half from 0.932 to 0.998 (32 kbp).
- **Moments only where windows are anchored** (`GenomeMoments_*.csv`): Legendre moments on the 500 bp genome ladder leave retrieval
  unchanged (0.830 / 0.794) but drop arbitrary-phase localization from 94.5% to 63.7%. x and y alone cost retrieval (0.606 / 0.709).
  The genome ladder keeps means.
- Triangle inequality: 0 of 12,180 violations on the 30 bacteria, but 2.1% of 41.4 M on 437 catalogue genomes.

## 2. Gene profile and protein search (Paper A, gene part)

How the profile was found:
- **R6 → R6d, no-go** (`R6b_GeneAligned_Summary.csv`, `R6d_Pooled_Summary.csv`). Back-translating a protein and localizing it in
  genome windows reaches ≤ 3.3% top-1 (bacteria). Per-gene profiles over codon positions 1 and 2 reach 75.2% in one genome, but
  pooled over 30 bacteria only 42.5% at 1.5 kbp and 13.2% at 999 bp.
- **R6e/R6f**: the protein turned straight into an expected profile (mean over synonymous codons, exact for 17 amino acids) reaches
  50.2% at 1.5 kbp, which extrapolates to 31% at catalogue scale. Means alone stop there.
- **CDS sweep, the profile adopted** (`CdsSweep_Homolog.csv`, `CdsFollowUp_Search.csv`): per CDS, codon positions 1 and 2, windows
  of 500 bp of gene from the start codon, Legendre moments P0–P2 of x, y, z (9 numbers). From the protein alone: top-1 91.5% among
  16k genes (means: 50.2%), extrapolated 75.7% at catalogue scale (means: 8.1%). Orthologs 97.8 / 94.8 / 88.7 / 87.3% in O157:H7,
  Shigella, Salmonella, Klebsiella. 300 random species (225,827 genes): 97.0% top-1.
- **Storage** (`CdsPrecision_Search.csv`): 16-bit halves change no number (largest difference 0.0001); 8-bit costs 1.8–5.9 points.
- **Short genes** (`ShortProtein_Windows.csv`, `ShortProtein_Homologs.csv`): 50–499 aa get 5 equal windows, compared only with genes
  within ±10% of the length. Orthologs 84–90%. Under 50 aa (`ShortProtein_Under50*.csv`): exact gene 97.2 / 97.0% at 20–34 / 35–49
  aa, orthologs 81.1 / 70.4%; wider bands change nothing. Passes, not built (needs a BV-BRC re-download).

On the catalogue (58,141 BV-BRC records plus 131 annotated eukaryotes, 64.3M genes ≥ 50 aa, 6.14 GB of profiles):
- `ProteinSearch_CatalogueV2.csv`: 301 proteins, from the protein alone. Right family (PGfam) first in 100% in every length band;
  exact gene first 95.0 / 95.2 / 90.0 / 88.1 / 71.7% at 50–99 / 100–166 / 167–332 / 333–499 / 500+ aa. Through the API, median
  5.9 s; hosted (2 vCPU) median 7.4 s with the same accuracy. The first check (`ProteinSearch_Catalogue.csv`) spent 19.6 s of
  20.7 s reading and decoding.
- `ProteinNull_Check.csv`: p and E per hit from 411 real proteins scored against genes of other families.
- `S7T_FullPass.csv`: eukaryotes from their CDS files, 2.25M genes, 209 MB; self-search 96.4% first or tied (yeast 100%).
- `SearchBreakdown_*.csv`: reading rows into pooled buffers and scoring in place cut the time above the transfer floor by 58%
  (end to end −12%, 2.8–7.8 s); PostgreSQL's transfer is now 74–92% of a query. Clustering rows by bucket: read −17–18%. Two exact
  CPU cuts (`_EarlyStop`) gave nothing and were reverted.
- `GeneMatch_Check.csv`: reciprocal best hits by gene score alone, for genes without a family. Same-name agreement 99.8% (Podospora)
  and 92.9% (Nelumbo) against 0.5% and 0.2% by chance.

## 3. The catalogue as a system (Paper B)

- **Bulk load** (`S25_BulkLoad.csv`): 21,946 BV-BRC genomes, 19,964 stored, 1,982 duplicates skipped (12.6% of viruses), 0 failed;
  8.6 h active (1.46 s per genome). Duplicates are found by provider md5 before download and SHA-256 after (identical grouping).
- **Find-similar** compares only lengths within [L/2, 2L] (without the band, scores measure containment: E. coli's top hits were
  270 kbp phages). Warm search 19–76 ms; the band also cuts the cache from 837 MB to 24 MB.
- **Significance is reported, not ranked by.** Empirical p and E from a stored score null; parametric tails were 2–14× off.
- **Shared-neighbour re-ranking** (`S34b_NeighbourRanking.csv`): same-genus precision@10 0.131 → 0.195 for bacteria (p = 0.002),
  0.462 → 0.491 for viruses (300 queries).
- **Origin rotation** (`S32_RotationProduction.csv`): circular chromosomes scored as stored and from their origin, best kept.
  Bacterial precision@10 0.131 → 0.169 by score, 0.195 → 0.231 with neighbours. Warm E. coli search 76 → 382 ms.
- **Wrap-around for viruses** (`CircularPermutation_Check.csv`, `WrapAround_Search.csv`): many phages are circularly permuted. Scoring
  a virus candidate as its windows followed by themselves raises same-genus vs other-family AUC 0.904 → 0.934 and precision@10
  0.597 → 0.670, at 1.15× the time. Skipping neighbour re-ranking for viruses was checked and not taken (`VirusRerank_Check.csv`).
- **Localization** over the whole catalogue at the coarsest level leaving ≥ 4 coarse windows, all query phases, both strands
  (`S35_LocalizeCheckMin4_Summary.csv`, `S35_LocalizeEval_Summary.csv`): source in the top 100 in 30/30 at 10–200 kbp, 23–26/30 at
  4 kbp. Cached 0.17–1.14 s per fragment. On the bundled genomes (`S37_LocalizationValidation.csv`): true start 40/40 at 10 and 50 kbp
  on both strands, 14/20 at 2 kbp inside a named record.
- **Platform equals the library**: find-similar through the database reproduces the library's AUCs exactly once its two designed
  differences (both strands, length band) are added.
- **API**: warm find-similar 150–650 ms POST to done, localization 0.9–2.3 s streamed; the end-to-end scenario passes 10/10 on an
  empty database; `NoDnaTests` finds no trace of a planted marker in any column, log or temp file. Hosted on 2 vCPU / 8 GB:
  find-similar 5 s, localization 19 s.
- **Storage**: database 7.7 GB before the short-gene layer (about 45% of the DNA as 2-bit), 12.2 GB after it. Genome levels as float16 halved the logical size but saved only
  54 MB on disk (TOAST already compressed float32). 85.8% of gene profiles are distinct, so deduplication doesn't pay.

## 4. Organizing organisms (Paper B section, or a note)

- **Whole-genome first-order Z-curve trees are GC trees.** They beat permuted labels, but NJ lays the bacteria out by GC.
- **Karlin ρ\* beats the Z-curve for trees and relatives**: genus AUC 0.966 vs 0.876 for window summaries; UPGMA recovers 76 / 37 / 8
  GTDB genera / families / phyla of 150 / 99 / 24 against 25 / 16 / 6. Within a family Spearman with GTDB distance 0.616 vs 0.543.
- **The complete second-order Z-curve equals Karlin exactly** (`S52e_GapDecomposition.csv`, `S52c_*.csv`): Z-curve means and lag-1
  covariances rebuild the dinucleotide table; log ρ\* in Z-curve axes is a rotation (same trees, map |r| = 1.000000). Dropping the 7
  neighbour coefficients is what made the incomplete form trail. Windowing it adds nothing.
- **Relatives service** (`S55_RelativesCheck.csv`): E. coli K-12's 10 relatives 9/10 Escherichia/Shigella, SARS-CoV-2's 10/10
  coronaviruses, 4–9 ms per query.
- **Map, not a periodic table**: 2-D PCA separates phyla somewhat better than GC + length (kNN purity 0.516 vs 0.465) but not genera.
- **Amino-acid CVTree** tracks GTDB far better than Karlin (within-family Spearman 0.842 vs 0.586, 320 bacteria). Hashed to 8,192
  dimensions it keeps that accuracy.
- **Eukaryote pilot**: on 297 assemblies (one per family) the signature keeps classes and orders together (5-nearest purity 0.53 and
  0.21 against shuffles 0.08 and 0.014) and beats GC + length at every rank.

## 5. What the profile can't do (P9; `docs/article/method-scope.md`)

| Question | Result | File |
|---|---|---|
| E1: flag mislabelled or contaminated genomes | **Fail.** AUC 0.476–0.545, no better than GC + length; flags are endosymbionts and extreme GC | — |
| E1b: the same on amino-acid CVTree, keys fixed (round 2) | **Pass, mostly thanks to the key.** GTDB family differs, minority families only (388 / 7,788): CVTree AUC 0.759 (0.730–0.785), log ρ\* 0.727, GC + length 0.667; top-1% precision 60% vs 19% vs 1% (enrichment 12.0 / 3.9 / 0.3). CheckM poor fails for all arms (0.536–0.553). Genus-isolated positives (102, context): 0.826 / 0.835 / 0.746, enrichment 37 / 7.9 / 0. 74 of the top 200 still have GC < 0.35 or > 0.65 | `p9_e1b_cvtree.cs` |
| E2: a phage's host genus | **Weak pass.** Top-1 0.057 vs GC 0.024, but always guessing Escherichia gives 0.110. Shared-sequence localization: 0/200 | `P9_E2_*.csv` |
| E2 follow-up: longer k-mers | 6-mers lift host top-1 from 0.150 to 0.256; 6-mer counts are now stored | `P9_KmerCheck.csv` |
| E3: genomic islands | **Fail.** Window detector F1 0.128 vs 0.153 for GC alone; x and y dilute GC | `P9_E3_Islands.csv` |
| Inserts between two genomes | 72.5% of inserts ≥ 10 kbp at 2.27 false calls per Mbp; same-GC donors missed | `P9_InsertCheck.csv` |
| E4: CpG depletion from stored counts | **Pass.** Vertebrates' ρ\*(CG) Q3 0.494 below insects' Q1 0.950 and fungi's 0.822 | `P9_E4_CpG_*.csv` |
| E5: habitat beyond genus and GC | **Pass.** 82 of 108 genera above the within-genus shuffle (p 3·10⁻⁸) | `P9_E5_Habitat.csv` |
| Somatic mutations (tumour vs normal) | Invisible: 100 per Mbp moves the score by 0.00028, 9–620× below the closest strain pair | `SomaticMutationRate.csv` |
| Protein search by back-translation | No-go (section 2) | `R6b_*`, `R6d_*` |

Summary: 500 bp windows give identity only over many windows (≥ 10 kbp, localization). Below that they give composition.

Round 2 corpus (S9.9, 2026-11-12): live database, records ≤ 61,841, fingerprint `9a0baea856543d97f121ba45900d82c5` (eight table
hashes, checked by every round-2 script, about 40 s; archive `sequence_platform_2026-11-12_P9b.dump`). 47,306 genomes, 170 excluded
as partial: cellular 15,314 (all with 6-mers, amino-acid 4-mers and gene profiles), viral 31,405 (all 6-mers; 28,570 with genes),
eukaryotic 587 (no 6-mers; 131 with genes). `p9_corpus.cs` writes `P9_S99_Corpus.csv` (not committed).

E1b checks (`p9_e1b_check.cs`, 2026-11-12, round-2 corpus, provider families with ≥ 5 genomes). **E1's GTDB key was mostly
renames**: 87.0% of its 3,030 positives sit in provider families where most genomes are positive (Lysobacteraceae vs Xanthomonadaceae,
Corynebacteriaceae vs Mycobacteriaceae, ...), which no neighbour score can flag; 388 positives in 68 families are a minority in
theirs, and 102 in 33 families are a minority in their genus too (genus-wide moves such as Bacillaceae splits drop out). CheckM poor:
3.9% wholesale, 1,161 minority positives in 203 families. Viral NCBI family differs: 322 positives, 202 minority but in 11 families
(Pestivirus out of Flaviviridae alone is 152), 33 isolated: not testable. Timing at 2,000 genomes (hashed CVTree 8,192, k = 10 cosine,
4 threads): cellular 0.9 s read / 9.8 s build / 9.8 s kNN, viral 0.5 / 8.3 / 9.8 s; extrapolated 10.9 min cellular, 35.3 min viral.

Group-contrast runner for visitors ("can log ρ\* tell A from B?", kNN share AUC with bootstrap CI, shuffle null, GC + length
baseline; `contrast_check.cs`, 2026-11-12): **stopped.** Statistics sound: E5 replay passes (46 genera, excess > 0 in 40, < 0 in 6,
p 1.6·10⁻⁷) and random halves of 121 genera never "separate" (0 %; "beyond GC" alone fires in 5 %, so it must require "separates").
But 2,000 vs 2,000 takes 4.4 s at 4 and at 2 threads (limit 2 s; the sequential bootstrap and shuffles dominate, not the kNN), and only
about 7 real contrasts reach 50 genomes per side: 3 within-genus habitat pairs, 3 phage host-genus pairs; the other 22 counted are
Microviridae sample animals and a common name ("Tomato"). Per-contrast power is low: 6 of 46 replay genera separate on their own.
