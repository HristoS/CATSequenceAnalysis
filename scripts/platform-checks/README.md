# Platform checks

Measurement scripts behind the platform rows of `results/README.md` (P3, 2026-10-21/22). They read the local catalogue
built by `SequencePlatform.Ingest`, so their numbers depend on what that database holds; the rows say when each was
measured. Analysis-library measurements live in `SequenceAnaliseTests` instead.

Run them from the repository root, with the connection string in the environment (the file-based apps don't see
`SequencePlatform.Ingest/appsettings.json`):

```powershell
$env:SEQUENCE_PLATFORM_DB = "Host=localhost;Port=5432;Username=...;Password=...;Database=sequence_platform"
dotnet run -c Release scripts/platform-checks/bandsearch.cs
```

| Script | Result | What it measures |
|---|---|---|
| `cacheload.cs` | `S33_ProfileCache.csv` | Load time and heap of one cached level; pass the window size (`-- 500`), one level per process |
| `scan.cs` | `S34_ScanTiming.csv` | One find-similar scan per query from in-memory levels, before the length band existed (it re-implements that unbanded rule) |
| `bandsearch.cs` | `S34_BandedSearch.csv` | `SimilaritySearch` timings, cached and streamed. Measured before the neighbour graph; with one built, searches also re-rank |
| `mutual.cs` | `Check_MutualNeighbours.csv` | E. coli, SARS-CoV-2 and HPV16's top 300 re-ranked by mutual rank and shared top-20 |
| `mutual_sample.cs` | `Check_MutualNeighboursSample.csv` | The same over 300 random catalogue queries, re-ranking computed in the script |
| `neighbour_eval.cs` | `S34b_NeighbourRanking.csv` | The same evaluation through `SearchRecordAsync`, with and without the stored graph (`build-neighbours` first) |
| `composition.cs` | `S34b_PlasmidComposition.csv` | Mean Z-curve axes of the Vibrio fischeri plasmid, coronaviruses and others of their length and GC |
| `localize_check.cs` | `S35_LocalizeCheck*.csv` | Coarse-to-fine localization of bundled-genome fragments over the whole catalogue (needs the bundled FASTA files in `data/genomes/`); pass the output CSV path, the summary goes to the console |
| `localize_eval.cs` | `S35_LocalizeEval*.csv` | The same fragments through `Localization.LocalizeAsync` (cached, and a few streamed), for comparison with the check and timing |
| `origin_catalogue.cs` | `S32_OriginCatalogue*.csv` | Replication origin of every stored bacterial chromosome from its profile's cumulative GC skew, against where the record starts |
| `rotated_twins.cs` | `S32_RotatedTwins.csv` | Records that are the same molecule from another start or strand (same length, dinucleotide counts within 2) |
| `rotation_eval.cs` | `S32_RotationRanking*.csv` | S3.4b's 150 bacterial queries by score, as stored and with circular chromosomes rotated to their origin at profile level |
| `timing.cs` | `S37_Timing.csv` | One query of each kind (find-similar by FASTA and by record, localize over the catalogue and inside a record): profiling, first, warm and streamed times |
| `worker_smoke.py` | `S43_WorkerTiming.csv` | Searches through the running API and its worker (P4 S4.3): POST time and POST-to-Done time per query; the worker's own time is in the API log. Needs the API running and an API key file |
| `e2e_scenario.py` | `S45_EndToEnd.csv` | The API's end-to-end scenario (P4 S4.5): `prepare` writes its FASTA bodies from `data/genomes/` into `SequencePlatform.Api/scenario/`; `run BASE_URL KEY_FILE` sends `Scenario.http`'s requests to an API on a fresh database and checks every answer. See the API README |
| `gc_stratum_check.cs` | `S51_GcStratumCheck.csv` | P5 S5.1 gate: AUC of 1 − `Compare` and of GC-only per GTDB rank and \|ΔGC\| stratum on 437 genomes, bootstrap CIs; pass the output path and the level's window size (`-- out.csv 32000`) |
| `r2_candidates.cs` | `S51_Candidates*.csv` | P5 S5.1 (R2): candidate tree distances (strand and origin, coverage, per-genome summaries) against GC and GC plus length on the same sample, with triangle violations; pass the output path and the factor over 500 bp (`-- out.csv 64`) |
| `r3_trees.cs` | `S52_Trees_Catalogue*.csv`, `results/trees/S52_catalogue_*.nwk` | P5 S5.2 (R3): NJ and UPGMA trees from Z-curve distances on the S5.1 sample at 32 kbp, clades recovered per GTDB and NCBI rank with a permuted-label null, and how far each distance follows GC; pass the output path, Newick files go to `trees/` beside it |
| `r1_baselines.cs` | `S52b_Baselines_Catalogue*.csv` | P5 (E1, E2): R1's composition baselines (Karlin δ\*, Natural Vector, Zheng-style curve, Z-curve mean, GC, GC + length) beside `Compare` and the window summary on the S5.1 sample, from stored data only: AUC by rank and GC stratum with bootstrap CIs, NJ / UPGMA clade recovery; pass the output path |
| `e3_combination.cs` | `S52b_E3_KarlinPlusZCurve*.csv` | P5 E3: Karlin δ\* combined with the window summary, against Karlin combined with GC + length as the control, on the S5.1 sample; AUC change with bootstrap CIs and NJ / UPGMA clades; pass the output path |
| `second_order_global.cs` | `S52c_SecondOrderGlobal*.csv` | P5: the global second-order Z-curve from the stored dinucleotide counts against Karlin δ\* on the S5.1 sample, same AUC harness and trees; pass the output path |
| `r4_map.cs` | `S53_Map_*.csv` | P5 S5.3 (R4): PCA map of the whole catalogue (cellular and viral, one genome per species) from the second-order Z-curve, with Karlin and GC + length maps; kNN purity, silhouette and the gap test per rank; pass the output path, the loadings and coordinates go beside it |
| `second_order_log.cs` | `S52d_SecondOrderLog*.csv` | P5: second-order Z-curve variants (independent values only; log ρ\* in Z-curve coordinates) against Karlin δ\* on the S5.1 sample: AUC by rank and GC stratum, NJ / UPGMA clades; pass the output path |
| `second_order_log_neighbours.cs` | `S52d_SecondOrderLog_Neighbours.csv` | P5: 10-nearest-neighbour purity per rank over the whole catalogue (r4_map.cs's corpora) for Karlin ρ\* and the second-order variants, paired bootstrap against Karlin; pass the output path |
| `lineage_tightness.cs` | `S54_LineageTightness*.csv` | P5 measure A: per-lineage purity and spread for the second-order Z-curve and Karlin ρ\*, sign tests by domain, GC class, viral realm and host type; pass the summary path |
| `gtdb_tree_distance.cs` | `S54_GtdbTreeDistance.csv` | P5 measure B: Spearman of signature distances with GTDB reference-tree patristic distance by depth band; pass the folder holding GTDB's `bac120.tree.gz`, `ar53.tree.gz` and `*_taxonomy.tsv.gz`, then the output path |
| `tail_calibration.py` | `S36_TailCalibration.csv` | Hold-out calibration of three tail models on the full `score_null` cells; reads `psql -A -t -c "select source, window_size, scores from score_null where cardinality(scores) >= 10000"` on stdin |
| `column_audit.py` | — | Lists migration columns no platform code references (the 2026-10-22 audit: only timestamps and the unused `search_run` / `search_hit`); pass the repository root |

Record ids (the plasmid is 7787, SARS-CoV-2 19634 in `composition.cs`) are those of the 2026-10 local catalogue.
