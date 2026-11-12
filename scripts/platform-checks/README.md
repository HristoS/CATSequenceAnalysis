# Platform checks

File-based C# apps (and three Python scripts) behind the platform results in `results/README.md`. Most read the local catalogue
built by `SequencePlatform.Ingest`, so their numbers depend on what that database held when they ran. Analysis-library measurements
live in `SequenceAnaliseTests` instead.

Kept here: scripts that produce a result a paper or `CLAUDE.md` cites, the P9 experiments (they publish to the UI), operational
checks, and current work. One-off checks whose decision is made were removed on 2026-11-12 and are at commit `f5cb655`
(`git show f5cb655:scripts/platform-checks/<name>`).

Run from the repository root with the connection string in the environment (file-based apps don't see
`SequencePlatform.Ingest/appsettings.json`):

```powershell
$env:SEQUENCE_PLATFORM_DB = "Host=localhost;Port=5432;Username=...;Password=...;Database=sequence_platform"
dotnet run -c Release scripts/platform-checks/localize_eval.cs
```

| Script | Result | What it does |
|---|---|---|
| `e2e_scenario.py` | — | The API's end-to-end scenario: `prepare` writes FASTA bodies into `SequencePlatform.Api/scenario/`; `run BASE_URL KEY_FILE` replays `Scenario.http` and checks every answer |
| `worker_smoke.py` | — | POST and POST-to-Done times through a running API and worker. Needs an API key file |
| `column_audit.py` | — | Lists migration columns no platform code references; pass the repository root |
| `localize_check.cs` | `S35_LocalizeCheckMin4_Summary.csv` | Coarse-to-fine localization of bundled-genome fragments over the catalogue (needs `data/genomes/`) |
| `localize_eval.cs` | `S35_LocalizeEval_Summary.csv` | The same fragments through `Localization.LocalizeAsync` |
| `neighbour_eval.cs` | `S34b_NeighbourRanking.csv` | Find-similar with and without the stored neighbour graph, 300 queries |
| `rotation_eval.cs` | `S32_RotationProduction.csv` | Find-similar ranking with circular chromosomes rotated to their origin |
| `circular_permutation_check.cs` | `CircularPermutation_Check.csv` | Wrap-around scoring for circularly permuted viruses, pairwise |
| `wrap_search_eval.cs` | `WrapAround_Search.csv` | Wrap-around through `SimilaritySearch`, 300 virus queries |
| `virus_rerank_check.cs` | `VirusRerank_Check.csv` | Whether to skip neighbour re-ranking for viruses |
| `gc_stratum_check.cs` | `S51_GcStratumCheck.csv` | AUC of 1 − `Compare` and GC per GTDB rank and \|ΔGC\| stratum, 437 genomes |
| `second_order_global.cs` | `S52c_SecondOrderGlobal*.csv` | Second-order Z-curve against Karlin δ\* on the same sample |
| `gap_decomposition.cs` | `S52e_GapDecomposition.csv` | Where the second-order Z-curve's gap to Karlin comes from |
| `relatives_check.cs` | `S55_RelativesCheck.csv` | `SignatureIndex`: load time, relatives, query and tree times |
| `protein_search_check.cs` | `ProteinSearch_Catalogue.csv` | First protein-search check at catalogue scale |
| `protein_search_check_v2.cs` | `ProteinSearch_CatalogueV2.csv` | Protein search on the length-bucketed gene layer, 301 proteins |
| `protein_null_check.cs` | `ProteinNull_Check.csv` | Scores of proteins against genes of other families, for p and E |
| `search_breakdown.cs` | `SearchBreakdown_*.csv` | Protein search and streamed localization stage by stage; a later run compares hits for equality |
| `gene_match_check.cs` | `GeneMatch_Check.csv` | Reciprocal best hits by gene score between two related eukaryotes |
| `fetch_eukaryotes.cs` | — | Downloads the eukaryote pairs the splice and gene checks use into `data/genomes/eukaryotes/` (git-ignored), checked against the stored SHA-256 |
| `splice_join_check.cs`, `splice_twist_check.cs` | — | Current work (2026-11-12): splice sites from Z-curve and DNA shape; no database |
| `p9_corpus.cs`, `p9_e1_outliers.cs` … `p9_e5_habitat.cs` | `P9_*.csv` | P9 experiments on the frozen corpus (`shared/ExperimentCorpus.cs`); `--publish` writes them to `ExperimentStore` |
| `p9_corpus.cs` (round 2, default) | `P9_S99_Corpus.csv` | Round-2 corpus on the live database, records ≤ 61,841, eight table hashes checked against `ExperimentCorpus` (D-P9d); `--fingerprint` prints them, `--round1` reruns S9.0 on the snapshot |
| `p9_e1b_check.cs`, `p9_e1b_cvtree.cs` | `P9_E1b_*.csv` | S9.10: E1's keys checked for wholesale renames and the run timed; then E1b on hashed amino-acid CVTree neighbours (round 2) |
| `contrast_check.cs` | `ContrastCheck_*.csv` | Group-contrast runner for visitors (E5 replay, random-halves control, timing, coverage); stopped, kept for its result line |
| `p9_kmer_check.cs`, `insert_check.cs` | `P9_KmerCheck.csv`, `P9_InsertCheck.csv` | P9 follow-ups: longer k-mers for phage host; inserts between two genomes |
