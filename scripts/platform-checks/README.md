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
| `tail_calibration.py` | `S36_TailCalibration.csv` | Hold-out calibration of three tail models on the full `score_null` cells; reads `psql -A -t -c "select source, window_size, scores from score_null where cardinality(scores) >= 10000"` on stdin |
| `column_audit.py` | — | Lists migration columns no platform code references (the 2026-10-22 audit: only timestamps and the unused `search_run` / `search_hit`); pass the repository root |

Record ids (the plasmid is 7787, SARS-CoV-2 19634 in `composition.cs`) are those of the 2026-10 local catalogue.
