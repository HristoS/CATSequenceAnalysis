# Results: performance and reliability evidence for ZCurveWindowedProfile

These CSVs are the merged/summary outputs of the benchmark suite in
`SequenceAnaliseTests/Trilateration*.cs`, `GenomeScaleStreamingTests.cs`,
`RealGenomeDistanceMatrixBenchmark.cs`, `PhageMotivatingExampleTests.cs`,
`WindowingStrategyBenchmark.cs`, and `ExternalBaselineComparisonBenchmark.cs`, covering construction cost, comparison cost, reliability
against known mutation rates from 20,000 bp up to a synthetic 3,100,000,000 bp genome-scale
reference, and (as of 2026-08-25) validation against real bacterial genomes. They're committed as a
snapshot so these figures can be checked without re-running the benchmarks (some of which take
minutes to nearly an hour at the largest length tested). The "Section" column below refers to a
manuscript currently in preparation and not yet public; it is included here for the authors' own
cross-referencing and can be disregarded by anyone without that document.

**2026-09-04, second pass (Mash head-to-head).** `ReportRetrievalClassification.csv` gained a `Method`
column and Mash rows; `ReportRetrievalOperatingPoint.csv` is new. The finding is not favorable to
`ZCurveWindowedProfile` and is reported as measured: at both relationship thresholds **Mash ranks
better** (AUC 0.969 vs. 0.835 and 0.802 vs. 0.789; best-operating-point F1 0.480 vs. 0.229 and 0.749 vs.
0.557), on identical pairs and labels. §8.7's "Mash saturates, this method stays graded" result is
unaffected and still reproduces — but it describes the shape of the score, not retrieval accuracy, and
§8.8 now says so explicitly rather than leaving the stronger reading available. Mash/FastANI raw-file
parsing moved to a shared `ExternalBaselineData.cs` (same rationale as `GenomeTaxonomy.cs`) so both
benchmarks read the committed TSVs through one code path; `ExternalBaseline_*.csv` re-ran bit-for-bit
identical after that refactor.

**2026-09-04 addition (retrieval framing, window-size sweep, off-floor correlation).**
`ReportRetrievalClassification.csv` and `WindowSizeSensitivity.csv` are new, addressing two gaps §9
previously listed as open: a threshold-based retrieval framing (ROC-AUC at two relationship
thresholds, computed by a new `RetrievalClassificationBenchmark.cs`) and a window-size sensitivity
sweep on real genomes (w ∈ {250, 500, 1000, 2000} bp, by a new `WindowSizeSensitivityBenchmark.cs`; a
100 bp point was tried and dropped on cost grounds — see that file's doc comment). `ExternalBaseline_Correlation.csv`
gained a third row, `ZCurve_vs_Mash_OffFloorOnly` (22 pairs, r=0.3919), backing a figure the manuscript's
§8.7 prose already stated but that no committed artifact previously supported
(`docs/zcurve-windowed-profile-factcheck.md`'s 2026-09-04 pass).

**2026-08-29 re-measurement (genome-scale search), superseded by the 2026-08-31 entry below.** `GenomeScalePlantedGeneLocalization_3100000000.csv`
was regenerated after the same optimization pass described below: `SearchMs` dropped from 590 to
371.3 (the fragment query here is only 4 windows regardless of reference size, so this was already
the cheap case, but the array-access fix still removes real per-call overhead across the ~25 million
window-distance computations the search performs). `BuildMs` came out higher this run (116,012 vs.
the previous 71,100) despite that code path being untouched by the optimization — read as ordinary
single-shared-machine run variance (see §8's Methodology note on hardware), not a regression, and
flagged as such in the manuscript rather than silently updated without comment.

**2026-08-29 addition (baselines).** `ExternalBaseline_*.csv` add Mash and FastANI as alignment-free
baselines on the exact same 30-genome corpus and taxonomy labels as `RealGenomeDistanceMatrix_*.csv`,
addressing §9's "no comparison against an established alignment-free baseline" gap. Mash and FastANI
themselves are not .NET code — they were run once, externally, via Docker; their raw output is
committed under `results/raw/` (see `results/raw/README.md` for exact versions/commands), and
`ExternalBaselineComparisonBenchmark.cs` parses that committed output, so the analysis CSVs here are
still regenerable via `dotnet test` without requiring Docker. The headline finding is not simple
agreement: Mash's k-mer/MinHash similarity saturates near zero for this corpus's more distant
(different-phylum) pairs, while `ZCurveWindowedProfile` stays graded across the same range — so
correlation between the two is much stronger when restricted to the pairs where Mash still carries
signal than across the full corpus (both figures are in `ExternalBaseline_Correlation.csv`).
FastANI's own coverage gap (it returns no result at all for genomes below its effective identity
threshold) turned out to be a citable finding in its own right, not a run-configuration problem: see
`ExternalBaseline_Summary.csv`'s `FastAniPairCount` column.

**2026-08-31 regeneration (vectorized comparison).** `ReportTrilaterationBuildTime.csv`,
`ReportTrilaterationCompareTime_WholeSequence.csv`, `ReportTrilaterationLocalizationTime.csv`, and
`ExternalBaseline_*.csv` were regenerated after `ZCurveWindowedProfile` switched its window store
from one interleaved `(X,Y,Z)` array of doubles to three per-axis float arrays, and its alignment
search to explicit SIMD over them. Comparison got substantially cheaper: whole-vs-whole at
10,240,000 bp went 0.3534 ms to 0.06358 ms, fragment localization 4.1809 ms to 2.378 ms, and the
full 435-pair real-genome sweep 133.2 s to 28.6 s, and genome-scale fragment search 371.3 ms to 141.7 ms. Storage halved as a side effect of the float element type: a genome-scale profile is now 74.4 MB rather than 148.8 MB, an 83:1 compression ratio against the raw sequence rather than 42:1 (`GenomeScaleStorageFootprint_3100000000.csv`, whose byte calculation was corrected from `sizeof(double)` to `sizeof(float)` to match). Scores are unchanged to the precision reported
anywhere (the four relationship-level means stay 0.8187 / 0.7914 / 0.7573 / 0.6063); float storage
shifts them by 1e-9 to 1.6e-8, which is far below any reported digit and is safe because each stored
value is a mean of `w` terms from {-1, 0, +1} and so takes one of only `2w+1` values in [-1, 1].

Construction figures in the same file are *not* comparable to the 2026-08-25 snapshot they replace,
and the difference is machine state rather than code. Untouched `CATProfile` construction measured
1.05-2.6x slower than its own published figures across several runs on 2026-08-30, and running the
build benchmark in complete isolation did not recover it, so the earlier session's machine was
simply faster. The numbers here come from a freshly restarted machine, where throughput was finally
stable across all eight lengths (CAT 7.2-8.5M bases/s, ZCurve 39.6-49.7M), and CAT is re-reported
alongside ZCurve so the two stay comparable within this snapshot. On that basis there is no evidence
the per-axis layout hurt construction: ZCurve is slightly slower than the old snapshot at short
lengths and slightly faster at the longest, and the CAT/ZCurve ratio (4.7-6.7x) sits in its usual
band.

Two genome-scale rows are known-stale and flagged as such in the prose: the 3,100,000,000 bp row of
the comparison-cost table and its reliability counterpart both come from
`GenomeScaleComparisonAndReliability`, which was not re-run because it costs the better part of an
hour. They are upper bounds now, not current measurements.

**2026-08-29 optimization pass.** `ExternalBaseline_Runtime.csv`'s `ZCurveWindowedProfile` row was
regenerated after a targeted performance pass on the measured bottleneck: `ZCurveWindowedProfile`'s
`Windows` storage is now backed by a concrete array instead of an interface-typed `List<T>`, and
`BestAlignment`/`Compare` (the hot loop identified above) read that array directly rather than
through the public `IReadOnlyList<T>` property, removing per-access interface-dispatch overhead from
a loop called up to `O(windowCount²)` times per real-genome pair. Same technique applied to
`ExternalBaselineComparisonBenchmark.cs` itself (`CollectionsMarshal.GetValueRefOrAddDefault` for the
FastANI parsing dictionary, pre-sized arrays instead of growable lists for the 30-profile/435-pair
collections, `for` loops in place of `foreach` where the collection is indexable). Every other CSV in
this addition is unchanged and bit-for-bit identical — this was a constant-factor implementation
change, not an algorithmic one; see the corresponding §9 bullet for what a real algorithmic fix would
still need to do. Public API is unaffected (`Windows` still returns `IReadOnlyList<T>`; every other
caller in the codebase is unchanged).

**2026-08-29 addition (motivating claims).** `ReportPhageMotivatingExample.csv` and
`ReportWindowingStrategy.csv` are new: §2.3's phiX174-vs-lambda-phage claim and §5.3's
three-construction drift table previously had no CSV or even test backing them
(`docs/zcurve-windowed-profile-factcheck.md`, discrepancies 3 and 4) — both are now computed by
dedicated tests and committed the same way §8.5's figures were fixed in the 2026-08-25 pass below.
The two phage genomes are committed under `data/genomes/` alongside (not part of) the 30-genome
bacterial corpus; see `data/genomes/README.md`.

**2026-08-25 regeneration.** Every file below (except the real-genome ones, which are new) was
regenerated from a single clean run against the commit that fixed a `Benchmark.Match` regression
(`30a5f65`, "fix the build") — the previous snapshot predated that fix and several of its numbers
(`CATProfile` construction time in particular, which changed substantially) no longer matched the
manuscript prose; see `docs/zcurve-windowed-profile-factcheck.md` for the discrepancies this
regeneration resolves. `GenomeScalePlantedGeneLocalization_*.csv` and
`GenomeScaleStorageFootprint_*.csv` are new in this pass: the tests behind them
(`GenomeScaleStreamingTests.PlantedGeneLocalizationAtGenomeScale` /
`.ProfileStorageFootprintAtGenomeScale`) previously only logged to console, so §8.5's build-time,
search-time, and storage-footprint figures were not actually reproducible from a committed snapshot
despite §11's claim that they were — they are now.

Regenerating them from scratch: run the relevant `[Test]` classes via `dotnet test` (see each
file's own doc comment for specifics — several are parameterized across the length sweep
20,000–10,240,000 bp and are individually cheap, but the full sweep plus the `Combine*` merge step
takes on the order of 10–15 minutes; the separate 3,100,000,000 bp genome-scale row
(`GenomeScaleComparisonAndReliability`) is a single ~20–60 minute run on its own, machine-dependent;
`RealGenomeDistanceMatrixBenchmark`'s all-pairs sweep across 30 real genomes takes a few minutes).
Output lands in `SequenceAnaliseTests/bin/Debug/net10.0/`; these are copies of the `Report*.csv`,
`GenomeScale*.csv`, and `RealGenomeDistanceMatrix*.csv` files from that folder, not regenerated
automatically.

| File | Section | Content |
|---|---|---|
| `ReportTrilaterationBuildTime.csv` | §8.1 | Profile construction time, CAT vs `ZCurveWindowedProfile`, 20,000–10,240,000 bp |
| `ReportTrilaterationCompareTime_WholeSequence.csv` | §8.2 | Whole-vs-whole comparison time, same length sweep |
| `ReportTrilaterationLocalizationTime.csv` | §8.2 | Fragment-localization (`BestAlignment`) search time, same sweep |
| `ReportReliabilityAtScale_Summary.csv` | §8.3 | Known-mutation-rate reliability grid, mean/stdev per condition, 20,000–10,240,000 bp |
| `ReportReliabilityAtScale_Raw.csv` | §8.3 | Same grid, per-seed raw scores (5 seeds/length) |
| `ReportTrilaterationAccuracy_Summary.csv` | §8.4 | Small-scale grid: `ZCurveWindowedProfile` vs Needleman–Wunsch, mean/stdev, 500–20,000 bp |
| `ReportTrilaterationAccuracy_Raw.csv` | §8.4 | Same grid, per-seed raw scores |
| `ReportNWCalibration.csv` | §8.4 | NW per-call timing and Go/SpotCheck/Skip verdict per length (the calibration step that decides §8.4's seed count) |
| `GenomeScaleReliability_3100000000.csv` | §8.3, §8.5 | The 3,100,000,000 bp reliability row (1 seed; see file for why) |
| `GenomeScaleCompareTime_3100000000.csv` | §8.2, §8.5 | The 3,100,000,000 bp whole-vs-whole compare and localization timing row |
| `GenomeScalePlantedGeneLocalization_3100000000.csv` | §8.5 | Streaming build time and planted-fragment localization result at 3,100,000,000 bp |
| `GenomeScaleStorageFootprint_3100000000.csv` | §8.5 | Profile size vs. raw-sequence size and compression ratio at 3,100,000,000 bp |
| `RealGenomeDistanceMatrix_Raw.csv` | §8.6 | All 435 pairwise `ZCurveWindowedProfile.Compare` scores among the 30 real genomes in `data/genomes/`, each labeled by relationship level |
| `RealGenomeDistanceMatrix_Summary.csv` | §8.6 | Same data, aggregated to mean/stdev/min/max per relationship level (same species / same genus / same phylum / different phylum) |
| `ReportPhageMotivatingExample.csv` | §2.3 | `CATProfile.Compare` between the real phiX174 and lambda phage genomes — the false-positive that motivates the whole manuscript |
| `ReportWindowingStrategy.csv` | §5.3 | Unrelated-pair drift (short length vs. long length) under all three profile constructions: raw endpoint, fixed checkpoint count, fixed window size |
| `ReportRegionalStructureDetection.csv` | §5.1 | AT-rich/GC-rich structured sequence vs. its shuffle: endpoint (blind to the shuffle) vs. k=20 trajectory (partially detects it) |
| `ExternalBaseline_Raw.csv` | §8.7 | All 435 pairs: `ZCurveWindowedProfile.Compare`, Mash distance/similarity, and FastANI % identity (where available) side by side |
| `ExternalBaseline_Summary.csv` | §8.7 | Same data, aggregated to mean/stdev per relationship level for all three methods, including each method's pair coverage |
| `ExternalBaseline_Correlation.csv` | §8.7 | Pearson r between `ZCurveWindowedProfile` and each baseline: all 435 pairs, excluding DifferentPhylum, off Mash's resolution floor only (22 pairs), and vs. FastANI (14 pairs) |
| `ExternalBaseline_Runtime.csv` | §8.7 | Wall-clock for building+comparing the full 30-genome corpus: `ZCurveWindowedProfile` (live), Mash and FastANI (external Docker run, see `results/raw/README.md`) |
| `ReportRetrievalClassification.csv` | §8.8 | ROC-AUC (Mann-Whitney rank-sum) for `ZCurveWindowedProfile` **and Mash**, at two binary relationship thresholds over the same 435 pairs: same-genus-or-closer vs. rest, and same-phylum-or-closer vs. different-phylum |
| `ReportRetrievalOperatingPoint.csv` | §8.8 | Precision/recall/F1 at the F1-maximizing cutoff, both methods, both thresholds |
| `WindowSizeSensitivity.csv` | §8.9 | Relationship-level mean/stdev at four window sizes (250/500/1000/2000 bp) on the same 30-genome corpus |

`ReportCATWithKMP.csv` and `GenomeScaleCeilingProbe.csv` (present in the build output folder but
not here) are excluded: the former belongs to the pre-existing `Tests.cs` experiment suite, unrelated
to this evaluation; the latter is a superseded probe (string-constructor ceiling at 1,070,000,000 bp)
whose result isn't currently referenced anywhere — see the note in
`TrilaterationPerformanceBenchmarks.cs` if reviving it.
