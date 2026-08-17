# Results: performance and reliability evidence for ZCurveWindowedProfile

These CSVs are the merged/summary outputs of the benchmark suite in
`SequenceAnaliseTests/Trilateration*.cs` and `GenomeScaleStreamingTests.cs`, covering construction
cost, comparison cost, and reliability against known mutation rates from 20,000 bp up to a
synthetic 3,100,000,000 bp genome-scale reference. They're committed as a snapshot so these figures
can be checked without re-running the benchmarks (some of which take minutes to nearly an hour at
the largest length tested). The "Section" column below refers to a manuscript currently in
preparation and not yet public; it is included here for the authors' own cross-referencing and can
be disregarded by anyone without that document.

Regenerating them from scratch: run the relevant `[Test]` classes via `dotnet test` (see each
file's own doc comment for specifics — several are parameterized across the length sweep
20,000–10,240,000 bp and are individually cheap, but the full sweep plus the `Combine*` merge step
takes on the order of 10–15 minutes; the separate 3,100,000,000 bp genome-scale row
(`GenomeScaleComparisonAndReliability`) is a single ~20–60 minute run on its own, machine-dependent).
Output lands in `SequenceAnaliseTests/bin/Debug/net10.0/`; these are copies of the `Report*.csv`
and `GenomeScale*.csv` files from that folder, not regenerated automatically.

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
| `GenomeScaleCompareTime_3100000000.csv` | §8.2, §8.5 | The 3,100,000,000 bp build/compare/localization timing row |

`ReportCATWithKMP.csv` and `GenomeScaleCeilingProbe.csv` (present in the build output folder but
not here) are excluded: the former belongs to the pre-existing `Tests.cs` experiment suite, unrelated
to this evaluation; the latter is a superseded probe (string-constructor ceiling at 1,070,000,000 bp)
whose result isn't currently referenced anywhere — see the note in
`TrilaterationPerformanceBenchmarks.cs` if reviving it.
