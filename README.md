# CATSequenceAnalysis

Research code for a PhD project introducing **CAT**, a GPS-trilateration-inspired method for
alignment-free DNA sequence similarity profiling (Gancheva & Stoev, 2019, 2024), and its extension
with a windowed Z-curve profile that fixes CAT's large-length discriminative collapse and adds
streaming construction and localization at genome scale.

A full write-up (theory, diagrams, and the empirical evidence behind the Z-curve extension —
performance and reliability from 20,000 bp up to a synthetic 3,100,000,000 bp genome-scale
reference, plus validation against 30 real bacterial genomes) is in preparation as a manuscript and
is not included in this repository ahead of submission. The raw results backing that document's
numbers are in [results/](results/) (see `results/README.md` for which file maps to which reported
result), so every reported figure remains independently reproducible from the code here even without
the write-up itself.

## Project layout

- `SequenceAnalise/` — the library:
  - `WindowedProfile.cs` — an alignment-free method using per-window k-mer frequency vectors.
  - `ZCurveProfile.cs` — `ZCurveProfile` (single-point Z-curve endpoint) and
    `ZCurveWindowedProfile` (the windowed extension this project is building toward publication,
    including `FromStream`/`FromPrecomputedWindows` for genome-scale, non-materialized
    construction).
  - `FastaReader.cs` — gzip-transparent FASTA parser used to load real reference genomes (see
    `data/genomes/`) into `ZCurveWindowedProfile`/`WindowedProfile`, alongside `Generator.cs`'s synthetic
    sequences.
  - `NeedlemanWunsch.cs` — global alignment, used as the accuracy ground truth.
  - `Generator.cs` — (optionally seeded) random DNA generation.
- `SequenceAnaliseTests/` — NUnit 4 tests, two kinds:
  - Correctness tests (`WindowedProfileTests.cs`, `ZCurveProfileTests.cs`,
    `ZCurveWindowedProfileTests.cs`, `GeneratorTests.cs`, `FastaReaderTests.cs`,
    `RealGenomeProfileTests.cs`) — fast, deterministic.
  - Long-running benchmark/experiment suites — `Trilateration*.cs` +
    `GenomeScaleStreamingTests.cs` (the performance and reliability evidence for
    `ZCurveWindowedProfile`; some individual cases run for tens of minutes at the largest lengths
    tested), and `RealGenomeDistanceMatrixBenchmark.cs` (all pairwise comparisons among the 30 real
    genomes in `data/genomes/`, a few minutes).
  - `GenomeCatalog.cs` — locates `data/genomes/` from the test binary's runtime location; not a test
    itself.
- `data/genomes/` — 30 real, complete RefSeq bacterial reference genomes (gzip FASTA), spanning six
  phyla and several deliberate same-species/same-genus pairs; see `data/genomes/README.md` for
  accessions, checksums, and why these 30 in particular.
- `results/` — CSV snapshots of the benchmark results reported above (not currently backed by a
  public write-up; see [How to run tests](#how-to-run-tests) to reproduce them directly).

## Tests

All benchmark reports are written to `SequenceAnaliseTests/bin/Debug/net10.0`; see
`results/README.md` for their column legend and what each file backs.

## How to run tests

Open the `*.sln` file and run the desired test from
[Test Explorer](https://learn.microsoft.com/en-us/visualstudio/test/run-unit-tests-with-test-explorer?view=vs-2022),
or use `dotnet test --filter "Name=<TestName>"`. Several benchmarks in `Trilateration*.cs`/
`GenomeScaleStreamingTests.cs` take tens of minutes at their largest length (3,100,000,000 bp) — see
each file's own doc comment, and `results/README.md`, before running the full sweep.
`RealGenomeDistanceMatrixBenchmark`'s all-pairs sweep across the 30 real genomes in `data/genomes/`
is much cheaper (a few minutes) since none of those genomes exceed ~8.7 Mbp.

> **Warning:**
> For locales like German, where the CSV separator is `;` and the decimal separator is `,` instead
> of `.`, opening these reports in Excel can look wrong —
> [see this thread](https://answers.microsoft.com/en-us/msoffice/forum/all/power-query-how-to-import-a-csv-file-that-does-not/e19fe942-80ad-4ac1-a164-061e7cac9486).

## From CAT to the windowed Z-curve profile

The original CAT method compressed a DNA sequence into a single point via trilateration against
three fixed reference cycles. That single-point construction collapses at large sequence lengths
(unrelated sequences of similar composition become indistinguishable), so it has been replaced by
`ZCurveWindowedProfile` (`SequenceAnalise/ZCurveProfile.cs`), and the CAT implementation itself has
been removed from this repository — it remains available in git history (last present in commit
`e152aa4`). See each class's own doc comments, and `results/README.md`, for how the reported
performance/reliability figures are reproduced. Full theoretical motivation and diagrams are in the
manuscript in preparation (see above).
