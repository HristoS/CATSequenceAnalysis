# CATSequenceAnalysis

Research code for a PhD project introducing **CAT**, a GPS-trilateration-inspired method for
alignment-free DNA sequence similarity profiling (Gancheva & Stoev, 2019, 2024), and its extension
with a windowed Z-curve profile that fixes CAT's large-length discriminative collapse and adds
streaming construction and localization at genome scale.

A full write-up (theory, diagrams, and the empirical evidence behind the Z-curve extension —
performance and reliability from 20,000 bp up to a synthetic 3,100,000,000 bp genome-scale
reference) is in preparation as a manuscript and is not included in this repository ahead of
submission. The raw results backing that document's numbers are in [results/](results/) (see
`results/README.md` for which file maps to which reported result), so every reported figure remains
independently reproducible from the code here even without the write-up itself.

## Project layout

- `SequenceAnalise/` — the library:
  - `DNA.cs` — `CATProfile`, `Benchmark`, `BenchmarkProfile`, `BenchmarkRepository`: the original
    trilateration method.
  - `WindowedProfile.cs` — an alignment-free method using per-window k-mer frequency vectors.
  - `ZCurveProfile.cs` — `ZCurveProfile` (single-point Z-curve endpoint) and
    `ZCurveWindowedProfile` (the windowed extension this project is building toward publication,
    including `FromStream`/`FromPrecomputedWindows` for genome-scale, non-materialized
    construction).
  - `NeedlemanWunsch.cs`, `PatternSearchingKMP.cs` — baseline alignment algorithms.
  - `Generator.cs`, `Search.cs` — sequence generation and the experiment harness.
- `SequenceAnaliseTests/` — NUnit 4 tests, two kinds:
  - Correctness tests (`WindowedProfileTests.cs`, `ZCurveProfileTests.cs`,
    `ZCurveWindowedProfileTests.cs`, `GeneratorTests.cs`) — fast, deterministic.
  - Long-running benchmark/experiment suites — `Tests.cs` (the original CAT vs. Needleman–Wunsch/KMP
    comparison harness; see the report legend below) and `Trilateration*.cs` +
    `GenomeScaleStreamingTests.cs` (the performance and reliability evidence for
    `ZCurveWindowedProfile`; some individual cases run for tens of minutes at the largest lengths
    tested).
- `results/` — CSV snapshots of the benchmark results reported above (not currently backed by a
  public write-up; see [How to run tests](#how-to-run-tests) to reproduce them directly).

## Tests

Three `TestCaseSource` groups in `Tests.cs` parameterize the original CAT experiment harness across
sequence lengths:
  - `SequenceLengthProvider` — source for `CanContributeForFormingTriangle` and
    `CollisionReportFor_1000`.
    - `CanContributeForFormingTriangle` — tests whether sequences of the given length satisfy the
      triangle-inequality precondition trilateration depends on.
    - `CollisionReportFor_1000` — reports how many collisions occur across all sequences of a given
      length.
  - `CompareCATwithNWProvider` — source for `CompareCATwithNW` and `CompareCATwithNWLong`.
    - `CompareCATwithNW`/`CompareCATWithKMP` compare CAT against Needleman–Wunsch or
      Knuth–Morris–Pratt on subsequences (first half, second half, middle) and random-length
      sequences taken from an input sequence.
  - `CompareCATwithNWLongProvider` — source for `CompareCATwithNWLong`, same as `CompareCATwithNW`
    but each experiment runs as its own task for faster execution.

- `CombineCompareCATWithKMP` — run after `CompareCATWithKMP` to merge all per-case reports into one CSV.
- `CombineCompareCATWithNW` — run after `CompareCATwithNW` to merge all per-case reports into one CSV.

All reports from this harness are generated in `SequenceAnaliseTests\bin\Debug\net10.0`. The newer
`Trilateration*.cs`/`GenomeScaleStreamingTests.cs` benchmarks write to the same folder; see
`results/README.md` for their column legend and what each file backs.

### Reports' header legend (`Tests.cs` harness: `CompareCATwithNW`/`CompareCATWithKMP`)

- Dna1 Length / Dna2 Length — length of the 1st/2nd DNA sequence
- Kind — whether the 2nd sequence is a subsequence of the 1st (first half, second half, middle) or
  randomly generated
- CAT Result / NW Result / KMP Result — similarity in [0, 1]
- CAT ms / NW ms / KMP ms — time in ms to perform that comparison
- dna string — raw sequence data
- totalRecords — total combinations for a sequence of the given length
- cat matches — total collisions found
- `[0.1 - 1.0]` — Needleman–Wunsch similarity for a sequence flagged as a CAT collision

## How to run tests

Open the `*.sln` file and run the desired test from
[Test Explorer](https://learn.microsoft.com/en-us/visualstudio/test/run-unit-tests-with-test-explorer?view=vs-2022),
or use `dotnet test --filter "Name=<TestName>"`. Some experiments in `Tests.cs` take up to a couple
of hours depending on the machine and sequence length; several benchmarks in `Trilateration*.cs`/
`GenomeScaleStreamingTests.cs` take tens of minutes at their largest length (3,100,000,000 bp) — see
each file's own doc comment, and `results/README.md`, before running the full sweep.

> **Warning:**
> For locales like German, where the CSV separator is `;` and the decimal separator is `,` instead
> of `.`, opening these reports in Excel can look wrong —
> [see this thread](https://answers.microsoft.com/en-us/msoffice/forum/all/power-query-how-to-import-a-csv-file-that-does-not/e19fe942-80ad-4ac1-a164-061e7cac9486).

## The CAT method, briefly

CAT compresses a DNA sequence into a single point via trilateration against three fixed reference
cycles (`SequenceAnalise/DNA.cs`: `CATProfile`, `Benchmark`). `ZCurveWindowedProfile`
(`SequenceAnalise/ZCurveProfile.cs`) extends this with a windowed Z-curve profile that avoids a
large-sequence-length ceiling present in the original single-point construction; see each class's
own doc comments, and `results/README.md`, for how the reported performance/reliability figures are
reproduced. Full theoretical motivation and diagrams are in the manuscript in preparation (see
above).
