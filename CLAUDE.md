# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Research code for a PhD project introducing **CAT** (a GPS-trilateration-inspired method for profiling DNA
sequences so similarity search is fast, roughly O(n), instead of pairwise O(n·m) alignment). Published as:
https://www.mdpi.com/2073-4425/15/3/341. The repo implements CAT plus two baseline algorithms
(Needleman–Wunsch global alignment, KMP exact substring search) and an experiment harness that generates CSV
reports comparing CAT's accuracy/speed against those baselines.

This is experiment/benchmark code supporting a paper, not a library with a stable public API — expect
dead code (`Program.cs` is an unused template, several methods in `Search.cs` and `DNA.cs` are commented out
alternates kept for reference), and NUnit `[Test]` methods that are actually long-running data-generation jobs
rather than correctness assertions.

## Commands

```
dotnet build SequenceAnalise.sln          # build both projects
dotnet test SequenceAnaliseTests          # run everything in Tests.cs (see warning below)
dotnet test --filter "Name=CanContributeForFormingTriangle"   # run a single test/experiment
```

- Open `SequenceAnalise.sln` in Visual Studio and use Test Explorer for the same effect with more visibility
  into per-case progress: https://learn.microsoft.com/en-us/visualstudio/test/run-unit-tests-with-test-explorer
- **Warning:** most "tests" in `SequenceAnaliseTests/Tests.cs` are randomized experiments, not fast
  deterministic checks — several `TestCaseSource` cases go up to sequence length 50,000 and can take minutes to
  hours depending on the machine. Don't run the full suite expecting quick CI-style feedback; target a specific
  test/case when iterating.
- All experiment output (`.csv`) is written to `SequenceAnaliseTests/bin/Debug/net10.0/` — see README.md for the
  column legend for each report type.
- Target framework is net10.0 for both projects (upgraded from net6.0).

## Architecture

**`SequenceAnalise/`** (library)
- `DNA.cs` — the CAT method itself.
  - `Benchmark` / `BenchmarkRepository`: three fixed 4-base reference cycles (A→"ACGT", T→"GTAC", C→"AGTC")
    used as comparison rulers. `Benchmark.Match(index, base)` scores how well a DNA base at a given index
    matches the cyclic benchmark, returning `(near, exact)` in one pass (near-match and exact-match both start
    by checking the same offset-0 position, so they're computed together instead of via two separate calls);
    `NearMatch`/`ExactMatch` remain as thin public wrappers over `Match` for API compatibility. A
    recency-weighted bonus (`prevMatches`) in `CATProfile`'s constructor rewards runs of consecutive matches.
  - `CATProfile`: walks a DNA string once against all three benchmarks to produce a single scalar
    `dnaDistance` per benchmark (A/T/C). The original approach turns these into 2D points via a
    triangulation/trilateration calculation (`BenchmarkProfile.BenchmarkCOS/D/H`, law-of-cosines against fixed
    benchmark separations `acDistance`/1.0), forming a triangle per sequence — but `a`/`t`/`c`'s (D,H)
    coordinates are each computed in a *different* local frame (origin at that benchmark, axis toward a
    different other benchmark), so treating them as one polygon (`CompareAreas`/`BuildPath`) doesn't have a
    clean geometric interpretation; see "known limitations" below. `CATProfile.Point` is a newer addition that
    fixes this via genuine 3-reference trilateration (one point in one shared frame, using all three derived
    benchmark separations including the previously-unused T-C one) — see `CompareTrilateration` below.
  - Two sequences are compared by `CATProfile.Compare` (currently wired to `CompareDistances`, a Euclidean
    distance between the two profiles' per-benchmark D/H values — this one is self-consistent even though the
    three frames don't align with each other, since it never mixes a/t/c together, only same-benchmark
    x-vs-y). `CompareAreas`, using `Clipper2Lib` polygon intersection over the same three points treated as one
    triangle, is the earlier/alternate metric — see the frame caveat above; `BuildPath` is computed lazily now,
    not on every construction.
  - `CanContribute()` checks the triangle-inequality constraints that make the three distances geometrically
    valid as a triangle — this must hold for trilateration to be meaningful, and is what
    `CanContributeForFormingTriangle` tests against (currently only exercised at sequence lengths 10-15; see
    below for why that's a gap, not just a runtime-cost tradeoff).
- `WindowedProfile.cs` — a separate, complementary alignment-free method (does NOT modify `CATProfile`). Splits
  a sequence into fixed-size windows (default 500 bp), represents each window by its k-mer frequency vector
  (default k=4, stored sparsely so memory is O(window size), not 4^k), and compares two profiles by sliding the
  shorter window-vector along the longer (`BestAlignment` → score + window offset; `Compare` → symmetric score).
  This exists because `CATProfile`'s single-point compression provably cannot separate unrelated sequences that
  share base composition — established empirically: real phiX174 vs lambda phage score ~0.99 via `Compare`, and
  every attempt to fix it by adding per-position signal to the 3-scalar profile (self-match, k-mer injections,
  position-tolerant matching) either introduced a length confound or washed toward a constant, because routing
  a sequence through the period-4 all-base benchmark caps the signal at composition. `WindowedProfile` sidesteps
  that by reading k-mer content directly: on real sequences it localizes exact fragments to the correct window
  offset (offset × windowSize ≈ base position), detects diverged fragments down to ~70% identity, and separates
  unrelated sequences (cosine ~0.5–0.67) from true fragments (~0.77–1.0). Tradeoff vs a single CAT point: O(n)
  still, but memory scales with the number of windows, so it suits reads / small-medium genomes / read-vs-genome
  localization rather than whole-large-genome-vs-genome at once. Tested in `WindowedProfileTests.cs`.
- `NeedlemanWunsch.cs` — classic global alignment (match +2 / mismatch −1 / gap −2). Rewritten to O(min(len))
  space via two rolling rows (was a full `int[,]` matrix — ~10GB at 50,000×50,000) and direct `char`
  comparison (was `string.Substring(i,1)` per cell). Verified byte-identical output against the original
  matrix/traceback implementation across thousands of randomized cases. Returns matches / min(len1, len2), not
  a standard identity/alignment-length ratio.
- `PatternSearchingKMP.cs` — textbook KMP exact substring search, used as the "fast but exact-only" baseline.
- `Search.cs` — the experiment harness:
  - `GetPermutations` / `GetPermutationsQ`: enumerate all 4^n DNA strings of a given length (treating a
    sequence as a base-4 number, A=0,C=1,G=2,T=3), used to exhaustively test collisions at small lengths.
    `IterateInParallel` partitions the full 4^n space into `Environment.ProcessorCount` (or an explicit
    `degreeOfParallelism` argument) contiguous shards via `ValueToSequence` (inverse of `GetPermutationsQ`'s
    value encoding) — generalized from an earlier hardcoded 3-way split (A→C, C→G, G→T), which only used 3
    cores regardless of the machine. Verified exact coverage/no duplicates/no gaps across many
    length×parallelism combinations, including edge cases where parallelism exceeds the permutation count.
  - `GenerateCollisionReport`: for random reference sequences (seeded via a `seed` parameter for
    reproducibility — each reference gets `new Random(HashCode.Combine(seed, i))`), finds every same-length
    permutation CAT scores as an exact match (`CATProfile.Compare == 1d`), then runs Needleman–Wunsch on each
    collision to report how genuinely similar those CAT-collisions are — this is the core accuracy evidence for
    CAT. Profiled: for a representative case (length 10, 1000 references), this loop is >99.9% of total runtime
    — CSV writing is a rounding error (~0.04%) by comparison, so don't look to file I/O for speed here.
  - `CompareCATWithNW` / `CompareCATWithKMP`: paired timing + score comparison for a single sequence pair.
  - `Subsequence{FirstHalf,SecondHalf,InTheMiddle}` / `ExactSame`: generate the test-case variants used
    throughout `Tests.cs`.
- `Generator.cs` — DNA string generation via `Random.Shared` by default (was `RandomNumberGenerator`, a crypto
  RNG that added ~10x overhead for no benefit here since generating test data isn't adversarial). Also exposes
  `GetRandomDNA(length, seed)` / `GetRandomDNA(length, Random)` for reproducible runs — wired through
  `GenerateCollisionReport` and all three `Tests.cs` comparison tests, each parallel experiment/task getting
  its own seeded `Random` instance (deterministic per experiment index, safe under concurrency since each
  `Random` is confined to one task). A `FastRng.ThreadSafe` (AES-NI based) package was evaluated as a
  replacement and measured 3-10x *slower* than `Random.Shared` for this access pattern (single `Next(4)` calls)
  — reverted; see git history if revisiting.

**`SequenceAnaliseTests/`** — NUnit 4 project (upgraded from NUnit 3; no test code changes were needed). See
README.md for the full report-column legend and the caveat about CSV locale issues (decimal/list separators)
when opening reports in Excel on non-US locales. `GeneratorTests.cs` has a `GeneratorSpeed` benchmark
(parameterized 1K-1M) that also asserts output length and valid-base-only as a cheap correctness check.

## Known limitations / open questions for the CAT methodology

These were found by direct empirical testing (not just code reading) and are worth resolving or at least
documenting in the paper before relying on CAT at genome scale — flagged for discussion, not fixed unilaterally.

- **Discriminative power collapsed at large N — root cause found, and a fix is now implemented.** Constructing
  `CATProfile` for sequences of increasing length, `dnaDistance` (A/T/C) converged toward a *small number of
  fixed points determined only by periodicity relative to the benchmark's period-4 cycle* — not by sequence
  content. Measured directly: the per-position outcome (which of the 4 benchmark offsets matches) landed at
  ~25%/25%/25%/25% for uniform-random, homopolymers of *any* base, AT-rich, and GC-rich compositions alike, and
  even for a deterministic period-3 repeat (`AAT AAT...`) — composition is invisible to the benchmark-offset
  check whenever the sequence's own period is coprime with 4. Only period-2/period-4-resonant patterns
  (`ATAT...`) diverged. This ruled out memory span as the cause: sliding-window and exponential-decay variants
  of the `prevMatchesX` recurrence were tested and made no difference to the spread between composition
  categories (only to noise), because they only change how a per-position outcome is averaged, not what that
  outcome depends on. **Fix:** added `selfMatchWeight` (`CATProfile.cs`, currently `0.3d`) — an "adjacent base in
  the DNA sequence itself matches the previous one" signal (`dnaString[i] == dnaString[i-1]`), added to the
  near-match score for all three benchmarks each step. This is a real, principled statistic (its long-run rate
  is `Σp_base²` over the base frequencies — the same quantity used as a homozygosity/repeat-content measure in
  population genetics) that the benchmark-offset check never captured. Verified at length 1,000,000: a random
  vs. 70%-GC-biased sequence went from indistinguishable (`Compare` ≈ 0.9998 either way) to a real, ~46-standard-
  deviation separation (`Compare(random, GC-rich)` ≈ 0.984 vs. `Compare(random, unrelated-same-composition)` ≈
  0.999); near-duplicate detection (97%-similar sequences) is undisturbed (still ≈ 0.9999); `CanContribute()`
  still holds for every composition tested, including exhaustively across all 4^10 length-10 sequences
  (`CanContributeForFormingTriangle(10)`). `lookback=1` (immediate neighbor) was clearly the strongest choice —
  lookback 2 or 3 gave much weaker separation in testing.
- **`CanContribute()` (the triangle-inequality precondition for trilateration) is only tested at length 10-15**
  (`Tests.cs`'s `SequenceLengthProvider`), and only against uniform-random sequences. It held in every
  configuration tested here up to length 10,000,000, including skewed compositions and repeats — but that's a
  different (looser) claim than "proven to hold in general," and real genomes have composition/structure well
  outside i.i.d. uniform.
- **Benchmark separations are derived, not arbitrary — confirmed and now used properly.** Comparing the
  benchmark strings position-by-position: A="ACGT" vs T="GTAC" share 0 of 4 positions -> distance 1.0; A vs
  C="AGTC" and T vs C each share exactly 1 of 4 -> distance 0.75. That's `1 - (shared positions)/4`, and A-T-C
  therefore form a fixed, valid triangle (0.75+0.75>1, 1+0.75>0.75). The original code only ever used the A-T
  and A-C separations (`c.Calculate(a)`, never `c.Calculate(t)`) — the derivable T-C=0.75 relationship sat
  unused, and `a`/`t`/`c` were each computed in their own local 2D frame (different origin, different reference
  axis) via independent 2-point law-of-cosines projections rather than one shared frame. `CompareDistances`
  sidesteps the frame mismatch (it only ever compares "x's `a`" to "y's `a`", same convention both sides), but
  `CompareAreas`/`BuildPath` combined the three into one polygon as if co-planar, which didn't have a clean
  geometric interpretation.
  **Added a real fix for this, alongside the original (not replacing it yet):** `CATProfile.Point` places A/T/C
  at fixed coordinates using all three derived separations (A=(0,0), T=(1,0), C solved from |A-C|=0.75 and
  |T-C|=0.75), then locates each sequence via standard GPS-style trilateration linearization from its three
  measured distances — one point in one shared frame instead of three separately-framed projections.
  `CompareTrilateration(x, y)` compares two profiles via straight Euclidean distance between their points.
  Verified: identical sequences produce identical points (`CompareTrilateration` = 1 exactly); 2000 randomized
  regression trials confirm this addition didn't disturb `Compare`/`CompareAreas`/`CanContribute`; and it's
  measurably *more* discriminative at moderate length (length 1000, 97%-similar vs. unrelated sequence:
  0.999086 vs. 0.971314 for trilateration, vs. 0.999004 vs. 0.977054 for the original `CompareDistances` — a
  wider separation between "similar" and "unrelated"). `Compare()` is still wired to `CompareDistances` pending
  a decision on whether to switch it — `CompareTrilateration` is available for A/B use now via
  `CATProfile.CompareTrilateration`/`.Point`.
- **Magic constants without documented derivation:** the near-match weights `baseDistance = { 0, 0.6, 0.4, 0.6 }`
  in `Benchmark` (unlike the benchmark separations above, no derivation for these has been found yet). May well
  be empirically tuned in the published paper — worth cross-referencing so future changes don't silently drift
  from whatever justified the original values.
- **Empty-string input produces NaN** throughout (`dnaDistance = 0/0`) since the near/exact-match accumulation
  loop never runs — not handled or validated against.

## Scale: can this process real DNA sequences?

Researched (not assumed) because the current implementation stores sequences as `string`/`char` and indexes
with `int` throughout (`CATProfile`, `NeedlemanWunsch`, KMP):

- **.NET string capacity:** capped at ~1.07 billion UTF-16 characters (2GB object limit ÷ 2 bytes/char).
  `gcAllowVeryLargeObjects` raises the array-length cap but explicitly does *not* raise the string cap.
- **Real sequence sizes, for scale:**
  - Typical virus: 7-20 kbp; smallest (circoviruses) ~1.8-4.5 kb; largest ("giant viruses" like Pandoravirus)
    up to ~2.5 Mbp.
  - Bacteria: ~0.5 Mbp to >10 Mbp typically (*E. coli* ~4.6 Mbp); largest known prokaryote genome
    (*Ktedonobacter racemifer*) ~13.7 Mbp.
  - Human: largest chromosome (chr1) ~249 Mbp; full genome ~3.05 Gbp gapless (T2T-CHM13) to ~3.24 Gbp
    (GRCh38.p10, with gaps) — **this already exceeds the .NET string capacity by ~3x**, so a whole human genome
    cannot be loaded into a single `string` with the current design, only a chromosome-at-a-time.
  - Largest known genome (2024 record): the fork fern *Tmesipteris oblanceolata*, ~160 Gbp — over 100x the
    .NET string cap.
  - Sequencer read lengths (the more immediate practical unit, since raw reads are what's actually processed
    before assembly): Illumina short reads 50-300 bp; PacBio HiFi ~10-25 kb; Oxford Nanopore from short up to
    ultra-long (>4 Mb reported), with typical N50 around 35 kb on PromethION.
- **Net takeaway:** chromosome-scale or smaller (up to ~250M bases) fits comfortably within current types with
  room to spare. Whole-large-genome-scale (multi-Gbp) needs a different representation entirely (packed 2-bit
  encoding instead of UTF-16 `char`, `long` indices, and likely a streaming/chunked construction of
  `CATProfile` rather than materializing the whole sequence at once) *before* the string-capacity ceiling is
  even the binding constraint — and per the discriminative-power finding above, resolving numerical behavior
  at scale matters at least as much as fixing storage.

Sources: [Human genome — Wikipedia](https://en.wikipedia.org/wiki/Human_genome),
[The complete sequence of a human genome (T2T) — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9186530/),
[Chromosome 1 — Jeen Health](https://jeen.health/library/chromosomes/1/),
[gcAllowVeryLargeObjects — Microsoft Learn](https://learn.microsoft.com/en-us/dotnet/framework/configure-apps/file-schema/runtime/gcallowverylargeobjects-element),
[Genome size — Wikipedia](https://en.wikipedia.org/wiki/Genome_size),
[Daily briefing: tiny fork fern has the biggest genome ever found — Nature](https://www.nature.com/articles/d41586-024-01674-5),
[Short, long, or ultra-long read length — Oxford Nanopore](https://nanoporetech.com/blog/short-long-or-ultra-long-which-read-length-is-right-for-you),
[PacBio HiFi vs Nanopore vs Illumina — CD Genomics](https://www.cd-genomics.com/resource-pacbio-hifi-vs-nanopore-vs-illumina.html),
[E. coli genome size — Biology Insights](https://biologyinsights.com/what-is-the-e-coli-genome-size/),
[Virus Genome — ScienceDirect](https://www.sciencedirect.com/topics/immunology-and-microbiology/virus-genome).
