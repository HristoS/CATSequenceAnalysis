# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Research code for a PhD project. It began with **CAT** (a GPS-trilateration-inspired method for profiling DNA
sequences so similarity search is roughly O(n) instead of pairwise O(n·m) alignment), published as
https://www.mdpi.com/2073-4425/15/3/341. CAT has since been superseded by a **windowed Z-curve profile**
(`ZCurveWindowedProfile`), and the CAT implementation (`CATProfile`, `Benchmark*`, `DNA`, the `Search`
collision/permutation harness, KMP baseline) has been **removed** from the codebase — it's in git history, last
present in commit `e152aa4`. Don't reintroduce it; check out that commit if old CAT figures need reproducing.

This is experiment/benchmark code supporting a manuscript in preparation, not a library with a stable public
API — expect NUnit `[Test]` methods that are actually long-running data-generation jobs rather than correctness
assertions.

## Commands

```
dotnet build SequenceAnalise.sln          # build every project
dotnet test SequenceAnaliseTests --filter "FullyQualifiedName~ZCurveWindowedProfileTests"   # one class
dotnet test SequenceAnaliseTests -- 'NUnit.Where=test =~ /BuildTime\(20000\)/'             # one parameterized case
dotnet build SequencePlatform.slnf        # platform projects only (library, Core, Api, Ingest, Tests)
dotnet test SequencePlatform.Tests        # platform tests, ~1 s; need a reachable PostgreSQL (see below)
dotnet run --project SequencePlatform.AppHost   # API, worker and dashboard under Aspire, on the local database
```

- Open `SequenceAnalise.sln` in Visual Studio and use Test Explorer for more visibility into per-case progress.
- **Warning:** the `*Benchmark*.cs`, `Trilateration*.cs`, and `GenomeScaleStreamingTests.cs` suites are
  experiments, not fast checks — some cases run tens of minutes (up to a synthetic 3.1 Gbp reference). Don't run
  the full suite expecting CI-style feedback. The fast deterministic correctness tests are
  `ZCurveProfileTests`, `ZCurveWindowedProfileTests`, `ZCurveGeneProfileTests`, `WindowedProfileTests`, `FastaReaderTests`,
  `GeneratorTests`, `PhylogenyTests`, `GenomeSignaturesTests`, `OrdinationTests`, and `RealGenomeProfileTests`.
- `--filter "Name=Foo(123)"` does not work for parameterized cases with this adapter; use `NUnit.Where` as above.
- All experiment output (`.csv`) is written to `SequenceAnaliseTests/bin/Debug/net10.0/`; committed snapshots live
  in `results/` — see `results/README.md` for which file backs which manuscript section.
- Target framework is net10.0 for every project. The library (`SequenceAnalise`) has no NuGet dependencies; the
  platform projects use Npgsql.

## Architecture

**`SequenceAnalise/`** (library)
- `ZCurveProfile.cs`
  - `ZCurveProfile`: maps each base onto three biochemical dichotomies (purine/pyrimidine, amino/keto,
    strong/weak H-bond) and samples the running average at `CheckpointCount` points. With one checkpoint it's
    pure composition, so it can't separate same-composition sequences — kept as the "endpoint" baseline.
  - `ZCurveWindowedProfile`: **the main method.** Fixed-size windows (default 500 bp), each summarized by its
    local Z-curve (x,y,z) average; compared by sliding the shorter profile along the longer (`BestAlignment` →
    score + window offset, `Compare` → symmetric score). Fixed window *size* (not count) keeps each window's
    noise floor constant, which is what avoids the large-N collapse that single-point methods (including CAT)
    hit. Storage is structure-of-arrays `float[]` per axis and `BestAlignment` is explicitly SIMD-vectorized
    (two shapes: wide query vectorizes across windows, narrow query across offsets) — see the in-code
    comments before touching the layout; scalar/interleaved alternatives were measured slower.
    `FromStream` builds from a `TextReader` without materializing a string (genome scale, past the .NET string
    cap); `FromPrecomputedWindows` rehydrates a cached profile. `Downsample(factor, phase)` derives the
    `WindowSize*factor` profile from an existing one (equal to building it from DNA up to float rounding), which is
    what `MultiResolutionBenchmark` uses. What matters is the *shorter* sequence's window count, not the window size:
    AUC holds while it keeps about 4–8+ windows, and at 0 windows `Compare` returns 0 (looks "unrelated").
    `ZCurveProfileLadder` (same file) stores 500 bp / 2 / 8 / 32 kbp levels and compares at the coarsest level keeping ≥ 8
    windows, returning null when not comparable. Scores from different levels aren't on a common scale, so a search uses
    `CompareForQuery` (one level chosen by the query; too-short candidates go to a separate short list). The viral test corpus is
    `data/genomes/viral/` (taxonomy: `SequenceAnaliseTests/ViralTaxonomy.cs`). int16-sum storage was
    prototyped and was *slower* (0.94x), so don't retry it without a new idea.
- `ZCurveGeneProfile.cs` — the per-gene profile for protein search (S7.R, S7.S): each CDS in coding orientation, codon positions 1
  and 2, Legendre moments P0–P2 of x, y, z per window (9 numbers, stored as 16-bit halves; 8 bits lost up to 6 points). Window rule
  by protein length: from 500 aa, windows of 333 such bases (500 bp of gene) from the start codon; 50–499 aa, 5 equal windows over
  the whole gene, compared only with genes within ±10% of the query's length (`results/ShortProtein_*.csv`); under 50 aa not kept.
  Binary format ZCGP v2 (varint-packed gene headers; the window count follows from the length); v3 (S7.T, NCBI eukaryotes) adds per gene
  the sequence it lies on and its introns, and is written only when some gene needs them, so BV-BRC rows stay v2. Queries come from a gene
  (`QueryFromGene`) or a protein (`QueryFromProtein`, expected axes over synonymous codons). Moments work only because windows are
  anchored at the start codon; on the unanchored genome ladder they hurt (`results/GenomeMoments_*.csv`), so `ZCurveWindowedProfile`
  keeps means.
- `WindowedProfile.cs` — same windowed/sliding-alignment structure, but each window is a sparse k-mer frequency
  vector (default k=4) compared by cosine. Richer signal, heavier memory; used as a comparison point.
- `FastaReader.cs` — gzip-transparent FASTA parser for the real genomes in `data/genomes/`. Passes IUPAC
  ambiguity codes through unchanged (both profile types treat them as zero-contribution).
- `Phylogeny.cs` — neighbour-joining (midpoint-rooted for display) and UPGMA trees, Newick, leaf order, and per-taxon
  clade recovery judged on the unrooted tree (P5). Validated on the textbook examples in `PhylogenyTests`.
- `GenomeSignatures.cs` — whole-genome composition signatures (P5 baselines): window summary, Z-curve mean, Zheng-style
  curve, Natural Vector, Karlin ρ\*/δ\*, and the second-order Z-curve (lag-1 axis correlations: dinucleotides in Z-curve form;
  `SecondOrderLog` is the interaction part of log ρ\* in Z-curve axes). `ToZCurveBasis` rotates any dinucleotide vector into Z-curve
  coordinates exactly; with all 16 coefficients the Z-curve form equals Karlin, while the interaction-only forms trail it (they drop
  neighbour information, not composition). The profile-derived ones solve each window's base counts from its Z-curve averages.
- `Ordination.cs` — PCA and a general symmetric eigen-decomposition (cyclic Jacobi), used by the P5 map and `GenomeSignatures`.
- `NeedlemanWunsch.cs` — global alignment (match +2 / mismatch −1 / gap −2), O(min(len)) space via two rolling
  rows. Returns matches / min(len1, len2), not a standard identity ratio. Used as accuracy ground truth at
  lengths where O(n·m) is affordable (`TrilaterationAccuracyBenchmarks` calibrates that).
- `Generator.cs` — random DNA via `Random.Shared`, or seeded (`GetRandomDNA(length, seed)` /
  `GetRandomDNA(length, Random)`) for reproducible experiments. Base order (A,T,C,G) is fixed — changing it
  changes every seeded sequence and therefore every committed result.

**Platform** (`docs/database-decision.md` rev 2). **The platform never stores DNA sequences.**
- `SequencePlatform.Core/` (library, Npgsql): `Database` (connection settings, create-if-missing, migrations), the schema
  as embedded `Migrations/*.sql` (append-only once applied; versions are the file names), `ProfileLevels` (ladder levels as
  ZCWP bytes inline up to 1 MB: insert in the caller's transaction, read one, stream one level of all sequences; written as ZCWP
  version 2, float16, since 2026-11-05, from `RecordStore.Ladder`, which rounds every level so a query's ladder equals the stored one;
  `convert-levels-half` rewrote the old rows. On disk that saved almost nothing, since TOAST already compressed float32 by half),
  `Catalogue` (keyset-paged record listing with kind/rank/length/provider/name filters), `ConnectionString`, `ProfileCache`
  (one level, optionally one length range, in memory, refreshed incrementally), `SimilaritySearch` (find-similar: only records
  of length [rL, L/r], default r = 0.5, at one level, both query strands; unbanded scores measure containment, not
  similarity; each hit carries an empirical p and E from `score_null`, built by `ScoreNullBuilder` / `build-null`; the top
  100 are re-ranked by shared top-20 neighbours from `neighbour_list` (`NeighbourGraph` / `build-neighbours`), else by score), `Localization` (localize a fragment ≥ 4 kbp over all records at least as long: coarsest level leaving it 4 windows,
  every query phase, both strands, top 100 refined at 500 bp; shorter fragments only inside a named record), `OriginRotation`
  (circular bacterial chromosomes ≥ 500 kbp laid out from their replication origin, found from the stored profile's cumulative GC
  skew; their 2/8/32 kbp levels stored under method version `1-origin` by `build-origins`; find-similar and the genome null score
  both versions of a record and of the query and keep the best; localization stays unrotated), `SearchLog` (a search as an asynchronous job, migration 0009: keeps both strands' fine query profile and
  SHA-256, never DNA; result lists ranked / identical / not_comparable; D7 expiry after 30 days or on DELETE keeps only hashes,
  length and summary), `GeneProfiles` / `ProteinKmers` (the gene layer: since S7.S, migration 0017, `ZCurveGeneProfile` rows in
  `gene_profile` by length bucket, 0 for genes of 500 aa and up and 1–25 for 50 · 1.1^b aa, replacing 0014's `zcurve-cds`
  profile_level row (0018 deleted the last ones); amino-acid 4-mer counts (3-mers dropped by 0016: `ProteinKmers.Shorter` derives them); BV-BRC protein
  families; all filled by `backfill-genes`, and for annotated NCBI eukaryotes by `backfill-genes-ncbi` from the CDS file, with the
  sequences genes lie on listed in `annotations.ncbi.gene_layer` and kept on protein hits as `search_hit.gene_sequence`, migration 0020),
  `ProteinSearch` (P7, migrations 0015 and 0017: a protein ≥ 50 aa, never stored, as expected Z-curve windows, compared with the long
  genes or with the genes within ±10% of its length, reading only those buckets from the database and scoring them in place (`ZCurveGeneProfile.ScoreRow`, no decode and no
  per-query hash: `verify-genes` checks the rows; since then PostgreSQL's transfer is most of a query, `results/SearchBreakdown_*.csv`, and `cluster-genes` keeps each bucket's rows together); each hit carries p and E against genes of other families, `ProteinScoreNull`, migration 0019, built by `build-protein-null`).
- `SequencePlatform.Api/`: ASP.NET Core minimal API on localhost:5080 (P4; see its README). Per-user API keys (`Users` in Core,
  only the SHA-256 stored); `CatalogueFilter.ViewerId` makes an upload visible to its owner only, and no viewer means public records
  only. The neighbour graph and score null are built from public records only. Uploads (`POST /api/sequences`) take the raw
  FASTA or ZCWP body, never multipart (`IFormFile` buffers to a temp file); FASTA streams through `FastaReader.ReadRecordStreams`.
  Searches (`/api/searches/similar`, `/localize`, `/protein`) are `SearchLog` jobs: 202 + `Location`, then GET pages the hit lists.
  `SignatureIndex` (Core) serves `GET /api/sequences/{id}/relatives`, `/signature?view=zcurve|karlin` and `GET /api/trees`
  synchronously from an in-memory table of dinucleotide signatures (log ρ\* = complete second-order Z-curve; P5 S5.5).
  `SearchWorker` (a hosted service in the API) claims them with `for update skip locked`, wakes on `LISTEN search_queued`, runs
  one at a time, and requeues jobs left running at startup (single worker until P6). Tests run the API in process
  (`ApiTests`, `UploadTests`, `SearchApiTests` with the worker off; `WorkerTests` with it on). `NoDnaTests` runs it on Kestrel
  and scans every column, log entry and temp file for a planted marker after upload and search (S4.4). `Scenario.http` is the
  end-to-end walk-through; `scripts/platform-checks/e2e_scenario.py` makes its bodies (`prepare`) and replays it with checks (`run`).
  `Landscape` (Core) serves `GET /api/landscape`: PCA of the signature, public records only, one genome per species per group, or `all`
  (exactly what `SignatureIndex` holds), one map per group because the Z-curve and Karlin views are a rotation. `ExperimentStore` (Core,
  migration 0012) holds P9 experiment results that `scripts/platform-checks/p9_*.cs --publish` computed on the frozen snapshot
  `sequence_platform_p9` (shared loader `scripts/platform-checks/shared/ExperimentCorpus.cs`, via `#:include`); `/api/experiments` reads them. The web UI is `SequencePlatform.Api/wwwroot/`: static files in front of the API, React 18.3.1
  and htm vendored with no build step and no Node (React 19 has no UMD builds), hash routes, `default-src 'self'`.
- `SequencePlatform.AppHost/` and `SequencePlatform.ServiceDefaults/` (P6): Aspire. Run mode starts `migrate`, then `api` and `worker` as two
  processes of the API project on the local database; publish mode writes Docker Compose with its own empty `postgres:17`, where the API and
  worker connect as the least-privilege role `sequence_app` (migration 0021: reads all, writes only uploads and search jobs; `migrate` sets its
  generated password; the platform tests run the API as it, via `-c role=`). No broker:
  PostgreSQL stays the queue. Exactly one worker. Generated output goes to the git-ignored `deploy/` (its `.env` holds a password).
- `SequencePlatform.Ingest/`: console app that fills the DB from BV-BRC: download → profile (`FromStream`, never a whole
  string) → store ladder levels and annotations → delete the file. `ingest-bvbrc --per-species` adds one genome per species not yet stored;
  `ingest-ncbi` is the eukaryote pilot: one record per assembly with summed nuclear dinucleotide counts and no profile levels, streamed
  with no temp file, human refused in code; annotated assemblies also get the gene layer (no 6-mers: those would need the genome again). Commands `migrate`, `ingest-bvbrc`, `ingest-ncbi`, `load-taxonomy`, `stats`, `list`,
  `profile`, `search-similar`, `localize`, `build-null`, `build-neighbours`, `build-origins`, `convert-levels-half`, `backfill-kmers`,
  `backfill-genes`, `backfill-genes-ncbi`, `build-protein-null`, `verify-genes`, `cluster-genes`, `add-user`, `disable-user`; see its README.
  6-mer counts per record (`record_kmer`, migration 0013, D8 widened 2026-11-04) are stored at ingest; `backfill-kmers` re-downloads older
  records and stores counts only where the SHA-256 matches. Measurement scripts for the platform's
  `results/` rows: `scripts/platform-checks/` (run from the repo root with `SEQUENCE_PLATFORM_DB` set).
- `SequencePlatform.Tests/` (NUnit 4): each run creates `sequence_platform_test_<random>`, migrates it, and drops it.
- Connection string from the `SEQUENCE_PLATFORM_DB` env var, else the git-ignored `SequencePlatform.Ingest/appsettings.json`
  (`ConnectionStrings:SequencePlatform`; the test project links the same file); never commit it. Local development DB:
  `sequence_platform` on the local PostgreSQL 17. A repo-level `nuget.config` restricts restore to nuget.org, because the
  machine-wide config has an unreachable private feed.
- Platform work is planned in phases P0–P7: tracker `docs/platform/README.md` (git-ignored, local only); start sessions with `/platform-session`.

**`SequenceAnaliseTests/`** — NUnit 4. Besides the correctness tests listed above: `Trilateration*` (build/compare
/localization timing, NW-ground-truth accuracy grid, known-mutation reliability at scale), `GenomeScaleStreaming`
(3.1 Gbp streaming build, planted-gene localization, storage footprint), `RealGenomeDistanceMatrix` /
`RetrievalClassification` / `WindowSizeSensitivity` / `ExternalBaselineComparison` (30-genome real-data
evaluation, incl. Mash/FastANI from committed TSVs in `results/raw/`), `WindowingStrategy` and
`DistanceMetricComparison` (design-choice evidence). Shared helpers: `BenchmarkTiming`, `AccuracyBenchmarkShared`,
`GenomeCatalog`, `GenomeTaxonomy`, `ExternalBaselineData`.

## Background: why CAT was replaced

Found empirically before the removal, kept here because the Z-curve design choices only make sense against it:
CAT walked the sequence against three period-4 benchmark strings and reduced it to three scalars → one point.
Any such single-point compression converges, as length grows, to a function of base composition only, so
unrelated sequences with similar composition become indistinguishable (real phiX174 vs lambda phage scored
~0.99). Per-position tweaks to the 3-scalar profile didn't fix it. Windowing — keeping a vector of local summaries
whose count grows with the sequence — is what does.

## Scale: can this process real DNA sequences?

Researched (not assumed) because the current implementation stores sequences as `string`/`char` and indexes
with `int` throughout (`ZCurveWindowedProfile`'s string constructor, `WindowedProfile`, `NeedlemanWunsch`):

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
- **Net takeaway:** chromosome-scale or smaller (up to ~250M bases) fits comfortably in a `string`.
  Whole-large-genome scale (multi-Gbp) is handled by `ZCurveWindowedProfile.FromStream`, which never
  materializes the sequence (exercised at a synthetic 3.1 Gbp in `GenomeScaleStreamingTests`); the
  string-based constructors, `WindowedProfile`, and `NeedlemanWunsch` remain bounded by the string cap.

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
