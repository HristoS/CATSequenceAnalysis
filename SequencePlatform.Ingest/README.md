# SequencePlatform.Ingest

A small test app that fills the platform database (`docs/database-decision.md` rev 2) from
[BV-BRC](https://www.bv-brc.org/). For each genome it downloads the FASTA, profiles it, stores the profile,
then deletes the download. **The DNA sequence is never stored.** What goes into the database:

- `sequence_source`: provider, API URL, BV-BRC genome id, retrieval time, SHA-256 and size of the downloaded file, licence note.
- `sequence_record` (one per FASTA record: chromosome, plasmid or viral segment): kind, length, GC fraction derived from the
  profile, NCBI taxid, lineage (superkingdom → species; realm for viruses), and BV-BRC annotations (host, disease,
  isolation source, assembly/RefSeq accessions, CheckM, …) in `annotations` jsonb. Also the two D8 opt-ins, computed in
  the same streaming pass as the profile (`SequenceFeatureReader` in Core): `sequence_sha256` (SHA-256 of the record's
  letters, uppercased, nothing else) and `dinucleotide_counts` (16 adjacent A/C/G/T pair counts, AA…TT). Plus BV-BRC's
  own `sequence_md5` as `provider_md5`, which is used for de-duplication.
- `record_kmer` (migration 0013, D8 widened 2026-11-04): the record's 6-mer counts, 4,096 on one strand (AAAAAA … TTTTTT, other
  letters break a run), counted in the same pass. Every shorter k sums out of them. The eukaryote pilot doesn't store them yet.
- `profile_level`: `ZCurveWindowedProfile` ladder levels 500 bp / 2 / 8 / 32 kbp in the ZCWP binary format
  (`ZCurveWindowedProfile.ToBinary`). Levels with zero windows aren't stored.

## Usage

The connection string comes from the `SEQUENCE_PLATFORM_DB` environment variable or, when that's unset, from
`SequencePlatform.Ingest/appsettings.json`. That file is git-ignored (it holds a password) and copied to the build output:

```json
{ "ConnectionStrings": { "SequencePlatform": "Host=localhost;Port=5432;Username=...;Password=...;Database=sequence_platform" } }
```

```powershell
$env:SEQUENCE_PLATFORM_DB = "Host=localhost;Port=5432;Username=...;Password=...;Database=sequence_platform"   # optional, overrides the file
dotnet run --project SequencePlatform.Ingest -- migrate
dotnet run --project SequencePlatform.Ingest -- ingest-bvbrc --group viruses --limit 100
dotnet run --project SequencePlatform.Ingest -- ingest-bvbrc --group bacteria --limit 30 --offset 30
dotnet run -c Release --project SequencePlatform.Ingest -- ingest-bvbrc --group viruses --per-species --limit 100
dotnet run -c Release --project SequencePlatform.Ingest -- ingest-ncbi --taxon fungi --limit 100 --dry-run
dotnet run --project SequencePlatform.Ingest -- stats
dotnet run --project SequencePlatform.Ingest -- list --kind bacterium --rank phylum --taxon Pseudomonadota --limit 20
dotnet run --project SequencePlatform.Ingest -- list --min-length 1000000 --name borrel --after 120
dotnet run -c Release --project SequencePlatform.Ingest -- search-similar --fasta query.fna.gz --top 20
dotnet run --project SequencePlatform.Ingest -- profile --fasta query.fna.gz --out query.zcwp
dotnet run -c Release --project SequencePlatform.Ingest -- search-similar --profile query.zcwp --kind virus
dotnet run -c Release --project SequencePlatform.Ingest -- localize --fasta fragment.fna --kind bacterium
dotnet run -c Release --project SequencePlatform.Ingest -- localize --fasta read.fna --record 1234
```

- `migrate` creates the database if it's missing and applies `SequencePlatform.Core/Migrations/*.sql` in order, recorded in
  `schema_migration`. Every command runs this first.
- `list` pages the catalogue by id: each page ends with the `--after` value for the next one. Filters: `--kind`,
  `--rank` with `--taxon` (superkingdom, realm, kingdom, phylum, class, order, family, genus, species), `--min-length`,
  `--max-length`, `--provider`, `--name` (case-insensitive substring). `--limit` is 1–1000, default 50.
- Database access, migrations, profile-level storage, the catalogue query and the provider-independent profile-and-store
  step (`RecordStore`) live in `SequencePlatform.Core`; this project is the BV-BRC crawler and CLI on top of it. Platform tests: `dotnet test SequencePlatform.Tests`.
- `--group bacteria` and `--group archaea` select complete genomes marked Reference or Representative (6,036 and 374 on
  2026-10-19). `--group viruses` selects complete Reference genomes (15,536).
- `--per-species` widens the catalogue instead of filling it: from all of a group's complete genomes (71,364 bacteria, 1,387 archaea, 10.3 M
  viruses on 2026-10-30), one genome for each species not yet stored, Reference or Representative first, in a fixed random order, up to
  `--limit`. Bacteria and archaea are listed once and picked from. Viruses can't be listed, so each missing species is one request as the
  run reaches it, and an interrupted run loses at most one batch of 50. BV-BRC refuses species names with an apostrophe; those are logged and
  skipped. Item groups are `<group>-species`, so failed items are retried with `--retry-failed --group viruses-species`.
- `--contact you@example.org` adds a contact address to the User-Agent. Recommended for large runs.
- `ingest-ncbi --taxon fungi|plants|animals --limit N [--max-bp B] [--parallel P] [--retry-failed] [--dry-run]` is the eukaryote pilot: reference assemblies at chromosome
  level or better, one per family in a fixed random order, up to `B` bases each (default 2 Gbp). `--parallel P` streams `P` genomes at once,
  default 1. A stream dropped midway marks that genome failed, and only `--retry-failed` fetches it again. `--dry-run` prints what would be
  fetched, by class. The FASTA streams from NCBI's FTP tree through hash, gzip and `SequenceFeatureReader` with no temp file, so the DNA never
  reaches the disk. A genome becomes one record, `kind = 'eukaryote'`, `data_class = 'public_eukaryote'`, holding its nuclear dinucleotide counts
  summed and no profile levels: it has relatives and is on the map and in trees, but find-similar and localize can't score it, and in a
  length band it appears in the Not comparable list. Since migration 0022 the same pass also stores each nuclear sequence of at least
  100 kbp as chromosome levels (see `backfill-chromosomes`). Organelles are left out by their FASTA description. Human assemblies are refused in
  code. `--limit` and `--max-bp` decide what is picked, so resume a run with the same values. NCBI's FTP gave about 0.5 MB/s
  (1.6 Mbp/s) per connection on 2026-10-30, so size a run by bases rather than genomes. Since S7.T, each genome stored also gets its
  gene layer, as `backfill-genes-ncbi` builds it.
- `backfill-genes-ncbi [--limit N] [--accession A] [--dry-run] [--check-search N]` (P7 S7.T, migration 0020) builds the gene layer for
  NCBI eukaryotes without downloading a genome. For each record with no 4-mer row it looks the assembly's annotation and nuclear genetic
  code up in NCBI Datasets, then streams the assembly report (which sequences are nuclear) and `*_cds_from_genomic.fna.gz` (each coding
  sequence already spliced, in coding orientation, with its location). Per gene (GeneID) it keeps the longest isoform and leaves out
  partial, pseudo and `exception`-flagged coding sequences, since the windows are anchored at the start codon. It stores the
  `gene_profile` rows (ZCGP version 3 where genes lie on several sequences or are spliced), amino-acid 4-mers, and under
  `annotations.ncbi.gene_layer` the CDS file's URL and SHA-256, the annotation's name and release, the counts left out and the names of
  the sequences the genes lie on, which protein hits refer to. An assembly with no annotation gets `no_cds` and is skipped by later runs.
  Genetic codes 1 and 12 (CUG-Ser yeasts) are read; any other is skipped and stays pending. Nothing of the DNA or the proteins is kept.
  `--dry-run` builds and reports without writing (the largest bucket row against gene rows' 4 MB limit, `GeneProfiles.MaxRowBytes`; genome levels keep 1 MB); `--check-search N` then searches
  N proteins of each assembly through `ProteinSearch` over the whole catalogue. On 2026-11-06, 131 of the 587 stored eukaryotes
  were annotated: 4–38 s per assembly; a 3 Gbp mammal kept 20,457 genes of 47,953 coding sequences (26,268 other isoforms).
- `backfill-chromosomes [--limit N] [--accession A] [--parallel P]` (migration 0022, 2026-11-11) stores each eukaryote's chromosome levels:
  every nuclear sequence of at least 100 kbp at 500 bp and 32 kbp windows (`chromosome_level`, ZCWP version 2), for drawing one
  chromosome's path. Since 0022 `ingest-ncbi` builds them in the same pass that counts dinucleotides; this re-downloads genomes stored
  before, smallest first. It keeps them only when the file gives the record's nuclear length and dinucleotide counts again, and records
  the outcome under `annotations.ncbi.chromosome_levels` (`stored` with sequences kept and the shorter ones left out, or `mismatch`), so a
  finished genome isn't picked again and a stopped run resumes. A failed download leaves the genome pending. The 587 genomes are 256 Gbp
  (about 77 GB of gzip from NCBI) and 3.1 GB of levels: the first four (1.7 Gbp, 2026-11-11) took 21 MB, and a 0.8 Gbp plant 20 min
  over one connection. Nothing reads these rows but the chromosome endpoints: a eukaryote is still not
  searched as one genome.
- `search-similar --fasta X | --profile X.zcwp` (P3 S3.4) ranks the records of comparable length, [r × L, L / r] for a query
  of length L, by `ZCurveWindowedProfile` similarity at one ladder level, the better of the query's two strands. `--ratio` is
  0.25–1, default 0.5; `list`'s filters narrow the candidates further; `--top` defaults to 20. Records with the same
  `sequence_sha256` on either strand are listed first as identical sequences. Each FASTA record is its own query. Very different
  lengths aren't compared, because there the score measures whether the shorter sequence is contained in the longer
  (`results/S34_LengthGate.csv`). The query file is read and kept in memory only; nothing of it is stored.
- Each hit shows `log10 E`: how many unrelated records of the same shape would be expected to score as high among this
  search's candidates, from the `score_null` table. `<=` marks an upper bound (the null sample has fewer than 10 scores
  that high); without a built null the column is empty. The ranking is by score, not by E. A query whose shortest
  candidates keep fewer than 8 windows (under about 8 kbp at ratio 0.5) gets a warning: its ranking mostly reflects GC.
- `localize --fasta X | --profile X.zcwp [--record ID] [list's filters] [--top N]` (P3 S3.5) says where each fragment comes from:
  record, offset on its stored strand (to about 500 bp), strand, score and `log10 E` from the fragment null. Candidates are all
  records at least as long as the fragment, scanned at the coarsest level that leaves the fragment 4 windows, then the best 100
  refined at 500 bp. Over the catalogue the fragment must be at least 4 kbp; under 8 kbp it gets a weak-signal warning. A shorter
  fragment can be placed inside one record with `--record ID`. When the top hit isn't the source it was mostly a close relative
  (another strain; `results/S35_LocalizeCheckMin4_Summary.csv`). Not found: a fragment across the origin of a circular record.
- `build-null [--per-cell N] [--genome-draws N] [--fragment-draws N] [--pool N] [--seed S] [--csv PATH]` (P3 S3.6,
  migration 0007) samples best-alignment scores between unrelated records (different superkingdom or phylum; viruses of
  different families) from the stored profiles, as find-similar and localization would score them (a virus candidate wrapped around
  its end, as find-similar scores it since 2026-11-12: rerun it and `build-neighbours` after deploying that), and replaces the
  `score_null` rows. Defaults take about 2 minutes on the 2026-10 catalogue; rerun it after large ingests. `--csv` also
  writes every cell's scores to a file.
- `build-neighbours [--ratio R]` (P3 S3.4b, migration 0008) runs find-similar for every stored record against the other
  assemblies and keeps each record's top 20 in `neighbour_list` (about 4 minutes on the 2026-10 catalogue). Once that graph
  exists for a search's ratio, `search-similar` re-orders its top 100 by `shared`: how many of the query's top 20 are also
  in the hit's own top 20. That raised same-genus precision@10 from 0.131 to 0.195 for bacteria. A composition
  look-alike can rise too (a plasmid with coronavirus-like GC and skew comes fourth for SARS-CoV-2); its E-value shows it
  isn't significant. Rebuild the graph after large ingests; hits missing from it are searched live.
- `build-origins` (S3.2's rotation) lays every circular bacterial chromosome of at least 500 kbp with a clear GC skew out from
  its replication origin, found from its stored 500 bp profile (the cumulative skew's minimum; no download), and stores its
  2, 8 and 32 kbp levels under method version `1-origin` (4,951 of 5,944 chromosomes, 34 s on the 2026-10 catalogue).
  `search-similar` then scores those records, and a query with a clear skew, both as stored and rotated, and keeps the
  better score; `strand*` marks a hit whose best score came from a rotated form. Rerun `build-null` and `build-neighbours`
  afterwards. Localization ignores the rotated levels, so its offsets stay in the deposit's coordinates.
- `convert-levels-half` (2026-11-05) rewrites every record whose stored levels are still ZCWP version 1 (float32) as version 2
  (float16): each record's 2, 8 and 32 kbp levels are rebuilt from its stored 500 bp level and rounded as ingest now stores them,
  in batches that commit together, so a stopped run resumes. It then reruns `build-origins`. No download. Afterwards rebuild the
  null and the neighbour graph and `VACUUM FULL profile_level`. On the local catalogue: 60,847 records in 408 s, 11,304 origins
  in 35 s more. The logical size halves, but on disk it saves almost nothing, because PostgreSQL's TOAST compression had already
  halved the float32 rows (`results/P7_StorageReduction.csv`).
- `backfill-kmers [--limit N] [--sample SEED]` gives BV-BRC records stored before migration 0013 their 6-mer counts: it downloads
  each genome again, reads it once through `SequenceFeatureReader`, deletes the file, and stores the counts only for a record whose
  SHA-256 and dinucleotide counts match what was stored (`done`); otherwise it records `mismatch` or `missing`, so a rerun skips it.
  Nothing else is recomputed. `--sample` takes genomes in a fixed random order (a trial); a genome that fails to download is left
  for the next run.
- `backfill-genes [--limit N] [--sample SEED] [--check-search N] [--parallel P]` (P7, migrations 0014–0018) builds the gene layer
  for BV-BRC records: it downloads each genome and its coding sequences (BV-BRC's own calls and PGfam numbers), checks the SHA-256 as
  `backfill-kmers` does, and stores per record its amino-acid 4-mer counts (`protein_kmer`) and its `ZCurveGeneProfile` as
  `gene_profile` rows by length bucket (genes of 500 aa and up in bucket 0, 50–499 aa in buckets 1–25), plus 6-mer counts where
  they are missing; then it deletes the file. A
  record is pending while it has no 4-mer row, or no `gene_profile` row. `--parallel` downloads that many genomes at once.
  `--check-search N` then searches N proteins per record and length band from the protein alone, as `ProteinSearch` compares them,
  and prints top-1 measured within the pass and extrapolated to the catalogue. On the local catalogue, the S7.S pass took 299.7 min
  at 4 in parallel (0.41 s per genome): 58,141 records, 64.3M genes, 6.14 GB of gene profiles (`results/S7S_FullPass.csv`).
- `build-protein-null [--per-bucket N] [--seed S]` (P7, migration 0019) builds p and E for protein hits: per query length bucket (the
  26 gene_profile buckets), N real BV-BRC proteins with a protein family (default 30, at most 2 per bucket from each stored genome, in a
  seeded order), each scored as a protein search scores it against the genes of other families; the cumulative counts go to
  `protein_score_null`. The proteins are fetched for the build and never stored. Rebuild after a large ingest. The check behind it:
  `results/ProteinNull_Check.csv`.
- `verify-genes` checks every `gene_profile` row against its stored SHA-256 and lists the rows that don't match (exit 2 if any). Protein
  search stopped hashing the rows it reads on 2026-11-09, since that cost as much as decoding them (`results/SearchBreakdown_Before.csv`);
  rows are still hashed when written, and `GeneProfiles.ReadAsync` (a record's gene view) still checks what it reads. Run it after a
  restore or a disk problem, or before a measurement that must trust the layer.
- `cluster-genes` rewrites `gene_profile` (and its TOAST table) in bucket order (`CLUSTER ... USING gene_profile_bucket`, then `VACUUM (ANALYZE)`).
  Rows are written record by record, so a bucket's rows end up spread over the whole table, and a protein search reading 2-3 buckets
  touches most of it; locally the read took 17-18% less afterwards (2026-11-10, `results/GeneProfileOrder_Check.csv`). It holds an
  exclusive lock for the rewrite (293 s locally for 6.6 GB, plus 12 s of vacuum: no protein search or gene write meanwhile) and needs
  free disk about the table's size. New records are appended unsorted, so run it again after a large ingest or gene backfill. Order changes no result.
- `profile --fasta X --out X.zcwp` writes the query-side profile without touching the database: `X.zcwp` and the reverse strand
  `X.rc.zcwp`, both from the DNA. `search-similar --profile X.zcwp` picks up `X.rc.zcwp` when it's there and otherwise searches the
  forward strand only. A file with several records writes `X.<record id>.zcwp` per record.
- `load-taxonomy [--refresh] [--skip-gtdb]` (S2.4, migration 0006) fills `taxon` with NCBI's current lineage for every taxid the
  catalogue references (NCBI Datasets API, 200 per request; a merged taxid resolves to the taxon it was merged into) and
  `gtdb_taxon` with GTDB's lineage for every referenced bacterial or archaeal assembly (latest release, matched on the
  accession without its GCA_/GCF_ prefix). A plain run fetches only taxids not yet in `taxon`, so run it after each ingest;
  `--refresh` re-fetches all. `sequence_record` keeps BV-BRC's lineage as ingested; the view `lineage_difference` lists, per
  record and rank, where it and NCBI's differ (`differs`, `ncbi_only`, `provider_only`).
- `add-user NAME` (P4, D3) creates an API user and prints their key once (`sp_...`); only its SHA-256 is stored. `disable-user NAME`
  revokes it; the user's uploads stay. See `SequencePlatform.Api/README.md`.

## Behaviour

- **Polite:** requests are sequential, except the genomes of `ingest-ncbi --parallel`, and HTTP 429/5xx are retried with exponential backoff
  that honours `Retry-After`. NCBI's requests start at least 1 s apart; BV-BRC's follow each other with no gap since 2026-11-04 (one
  in flight at a time; the gap had become half of each species' time). BV-BRC publishes no rate limit.
- **Resumable:** every listed genome gets a `crawl_item` row (`pending`, then `done`, `duplicate`, `too_large` or
  `failed`, with attempt count and last error). Items in a final state are skipped without any request, so re-running a
  batch downloads nothing it already handled. One genome is one transaction, and its `done` mark commits in it, so a
  failure leaves nothing partial behind. Failed items are skipped by normal runs and retried with
  `ingest-bvbrc --retry-failed [--group G] [--limit N]`; `stats` shows the counts per group and status.
- `too_large` means a record's 500 bp level is over 1 MB and wasn't stored (it needs object storage); the coarser levels
  that fit are kept. It doesn't occur for viruses or bacteria.
- **De-duplicating:** a genome is skipped as `duplicate` when every one of its sequences is already stored. The first copy
  crawled wins, and `crawl_item.duplicate_of` points at it. BV-BRC lists one assembly under two genome ids (same contig
  accessions) and a GenBank record next to its RefSeq copy (HM119410 / NC_055339, different accessions). Both are
  recognised before downloading: either by accession, or by BV-BRC's per-contig `sequence_md5`, stored as
  `sequence_record.provider_md5`. A contig without an MD5 is downloaded and profiled, then skipped if its `sequence_sha256`
  is stored. A genome that is only partly stored is kept whole, and the run says how many of its records were already
  stored. The assembly accession isn't used as a key, because the segments of one virus share it.
- FASTA comes from the BV-BRC API (`Accept: application/dna+fasta`); its FTP host isn't reachable from every network.

## Known limitations

- **No object storage yet.** A level over 1 MB (chromosomes above about 43 Mbp at 500 bp) is skipped with a message. That
  doesn't happen for viruses or bacteria.
- **Segmented viruses:** BV-BRC lists each segment as its own genome, so each becomes its own source and record. Grouping
  segments per virus is a search-time concern.
