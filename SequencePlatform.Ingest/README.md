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
- `--contact you@example.org` adds a contact address to the User-Agent. Recommended for large runs.
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
  different families) from the stored profiles, as find-similar and localization would score them, and replaces the
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
- `profile --fasta X --out X.zcwp` writes the query-side profile without touching the database: `X.zcwp` and the reverse strand
  `X.rc.zcwp`, both from the DNA. `search-similar --profile X.zcwp` picks up `X.rc.zcwp` when it's there and otherwise searches the
  forward strand only. A file with several records writes `X.<record id>.zcwp` per record.
- `load-taxonomy [--refresh] [--skip-gtdb]` (S2.4, migration 0006) fills `taxon` with NCBI's current lineage for every taxid the
  catalogue references (NCBI Datasets API, 200 per request; a merged taxid resolves to the taxon it was merged into) and
  `gtdb_taxon` with GTDB's lineage for every referenced bacterial or archaeal assembly (latest release, matched on the
  accession without its GCA_/GCF_ prefix). A plain run fetches only taxids not yet in `taxon`, so run it after each ingest;
  `--refresh` re-fetches all. `sequence_record` keeps BV-BRC's lineage as ingested; the view `lineage_difference` lists, per
  record and rank, where it and NCBI's differ (`differs`, `ncbi_only`, `provider_only`).

## Behaviour

- **Polite:** requests are sequential and at least 1 s apart, and HTTP 429/5xx are retried with exponential backoff
  that honours `Retry-After`. BV-BRC publishes no rate limit.
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
