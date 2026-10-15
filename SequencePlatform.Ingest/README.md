# SequencePlatform.Ingest

A small test app that fills the platform database (`docs/database-decision.md` rev 2) from
[BV-BRC](https://www.bv-brc.org/). For each genome it downloads the FASTA, profiles it, stores the profile,
then deletes the download. **The DNA sequence is never stored.** What goes into the database:

- `sequence_source`: provider, API URL, BV-BRC genome id, retrieval time, SHA-256 and size of the downloaded file, licence note.
- `sequence_record` (one per FASTA record: chromosome, plasmid or viral segment): kind, length, GC fraction derived from the
  profile, NCBI taxid, lineage (superkingdom → species; realm for viruses), and BV-BRC annotations (host, disease,
  isolation source, assembly/RefSeq accessions, CheckM, …) in `annotations` jsonb.
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
```

- `migrate` creates the database if it's missing and applies `SequencePlatform.Core/Migrations/*.sql` in order, recorded in
  `schema_migration`. Every command runs this first.
- `list` pages the catalogue by id: each page ends with the `--after` value for the next one. Filters: `--kind`,
  `--rank` with `--taxon` (superkingdom, realm, kingdom, phylum, class, order, family, genus, species), `--min-length`,
  `--max-length`, `--provider`, `--name` (case-insensitive substring). `--limit` is 1–1000, default 50.
- Database access, migrations, profile-level storage and the catalogue query live in `SequencePlatform.Core`; this
  project is the BV-BRC crawler and CLI on top of it. Platform tests: `dotnet test SequencePlatform.Tests`.
- `--group bacteria` selects complete genomes marked Reference or Representative. `--group viruses` selects complete Reference genomes.
- `--contact you@example.org` adds a contact address to the User-Agent. Recommended for large runs.

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
- **De-duplicating:** BV-BRC sometimes lists one assembly under two genome ids. A genome whose contig accessions are
  all stored already is skipped before downloading.
- FASTA comes from the BV-BRC API (`Accept: application/dna+fasta`); its FTP host isn't reachable from every network.

## Known limitations

- **No object storage yet.** A level over 1 MB (chromosomes above about 43 Mbp at 500 bp) is skipped with a message. That
  doesn't happen for viruses or bacteria.
- **Segmented viruses:** BV-BRC lists each segment as its own genome, so each becomes its own source and record. Grouping
  segments per virus is a search-time concern.
