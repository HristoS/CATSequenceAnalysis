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

The connection string comes only from an environment variable, so keep it out of the repo:

```powershell
$env:SEQUENCE_PLATFORM_DB = "Host=localhost;Port=5432;Username=...;Password=...;Database=sequence_platform"
dotnet run --project SequencePlatform.Ingest -- migrate
dotnet run --project SequencePlatform.Ingest -- ingest-bvbrc --group viruses --limit 100
dotnet run --project SequencePlatform.Ingest -- ingest-bvbrc --group bacteria --limit 30 --offset 30
dotnet run --project SequencePlatform.Ingest -- stats
```

- `migrate` creates the database if it's missing and applies `Migrations/*.sql` in order, recorded in `schema_migration`.
  Every command runs this first.
- `--group bacteria` selects complete genomes marked Reference or Representative. `--group viruses` selects complete Reference genomes.
- `--contact you@example.org` adds a contact address to the User-Agent. Recommended for large runs.

## Behaviour

- **Polite:** requests are sequential and at least 1 s apart, and HTTP 429/5xx are retried with exponential backoff
  that honours `Retry-After`. BV-BRC publishes no rate limit.
- **Resumable:** genomes already stored are skipped without any request, so re-running a batch continues where it
  stopped. One genome is one transaction, so a failure leaves nothing partial behind.
- **De-duplicating:** BV-BRC sometimes lists one assembly under two genome ids. A genome whose contig accessions are
  all stored already is skipped before downloading.
- FASTA comes from the BV-BRC API (`Accept: application/dna+fasta`); its FTP host isn't reachable from every network.

## Known limitations

- **No object storage yet.** A level over 1 MB (chromosomes above about 43 Mbp at 500 bp) is skipped with a message. That
  doesn't happen for viruses or bacteria.
- **Segmented viruses:** BV-BRC lists each segment as its own genome, so each becomes its own source and record. Grouping
  segments per virus is a search-time concern.
