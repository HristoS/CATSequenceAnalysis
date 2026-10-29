# SequencePlatform.Api

The platform's HTTP API (P4): the catalogue, uploads, and find-similar and localize as asynchronous jobs, run by a worker in
the same process.
It stores no DNA: uploads and queries are profiled and discarded, as in the Ingest. `NoDnaTests` checks this on Kestrel: a genome
carrying a unique 1 kbp marker goes through upload, find-similar, localize and rejected bodies, and afterwards no column of any
table (bytea included), no log entry at any level and no temp file holds any 24 bp piece of the marker or of its reverse
complement; a body over 32 MiB gets 413 and stores nothing.

## Run it

1. Migrate the database and create a user. The key is printed once; only its hash is stored.

   ```
   dotnet run --project SequencePlatform.Ingest -- migrate
   dotnet run --project SequencePlatform.Ingest -- add-user alice
   ```

2. Give the API the connection string, in any of these ways (the first one found wins):
   - `dotnet user-secrets set ConnectionStrings:SequencePlatform "Host=localhost;..." --project SequencePlatform.Api`
   - the environment variable `ConnectionStrings__SequencePlatform`
   - the environment variable `SEQUENCE_PLATFORM_DB`, as the Ingest reads it

3. `dotnet run --project SequencePlatform.Api` listens on http://localhost:5080. The OpenAPI document is at `/openapi/v1.json`
   (Development only, no key needed).

Every other request needs `Authorization: Bearer <key>`. **Localhost only** until S6.4: the key is a plain bearer token over HTTP.
Revoke a key with `dotnet run --project SequencePlatform.Ingest -- disable-user alice`.

## Visibility (D3)

Downloaded records are public. An upload is visible only to its owner: someone else's upload is a 404 everywhere, the same as a
record that doesn't exist, and it never becomes a search candidate for anyone else. The neighbour graph and the score null are
built from public records only, so private uploads don't shape anyone else's ranking or p-values either.

## Endpoints

| Method | Route | |
|---|---|---|
| GET | `/api/sequences` | Keyset pages (`after`, `limit` ≤ 1000); filters `kind`, `rank` + `taxon`, `minLength`, `maxLength`, `provider`, `name` |
| GET | `/api/sequences/{id}` | Metadata, provider / NCBI / GTDB lineage, annotations, stored levels, provenance |
| GET | `/api/sequences/{id}/profile?windowSize=` | One level as ZCWP bytes (`application/x-zcwp`), 500 / 2000 / 8000 / 32000 |
| POST | `/api/sequences` | Upload `text/x-fasta` (any number of records) or `application/x-zcwp` (the 500 bp level, with `lengthBp`); `kind` required, `name`, `taxid`, `dataClass` (default `restricted`), `recordId` for a profile. 201 with the new records |
| POST | `/api/searches/similar` | Queue find-similar for one FASTA record or a ZCWP query (forward strand, optionally the reverse right after it); catalogue filters, `ratio` (0.25–1, default 0.5), `top` (1–1000, default 50), `lengthBp` for a profile. 202 + `Location` |
| POST | `/api/searches/localize` | The same for a fragment: ≥ 4 kbp over the catalogue, or any length with `record`; `top` 1–100 (default 20) |
| GET | `/api/searches/{id}` | Status (Queued, Running, Done, Failed, Expired), parameters, summary, and once done one page of one hit list: `list` = Ranked / Identical / NotComparable, `after`, `limit` |
| DELETE | `/api/searches/{id}` | Drop the query profile and results now (D7); 204. The search stays as Expired with its summary |

Requests for each are in [Sequences.http](Sequences.http) and [Searches.http](Searches.http), and a whole walk-through in
[Scenario.http](Scenario.http).

## End-to-end scenario

[Scenario.http](Scenario.http) goes from an empty database to answers, using only the bundled genomes in `data/genomes/`:
upload, list, find similar, poll, read both lists, localize a planted gene, delete. From the repository root:

1. Make the bodies (plain FASTA into `SequencePlatform.Api/scenario/`, git-ignored):
   `python -I scripts/platform-checks/e2e_scenario.py prepare`. They are the 29 viral genomes as one file, SARS-CoV-2 alone,
   its S gene (3,822 bp), and E. coli K-12 MG1655 with that S gene inserted at 2,000,000.
2. Point `SEQUENCE_PLATFORM_DB` at a database that doesn't exist yet, run the Ingest's `migrate` (it creates it) and `add-user`,
   and start the API on it ("Run it" above).
3. Send the requests in order, polling steps 6 and 9 until `Done`. Or let the script send them and check every answer:
   `python -I scripts/platform-checks/e2e_scenario.py run http://localhost:5080 KEY_FILE`, with the key in `KEY_FILE`.

What comes back (2026-10-29, Release, `results/S45_EndToEnd.csv`):

| Step | Answer |
|---|---|
| Upload the viruses, then the planted K-12 | 201: 29 records; 1 record of 4,645,474 bp, 9,290 windows, levels 500 bp to 32 kbp |
| List your uploads | the 30 records, owned and `restricted` |
| Find similar to SARS-CoV-2 | Done 0.28 s after POST. Band 14,952–59,806 bp, 8 ranked: SARS-CoV-2 1.000, MERS 0.879, HCoV-229E 0.877, HCoV-OC43 0.851, then human adenovirus A 0.769 |
| The identical list | the uploaded SARS-CoV-2, by SHA-256 |
| Localize the S gene in the planted record | Done 0.11 s after POST: offset 2,000,000, forward strand, score 1.000 |
| Delete the search | 204; it reads back as Expired with its summary and no hits |

A fresh database has no score null or neighbour graph (the Ingest's `build-null` and `build-neighbours` make them from public
records), so the hits carry no p or E and are ordered by score alone.

## Uploads

The body is the file itself and the metadata goes in the query string. A multipart form would go through `IFormFile`, which
buffers anything over 64 KB to a temp file. FASTA streams in one pass: each record is profiled as it arrives and never exists as a string, and the
body is hashed on the way, so only the profile, the hashes and the metadata are stored, exactly as for a download. A request
body is capped at 32 MiB (413 beyond), which holds the largest known prokaryote as FASTA. E. coli K-12 (4.7 MB) took 0.75–1.4 s
on the local machine (two single uploads, not a benchmark).

## Searches

A search is a job: POST checks everything it can at once (filters, ratio, top, the fragment's length, that a named record is
visible to you), profiles the query, and answers 202 with the job's address. The query's DNA lives only for that request: the
job keeps the fine profile and SHA-256 of both strands. A ZCWP query has no hash, so it finds no Identical records, and without a
reverse strand it is searched on the forward strand only. Searches are private: anyone else's id is a 404. Results are kept 30
days (D7), then the search is Expired and keeps only its parameters, the query's length and hash, and the summary.

## The worker

`SearchWorker` runs queued searches one at a time in the API process: each search already uses every core. It wakes on a
PostgreSQL notification (`search_queued`) and also polls every 30 s, so a lost listener only delays jobs. It assumes it is the
only worker (one API process until P6): at startup a job still `Running` lost its worker and is requeued, and after three claims
it is failed instead. Every 10 minutes it applies D7's expiry and reloads the score null and neighbour graphs, which the Ingest
rebuilds; the cached profiles are refreshed before every search, so a record stored a moment ago is a candidate.

Find-similar keeps its band levels in memory (24 MB at ratio 0.5). Localization streams from the database by default
(`SearchWorker:CacheLocalization=false`): about 1-2.3 s per fragment instead of 0.2-1 s, without 823 MB in the process. Other
settings in the `SearchWorker` section: `Enabled`, `PollInterval`, `StuckAfter`, `MaintenanceInterval`.

On the local catalogue (`results/S43_WorkerTiming.csv`, Release): warm find-similar takes 150-650 ms from POST to Done (the
upper end is profiling a 4.6 Mbp FASTA query in the POST), localization 0.9-2.3 s; the first search after a start takes about
3 s while the caches, null and graph load. The process used 451 MB.
