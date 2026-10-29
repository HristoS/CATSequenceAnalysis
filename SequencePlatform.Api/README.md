# SequencePlatform.Api

The platform's HTTP API (P4): the catalogue, uploads, and find-similar and localize as asynchronous jobs. Queued searches run once
the worker exists (S4.3).
It stores no DNA: uploads and queries are profiled and discarded, as in the Ingest.

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

Requests for each are in [Sequences.http](Sequences.http) and [Searches.http](Searches.http).

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
