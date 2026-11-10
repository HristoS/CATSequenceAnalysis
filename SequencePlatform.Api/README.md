# SequencePlatform.Api

The platform's HTTP API (P4) and web UI: the catalogue, uploads, find-similar, localize and protein search (P7) as asynchronous
jobs, and the map and trees of organisms. A worker runs the jobs, in the same process by default and in its own under Aspire.
It stores no DNA: uploads and queries are profiled and discarded, as in the Ingest. `NoDnaTests` checks this on Kestrel: a genome
carrying a unique 1 kbp marker goes through upload, find-similar, localize and rejected bodies, and afterwards no column of any
table (bytea included), no log entry at any level and no temp file holds any 24 bp piece of the marker or of its reverse
complement; a body over 32 MiB gets 413 and stores nothing. A marked 1,000 aa protein goes through protein search the same way,
and no 24-letter piece of it is stored or logged either.

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

Every other request needs `Authorization: Bearer <key>`. Run locally, the key is a plain bearer token over HTTP. Hosted, Caddy holds
the TLS in front of Kestrel, which listens on localhost only, and the API connects to PostgreSQL as the least-privilege role `sequence_app`
(migration 0021; `hosting/tencent/README.md`). The API tests connect as that role too, so a write it isn't granted fails a test.
Revoke a key with `dotnet run --project SequencePlatform.Ingest -- disable-user alice`.

## Run it with Aspire

`dotnet run --project SequencePlatform.AppHost` runs the database migration, then the API and the worker as two processes of this
project (`SearchWorker__Enabled=false` on the API), with logs, traces and metrics on the dashboard whose login URL it prints. It uses the
database in the git-ignored `SequencePlatform.Ingest/appsettings.json`, and the API stays on http://localhost:5080. Aspire gives no
cluster: it starts the services here and generates their deployment files. Keep one worker, because at startup it requeues every
running job as its own.

The same services as Docker Compose, on their own empty `postgres:17`, with the API on 127.0.0.1:5080 only:

```
dotnet run --project SequencePlatform.AppHost --no-launch-profile -- --operation publish --step deploy --output-path deploy/compose
```

`--step publish` writes the files without building or starting anything; leaving `--step` out runs the whole pipeline, teardown
included. The output is git-ignored because its `.env` holds the database password. Stop it with
`docker compose -p <project> down -v`, where `docker ps` shows the project name. Use `http://127.0.0.1:5080` with this stack: on Windows
`localhost` tries IPv6 first and each request loses 2 s. Create a user by running the `migrate` image with `add-user NAME` and
`SEQUENCE_PLATFORM_DB` on the stack's network.

## The web UI

http://localhost:5080 serves the UI from `wwwroot/`: sign in with an API key, kept in the tab's `sessionStorage`. It is React 18.3.1
and htm with no build step and no Node: `vendor/` holds their pinned files, each page is an ES module in `app/`, and components use
`htm` templates rather than JSX. Pin React to 18: version 19 has no UMD builds. The page's Content-Security-Policy is
`default-src 'self'`, so it loads nothing from elsewhere. Static files are served before authentication because they hold no data;
every `/api` call still needs the key. A FASTA file goes from the browser straight to the API, so the no-DNA guarantee below covers it.

A record's page shows its lineages, nearest relatives, dinucleotide signature (Karlin, Z-curve coefficients or the gene layer's mean
window), Z-curve path and, for BV-BRC genomes, its genes: pick one to see its Z-curve over codon positions 1 and 2, redrawn from the
stored windows; the path panel can also build the genome's path from its genes alone (Source: Genes), and highlights the picked gene. `#/compare?ids=a,b,c` (from the catalogue's checkboxes, a record's nearest relatives, or the top hits of a search) puts 2 or 3 records on shared plots, each in its own colour that the legend switches on and off
everywhere: their Z-curve paths, from the genome windows or from the genes, with the picked protein's genes highlighted (circular chromosomes from the replication origin, `results/CompareOverlay_Check.csv`; scaled to the
same genome length or in bases), their signatures as grouped bars, and their copies of a protein family they all carry (or any gene
picked per record) as overlaid gene curves and heat strips.

## Visibility (D3)

Downloaded records are public. An upload is visible only to its owner: someone else's upload is a 404 everywhere, the same as a
record that doesn't exist, and it never becomes a search candidate for anyone else. The neighbour graph and the score null are
built from public records only, so private uploads don't shape anyone else's ranking or p-values either.

## Endpoints

| Method | Route | |
|---|---|---|
| GET | `/api/sequences` | Keyset pages (`after`, `limit` ≤ 1000); filters `kind` (a stored kind, or `prokaryote` for bacterium and archaeon), `rank` + `taxon`, `minLength`, `maxLength`, `provider`, `name`. Searches take the same filters |
| GET | `/api/sequences/facets` | Filter choices with record counts over what the caller can see: `kinds`, `kindGroups` (with their member kinds) and `providers`, largest first, plus `total` |
| GET | `/api/sequences/{id}` | Metadata, provider / NCBI / GTDB lineage, annotations, stored levels, provenance, and `originBp`: where its `1-origin` levels start (its replication origin), when it has them |
| GET | `/api/sequences/{id}/profile?windowSize=&origin=` | One level as ZCWP bytes (`application/x-zcwp`), 500 / 2000 / 8000 / 32000; `origin=true` gives a circular chromosome's level laid out from its replication origin (2000 bp and up, where stored) |
| POST | `/api/sequences` | Upload `text/x-fasta` (any number of records) or `application/x-zcwp` (the 500 bp level, with `lengthBp`); `kind` required, `name`, `taxid`, `dataClass` (default `restricted`), `recordId` for a profile. 201 with the new records |
| POST | `/api/searches/similar` | Queue find-similar for one FASTA record or a ZCWP query (forward strand, optionally the reverse right after it); catalogue filters, `ratio` (0.25–1, default 0.5), `top` (1–1000, default 50), `lengthBp` for a profile. 202 + `Location` |
| POST | `/api/searches/localize` | The same for a fragment: ≥ 4 kbp over the catalogue, or any length with `record`; `top` 1–100 (default 20) |
| POST | `/api/searches/protein` | Queue a protein search: one protein as a FASTA record or bare letters (`text/x-fasta` or `text/plain`), ≥ 333 aa; catalogue filters, `code` (11, default, or 4), `top` (1–1000, default 50). 202 + `Location` |
| GET | `/api/searches/{id}` | Status (Queued, Running, Done, Failed, Expired), parameters, summary, and once done one page of one hit list: `list` = Ranked / Identical / NotComparable, `after`, `limit` |
| GET | `/api/searches/{id}/query` | A protein search's query as the windows it was compared with (the expected Z-curve, never the letters), in the gene windows' format, so the page can draw it over a hit gene; 404 for someone else's search, a DNA search, or once expired |
| DELETE | `/api/searches/{id}` | Drop the query profile and results now (D7); 204. The search stays as Expired with its summary |
| GET | `/api/sequences/{id}/relatives?k=` | The `k` (1–100, default 10) genomes nearest by dinucleotide signature, in the record's group (bacteria and archaea, or viruses), nearest first. 400 for a plasmid or a profile-only upload |
| GET | `/api/sequences/{id}/signature?view=` | The record's dinucleotide signature: `zcurve` (default, the 16 Z-curve coefficients) or `karlin` (log ρ\* and ρ\* per dinucleotide) |
| GET | `/api/sequences/{id}/genes` | The record's stored genes (BV-BRC, 50 aa and up), in order of start, with strand, length in aa, family and product, plus `meanWindow`: the 9 window features averaged over all of them, in the signature's format. No genes for uploads and NCBI eukaryotes |
| GET | `/api/sequences/{id}/genes/path?origin=` | The record's Z-curve path built from its genes alone, in chromosome order (from the replication origin with `origin=true`, when it has one): each gene window one step of its means over codon positions 1 and 2, reverse-strand genes flipped onto the forward strand, intergenic DNA left out; flat `coordinates`, and each gene's `firstPoint` (`results/GenePath_Check.csv`) |
| GET | `/api/sequences/{id}/genes/{index}` | One gene's stored windows: mean, trend and curvature of x, y and z over codon positions 1 and 2, from which the page redraws the gene's Z-curve (`results/GeneCurve_Check.csv`) |
| GET | `/api/trees?ids=` or `?rank=&taxon=` | A UPGMA (default) or `method=nj` tree of 3–500 genomes from their signatures: Newick with record ids as labels, and the leaf records |
| GET | `/api/landscape?group=` | The public catalogue on two axes: `cellular` (default), `viral` or `eukaryotic`, one genome per species, or `all`, every public genome on one set of axes. PCA of the signature as parallel arrays with domains, lineages and each axis's heaviest dinucleotides |
| GET | `/api/experiments` | The published P9 experiments, newest first: question, verdict (`pass`, `fail`, `context`), snapshot, script, summary numbers |
| GET | `/api/experiments/{name}` | One experiment and a page of its findings: `kind`, `order` = `score` (default) or `-score`, `offset`, `limit` (1–1000, default 100) |
| GET | `/api/sequences/{id}/findings` | What the published experiments say about one record; empty when nothing |
| GET | `/health`, `/alive` | Open without a key: whether the process is up |

Requests for each are in [Sequences.http](Sequences.http) and [Searches.http](Searches.http), and a whole walk-through in
[Scenario.http](Scenario.http).

## Experiments (P9)

Experiments run offline: `scripts/platform-checks/p9_*.cs` compute on the frozen snapshot database `sequence_platform_p9` and, with
`--publish`, write their run and findings here (`ExperimentStore`, migration 0012). A publish replaces the previous run of the same
name. Findings about records that have since been deleted, or that aren't public, are skipped. The API only reads them; the UI shows
them under Experiments.

## Relatives, signatures and trees (P5 S5.5)

These answer at once from `SignatureIndex`, an in-memory table of every genome's dinucleotide signature: log ρ\* over both strands from
the stored dinucleotide counts, compared by Euclidean distance. It is Karlin's genomic signature and also the complete second-order
Z-curve (its Z-curve coefficients are a rotation of it), so both views describe one vector. One record stands for each catalogue genome
(the largest chromosome of a bacterium or archaeon; every virus record; every record of an upload); plasmids are left out. The index
loads in the background when the API starts (about 3 s for the local catalogue) and reloads when records are added or deleted. On the
local catalogue (2026-10-30, `results/S55_RelativesCheck.csv`): a relatives query takes 4–9 ms (median), E. coli K-12's 10 relatives
are 9 Escherichia or Shigella and one Citrobacter, SARS-CoV-2's are 10 coronaviruses, and a 500-genome tree takes 0.4 s (UPGMA) or 1.0 s
(neighbour-joining).

The map is built on first request from public records only, so nobody's upload moves it, and is kept until the record count or highest
id changes: after `load-taxonomy` alone, restart the API. On the per-group maps cellular genomes need a GTDB assignment and viruses of at
least 1 kbp carry NCBI's lineage; the `all` map holds exactly what `SignatureIndex` holds, labelled by GTDB where assigned and NCBI otherwise.
Eukaryotes are in the pilot's scope only: one record per assembly, from nuclear counts. Responses are compressed for clients that ask.
It uses log ρ\*, and the Z-curve and Karlin
views are a rotation of one vector, so there is one map. Two axes keep phyla apart but not genera (`results/S53_Map_Metrics.csv`), so
use relatives to find what a genome resembles. On the local catalogue (`results/P6_LocalStack.csv`, Debug): 5,480 cellular and 9,192 viral
genomes, 1.4 s and 1.2 s cold, 70 ms warm, 0.5 and 0.8 MB.

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

## Protein search

Which stored genes encode a protein (P7): the protein is turned into the Z-curve windows its gene would have on codon positions 1
and 2 (synonymous codons weighted equally, genetic code 11 or 4) and compared with the genes BV-BRC annotates on the matching
records (`gene_profile`, migration 0017, filled by `backfill-genes`). A protein of 500 aa or more has windows of 333 positions from
the start codon and meets the genes of 500 aa and up. A protein of 50–499 aa is split into 5 equal windows and meets only the genes
within 10% of its length, so the search reads only those length buckets. Proteins under 50 aa are refused with 400. Uploads and NCBI
eukaryotes have no gene rows, so only BV-BRC records are searched. The job keeps the windows, the length and the SHA-256 of the letters, never the letters (migration
0015). A hit is one gene: its record, `geneStartBp` (0-based) and `geneEndBp` (exclusive, stop codon included) on the stored
strand, its `strand`, its BV-BRC protein family (`familyId`, PGF_…) and that family's `product`. There is one list, Ranked, and no p
or E yet.

Measured over the local catalogue (58,141 records, 64.3M genes), from BV-BRC's own proteins alone
(`results/ProteinSearch_CatalogueV2.csv`, 301 proteins from 13 genomes): the right family comes first for every protein that has
one, in every length band. The exact gene comes first for 95.0% at 50–99 aa, 95.2% at 100–166, 90.0% at 167–332, 88.1% at 333–499
and 71.7% at 500 aa and up, and is in the top 10 for 95–100%. When the exact gene isn't first, the gene that is belongs to the
same family. A search takes 2–8 s. `weakSignal` stays false for proteins.

Each hit carries p and E (`log10P`, `pKind`, `log10E`) once `build-protein-null` has run: E is how many genes of other protein
families would score this high in a search of this size, from real proteins scored against the gene layer. `pKind` Bound means fewer
than 10 null pairs scored this high, so E is an upper bound (the page shows "< E"). On the check behind it (`results/ProteinNull_Check.csv`),
60–66% of proteins have E < 0.01 at their best same-family hit and 75–81% have E < 1; E near 1 or above means the hit is unconfirmed.

## The worker

`SearchWorker` runs queued searches one at a time: each search already uses every core. It runs in the API process by default, and
in a process of its own when the API sets `SearchWorker:Enabled=false` and a second instance keeps it on, as Aspire does. PostgreSQL
is the queue, so nothing else changes: across containers a job waits 2–3 ms, and warm find-similar takes 34–40 ms from POST to Done
(`results/P6_LocalStack.csv`). It wakes on a
PostgreSQL notification (`search_queued`) and also polls every 30 s, so a lost listener only delays jobs. It assumes it is the
only worker: at startup a job still `Running` lost its worker and is requeued, and after three claims
it is failed instead. Every 10 minutes it applies D7's expiry and reloads the score null and neighbour graphs, which the Ingest
rebuilds; the cached profiles are refreshed before every search, so a record stored a moment ago is a candidate.

Find-similar keeps its band levels in memory (24 MB at ratio 0.5). Localization streams from the database by default
(`SearchWorker:CacheLocalization=false`, the server's setting): 0.85-9.6 s per fragment from 128 down to 4 kbp. Local runs (the `http`
launch profile and the AppHost in run mode) set it to true: every level in memory, 0.5-1.6 s, for about 0.8 GB
(`results/SearchBreakdown_*.csv`). Protein search streams the gene rows its length band needs (0.25-1.05 GB) from the database for
every query, 2.8-7.8 s locally, so it holds the worker that long and other searches queue behind it. It scores them in place without
hashing them (`verify-genes` in the Ingest checks the hashes), so most of that time is PostgreSQL's transfer. Other
settings in the `SearchWorker` section: `Enabled`, `PollInterval`, `StuckAfter`, `MaintenanceInterval`.

On the local catalogue (`results/S43_WorkerTiming.csv`, Release): warm find-similar takes 150-650 ms from POST to Done (the
upper end is profiling a 4.6 Mbp FASTA query in the POST), localization 0.9-2.3 s; the first search after a start takes about
3 s while the caches, null and graph load. The process used 451 MB.
