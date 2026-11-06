# Seeding the hosted database from a local catalogue

Re-crawling BV-BRC on the server takes a day or more (the 2026-10 bulk load was 8.6 h of active time on a workstation), so the hosted
catalogue is restored from a dump of a local `sequence_platform` database instead. The dump holds profiles and metadata, never DNA.

## On the machine with the catalogue

**1. Pick the commit.** Push `feature/localPlatfromBuild`, then note the commit the server should run with `git rev-parse HEAD`.

**2. Bring the database up to that commit**, from the repository root:

```powershell
dotnet run --project SequencePlatform.Ingest -- migrate
```

If you've ingested since the derived tables were last built, rebuild them, in this order:

```powershell
dotnet run -c Release --project SequencePlatform.Ingest -- load-taxonomy
dotnet run -c Release --project SequencePlatform.Ingest -- build-origins
dotnet run -c Release --project SequencePlatform.Ingest -- build-null
dotnet run -c Release --project SequencePlatform.Ingest -- build-neighbours
```

Then save the counts, so the restore can be checked against them:

```powershell
dotnet run --project SequencePlatform.Ingest -- stats > stats.txt
```

**3. Dump** with PostgreSQL 17's `pg_dump`, the server's major version. The username is the one in the git-ignored
`SequencePlatform.Ingest/appsettings.json`. `pg_dump` asks for the password, so nothing goes on the command line.

**Windows (PowerShell):**

```powershell
$pg = "C:\Program Files\PostgreSQL\17\bin"
& "$pg\pg_dump.exe" --host localhost --port 5432 --username YOUR_DB_USER --dbname sequence_platform `
  --format custom --no-owner --no-privileges `
  --exclude-table-data search_run --exclude-table-data search_hit `
  --file sequence_platform.dump
```

**macOS / Linux:**

```bash
pg_dump --host localhost --port 5432 --username YOUR_DB_USER --dbname sequence_platform \
  --format custom --no-owner --no-privileges \
  --exclude-table-data search_run --exclude-table-data search_hit \
  --file sequence_platform.dump
```

What the options do:
- `--format custom`: compressed, and restorable with `pg_restore`.
- `--no-owner --no-privileges`: the server's database has its own role, so the local role names don't matter.
- `--exclude-table-data search_run / search_hit`: leaves out the local search history. The tables are kept, just empty.

**4. Check it:**

```powershell
& "$pg\pg_restore.exe" --list sequence_platform.dump | Select-Object -First 15
Get-FileHash sequence_platform.dump -Algorithm SHA256
```

On macOS/Linux, use `pg_restore --list … | head` and `shasum -a 256 sequence_platform.dump`. Expect a file of roughly the
database's size on disk, about 0.5 GB for the 2026-10 catalogue. The 2026-11-06 catalogue (61,060 records, with the gene layer) dumped
to 10.0 GB in 28 minutes from a 12.2 GB database: the 16-bit gene and level payloads barely compress.

**5. Upload** to the server's private upload folder. Admins get a Linux user and an SSH key from whoever runs the server; see the
[hosting README](README.md).

```powershell
scp -i ~/.ssh/cat_platform sequence_platform.dump stats.txt YOUR_USER@SERVER_IP:/opt/cat/incoming/
```

Then send the commit hash and the SHA-256 to the person doing the restore. **Never put the dump on GitHub or anywhere public**, and
never send `appsettings.json`, a connection string or a password: the server has its own.

## What the dump carries over

- **The catalogue:**
  - sources and records, with annotations
  - profile levels
  - taxonomy, the score null, neighbour lists, origins and the crawl state
- **The local API users** (`app_user`: names and SHA-256 hashes of their keys).
  - Disable the ones that shouldn't work on the public server after the restore. Keep the ones that should (on 2026-11-06,
    `hristo` stayed active).
  - Then create fresh users there.
- **Private uploads, if any.**
  - They stay owned by their (now disabled) local user, so they're invisible to everyone else.
  - Delete them if they shouldn't be on the server.

## On the server

From `SequencePlatform.AppHost/deploy/compose/`, with `COMPOSE_PROJECT` from `docker compose ls`.

> First run end to end on 2026-11-06 (commit `0bd3557`): upload 10.0 GB in about 37 min at about 4.5 MB/s, `pg_restore` 8 min,
> `vacuum analyze` 1 min, then `stats` matched the local counts exactly. Take a Lighthouse snapshot first if the server holds anything
> worth keeping, because this replaces the whole database, **including the server's own API users**: save the rows of the ones to
> keep beforehand (`copy (select * from app_user where name = '...') to stdout with (format csv, header)`) and insert them again
> afterwards without the `id` column, so their keys keep working.

```bash
c="docker compose -p COMPOSE_PROJECT -f docker-compose.yaml --env-file .env.Production"
sha256sum /opt/cat/incoming/sequence_platform.dump        # compare with the uploader's SHA-256

$c stop api worker
$c exec -T postgres sh -c 'PGPASSWORD="$POSTGRES_PASSWORD" psql -U postgres -v ON_ERROR_STOP=1 \
  -c "drop database if exists sequence_platform with (force)" -c "create database sequence_platform"'
$c exec -T postgres sh -c 'PGPASSWORD="$POSTGRES_PASSWORD" pg_restore -U postgres -d sequence_platform --no-owner --no-privileges' \
  < /opt/cat/incoming/sequence_platform.dump

$c exec -T postgres sh -c 'PGPASSWORD="$POSTGRES_PASSWORD" psql -U postgres -d sequence_platform -c "vacuum analyze"'   # pg_restore doesn't analyze
$c run --rm --no-deps -T migrate migrate                 # schema up to the deployed commit
$c run --rm --no-deps -T migrate stats                   # compare with stats.txt
$c exec -T postgres sh -c 'PGPASSWORD="$POSTGRES_PASSWORD" psql -U postgres -d sequence_platform -Atc "select name from app_user"'
$c run --rm --no-deps -T migrate disable-user NAME       # for each imported user
sudo systemctl restart cat-platform                      # api and worker reload the catalogue, null and graphs
```

Then create the real users with `add-user`, and check that relatives, the landscape and a find-similar answer against the restored
catalogue.
