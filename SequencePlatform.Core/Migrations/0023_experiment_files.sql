-- 0023_experiment_files.sql: every experiment an article states gets a page with its result files (2026-11-12). Until now only the P9
-- scripts published here, with no files. Now:
--
-- * experiment_file holds a run's result tables (CSV, gzip-compressed) for download, so the page's numbers can be checked against the
--   rows behind them. Never sequences: an experiment's files are tables of numbers.
-- * A run says which article section it backs, the commit and machine that computed it, the command that reproduces it, and whether
--   it can be rerun from the database alone (runnable). Experiments that need DNA (the library benchmarks on bundled genomes) are
--   published from their CSVs and are never runnable here: the platform doesn't download DNA for experiments.
-- * origin separates the run an article cites ('article', computed offline on a frozen copy) from a later rerun on this database
--   ('server'), so a page can show both; a name is unique per origin.
--
-- Written only by the Ingest and the scripts as the owner; sequence_app reads it by 0021's default privileges.
alter table experiment_run
  add column article       text,
  add column origin        text not null default 'article' check (origin in ('article', 'server')),
  add column source_commit text,
  add column machine       text,
  add column command       text,
  add column runnable      boolean not null default false;

alter table experiment_run drop constraint experiment_run_name_key;
alter table experiment_run add constraint experiment_run_name_origin_key unique (name, origin);

create table experiment_file (
  run_id     bigint not null references experiment_run(id) on delete cascade,
  name       text not null check (name ~ '^[A-Za-z0-9_.-]+$'),
  byte_size  bigint not null check (byte_size >= 0),        -- uncompressed
  content    bytea not null,                                -- gzip
  primary key (run_id, name)
);
