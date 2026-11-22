-- 0024_experiment_jobs.sql: rerunning a database experiment from the experiments page (Phase 1, 2026-11-21).
--
-- * app_user.is_admin: only an admin may queue a rerun. Set by the Ingest (add-user --admin, set-admin); the API only reads it.
-- * experiment_job: a rerun as an asynchronous job, queued by the API and run by the experiment worker, a third process beside the API
--   and the search worker. PostgreSQL stays the queue, as for searches (0009): claimed with `for update skip locked`, announced on the
--   channel experiment_queued, requeued when its worker stopped. One queued or running job per experiment at a time.
-- * experiment_query: a stored query set, so an experiment whose queries came from outside the database (protein-homolog's proteins,
--   from BV-BRC) can be rerun from the database alone. Never sequences: per query, the record and the fields the experiment needs, and
--   the query profile (protein-homolog: QueryFromProtein's expected axes). Filled by the Ingest's seed-experiment-queries, as the owner.
--
-- The worker connects as sequence_app (0021) and publishes its result as an experiment_run of origin 'server'. Row-level security keeps
-- that role to server runs: it reads every run but inserts or deletes only server runs, and adds files and findings only to them, so a
-- flaw in the API can't replace a run an article cites. Table owners (migrate, the Ingest, the scripts) bypass the policies.
--
-- Like 0021, the grants below are lost by a restore with --no-privileges (hosting/tencent/seed-database.md): rerun the grant lines of
-- this file by hand after 0021's. The policies are part of the schema and survive.

alter table app_user add column is_admin boolean not null default false;

create table experiment_job (
  id            bigint generated always as identity primary key,
  public_id     uuid not null default gen_random_uuid() unique,
  experiment    text not null,
  requested_by  bigint not null references app_user(id),
  status        text not null default 'queued' check (status in ('queued', 'running', 'done', 'failed')),
  error         text,
  attempts      smallint not null default 0,
  created_at    timestamptz not null default now(),
  started_at    timestamptz,
  completed_at  timestamptz,
  run_id        bigint references experiment_run(id) on delete set null   -- the server run it published
);

create index experiment_job_queued on experiment_job (created_at) where status = 'queued';
create index experiment_job_experiment on experiment_job (experiment, created_at desc);
create unique index experiment_job_active on experiment_job (experiment) where status in ('queued', 'running');

create table experiment_query (
  experiment   text not null,
  position     integer not null,
  sequence_id  bigint not null references sequence_record(id) on delete cascade,
  fields       jsonb not null,
  profile      double precision[] not null,
  primary key (experiment, position)
);

grant insert, update on experiment_job to sequence_app;
grant insert, delete on experiment_run to sequence_app;
grant insert on experiment_file, experiment_finding to sequence_app;

alter table experiment_run enable row level security;
create policy experiment_run_read on experiment_run for select to sequence_app using (true);
create policy experiment_run_insert on experiment_run for insert to sequence_app with check (origin = 'server');
create policy experiment_run_delete on experiment_run for delete to sequence_app using (origin = 'server');

alter table experiment_file enable row level security;
create policy experiment_file_read on experiment_file for select to sequence_app using (true);
create policy experiment_file_insert on experiment_file for insert to sequence_app
  with check (exists (select 1 from experiment_run r where r.id = run_id and r.origin = 'server'));

alter table experiment_finding enable row level security;
create policy experiment_finding_read on experiment_finding for select to sequence_app using (true);
create policy experiment_finding_insert on experiment_finding for insert to sequence_app
  with check (exists (select 1 from experiment_run r where r.id = run_id and r.origin = 'server'));
