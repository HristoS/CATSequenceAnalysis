-- 0012_experiments.sql: P9 U0. Results of the P9 experiments, computed offline on a frozen snapshot and published here by their
-- scripts (scripts/platform-checks/p9_*.cs --publish) so the web UI can show them. One row per experiment: a publish replaces the
-- previous run of the same name, findings included. A finding is about one public record; it goes when the record goes.
create table experiment_run (
  id           bigint generated always as identity primary key,
  name         text not null unique,
  title        text not null,
  question     text not null,
  verdict      text not null check (verdict in ('pass', 'fail', 'context')),
  snapshot     text not null,
  script       text not null,
  summary      jsonb not null,
  published_at timestamptz not null default now()
);

create table experiment_finding (
  run_id       bigint not null references experiment_run(id) on delete cascade,
  sequence_id  bigint not null references sequence_record(id) on delete cascade,
  kind         text not null,
  score        double precision,
  payload      jsonb not null,
  primary key (run_id, kind, sequence_id)
);
create index experiment_finding_sequence_idx on experiment_finding (sequence_id);
