-- 0009_search_jobs.sql: search_run as an asynchronous job (P4 S4.1), and search_hit reshaped to hold what P3's searches return.
--
-- search_hit is dropped and recreated rather than altered: it was never written (0 rows on every database with 0002,
-- 2026-10-26), and its refined_score / offset_window described a 500 bp re-rank and a single hit list that P3 dropped.
--
-- The worker runs after an upload's DNA is gone, so a job carries everything the search needs: the fine profile and SHA-256
-- of both strands (the reverse strand is exact only when profiled from DNA, results/C3_ProfileFlip.csv) and the exact length.
-- Coarse levels and the origin rotation are derived from the fine profile. D7: after 30 days, or on DELETE, a run keeps only
-- its hashes, length, parameters and summary; profiles and hits go, and status becomes 'expired'.

drop table search_hit;

alter table search_run
  alter column query_sha256 drop not null,                -- null for a profile upload: there was no DNA to hash
  add column public_id             uuid not null default gen_random_uuid() unique,
  add column status                text not null default 'queued'
                                   check (status in ('queued','running','done','failed','expired')),
  add column error                 text,
  add column started_at            timestamptz,
  add column completed_at          timestamptz,
  add column expires_at            timestamptz not null default now() + interval '30 days',
  add column query_reverse_sha256  bytea,
  add column query_reverse_profile bytea,                 -- fine level, ZCWP; null when a profile upload gave one strand
  add column length_ratio          double precision,      -- find-similar's band ratio
  add column within_sequence_id    bigint,                -- localize inside this record only; no foreign key, so deleting
                                                          -- the record doesn't silently turn a queued job into a catalogue scan
  add column band_min_bp           bigint,
  add column band_max_bp           bigint,                -- null: unbounded
  add column ranked_total          integer,               -- find-similar: ranked candidates before top K
  add column candidates            integer,               -- localize: records scanned
  add column weak_signal           boolean,               -- find-similar FewWindows, localize WeakSignal
  add column neighbour_ranked      boolean;

comment on column search_run.level_window_size is
  'find-similar: the level every hit is scored at; localize: the coarse scan level';

create index search_run_queued on search_run (created_at) where status = 'queued';
create index search_run_expiry on search_run (expires_at) where status <> 'expired';

create table search_hit (
  search_run_id      bigint not null references search_run(id) on delete cascade,
  list_kind          text not null check (list_kind in ('ranked','identical','not_comparable')),
  position           integer not null,     -- 1-based order within the list; for 'ranked', the rank
  sequence_id        bigint not null references sequence_record(id) on delete cascade,
  strand             text check (strand in ('forward','reverse')),   -- null for 'not_comparable'
  score              double precision,     -- 'ranked' only
  coarse_score       double precision,     -- localize: the score that brought the record into the refined top
  offset_bp          bigint,               -- localize: where the fragment starts on the stored strand
  log10_p            double precision,
  p_kind             text check (p_kind in ('empirical','bound')),
  log10_e            double precision,
  shared_neighbours  smallint,             -- find-similar, neighbour-ranked top hits
  rotated            boolean not null default false,
  primary key (search_run_id, list_kind, position),
  check ((log10_p is null) = (p_kind is null))
);

create index search_hit_sequence on search_hit (sequence_id);
