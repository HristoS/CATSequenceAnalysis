-- 0002_search_log.sql: optional search log for reproducibility and result caching (docs/database-decision.md rev 2, section 5.1).

create table search_run (
  id                 bigint generated always as identity primary key,
  created_at         timestamptz not null default now(),
  search_kind        text not null check (search_kind in ('similar','localize')),
  query_sha256       bytea not null,      -- the query's bases are never stored
  query_length_bp    bigint not null,
  query_profile      bytea,               -- fine level in the §4.3 format (query DNA is discarded)
  level_window_size  integer,             -- the single level this query ranked at
  filters            jsonb not null default '{}'::jsonb,
  top_k              integer not null
);

create table search_hit (
  search_run_id      bigint not null references search_run(id) on delete cascade,
  sequence_id        bigint not null references sequence_record(id),
  rank               integer,
  coarse_score       double precision,    -- null = not comparable, never 0
  refined_score      double precision,    -- 500 bp re-rank for the top K
  offset_window      integer,             -- localization only
  primary key (search_run_id, sequence_id)
);
