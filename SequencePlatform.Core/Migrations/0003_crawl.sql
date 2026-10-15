-- 0003_crawl.sql: per-item crawl state (docs/data-sources.md section 4, rule 5: resume, don't restart).
-- One row per provider item ever listed. A run skips items that reached a final state, so a re-run
-- downloads nothing it already handled; 'failed' items are retried only on request.

create table crawl_item (
  id              bigint generated always as identity primary key,
  provider        text not null,
  external_id     text not null,          -- the provider's id, e.g. BV-BRC genome id 1002689.3
  item_group      text,                   -- the listing it came from: 'viruses', 'bacteria', 'archaea', ...
  status          text not null default 'pending'
                  check (status in ('pending','done','failed','duplicate','too_large')),
  attempts        integer not null default 0,
  last_error      text,
  source_id       bigint references sequence_source(id),
  updated_at      timestamptz not null default now(),
  unique (provider, external_id),
  check ((status in ('done','too_large')) = (source_id is not null))
);
create index crawl_item_status_idx on crawl_item (provider, status);

-- Sources ingested before crawl state existed count as done.
insert into crawl_item (provider, external_id, status, attempts, source_id)
select provider, accession, 'done', 1, id from sequence_source where accession is not null and provider <> 'user upload';
