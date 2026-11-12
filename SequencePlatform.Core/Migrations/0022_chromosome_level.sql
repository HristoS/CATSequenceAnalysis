-- 0022_chromosome_level.sql: windowed Z-curve levels per eukaryote chromosome (2026-11-11). An NCBI eukaryote is one record per assembly
-- with no profile_level rows (the pilot), so its Z-curve path couldn't be drawn. Each nuclear sequence of at least 100 kbp now keeps two
-- levels of its own, 500 bp and 32 kbp windows (ZCWP version 2, float16), about 6 B per 500 bp, so 3.1 GB over the 256 Gbp catalogue.
--
-- A table of its own rather than profile_level rows: everything that reads profile_level (the worker's cache, find-similar, localization,
-- score_null, the neighbour graph) reads every record's level of one window size, and a eukaryote's chromosomes are neither one record
-- nor comparable as one. There's no inline limit: the largest chromosome's 500 bp level is a few MB. Written only by the Ingest
-- (ingest-ncbi, backfill-chromosomes) as the owner; sequence_app reads it by 0021's default privileges.
create table chromosome_level (
  id               bigint generated always as identity primary key,
  sequence_id      bigint not null references sequence_record(id) on delete cascade,
  sequence_index   integer not null check (sequence_index >= 0),   -- the kept nuclear sequences, in the genome FASTA's order
  sequence_name    text not null,                                   -- its FASTA id (an accession.version)
  length_bp        bigint not null check (length_bp >= 1),
  method           text not null default 'zcurve-windowed',
  method_version   text not null,
  window_size      integer not null check (window_size >= 1),
  window_count     integer not null check (window_count >= 1),
  format_version   smallint not null,
  payload          bytea not null,
  byte_size        bigint not null,
  sha256           bytea not null,
  created_at       timestamptz not null default now(),
  unique (sequence_id, sequence_index, method, method_version, window_size)
);
