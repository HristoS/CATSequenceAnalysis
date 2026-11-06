-- 0017_gene_profile.sql: platform step S7.S (2026-11-05). A record's gene profiles (ZCGP version 2) move out of profile_level into their own
-- table, split by protein length so a search reads only what it compares. Bucket 0 holds the record's genes of 500 aa and more, in
-- 333-position windows; bucket b (1-25) its genes of 50 * 1.1^(b-1) to 50 * 1.1^b aa (capped at 499), each in 5 equal windows, which a
-- query is compared with only within 10% of its own length (results/ShortProtein_*.csv, GeneLayerV2_Sizes.csv). Genes under 50 aa aren't
-- kept. Every record the gene pass finishes has a bucket 0 row, possibly with no gene, so a rerun skips it. Each record's old 'zcurve-cds'
-- profile_level row is deleted when its new rows are written.

create table gene_profile (
  sequence_id    bigint not null references sequence_record(id) on delete cascade,
  bucket         smallint not null check (bucket between 0 and 25),
  format_version smallint not null,
  gene_count     integer not null check (gene_count >= 0),
  window_count   integer not null check (window_count >= 0),
  payload        bytea not null,
  byte_size      integer not null,
  sha256         bytea not null,
  primary key (sequence_id, bucket)
);

create index gene_profile_bucket on gene_profile (bucket);
