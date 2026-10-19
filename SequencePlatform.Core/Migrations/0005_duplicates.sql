-- 0005_duplicates.sql: exact-duplicate sequences are skipped at ingest (decided 2026-10-19). BV-BRC lists a GenBank
-- record and its RefSeq copy (HM119410 / NC_055339) under two genome ids; about 21% of the reference viruses so far.
-- A genome is a duplicate when every record is already stored: recognised before download by the provider's MD5,
-- or after profiling by sequence_sha256 (0004) when the provider gives no MD5.

alter table sequence_record
  -- The provider's MD5 of the record sequence as listed (BV-BRC genome_sequence.sequence_md5), so a duplicate is
  -- found without downloading it. Only compared within one provider; null when the provider gives none.
  add column provider_md5 text;

create index sequence_record_provider_md5_idx on sequence_record (provider_md5);

alter table crawl_item
  -- For status 'duplicate': the stored source holding the same sequences.
  add column duplicate_of bigint references sequence_source(id) on delete set null,
  add check (duplicate_of is null or status = 'duplicate');
