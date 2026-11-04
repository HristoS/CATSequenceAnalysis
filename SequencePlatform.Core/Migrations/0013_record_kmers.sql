-- 0013_record_kmers.sql: D8 widened (decided 2026-11-04): 6-mer counts per record, computed in the same streaming pass as the
-- profile (SequenceFeatureReader) for new records, and for records already stored by re-downloading them (Ingest backfill-kmers).
-- Still no bases: 4,096 counts, one strand, AAAAAA ... TTTTTT with the first base most significant; other letters break a run.
-- A backfilled record is stored only when the download is the sequence that was profiled (same SHA-256 and dinucleotide counts);
-- otherwise its row records why not, so a rerun doesn't download it again.
create table record_kmer (
  sequence_id   bigint not null references sequence_record(id) on delete cascade,
  k             smallint not null check (k = 6),
  counts        integer[] check (counts is null or array_length(counts, 1) = 4096),
  -- done: counts stored; mismatch: the download differs from the stored record; missing: the download no longer has the record.
  status        text not null check (status in ('done', 'mismatch', 'missing')),
  checked_at    timestamptz not null default now(),
  primary key (sequence_id, k),
  check ((status = 'done') = (counts is not null))
);
