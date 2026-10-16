-- 0004_record_features.sql: the D8 opt-ins (decided 2026-10-15), per record, computed at ingest in the
-- same streaming pass as the profile (SequenceFeatureReader). Still no bases: a hash and 16 counts.
-- Nullable, because a precomputed ZCWP profile upload (D2) arrives without its DNA.

alter table sequence_record
  -- SHA-256 of the record's ASCII letters, uppercased, nothing else: "identical" means this matches.
  add column sequence_sha256 bytea check (sequence_sha256 is null or length(sequence_sha256) = 32),
  -- Adjacent A/C/G/T pairs, in AA, AC, AG, AT, CA, ..., TT order; other letters break a pair.
  add column dinucleotide_counts bigint[] check (dinucleotide_counts is null or array_length(dinucleotide_counts, 1) = 16);

create index sequence_record_sha256_idx on sequence_record (sequence_sha256);
