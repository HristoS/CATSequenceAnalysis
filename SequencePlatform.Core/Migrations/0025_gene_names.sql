-- 0025_gene_names.sql: each stored gene's name as its provider gives it (2026-11-22, D8 widened by the author's decision). BV-BRC's
-- product for every coding sequence of a BV-BRC record, NCBI's [protein=...] for every gene of an annotated eukaryote. Until now a gene
-- had a name only through its protein family (protein_family.product), so genes without one, most virus and every eukaryote gene,
-- had none. The 40% rule set before the check (results/README.md section 2) passed only for eukaryotes (54.7% of genes informatively
-- named); BV-BRC's genes without a family are informatively named in 9-22% of cases, and the author chose to store every group anyway.
-- Names are annotation text, never sequence.
--
-- Names are kept once in gene_name_text; a record's genes point at them in one integer array, in the order of the record's stored gene
-- profile (by sequence, then start: ZCurveGeneProfile.Combine over its gene_profile rows), 0 where the provider gave no name. A record
-- with a gene_name row is done, so the backfill resumes where it stopped; a record whose gene profile is rewritten gets its names
-- rewritten in the same transaction. Written only by the Ingest; sequence_app reads both tables by 0021's default privileges.
create table gene_name_text (
  id   integer generated always as identity primary key,
  name text not null unique check (length(name) between 1 and 1000)
);

create table gene_name (
  sequence_id bigint primary key references sequence_record(id) on delete cascade,
  source      text not null check (source in ('bvbrc', 'ncbi')),
  names       integer[] not null,
  named       integer not null check (named >= 0),   -- how many genes have a name
  unmatched   integer not null default 0 check (unmatched >= 0),   -- stored genes the provider's list didn't have (no name)
  checked_at  timestamptz not null default now()
);
