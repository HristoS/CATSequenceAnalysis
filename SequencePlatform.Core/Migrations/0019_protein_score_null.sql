-- 0019_protein_score_null.sql: P7, p and E for protein hits (results/ProteinNull_Check.csv). Per query length bucket (the gene_profile
-- buckets: 0 from 500 aa, 1-25 for 50 * 1.1^(b-1) aa), how often genes of other BV-BRC families reach each score against real proteins
-- compared as a protein search compares them. Built by build-protein-null from proteins fetched at build time; no protein is stored.
-- exceed[i + 1] (1-based) = pairs scoring at least i / (cardinality(exceed) - 1); exceed[1] = pairs.
create table protein_score_null (
  method_version text not null,
  bucket         smallint not null check (bucket between 0 and 25),
  queries        integer not null check (queries > 0),
  pairs          bigint not null check (pairs > 0),
  exceed         bigint[] not null check (cardinality(exceed) > 1),
  built_at       timestamptz not null default now(),
  primary key (method_version, bucket)
);
