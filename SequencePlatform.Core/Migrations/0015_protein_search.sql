-- 0015_protein_search.sql: protein search (P7, 2026-11-05), a third search_run kind. The query is an amino-acid sequence, scored
-- against every stored gene row (0014) from the protein alone (results/ProteinSearch_Catalogue*.csv). As with DNA, the protein lives
-- only for the POST: a job keeps its SHA-256 (of the upper-case letters, a final stop dropped), its length and its expected Z-curve
-- windows (ProteinQuery), never the letters. For a protein job, query_profile holds those windows as little-endian float64, 9 per
-- window, and query_length_bp is the coding length the protein implies, 3 x query_length_aa.

alter table search_run drop constraint search_run_search_kind_check;
alter table search_run add constraint search_run_search_kind_check check (search_kind in ('similar', 'localize', 'protein'));

alter table search_run
  add column query_length_aa integer check (query_length_aa >= 1),
  add column genetic_code    smallint check (genetic_code in (4, 11)),   -- NCBI translation table the query was read with
  add constraint search_run_protein_check check ((search_kind = 'protein') = (query_length_aa is not null and genetic_code is not null));

-- A protein hit is one gene of a record: its coordinates on the stored record (0-based start, exclusive end, stop codon included, as
-- BV-BRC lists it), strand in search_hit.strand, and its BV-BRC PGfam number (protein_family.id; null when it has none).
alter table search_hit
  add column gene_start integer,
  add column gene_end   integer,
  add column family_id  integer;
