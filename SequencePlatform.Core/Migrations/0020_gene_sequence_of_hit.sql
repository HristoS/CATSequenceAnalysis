-- 0020_gene_sequence_of_hit.sql: P7 S7.T (2026-11-06), the gene layer for NCBI eukaryotes. A eukaryote record is a whole assembly, so its
-- genes (ZCGP version 3) lie on one of its sequences, numbered in the order the record's annotations list them (ncbi.gene_layer.sequences), and
-- gene_start/gene_end are on that sequence. A protein hit keeps the number; 0 for a BV-BRC record, whose genes lie on the record itself.
alter table search_hit add column gene_sequence integer check (gene_sequence >= 0);
