-- 0014_gene_layer.sql: platform step S7.R (decided 2026-11-04). The gene profile of a record is a fifth profile_level row, method
-- 'zcurve-cds', payload ZCGP v1 (ZCurveGeneProfile): per coding sequence its coordinates, strand and BV-BRC PGfam number, and 16-bit
-- Legendre-moment windows over codon positions 1 and 2. Still no bases. Its window_count is the total over the record's genes, which is
-- 0 when every gene is shorter than one window (500 bp), so that row alone may have no window.
alter table profile_level drop constraint profile_level_window_count_check;
alter table profile_level add constraint profile_level_window_count_check check (window_count >= 1 or (method = 'zcurve-cds' and window_count >= 0));

-- Amino-acid k-mer counts of a record's translated coding sequences (relatives and trees: results/T2_CvTreeGtdb.csv). Payload by
-- format_version; v1: 20^k uint16 counts, little-endian, letters in the order ACDEFGHIKLMNPQRSTVWY with the first residue most
-- significant, deflate-compressed; a count above 65,535 is stored as 65,535. Words don't cross gene boundaries. Like record_kmer, a
-- row records why a record has no counts, so a rerun of the pass doesn't download it again.
create table protein_kmer (
  sequence_id    bigint not null references sequence_record(id) on delete cascade,
  k              smallint not null check (k in (3, 4)),
  format_version smallint not null,
  gene_count     integer check (gene_count >= 0),
  payload        bytea,
  -- done: counts stored; mismatch: the download differs from the stored record; missing: the download no longer has the record;
  -- no_cds: BV-BRC lists no coding sequence on it.
  status         text not null check (status in ('done', 'mismatch', 'missing', 'no_cds')),
  checked_at     timestamptz not null default now(),
  primary key (sequence_id, k),
  check ((status = 'done') = (payload is not null))
);

-- BV-BRC protein families named once, not per gene: the number of a PGF_ id (PGF_00421071 -> 421071) and its product.
create table protein_family (
  id         integer primary key,
  product    text not null,
  fetched_at timestamptz not null default now()
);
