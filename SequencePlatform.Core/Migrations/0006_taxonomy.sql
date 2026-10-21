-- 0006_taxonomy.sql: S2.4. sequence_record keeps each provider's lineage as ingested (BV-BRC's for the bulk load); these
-- tables hold current reference lineages beside it. The 2026-10-20 check found BV-BRC's names stale on about 3% of
-- families (results/S24_LineageCheck.csv), so evaluations against taxonomy should read these, and lineage_difference shows
-- where the two disagree.

-- NCBI Taxonomy's current lineage for each taxid sequence_record references, from the NCBI Datasets API.
create table taxon (
  taxid            integer primary key,   -- as referenced by sequence_record.ncbi_taxid
  status           text not null check (status in ('current', 'merged', 'unknown')),
  current_taxid    integer,               -- NCBI's id today: taxid itself unless merged into another; null when unknown
  name             text,                  -- current scientific name
  rank             text,                  -- lowercase NCBI rank of current_taxid: species, strain, no rank, ...
  tax_superkingdom text,                  -- the domain for cellular life (Bacteria, Archaea); Viruses for viruses
  tax_realm        text,
  tax_kingdom      text,
  tax_phylum       text,
  tax_class        text,
  tax_order        text,
  tax_family       text,
  tax_genus        text,
  tax_species      text,
  fetched_at       timestamptz not null default now(),
  check ((status = 'unknown') = (current_taxid is null))
);

-- GTDB taxonomy (built from a genome phylogeny, so the fairer ground truth for R3) for the bacterial and archaeal
-- assemblies sequence_record references. Keyed by the assembly accession without its GCA_/GCF_ prefix, because a
-- GenBank assembly and its RefSeq copy share the digits and version, and BV-BRC mostly gives the GCA_ form.
create table gtdb_taxon (
  assembly_key     text primary key,      -- '000005845.2'
  gtdb_accession   text not null,         -- as GTDB lists it: RS_GCF_000005845.2 or GB_GCA_...
  release          text not null,         -- 'v232'
  tax_domain       text,
  tax_phylum       text,
  tax_class        text,
  tax_order        text,
  tax_family       text,
  tax_genus        text,
  tax_species      text,
  fetched_at       timestamptz not null default now()
);

-- One row per record and rank where the provider's lineage and NCBI's current one differ: 'differs' when both have a
-- name, 'ncbi_only' or 'provider_only' when one of them is empty.
create view lineage_difference as
select r.id as sequence_id, r.ncbi_taxid, v.rank, v.provider_value, v.ncbi_value,
       case when v.provider_value is null then 'ncbi_only'
            when v.ncbi_value is null then 'provider_only'
            else 'differs' end as outcome
from sequence_record r
join taxon t on t.taxid = r.ncbi_taxid
cross join lateral (values
  ('superkingdom', r.tax_superkingdom, t.tax_superkingdom),
  ('realm', r.tax_realm, t.tax_realm),
  ('kingdom', r.tax_kingdom, t.tax_kingdom),
  ('phylum', r.tax_phylum, t.tax_phylum),
  ('class', r.tax_class, t.tax_class),
  ('order', r.tax_order, t.tax_order),
  ('family', r.tax_family, t.tax_family),
  ('genus', r.tax_genus, t.tax_genus),
  ('species', r.tax_species, t.tax_species)) as v(rank, provider_value, ncbi_value)
where v.provider_value is distinct from v.ncbi_value;
