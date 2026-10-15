-- 0001_catalogue.sql: catalogue and profile levels (docs/database-decision.md rev 2, section 5.1).
-- The DNA sequence itself is never stored anywhere in this schema.

create table sequence_source (
  id              bigint generated always as identity primary key,
  provider        text not null,          -- 'NCBI RefSeq', 'ENA', 'user upload', ...
  url             text,                   -- where it was downloaded from; null only for user uploads
  accession       text,                   -- versioned: NC_045512.2, GCF_000005845.2
  retrieved_at    timestamptz not null,
  raw_sha256      bytea not null,         -- of the downloaded file; the file itself is discarded
  raw_byte_size   bigint not null,
  license         text,                   -- provider licence/terms, for attribution
  check (url is not null or provider = 'user upload'),
  unique (provider, accession)
);

-- One catalogue entry per FASTA record (per chromosome for eukaryotes). Holds no bases.
create table sequence_record (
  id              bigint generated always as identity primary key,
  source_id       bigint not null references sequence_source(id),
  record_id       text not null,          -- FASTA record id within the source file
  name            text,
  kind            text not null check (kind in ('virus','bacterium','archaeon','eukaryote','plasmid','organelle','other')),
  data_class      text not null check (data_class in ('public_microbe','public_eukaryote','human','restricted')),
  length_bp       bigint not null check (length_bp > 0),
  gc_fraction     real,                   -- derived from the profile: (mean z + 1) / 2
  ncbi_taxid      integer,
  tax_superkingdom text,                  -- NCBI superkingdom/domain: Bacteria, Archaea, Eukaryota, Viruses
  tax_realm       text,
  tax_kingdom     text,
  tax_phylum      text,
  tax_class       text,
  tax_order       text,
  tax_family      text,
  tax_genus       text,
  tax_species     text,
  annotations     jsonb not null default '{}'::jsonb,   -- free-form: host, isolation source, notes
  created_at      timestamptz not null default now(),
  unique (source_id, record_id)
);
create index sequence_record_family_idx on sequence_record (kind, tax_family, tax_genus);
create index sequence_record_phylum_idx on sequence_record (kind, tax_phylum, tax_genus);
create index sequence_record_length_idx on sequence_record (length_bp);
create index sequence_record_annotations_idx on sequence_record using gin (annotations);

-- One row per ladder level that has at least one window.
create table profile_level (
  id               bigint generated always as identity primary key,
  sequence_id      bigint not null references sequence_record(id) on delete cascade,
  method           text not null default 'zcurve-windowed',
  method_version   text not null,
  fine_window_size integer not null,      -- 500
  factor           integer not null check (factor >= 1),
  window_size      integer generated always as (fine_window_size * factor) stored,
  window_count     integer not null check (window_count >= 1),
  format_version   smallint not null,
  payload          bytea,                 -- inline when <= 1 MB
  object_key       text,                  -- otherwise: key in object storage
  byte_size        bigint not null,
  sha256           bytea not null,
  created_at       timestamptz not null default now(),
  check ((payload is null) <> (object_key is null)),
  unique (sequence_id, method, method_version, window_size)
);
-- Worker cache load: every profile of one method/version at one level.
create index profile_level_load_idx on profile_level (method, method_version, window_size);
