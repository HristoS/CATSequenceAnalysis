-- 0008_neighbour_list.sql: P3 step S3.4b. Each record's own find-similar top K (other assemblies only, best first) at
-- one length ratio, so a search can re-rank its top hits by how many neighbours the query and each hit share
-- (results/Check_MutualNeighboursSample.csv). Built offline from stored profiles (NeighbourGraph, CLI build-neighbours);
-- a rebuild replaces the rows of its method, version and ratio. neighbour_ids holds no foreign key, so a deleted record
-- can linger in other lists until the next rebuild; a deleted record's own row goes with it.
create table neighbour_list (
  method           text not null,
  method_version   text not null,
  length_ratio     real not null,
  sequence_id      bigint not null references sequence_record(id) on delete cascade,
  neighbour_ids    bigint[] not null,
  built_at         timestamptz not null default now(),
  primary key (method, method_version, length_ratio, sequence_id)
);
