-- 0010_users.sql: users and owners for a small group with private uploads (D3, P4 S4.2a).
--
-- A user holds an API key; only its SHA-256 is stored. The key is 32 random bytes, so a plain hash is enough: there is no
-- low-entropy password for a slow hash to protect. An upload belongs to its user and is visible only to them; every
-- downloaded source (owner null) is public. A search belongs to whoever ran it; runs from the Ingest have no owner.

create table app_user (
  id              bigint generated always as identity primary key,
  name            text not null unique,
  api_key_sha256  bytea not null unique check (length(api_key_sha256) = 32),
  created_at      timestamptz not null default now(),
  disabled_at     timestamptz                     -- set to revoke the key; the user's uploads stay
);

alter table sequence_source
  add column owner_id bigint references app_user(id),
  add constraint sequence_source_owner_upload check ((owner_id is not null) = (provider = 'user upload'));

create index sequence_source_owner on sequence_source (owner_id) where owner_id is not null;

alter table search_run
  add column owner_id bigint references app_user(id);

create index search_run_owner on search_run (owner_id, created_at);
