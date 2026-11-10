-- 0021_sequence_app.sql: the least-privilege role the API and the search worker connect as (P6 S6.4, database-decision.md §4.2 rule 3),
-- so a flaw in the internet-facing API can't use a superuser's powers (COPY ... PROGRAM, DDL, other databases, roles).
--
-- The role is created without login and without a password: credentials stay out of the repo. `migrate` gives it a login and the
-- password from SEQUENCE_APP_PASSWORD when that is set (the hosted Compose stack); migrations keep running as the owner. A role belongs
-- to the whole server, so a second database migrated on it (the test databases) reuses it and only adds its grants.
--
-- It reads every table. It writes only what the API and the worker write (found 2026-11-10): an upload (RecordStore.StoreAsync) inserts
-- its source, record, ladder levels and 6-mer counts, the counts with `on conflict do update`; a search job is inserted, claimed
-- (`for update skip locked` needs update), finished or expired in search_run, and its hits inserted, or deleted on expiry. LISTEN and
-- pg_notify need no grant. A later migration adding a table the API writes must grant that write itself; reads come by default below.
--
-- Keep this file safe to run again: a dump carries no roles, and hosting/tencent/seed-database.md restores with --no-privileges, so after
-- a restore it is rerun by hand (the restored schema_migration already lists it).

do $$
begin
  if not exists (select from pg_roles where rolname = 'sequence_app') then
    create role sequence_app nologin;
  end if;
end
$$;

do $$
begin
  execute format('grant connect on database %I to sequence_app', current_database());
end
$$;

grant usage on schema public to sequence_app;
grant select on all tables in schema public to sequence_app;
grant usage on all sequences in schema public to sequence_app;

grant insert on sequence_source, sequence_record, profile_level to sequence_app;
grant insert, update on record_kmer to sequence_app;
grant insert, update on search_run to sequence_app;
grant insert, delete on search_hit to sequence_app;

-- Tables and sequences that later migrations create (as the role running them) are readable without another grant.
alter default privileges in schema public grant select on tables to sequence_app;
alter default privileges in schema public grant usage on sequences to sequence_app;
