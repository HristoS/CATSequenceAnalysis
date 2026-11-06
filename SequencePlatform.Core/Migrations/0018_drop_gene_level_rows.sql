-- 0018_drop_gene_level_rows.sql: platform step S7.S (2026-11-05). Gene profiles live in gene_profile (0017); the S7.S pass replaced
-- every 'zcurve-cds' profile_level row, so the rows and 0014's exception for them go: every profile_level row has a window again.
delete from profile_level where method = 'zcurve-cds';
alter table profile_level drop constraint profile_level_window_count_check;
alter table profile_level add constraint profile_level_window_count_check check (window_count >= 1);
