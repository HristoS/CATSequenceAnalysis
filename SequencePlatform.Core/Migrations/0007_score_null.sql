-- 0007_score_null.sql: P3 step S3.6 (R5). Empirical null distributions of best-alignment scores between unrelated
-- records, one row per cell (SequenceAnalyses.ScoreNullCell): level, bucket of the shorter profile's window count,
-- bucket of the offsets tried. Built offline from stored profiles only (ScoreNullBuilder); a rebuild replaces the
-- rows of its method, version and source. p-values are empirical only; no tail model is stored, because none fitted
-- (results/S36_TailCalibration.csv).
--   genome:   whole records against unrelated records in a find-similar length band (the p and E SimilaritySearch reports)
--   fragment: slices of 500 bp windows against longer unrelated records (what localization, S3.5, will use)
--   permuted: whole records against window-permuted records of any taxonomy, a control only
create table score_null (
  method           text not null,
  method_version   text not null,
  source           text not null check (source in ('genome', 'fragment', 'permuted')),
  window_size      integer not null,
  windows_bucket   smallint not null,
  offsets_bucket   smallint not null,
  scores           real[] not null,             -- sorted ascending; both query strands, the better one
  built_at         timestamptz not null default now(),
  primary key (method, method_version, source, window_size, windows_bucket, offsets_bucket)
);
