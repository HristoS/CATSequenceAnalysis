-- 0011_search_attempts.sql: how many times the worker has claimed a search (P4 S4.3). A job still running when its worker
-- stopped is requeued; after three claims it is marked failed instead, so a query that brings the worker down can't loop.

alter table search_run
  add column attempts smallint not null default 0;
