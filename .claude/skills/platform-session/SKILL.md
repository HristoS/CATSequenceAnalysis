---
name: platform-session
description: Start or continue a work session on the sequence platform (SequencePlatform.* projects) using the plan in docs/platform/README.md. Use when the user says /platform-session, "continue the platform", "next platform step", or names a platform phase (P0–P9) or open item.
---

# Platform session

The plan is one file, `docs/platform/README.md` (git-ignored by choice, so it exists only on the author's machine).
If it is missing, stop and say so. Don't reconstruct the plan from memory.

## Start
1. Read `docs/platform/README.md`. It has the status and next-step candidates, the rules, what each phase built, the open work,
   and the decisions in force.
2. Pick the item: the argument if one was given, otherwise ask which of the "Next step" candidates to take.
3. Read only the code and `results/README.md` sections the item needs. Background that the plan no longer carries is in
   `docs/_archive/docs-before-cleanup-2026-11-12.zip` (old phase documents and the full decision log); open it only when a question
   needs it.
4. Before writing code, write the item's card (verdict, gain, cost, assumption and check, cheaper alternative, stop criterion)
   with current numbers, and wait for a go-ahead. This is the global CLAUDE.md rule. Skip it only for checks, documentation,
   or when the user says "just do it".

## Work
- Keep to the item's scope; a new idea goes into "Open work" as a line, not into the code.
- Apply the plan's rules: profile only and never DNA; migrations are append-only; no NuGet in `SequenceAnalise`; a measurement
  gets a summary line in `results/README.md` (and a CSV only if a paper or `CLAUDE.md` cites it).
- Load `cat-code-style` before writing C#.
- If a stop criterion triggers, stop and report it with the numbers. Don't work around it.

## End
1. Update "Status" and remove or tick the finished item in "Open work".
2. Add a row to "Decisions in force" for any decision (replace the row it supersedes).
3. Close with a short summary: what was done, numbers, what's next.
