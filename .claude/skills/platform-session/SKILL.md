---
name: platform-session
description: Start or continue a work session on the sequence platform (SequencePlatform.* projects) using the phase tracker in docs/platform/. Use when the user says /platform-session, "continue the platform", "next platform step", or names a platform phase (P0–P7). Loads only the tracker and one phase document to keep context small.
---

# Platform session

The plan lives in `docs/platform/` (git-ignored by choice, so it exists only on the author's machine).
If `docs/platform/README.md` is missing, stop and say so. Don't reconstruct the plan from memory.

## Start
1. Read `docs/platform/README.md`. It names the current phase and next step, the rules, and the decision log.
2. Pick the phase: the argument if one was given (`P0`–`P7`), otherwise the README's "Current phase".
3. Read only that phase document in `docs/platform/phases/`, and then only the files in its "Read only these" list.
   Open `docs/platform/overview.md` only when a question needs background that the phase document doesn't carry.
4. Before writing code for a step, restate the phase's **Card** (verdict, gain, cost, assumption and check, cheaper alternative,
   stop criterion) with current numbers, and wait for a go-ahead. This is the global CLAUDE.md rule. Skip it only for checks,
   documentation, or when the user says "just do it".

## Work
- Do the next unchecked step. Keep to its scope; a new idea goes into the phase's session log as a note, not into the code.
- Apply the README rules: profile only and never DNA; migrations are append-only; no NuGet in `SequenceAnalise`;
  numbers go to `results/` with a README row.
- Load `cat-code-style` before writing C#. Use the skills the phase document names.
- If a stop criterion triggers, stop and report it with the numbers. Don't work around it.

## End
1. Tick the finished steps (`- [x]`) in the phase document.
2. Append one row to the phase's session log: date, what was done, measurements, next step.
3. Update the phase row in the README table, the "Current phase / Next step" line, and the decision log for any decision or measurement.
4. Close with a short summary: steps done, numbers, what's next.
