---
name: cat-code-style
description: Apply this repo's C# naming and member-ordering conventions (extracted from DNA.cs, PatternSearchingKMP.cs, NeedlemanWunsch.cs) when writing or editing code in CATSequenceAnalysis. Use whenever adding new C# code, reviewing a diff for style, or doing a naming/ordering cleanup pass — never for changing algorithmic/optimization logic.
---

# CAT code style

Full derivation and evidence: [docs/code-conventions.md](../../../docs/code-conventions.md) in
this repo. Read it once if you need the "why"; this file is the quick-reference checklist for
applying it while writing or editing code.

The baseline was `SequenceAnalise/DNA.cs`, `SequenceAnalise/PatternSearchingKMP.cs`, and
`SequenceAnalise/NeedlemanWunsch.cs`. The first two were deleted with the CAT removal (view them with
`git show e152aa4:SequenceAnalise/DNA.cs`); `NeedlemanWunsch.cs` is still in the tree. When in doubt,
open it and match what it does.

## Naming checklist

- Types, public methods, public properties → **PascalCase**.
- Private/instance fields → **camelCase**, no `_` or `m_` prefix. (`dnaString`, not
  `_dnaString` or `DnaString`.)
- Locals and parameters → **camelCase**.
- `const` fields (public or private) → **camelCase**, not PascalCase, not SCREAMING_SNAKE_CASE.
  This is the one place this repo deliberately departs from typical .NET style — don't "fix"
  it toward PascalCase.
- Named tuple elements in a public signature → PascalCase (`(double X, double Y)`).
- `this.` — only required to disambiguate a same-named constructor/method parameter from a
  field. Don't chase full consistency beyond that; the baseline files don't either.

## Ordering checklist (per class)

1. `const` / derived `static readonly` fields
2. Instance fields
3. Public properties
4. Constructors
5. Methods, grouped by relationship — keep a private helper next to the one public method
   that calls it rather than sorting all `public` before all `private`. A helper used by
   exactly one caller belongs beside that caller, not in a separate "private methods" block.

Static utility classes with a single real entry point: public method first, private helpers
after it, in call order.

## Comments

- Default to no comment. Add one only to explain *why* (a derivation, an invariant, a
  workaround) — never to narrate *what* the next line does.
- No XML `///` doc comments on core algorithm types in this repo's style (see
  `docs/code-conventions.md` for the one dead-code exception).
- Keep an attribution comment (URL/citation) when a method is ported from a known external
  algorithm reference.
- Remove unused `using` directives.

## Hard boundary — do not cross

This skill governs naming and member ordering **only**. When applying it as a cleanup pass:

- Never change algorithmic behavior, computed values, or optimization-motivated structures
  (rolling DP arrays, sharding math, sparse vectors, tie-break logic, etc.).
- Renaming a `public`/`internal` member that other files reference requires updating every
  call site — grep for the symbol across both `SequenceAnalise/` and `SequenceAnaliseTests/`
  before renaming it, and treat a rename you can't fully trace as out of scope rather than
  guessing.
- If a style rule and a correctness/performance concern conflict, leave the code as-is and
  flag it rather than resolving the conflict unilaterally.
