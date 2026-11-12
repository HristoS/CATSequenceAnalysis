# Web UI layout guide

How every page of the web UI (`wwwroot/app/`) is put together. A new page or view follows these rules; the shared components that
implement them are in `wwwroot/app/html.js`, their styles in `wwwroot/app/styles.css`.

## G1. Page shell

Every page has the same three parts:

1. The top bar.
2. `PageHeader`: the title, a lede of one sentence on what the page is for, and the page's actions (buttons) on the right.
3. The body.

Sizes:
- Desktop only: pages are designed for screens 1280 px wide and up. There is no phone layout, and none is checked.
- The page is at most 1440 px wide (`--page-max`), with a side gutter of 24 px (`--gutter`).
- The top bar lines up with the page.
- The lede is at most 80 characters a line. Anything longer goes in "How to read this" or in a guide panel.

## G2. Five page templates

Every page uses exactly one:

| Template | Body | Pages |
|---|---|---|
| List | filter bar, table, pager | Catalogue, Experiments |
| Detail | summary panels in a grid (`.summarygrid`), then workspaces | Record, Compare, one experiment |
| Workspace | `PlotLayout`: a plot with a 320 px sidebar | Map, Tree |
| Form | `FormLayout`: the form (600 px) with a guide column (what you get, what it accepts, what is kept) | Search, Upload |
| Results | a summary, then the list; a `Workspace` when rows open into a plot | Search results |

```
List                         Detail                          Workspace / Form
┌ header ───────────────┐    ┌ header ───────────────┐       ┌ header ───────────────┐
│ filter bar            │    │ [panel][panel][panel] │       │ ┌ plot/form ┐┌ side ┐ │
│ table                 │    │ ┌ list ─┐┌ views ───┐ │       │ │           ││      │ │
│ pager                 │    │ │       ││ (sticky) │ │       │ └───────────┘└──────┘ │
└───────────────────────┘    └───────────────────────┘       └───────────────────────┘
```

## G3. Controls sit with what they control

- A control that changes one plot or table goes in that panel's toolbar, directly above it.
- A control for the whole page goes in the one bar under the header.
- Every panel starts with a `PanelBar`: its title, its "?", and actions on the whole panel on the right.
- There are no control rows floating between panels.

## G4. Controls keep their natural width

In filter bars (`.filterbar`) and toolbars (`.toolbar`), nothing stretches to fill the row. Widths follow the content:

- text and search inputs: 16 em;
- selects: as wide as their longest option;
- numbers: 8 em;
- a min–max pair: one field (`RangeField`, "Length (bp) [min] – [max]");
- fields that only make sense together (Rank and Taxon): side by side in a `.pair`.

Only a form column (`FormLayout`) stacks its fields at full width. Every control is 34 px tall (`--control-h`), so mixed rows line up.

## G5. Linked views are adjacent

When choosing a row draws or highlights something, the list and what it changes are on screen together (`Workspace`):

- The list is on the left, a narrow chooser (~360 px), and scrolls with the page.
- The views are on the right and stay in view (`position: sticky`).

Linked plots share one card (`PlotPair`):

- One toolbar for the card: View and Drift, shared by both plots (`usePlotSettings`). A View preset turns both; dragging turns one.
- Each plot has a one-line caption (`PlotCaption`). Options that change only that plot and rarely (the genome's source, layout,
  scale, sequence, window size) sit behind a button in its caption (`Popover`), labelled with the current choice.
- The two canvases have the same width and the same top, so they line up. Notes go under a canvas, never between caption and canvas.
- Per-record strips (heat strips) go in one row beneath both plots, side by side.
- A card never scrolls inside itself: plot heights are capped so the card fits the screen.

The selection shows everywhere:

- The picked row is marked (`tr.row.picked`).
- Each view names what it draws, with a Clear button.
- When the picked row is paged or filtered away, a `.showing` line offers "Show it in the list".

Before any row is picked, the views show a `.placeholder` saying what to click.

Current workspaces:

| Page | List | Views |
|---|---|---|
| Record | genes | one card: the Z-curve path and the picked gene, side by side |
| Compare | families and per-record gene pickers | one card: the whole genomes and the picked genes side by side, heat strips beneath |
| Protein search results | hits | the picked hit against the query |

## G6. Text

- **Lede:** one sentence.
- **"?" help (`Help`):** on every control and on every column header whose meaning isn't obvious. Help text says what the control does,
  not how it's built.
- **"How to read this" (`HowTo`):** under a plot, for the longer explanation.
- **Status messages:** loading, empty and error messages (`Message`) stand in the place of the content they replace, never in a banner
  elsewhere.

## G7. Spacing and sizes

- **Spacing:** steps of 4, 8, 12, 16 and 24 px.
- **Panels:** padding 14–16 px, and 16 px between panels.
- **Plots:**
  - never taller than the screen (`.canvaswrap`, `.zpath` max-height);
  - inside a plots card, 4:3 and at most 46% of the screen height (460 px).
- **Pagers:** First, Previous, a page choice, Next and Last, plus the row count (`Pager`).
