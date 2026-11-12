// React without a build step: the UMD builds in vendor/ set window.React, and htm stands in for JSX.
import htm from "../vendor/htm.mjs";

export const { useEffect, useMemo, useRef, useState } = window.React;
export const html = htm.bind(window.React.createElement);

export const formatBp = bp => (bp >= 1e6 ? `${(bp / 1e6).toFixed(2)} Mbp` : bp >= 1e3 ? `${(bp / 1e3).toFixed(1)} kbp` : `${bp} bp`);
export const formatGc = fraction => (fraction == null ? "–" : `${(fraction * 100).toFixed(1)} %`);

// The page shell every page follows: a header (title, one line on what the page is for, actions on the right), then one of three
// bodies: data panels stacked at full width, FormLayout (a form and its guide), or PlotLayout (a plot and its sidebar).
export function PageHeader({ title, lede = null, actions = null }) {
  return html`<header class="pagehead">
    <div class="pagetitle"><h1>${title}</h1>${lede && html`<p class="lede">${lede}</p>`}</div>
    ${actions && html`<div class="pageactions">${actions}</div>`}
  </header>`;
}

// A form beside its guide: what the form does, its limits and what is kept. One column on narrow screens, the guide below.
export function FormLayout({ form, guide }) {
  return html`<div class="formlayout"><div class="formcol">${form}</div><aside class="guide">${guide}</aside></div>`;
}

// A plot beside its sidebar of controls, selection and legend; the sidebar stays in view while the page scrolls.
export function PlotLayout({ plot, side }) {
  return html`<div class="plotlayout"><div class="plotcol">${plot}</div><aside class="sidecol">${side}</aside></div>`;
}

// A list beside the views it drives (docs/ui-layout.md G5): the list scrolls with the page, the views stay in view beside it, so choosing a
// row never needs a scroll to see what it changed. split "views" gives the views two thirds (a gene list beside plots), "even" half each
// (a wide hit table beside one plot). One column below 1100 px, views first.
export function Workspace({ list, views, split = "views" }) {
  return html`<div class=${`workspace ${split}`}><div class="wslist">${list}</div><div class="wsviews">${views}</div></div>`;
}

// A panel's title row: the title, its "?", and anything that acts on the whole panel on the right.
export function PanelBar({ title, help = null, children = null }) {
  return html`<div class="panelbar"><h2>${title}${help && html`<${Help} label=${`About ${title}`}>${help}<//>`}</h2>
    ${children && html`<div class="panelactions">${children}</div>`}</div>`;
}

// The View and Drift that a page's linked plots share (UI-LAYOUT.md G5). viewSet is a fresh object on every pick, so a plot dragged into
// a free angle returns to the preset even when the same preset is picked again; toolbar is the two switches for the plots' card.
export function usePlotSettings(defaultDrift = "removed") {
  const [viewSet, setViewSet] = useState({ view: "3d", n: 0 });
  const [drift, setDrift] = useState(defaultDrift);
  const pickView = view => setViewSet(old => ({ view, n: old.n + 1 }));
  const toolbar = html`<${Segmented} label="View" value=${viewSet.view} onChange=${pickView} help=${`${helpText.view} Sets both plots; dragging turns one.`}
      options=${[["3d", "3-D"], ["xy", "x–y"], ["xz", "x–z"], ["yz", "y–z"]]} />
    <${Segmented} label="Drift" value=${drift} onChange=${setDrift} help=${`${helpText.drift} Kept, a curve's overall direction is its composition, so genes of different GC point different ways; ` +
      "Removed compares the shape of the bends. On a single gene, read the drift-removed bends as a sketch."}
      options=${[["removed", "Removed"], ["kept", "Kept"]]} />`;
  return { viewSet, drift, toolbar };
}

// A button that opens a small panel of options used less often; it says the current choice in its label and closes on Escape or a click
// outside it.
export function Popover({ label, title, children }) {
  const [open, setOpen] = useState(false);
  const box = useRef(null);
  useEffect(() => {
    if (!open) return undefined;
    const away = e => { if (box.current && !box.current.contains(e.target)) setOpen(false); };
    const escape = e => { if (e.key === "Escape") setOpen(false); };
    document.addEventListener("pointerdown", away);
    document.addEventListener("keydown", escape);
    return () => { document.removeEventListener("pointerdown", away); document.removeEventListener("keydown", escape); };
  }, [open]);
  return html`<div class="popover" ref=${box}>
    <button type="button" class="popbutton" aria-expanded=${open} title=${title} onClick=${() => setOpen(!open)}>${label}<span aria-hidden="true">▾</span></button>
    ${open && html`<div class="poppanel" role="dialog" aria-label=${title}>${children}</div>`}
  </div>`;
}

// A plot's title line inside a PlotPair: one line of fixed height, so the canvases under the captions start at the same place.
export function PlotCaption({ title, children = null }) {
  return html`<div class="plotcap"><h3 title=${typeof title === "string" ? title : undefined}>${title}</h3>${children}</div>`;
}

// Linked plots in one card (UI-LAYOUT.md G5): the shared toolbar, two plot slots of equal width side by side (each a PlotCaption, then
// its canvas, then its notes), and an optional row of strips beneath both. Without right, the one plot takes the whole width.
export function PlotPair({ title = "Plots", toolbar, left, right = null, strips = null }) {
  return html`<section class="panel plotpair">
    <${PanelBar} title=${title}>${toolbar}<//>
    <div class=${right ? "plotrow" : "plotrow single"}><div class="plotslot">${left}</div>${right && html`<div class="plotslot">${right}</div>`}</div>
    ${strips && html`<div class="striprow">${strips}</div>`}
    <p class="hint">Drag a plot to rotate it, scroll to zoom, shift-drag to pan, double-click to reset; hover a point for what it is.</p>
  </section>`;
}

// A min–max pair as one field (G4): "Length (bp) [ min ] – [ max ]".
export function RangeField({ label, min, max, onMin, onMax, help = null, step = null }) {
  return html`<div class="field"><span>${label}${help && html`<${Help}>${help}<//>`}</span>
    <div class="range">
      <input type="number" min="0" step=${step} value=${min} onInput=${e => onMin(e.target.value)} placeholder="min" aria-label=${`${label}, at least`} />
      <span aria-hidden="true">–</span>
      <input type="number" min="0" step=${step} value=${max} onInput=${e => onMax(e.target.value)} placeholder="max" aria-label=${`${label}, at most`} />
    </div></div>`;
}

export function Message({ error, children }) {
  return html`<p class=${error ? "message error" : "message"} role=${error ? "alert" : "status"}>${children}</p>`;
}

// Paging and filtering of rows already loaded in the page, with no server call: `text(row)` is what the filter box matches (case
// ignored), `size` rows per page. A change in the row count (new data) or the filter goes back to the first page. locate(predicate)
// turns to the page holding the first row that satisfies it, clearing the filter when the filter hides it.
export function usePaged(rows, { size = 50, text = null } = {}) {
  const [page, setPage] = useState(0);
  const [needle, setNeedle] = useState("");
  useEffect(() => setPage(0), [rows.length]);
  const wanted = needle.trim().toLowerCase();
  const filtered = useMemo(() => (wanted && text ? rows.filter(r => (text(r) ?? "").toLowerCase().includes(wanted)) : rows), [rows, wanted]);
  const pages = Math.max(1, Math.ceil(filtered.length / size)), current = Math.min(page, pages - 1);
  return {
    rows: filtered.slice(current * size, (current + 1) * size), filtered, total: rows.length, needle, page: current, pages, setPage,
    setNeedle: value => { setNeedle(value); setPage(0); },
    locate: predicate => {
      let at = filtered.findIndex(predicate);
      if (at < 0) {
        at = rows.findIndex(predicate);
        if (at < 0) return;
        setNeedle("");
      }
      setPage(Math.floor(at / size));
    },
  };
}

// A page number to jump to: a list of pages while that stays short, else a number box applied on Enter or on leaving it.
const listedPages = 200;
function PageChoice({ page, pages, onChange }) {
  const [draft, setDraft] = useState(String(page + 1));
  useEffect(() => setDraft(String(page + 1)), [page]);
  if (pages <= listedPages)
    return html`<select class="pagechoice" aria-label="Page" value=${page} onChange=${e => onChange(Number(e.target.value))}>
      ${Array.from({ length: pages }, (_, i) => html`<option key=${i} value=${i}>${i + 1}</option>`)}
    </select>`;
  const commit = () => {
    const n = Math.round(Number(draft));
    if (Number.isFinite(n) && n >= 1 && n <= pages) onChange(n - 1); else setDraft(String(page + 1));
  };
  return html`<input class="pagechoice" type="number" min="1" max=${pages} aria-label="Page" value=${draft}
    onInput=${e => setDraft(e.target.value)} onBlur=${commit} onKeyDown=${e => { if (e.key === "Enter") { e.preventDefault(); commit(); } }} />`;
}

// First / Previous / page choice / Next / Last for usePaged, with how many rows match; the buttons only when there's more than one page.
export function Pager({ paged, noun = "rows" }) {
  if (paged.total === 0) return null;
  const count = paged.filtered.length === paged.total ? paged.total.toLocaleString() : `${paged.filtered.length.toLocaleString()} of ${paged.total.toLocaleString()}`;
  const first = paged.page === 0, last = paged.page >= paged.pages - 1;
  return html`<div class="pager">
    ${paged.pages > 1 && html`
      <button type="button" disabled=${first} onClick=${() => paged.setPage(0)} title="First page" aria-label="First page">«</button>
      <button type="button" disabled=${first} onClick=${() => paged.setPage(paged.page - 1)}>‹ Previous</button>
      <span class="pagelabel">Page <${PageChoice} page=${paged.page} pages=${paged.pages} onChange=${paged.setPage} /> of ${paged.pages.toLocaleString()}</span>
      <button type="button" disabled=${last} onClick=${() => paged.setPage(paged.page + 1)}>Next ›</button>
      <button type="button" disabled=${last} onClick=${() => paged.setPage(paged.pages - 1)} title="Last page" aria-label="Last page">»</button>`}
    <span class="hint">${count} ${noun}</span>
  </div>`;
}

// The filter box for usePaged: it filters only what's already loaded.
export function FilterBox({ paged, label = "Filter", placeholder = "", help = null }) {
  return html`<label class="field filterbox"><span>${label}${help && html`<${Help}>${help}<//>`}</span>
    <input type="search" value=${paged.needle} placeholder=${placeholder} onInput=${e => paged.setNeedle(e.target.value)} /></label>`;
}

// A "?" beside a label or column header: its text shows on hover, on keyboard focus, or on a tap. It's a span, not a button, so a <label>
// around it still names its own control, and its click is cancelled so the label doesn't take it. The popover is fixed to the viewport,
// so a scrolling table around it doesn't clip it.
export function Help({ children, label = "What is this?" }) {
  const [at, setAt] = useState(null);
  const show = e => {
    const r = e.currentTarget.getBoundingClientRect(), above = r.bottom + 180 > window.innerHeight;
    setAt({ left: Math.max(8, Math.min(r.left - 8, window.innerWidth - 336)), top: above ? r.top - 6 : r.bottom + 6, above });
  };
  const hide = () => setAt(null);
  // Fixed to where it opened, so it closes when anything scrolls rather than drift away from its "?".
  useEffect(() => {
    if (!at) return undefined;
    window.addEventListener("scroll", hide, { capture: true, passive: true });
    return () => window.removeEventListener("scroll", hide, { capture: true });
  }, [!!at]);
  return html`<span class="help" role="button" tabIndex="0" aria-label=${label} aria-expanded=${!!at}
      onMouseEnter=${show} onMouseLeave=${hide} onFocus=${show} onBlur=${hide}
      onClick=${e => { e.preventDefault(); e.stopPropagation(); if (at) hide(); else show(e); }}
      onKeyDown=${e => { if (e.key === "Escape") hide(); }}>?${at && html`<span class=${at.above ? "tip helptip above" : "tip helptip"} role="tooltip"
      style=${{ left: `${at.left}px`, top: `${at.top}px` }}>${children}</span>`}</span>`;
}

// The longer explanation of a plot or table, folded away under its own summary line.
export function HowTo({ title = "How to read this", children }) {
  return html`<details class="howto"><summary>${title}</summary><div class="hint">${children}</div></details>`;
}

// options are [value, text] or [value, text, tooltip]; help is the "?" text for the whole control.
export function Segmented({ label, value, options, onChange, help = null }) {
  return html`<div class="field">
    <span>${label}${help && html`<${Help} label=${`About ${label}`}>${help}<//>`}</span>
    <div class="segmented" role="group" aria-label=${label}>
      ${options.map(([v, text, tip]) => html`<button key=${v} type="button" title=${tip} aria-pressed=${v === value} onClick=${() => onChange(v)}>${text}</button>`)}
    </div>
  </div>`;
}

// The kinds a record can be, as the upload page offers them and the catalogue filters them (value, label).
export const recordKinds = [["bacterium", "Bacteria"], ["archaeon", "Archaea"], ["virus", "Viruses"], ["plasmid", "Plasmids"], ["eukaryote", "Eukaryotes"]];

// The texts of the controls that recur on several pages, so each reads the same everywhere.
export const helpText = {
  drift: "A sequence's overall composition makes its path lean steadily one way (each axis's best-fit line). Removed subtracts that lean, so " +
    "the local bends stand out; Kept draws the raw running sums.",
  view: "3-D can be dragged to rotate. x–y, x–z and y–z look straight down the third axis. x is purine − pyrimidine, y amino − keto, " +
    "z strong − weak (GC).",
  scale: "Same length divides each path by its genome's length, so a phage and a bacterium can share the axes. Bases keeps the raw counts.",
  layout: "Where each path starts. From origin: circular chromosomes start at their replication origin, which lines relatives up. Aligned to " +
    "first: each path is started (and turned onto the other strand if that fits better) where it best follows the first record. As " +
    "deposited: wherever each record's file starts.",
  source: "Windows: running sums of the stored genome windows. Genes: a path built from the stored genes alone (codon positions 1 and 2), " +
    "leaving out the DNA between genes.",
  geneDrift: "A gene's coding composition makes its curve lean steadily one way. Removed subtracts that lean to show the shape along the " +
    "gene; read it as a sketch, since the redrawn curve is only exact at window ends.",
  signature: "log ρ*: each dinucleotide's frequency against what its two bases' frequencies predict; 0 is as expected, below 0 avoided. " +
    "Z-curve coefficients: the same 16 numbers rotated into Z-curve axes. Gene layer: the mean gene window over codon positions 1 and 2.",
};
