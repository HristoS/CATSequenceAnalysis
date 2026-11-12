// React without a build step: the UMD builds in vendor/ set window.React, and htm stands in for JSX.
import htm from "../vendor/htm.mjs";

export const { useEffect, useMemo, useRef, useState } = window.React;
export const html = htm.bind(window.React.createElement);

export const formatBp = bp => (bp >= 1e6 ? `${(bp / 1e6).toFixed(2)} Mbp` : bp >= 1e3 ? `${(bp / 1e3).toFixed(1)} kbp` : `${bp} bp`);
export const formatGc = fraction => (fraction == null ? "–" : `${(fraction * 100).toFixed(1)} %`);

export function Message({ error, children }) {
  return html`<p class=${error ? "message error" : "message"} role=${error ? "alert" : "status"}>${children}</p>`;
}

// Paging and filtering of rows already loaded in the page, with no server call: `text(row)` is what the filter box matches (case
// ignored), `size` rows per page. A change in the row count (new data) or the filter goes back to the first page.
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
  };
}

// Previous / Next for usePaged, with how many rows match and which page this is; the buttons only when there's more than one page.
export function Pager({ paged, noun = "rows" }) {
  if (paged.total === 0) return null;
  const count = paged.filtered.length === paged.total ? paged.total.toLocaleString() : `${paged.filtered.length.toLocaleString()} of ${paged.total.toLocaleString()}`;
  return html`<div class="pager">
    ${paged.pages > 1 && html`<button type="button" disabled=${paged.page === 0} onClick=${() => paged.setPage(paged.page - 1)}>Previous</button>`}
    <span class="hint">${count} ${noun}${paged.pages > 1 ? ` · page ${paged.page + 1} of ${paged.pages}` : ""}</span>
    ${paged.pages > 1 && html`<button type="button" disabled=${paged.page >= paged.pages - 1} onClick=${() => paged.setPage(paged.page + 1)}>Next</button>`}
  </div>`;
}

// The filter box for usePaged: it filters only what's already loaded.
export function FilterBox({ paged, label = "Filter", placeholder = "" }) {
  return html`<label class="field filterbox"><span>${label}</span>
    <input type="search" value=${paged.needle} placeholder=${placeholder} onInput=${e => paged.setNeedle(e.target.value)} /></label>`;
}

export function Segmented({ label, value, options, onChange }) {
  return html`<div class="field">
    <span>${label}</span>
    <div class="segmented" role="group" aria-label=${label}>
      ${options.map(([v, text]) => html`<button key=${v} type="button" aria-pressed=${v === value} onClick=${() => onChange(v)}>${text}</button>`)}
    </div>
  </div>`;
}
