// React without a build step: the UMD builds in vendor/ set window.React, and htm stands in for JSX.
import htm from "../vendor/htm.mjs";

export const { useEffect, useMemo, useRef, useState } = window.React;
export const html = htm.bind(window.React.createElement);

export const formatBp = bp => (bp >= 1e6 ? `${(bp / 1e6).toFixed(2)} Mbp` : bp >= 1e3 ? `${(bp / 1e3).toFixed(1)} kbp` : `${bp} bp`);
export const formatGc = fraction => (fraction == null ? "–" : `${(fraction * 100).toFixed(1)} %`);

export function Message({ error, children }) {
  return html`<p class=${error ? "message error" : "message"} role=${error ? "alert" : "status"}>${children}</p>`;
}

export function Segmented({ label, value, options, onChange }) {
  return html`<div class="field">
    <span>${label}</span>
    <div class="segmented" role="group" aria-label=${label}>
      ${options.map(([v, text]) => html`<button key=${v} type="button" aria-pressed=${v === value} onClick=${() => onChange(v)}>${text}</button>`)}
    </div>
  </div>`;
}
