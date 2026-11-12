import { get } from "./api.js";
import { Message, Segmented, formatBp, formatGc, html, useEffect, useMemo, useRef, useState } from "./html.js";

const maxTreeLeaves = 500;
const minTreeLeaves = 3;
const gridCell = 12;
const hoverRadiusSquared = 64;
const loaded = new Map();
const groups = [["cellular", "Bacteria and archaea"], ["viral", "Viruses"], ["eukaryotic", "Eukaryotes"], ["all", "Everything"]];

const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const quantile = (values, q) => {
  const sorted = Float64Array.from(values).sort();
  return sorted[Math.min(sorted.length - 1, Math.max(0, Math.round(q * (sorted.length - 1))))];
};

// Every genome's colour class: domains on the map of everything, phyla on the others.
const classesOf = map => (map.group === "all" ? { indexes: map.domain, names: map.domains, shown: 4 } : { indexes: map.phylum, names: map.phyla, shown: 3 });

function topClasses(map) {
  const { indexes, names, shown } = classesOf(map), counts = new Map();
  indexes.forEach(i => { if (i >= 0) counts.set(i, (counts.get(i) || 0) + 1); });
  return [...counts.entries()].sort((a, b) => b[1] - a[1]).slice(0, shown).map(([i, n]) => ({ i, n, name: names[i] }));
}

function taxonIndex(map) {
  const entries = [];
  for (const [rank, list, indexes] of [["phylum", map.phyla, map.phylum], ["family", map.families, map.family], ["genus", map.genera, map.genus]]) {
    const counts = new Array(list.length).fill(0);
    indexes.forEach(k => { if (k >= 0) counts[k]++; });
    list.forEach((name, index) => { if (counts[index] >= 2) entries.push({ rank, index, name, count: counts[index] }); });
  }
  return entries;
}

const indexesOf = (map, rank) => ({ phylum: map.phylum, family: map.family, genus: map.genus })[rank];

// Draws every genome once on the base canvas and returns screen positions plus a hover grid; the marker has its own canvas so hovering never redraws thousands of points.
function drawPoints(canvas, marker, map, highlight) {
  const dpr = window.devicePixelRatio || 1, w = canvas.clientWidth, h = canvas.clientHeight, n = map.x.length;
  for (const c of [canvas, marker]) { c.width = Math.round(w * dpr); c.height = Math.round(h * dpr); c.getContext("2d").setTransform(dpr, 0, 0, dpr, 0, 0); }
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, w, h);

  const pad = 10;
  const x0 = quantile(map.x, 0.002), x1 = quantile(map.x, 0.998), y0 = quantile(map.y, 0.002), y1 = quantile(map.y, 0.998);
  const sx = v => pad + (v - x0) / (x1 - x0 || 1) * (w - 2 * pad), sy = v => h - pad - (v - y0) / (y1 - y0 || 1) * (h - 2 * pad);
  ctx.strokeStyle = css("--rule");
  ctx.lineWidth = 1;
  if (x0 < 0 && x1 > 0) { ctx.beginPath(); ctx.moveTo(sx(0), pad); ctx.lineTo(sx(0), h - pad); ctx.stroke(); }
  if (y0 < 0 && y1 > 0) { ctx.beginPath(); ctx.moveTo(pad, sy(0)); ctx.lineTo(w - pad, sy(0)); ctx.stroke(); }

  const classes = classesOf(map).indexes;
  const colorOf = new Map(topClasses(map).map((t, k) => [t.i, css(`--p${k + 1}`)])), grey = css("--p0"), hi = css("--hi");
  const taxa = highlight && indexesOf(map, highlight.rank);
  const matches = i => taxa[i] === highlight.index;
  const screen = new Array(n);
  const layers = highlight ? [i => !matches(i), matches] : [i => !colorOf.has(classes[i]), i => colorOf.has(classes[i])];
  layers.forEach((include, layer) => {
    for (let i = 0; i < n; i++) {
      if (!include(i)) continue;
      const x = sx(map.x[i]), y = sy(map.y[i]);
      screen[i] = [x, y];
      ctx.globalAlpha = highlight ? (layer ? 0.95 : 0.22) : (layer ? 0.7 : 0.45);
      ctx.fillStyle = highlight ? (layer ? hi : grey) : (layer ? colorOf.get(classes[i]) : grey);
      ctx.beginPath();
      ctx.arc(x, y, highlight ? (layer ? 3.4 : 1.7) : (layer ? 2.2 : 1.8), 0, 7);
      ctx.fill();
    }
  });
  ctx.globalAlpha = 1;

  const grid = new Map();
  screen.forEach((point, i) => {
    const key = `${Math.floor(point[0] / gridCell)},${Math.floor(point[1] / gridCell)}`;
    (grid.get(key) || grid.set(key, []).get(key)).push(i);
  });
  return { screen, grid };
}

function drawMarkers(marker, plot, indexes) {
  const ctx = marker.getContext("2d");
  ctx.clearRect(0, 0, marker.clientWidth, marker.clientHeight);
  for (const i of indexes) {
    if (i < 0 || !plot?.screen[i]) continue;
    const [x, y] = plot.screen[i];
    ctx.strokeStyle = css("--ink"); ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(x, y, 7, 0, 7); ctx.stroke();
    ctx.strokeStyle = css("--panel"); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(x, y, 9, 0, 7); ctx.stroke();
  }
}

function nearest(marker, plot, clientX, clientY) {
  if (!plot) return -1;
  const r = marker.getBoundingClientRect(), x = clientX - r.left, y = clientY - r.top, gx = Math.floor(x / gridCell), gy = Math.floor(y / gridCell);
  let best = -1, bestDistance = hoverRadiusSquared;
  for (let dx = -1; dx <= 1; dx++)
    for (let dy = -1; dy <= 1; dy++)
      for (const i of plot.grid.get(`${gx + dx},${gy + dy}`) || []) {
        const q = plot.screen[i], d = (q[0] - x) ** 2 + (q[1] - y) ** 2;
        if (d < bestDistance) { bestDistance = d; best = i; }
      }
  return best;
}

const taxonName = (list, k) => (k >= 0 ? list[k] : "unclassified");
// Which two of the three principal components the map shows: the server sends x, y and z for PC1, PC2 and PC3.
const pairs = [["12", "PC1 × PC2"], ["13", "PC1 × PC3"], ["23", "PC2 × PC3"]];
const component = (map, c) => [map.x, map.y, map.z ?? map.y][c];
const shownAs = (map, pair) => ({ ...map, x: component(map, Number(pair[0]) - 1), y: component(map, Number(pair[1]) - 1) });
const axisNote = (axis, i) =>
  `PC${i + 1} (${Math.round(axis.explainedVariance * 100)} %): ${axis.topLoadings.map(l => `${l.dinucleotide} ${l.weight < 0 ? "−" : "+"}`).join(", ")}`;

export function MapPage({ params }) {
  const [group, setGroup] = useState(groups.some(([g]) => g === params.get("group")) ? params.get("group") : "cellular");
  const [map, setMap] = useState(null);
  const [error, setError] = useState(null);
  const [highlight, setHighlight] = useState(null);
  const [selected, setSelected] = useState(-1);
  const [text, setText] = useState("");
  const [notice, setNotice] = useState(null);
  const [pair, setPair] = useState(pairs.some(([p]) => p === params.get("axes")) ? params.get("axes") : "12");
  const base = useRef(null), marker = useRef(null), tip = useRef(null), plot = useRef(null), hovered = useRef(-1), selectedRef = useRef(-1);
  const wanted = params.get("id");

  selectedRef.current = selected;

  useEffect(() => {
    let live = true;
    setMap(null); setError(null); setHighlight(null); setSelected(-1); setNotice(null);
    const pending = loaded.get(group) ?? get(`/api/landscape?group=${group}`);
    Promise.resolve(pending).then(result => {
      if (!live) return;
      loaded.set(group, result);
      setMap(result);
      if (wanted) {
        const at = result.ids.indexOf(Number(wanted));
        if (at >= 0) setSelected(at);
        else setNotice("That genome isn't on this map: the map holds one genome per species, with a GTDB assignment for bacteria and archaea.");
      }
    }).catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [group, wanted]);

  useEffect(() => {
    if (!map) return undefined;
    const shown = shownAs(map, pair);
    const redraw = () => {
      plot.current = drawPoints(base.current, marker.current, shown, highlight);
      drawMarkers(marker.current, plot.current, [hovered.current, selectedRef.current]);
    };
    redraw();
    const resize = new ResizeObserver(redraw);
    resize.observe(base.current);
    const theme = window.matchMedia("(prefers-color-scheme: dark)");
    theme.addEventListener("change", redraw);
    return () => { resize.disconnect(); theme.removeEventListener("change", redraw); };
  }, [map, highlight, pair]);

  useEffect(() => {
    if (map) drawMarkers(marker.current, plot.current, [hovered.current, selected]);
  }, [map, selected]);

  const taxa = useMemo(() => (map ? taxonIndex(map) : []), [map]);
  const highlightedIds = useMemo(() => {
    if (!map || !highlight) return [];
    const indexes = indexesOf(map, highlight.rank);
    return map.ids.filter((_, i) => indexes[i] === highlight.index);
  }, [map, highlight]);

  function hover(event) {
    const i = nearest(marker.current, plot.current, event.clientX, event.clientY);
    hovered.current = i;
    drawMarkers(marker.current, plot.current, [i, selected]);
    const box = tip.current;
    if (i < 0) { box.hidden = true; return; }
    box.hidden = false;
    box.firstChild.textContent = map.names[i] || "(no name)";
    box.lastChild.textContent = `${taxonName(map.genera, map.genus[i])} · ${taxonName(map.families, map.family[i])}\n${taxonName(map.phyla, map.phylum[i])}\nGC ${formatGc(map.gc[i])} · ${formatBp(map.lengthBp[i])}`;
    box.style.left = `${Math.max(8, Math.min(event.clientX + 14, window.innerWidth - 310))}px`;
    box.style.top = `${Math.min(event.clientY + 14, window.innerHeight - 120)}px`;
  }

  function leave() {
    hovered.current = -1;
    tip.current.hidden = true;
    drawMarkers(marker.current, plot.current, [selected]);
  }

  function find(event) {
    event.preventDefault();
    const wantedText = text.trim().toLowerCase();
    if (!wantedText) return;
    const found = taxa.find(t => t.name.toLowerCase() === wantedText) || taxa.find(t => t.name.toLowerCase().startsWith(wantedText));
    setNotice(found ? null : `No phylum, family or genus with two or more genomes matches "${text.trim()}".`);
    if (found) setHighlight(found);
  }

  const tops = map ? topClasses(map) : [];
  const treeHref = `#/tree?ids=${highlightedIds.join(",")}`;
  const treeable = highlightedIds.length >= minTreeLeaves && highlightedIds.length <= maxTreeLeaves;

  return html`<h1>Map of the catalogue</h1>
    <p class="lede">${group === "all" ? "Every genome" : "Every species once"}, placed by two of the first three principal components of its dinucleotide signature (log ρ*). This is a picture of composition:
      two axes keep phyla apart but not genera, so use a genome's relatives, not its neighbours here, to find what it resembles.</p>
    <div class="controls">
      <${Segmented} label="Genomes" value=${group} onChange=${setGroup} options=${groups} />
      ${map?.axes.length > 2 && html`<${Segmented} label="Axes" value=${pair} onChange=${setPair} options=${pairs} />`}
      <form class="find field" onSubmit=${find}>
        <label htmlFor="taxon">Highlight a phylum, family or genus</label>
        <div>
          <input id="taxon" list="taxa" autoComplete="off" placeholder="e.g. Streptomyces" value=${text} onInput=${e => setText(e.target.value)} />
          <button type="submit">Show</button>
          <button type="button" onClick=${() => { setText(""); setHighlight(null); setNotice(null); }}>Clear</button>
        </div>
        <datalist id="taxa">${taxa.map(t => html`<option key=${`${t.rank}${t.index}`} value=${t.name}>${t.rank}, ${t.count} genomes</option>`)}</datalist>
      </form>
    </div>
    ${error && html`<${Message} error>${error}<//>`}
    ${notice && html`<${Message}>${notice}<//>`}
    ${!map && !error && html`<${Message}>Loading the map…<//>`}
    ${map && html`<p class="status">${map.ids.length.toLocaleString()} genomes.${highlight ? ` ${highlight.name} (${highlight.rank}): ${highlightedIds.length.toLocaleString()} genomes.` : group === "all" ? " Colours show domains." : " Colours show the three largest phyla."}
      ${highlight && (treeable ? html` <a href=${treeHref}>Tree of these genomes</a>` : ` A tree needs ${minTreeLeaves} to ${maxTreeLeaves} genomes.`)}</p>
    <div class="mapgrid">
      <figure class="canvaswrap">
        <canvas ref=${base}></canvas>
        <canvas ref=${marker} class="over" role="img" aria-label="Map of genomes" onMouseMove=${hover} onMouseLeave=${leave}
          onClick=${e => { const i = nearest(marker.current, plot.current, e.clientX, e.clientY); if (i >= 0) setSelected(i); }}></canvas>
      </figure>
      <aside class="stack">
        ${selected >= 0 ? html`<section class="panel">
          <h2>${map.names[selected] || "(no name)"}</h2>
          <p>${taxonName(map.genera, map.genus[selected])} · ${taxonName(map.families, map.family[selected])} · ${taxonName(map.phyla, map.phylum[selected])}</p>
          <p>GC ${formatGc(map.gc[selected])} · ${formatBp(map.lengthBp[selected])}</p>
          <p class="actions">
            <a class="button" href="#/sequences/${map.ids[selected]}">Open record</a>
            ${map.genus[selected] >= 0 && html`<button type="button"
              onClick=${() => setHighlight({ rank: "genus", index: map.genus[selected], name: map.genera[map.genus[selected]] })}>Highlight genus</button>`}
          </p>
        </section>` : html`<p class="hint">Hover a point to see the genome; click to select it.</p>`}
        <div class="legend">
          ${highlight
            ? html`<span><i class="swatch hi"></i>${highlight.name}</span><span><i class="swatch p0"></i>all other genomes</span>`
            : html`${tops.map((t, k) => html`<span key=${t.i}><i class="swatch p${k + 1}"></i>${t.name} (${t.n.toLocaleString()})</span>`)}
                ${group !== "all" && html`<span><i class="swatch p0"></i>other phyla${group === "viral" ? " or none" : ""}</span>`}`}
        </div>
        <p class="hint">${map.axes.map(axisNote).join(" · ")}. Shown: ${pairs.find(([p]) => p === pair)[1]}, horizontal × vertical.</p>
        <p class="hint">Each axis lists its heaviest dinucleotides; + and − show which way they pull.</p>
      </aside>
    </div>`}
    <div class="tip" ref=${tip} hidden><b></b><span></span></div>`;
}
