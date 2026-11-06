import { get } from "./api.js";
import { Message, Segmented, formatBp, html, useEffect, useMemo, useRef, useState } from "./html.js";
import { drawCurve, views } from "./zcurve3d.js";

const rowsPerPage = 25;
const samplesPerWindow = 24;
const featureNames = ["x mean", "y mean", "z mean", "x trend", "y trend", "z trend", "x curvature", "y curvature", "z curvature"];

// One request per record however many panels ask: the gene list and the signature's gene-layer view share it.
const lists = new Map();
export function loadGenes(id) {
  if (!lists.has(id)) lists.set(id, get(`/api/sequences/${id}/genes`).catch(e => { lists.delete(id); throw e; }));
  return lists.get(id);
}

function useGenes(id) {
  const [state, setState] = useState({ list: null, error: null });
  useEffect(() => {
    let live = true;
    setState({ list: null, error: null });
    loadGenes(id).then(list => { if (live) setState({ list, error: null }); }).catch(e => { if (live) setState({ list: null, error: e.message }); });
    return () => { live = false; };
  }, [id]);
  return state;
}

// The record's mean gene window as signature bars: the third view of the signature panel.
export function GeneLayerSignature({ id }) {
  const { list, error } = useGenes(id);
  if (error) return html`<${Message}>${error}<//>`;
  if (!list) return html`<${Message}>Loading the gene layer…<//>`;
  if (list.geneCount === 0) return html`<${Message}>No genes are stored for this record: the gene layer covers BV-BRC genomes only.<//>`;
  const peak = Math.max(...list.meanWindow.map(v => Math.abs(v.value)), 1e-9);
  return html`<div class="bars">${list.meanWindow.map(v => html`<div class="bar" key=${v.name} title=${v.meaning}>
      <span class="mono">${v.name}</span>
      <div class="track"><i class=${v.value < 0 ? "negative" : "positive"} style=${{ width: `${(Math.abs(v.value) / peak) * 50}%` }}></i></div>
      <span class="mono num">${v.value.toFixed(3)}</span>
    </div>`)}</div>
    <p class="hint">The mean window of this record's ${list.geneCount.toLocaleString()} genes (${list.windowCount.toLocaleString()} windows): each axis over
      codon positions 1 and 2, the bases a protein fixes. Means are coding composition; trend and curvature say how it changes along a window,
      which is what lets a protein find its gene.</p>`;
}

export function GenesPanel({ id }) {
  const { list, error } = useGenes(id);
  const [text, setText] = useState("");
  const [page, setPage] = useState(0);
  const [selected, setSelected] = useState(null);
  useEffect(() => { setSelected(null); setPage(0); setText(""); }, [id]);

  const shown = useMemo(() => {
    if (!list) return [];
    const needle = text.trim().toLowerCase();
    return needle ? list.genes.filter(g => (g.product ?? "").toLowerCase().includes(needle) || (g.familyId ?? "").toLowerCase().includes(needle)) : list.genes;
  }, [list, text]);

  if (error) return html`<${Message}>${error}<//>`;
  if (!list) return html`<${Message}>Loading the genes…<//>`;
  if (list.geneCount === 0) return html`<${Message}>No genes are stored for this record: the gene layer covers BV-BRC genomes, genes of 50 aa and up.<//>`;
  const pages = Math.max(1, Math.ceil(shown.length / rowsPerPage)), current = Math.min(page, pages - 1);

  return html`<label class="field"><span>Product or family contains</span>
      <input value=${text} onInput=${e => { setText(e.target.value); setPage(0); }} placeholder="e.g. gyrase" /></label>
    <div class="tablewrap"><table>
      <thead><tr><th class="num">#</th><th class="num">Position</th><th>Strand</th><th class="num">aa</th><th>Product</th></tr></thead>
      <tbody>${shown.slice(current * rowsPerPage, (current + 1) * rowsPerPage).map(g => html`<tr key=${g.index}
          class=${selected === g.index ? "row picked" : "row"} tabIndex="0" onClick=${() => setSelected(g.index)}
          onKeyDown=${e => { if (e.key === "Enter") setSelected(g.index); }}>
        <td class="num">${g.index + 1}</td>
        <td class="num">${formatBp(g.startBp)}</td>
        <td>${g.strand === "Reverse" ? "−" : "+"}</td>
        <td class="num">${g.residues.toLocaleString()}</td>
        <td title=${g.familyId ?? ""}>${g.product ?? g.familyId ?? "–"}</td>
      </tr>`)}</tbody>
    </table></div>
    ${shown.length === 0 && html`<${Message}>No gene matches.<//>`}
    <div class="pager">
      <button type="button" disabled=${current === 0} onClick=${() => setPage(current - 1)}>Previous</button>
      <span>${shown.length.toLocaleString()} of ${list.geneCount.toLocaleString()} genes · page ${current + 1} of ${pages}</span>
      <button type="button" disabled=${current >= pages - 1} onClick=${() => setPage(current + 1)}>Next</button>
    </div>
    ${selected == null ? html`<p class="hint">Pick a gene to draw its Z-curve.</p>` : html`<${GeneView} id=${id} index=${selected} />`}`;
}

// The running sums inside window k of W positions at t in [-1, 1], from its scaled Legendre moments m0, m1, m2 (one axis):
// W / 2 * (m0 (t + 1) + sqrt(3) m1 (t^2 - 1) / 2 + sqrt(5) m2 (t^3 - t) / 2). Exact at the window's ends.
function genePath(windows, positions) {
  const points = [[0, 0, 0]], at = [0], marks = [];
  for (const m of windows) {
    const start = points[points.length - 1], from = at[at.length - 1];
    for (let s = 1; s <= samplesPerWindow; s++) {
      const t = -1 + 2 * s / samplesPerWindow;
      points.push(start.map((v, a) => v + positions / 2 *
        (m[a] * (t + 1) + Math.sqrt(3) * m[3 + a] * (t * t - 1) / 2 + Math.sqrt(5) * m[6 + a] * (t * t * t - t) / 2)));
      at.push(from + positions * s / samplesPerWindow);
    }
    marks.push(points.length - 1);
  }
  // Each axis's least-squares line against position: the gene's overall coding composition, which "drift removed" subtracts.
  const n = at.length, meanAt = at.reduce((a, b) => a + b, 0) / n, spread = at.reduce((a, p) => a + (p - meanAt) ** 2, 0) || 1;
  const fits = [0, 1, 2].map(a => {
    const mean = points.reduce((s, p) => s + p[a], 0) / n;
    const slope = points.reduce((s, p, i) => s + (at[i] - meanAt) * (p[a] - mean), 0) / spread;
    return [mean - slope * meanAt, slope];
  });
  return { points, at, marks, fits };
}

function GeneView({ id, index }) {
  const [gene, setGene] = useState(null);
  const [error, setError] = useState(null);
  // Kept by default: removing the drift leaves a shape about as large as the rebuilt curve's error (gene_curve_check.cs).
  const [drift, setDrift] = useState("kept");
  const [view, setView] = useState("3d");
  const canvas = useRef(null), tip = useRef(null), screen = useRef(null), angles = useRef(views["3d"]), drag = useRef(null);

  useEffect(() => {
    let live = true;
    setGene(null); setError(null);
    get(`/api/sequences/${id}/genes/${index}`).then(g => { if (live) setGene(g); }).catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [id, index]);
  const shape = useMemo(() => (gene ? genePath(gene.windows, gene.windowPositions) : null), [gene]);

  const redraw = () => {
    if (!shape || !canvas.current) return;
    const points = drift === "kept" ? shape.points : shape.points.map((p, i) => p.map((v, a) => v - (shape.fits[a][0] + shape.fits[a][1] * shape.at[i])));
    screen.current = drawCurve(canvas.current, points, { yaw: angles.current[0], pitch: angles.current[1], marks: shape.marks }).screen;
  };
  useEffect(() => { if (view !== "free") angles.current = views[view]; redraw(); }, [shape, drift, view]);
  useEffect(() => {
    if (!canvas.current) return;
    const observer = new ResizeObserver(redraw);
    observer.observe(canvas.current);
    return () => observer.disconnect();
  });

  const down = e => { drag.current = { x: e.clientX, y: e.clientY, angles: angles.current }; e.currentTarget.setPointerCapture(e.pointerId); tip.current.hidden = true; };
  const up = () => { drag.current = null; };
  const move = e => {
    if (drag.current) {
      const d = drag.current;
      angles.current = [d.angles[0] + (e.clientX - d.x) * 0.01, Math.max(-1.55, Math.min(1.55, d.angles[1] + (e.clientY - d.y) * 0.01))];
      if (view !== "free") setView("free"); else redraw();
      return;
    }
    const r = canvas.current.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
    let best = -1, bestDistance = 14 * 14;
    screen.current?.forEach((p, k) => { const d = (p[0] - mx) ** 2 + (p[1] - my) ** 2; if (d < bestDistance) { bestDistance = d; best = k; } });
    if (best < 0) { tip.current.hidden = true; return; }
    tip.current.textContent = `about amino acid ${Math.round(shape.at[best] / 2).toLocaleString()} · window ${Math.max(1, Math.ceil(best / samplesPerWindow))}`;
    tip.current.hidden = false;
    tip.current.style.left = `${Math.min(screen.current[best][0] + 12, r.width - tip.current.offsetWidth - 6)}px`;
    tip.current.style.top = `${Math.max(6, screen.current[best][1] - 14)}px`;
  };

  if (error) return html`<${Message}>${error}<//>`;
  if (!gene) return html`<${Message}>Loading the gene's windows…<//>`;
  const g = gene.gene;
  const peak = Math.max(...gene.windows.flat().map(Math.abs), 1e-9);
  return html`<h3 class="genetitle">${g.product ?? "Gene"} <span class="mono">${g.familyId ?? ""}</span></h3>
    <p class="hint">${`${formatBp(g.startBp)}–${formatBp(g.endBp)} on the ${g.strand === "Reverse" ? "reverse" : "forward"} strand · ` +
      `${g.residues.toLocaleString()} aa · ${g.windows} windows of ${gene.windowPositions} codon positions 1 and 2`}</p>
    <div class="toolbar">
      <${Segmented} label="Drift" value=${drift} onChange=${setDrift} options=${[["removed", "Removed"], ["kept", "Kept"]]} />
      <${Segmented} label="View" value=${view} onChange=${setView} options=${[["3d", "3-D"], ["xy", "x–y"], ["xz", "x–z"], ["yz", "y–z"]]} />
    </div>
    <div class="zpath">
      <canvas ref=${canvas} role="img" aria-label=${`Z-curve of gene ${index + 1} over codon positions 1 and 2; drag to rotate`}
        onPointerDown=${down} onPointerMove=${move} onPointerUp=${up} onPointerCancel=${up}
        onPointerLeave=${() => { if (!drag.current) tip.current.hidden = true; }}></canvas>
      <div class="ztip" ref=${tip} hidden></div>
    </div>
    <p class="hint">The gene's Z-curve over codon positions 1 and 2, redrawn from its stored windows alone (no bases): each window's mean,
      trend and curvature give a smooth piece that ends exactly where the true curve does; small rings mark the window ends. Between them
      it's an approximation: on E. coli genes the worst gap is about 11–13% of the whole path, but about 40% of what's left once the
      drift (the gene's overall coding composition) is removed, so read the drift-removed bends as a sketch.</p>
    <div class="heat" style=${{ gridTemplateColumns: `110px repeat(${gene.windows.length}, minmax(18px, 1fr))` }}>
      <span></span>${gene.windows.map((_, w) => html`<span class="mono heathead" key=${w}>${w + 1}</span>`)}
      ${featureNames.map((name, f) => html`<span class="mono" key=${name}>${name}</span>
        ${gene.windows.map((m, w) => html`<i key=${`${name}${w}`} title=${`${name}, window ${w + 1}: ${m[f].toFixed(3)}`}
          style=${{ background: `color-mix(in srgb, var(${m[f] < 0 ? "--p2" : "--p1"}) ${Math.round(Math.abs(m[f]) / peak * 100)}%, transparent)` }}></i>`)}`)}
    </div>
    <p class="hint">The 9 stored numbers per window, blue positive and orange negative, darker further from zero. x purine − pyrimidine,
      y amino − keto, z strong − weak (GC).</p>`;
}
