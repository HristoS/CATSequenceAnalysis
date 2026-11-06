import { get, getBytes } from "./api.js";
import { Message, Segmented, formatBp, html, useEffect, useMemo, useRef, useState } from "./html.js";
import { HeatStrip, genePath, loadGenes } from "./genes.js";
import { drawCurves, parse, path, views } from "./zcurve3d.js";

// Two or three records side by side on one set of plots: Z-curve paths, matching genes, signatures. Each record keeps its colour on
// every plot, and the legend's switches show or hide it everywhere.
const colours = ["--p1", "--p2", "--p3"];
export const maxCompared = 3;
const familiesShown = 30;
const genesListed = 200;

const idsOf = params => [...new Set((params.get("ids") ?? "").split(",").map(Number).filter(n => Number.isInteger(n) && n > 0))].slice(0, maxCompared);
export const compareHref = ids => `#/compare?ids=${ids.slice(0, maxCompared).join(",")}`;

export function ComparePage({ params }) {
  const ids = idsOf(params);
  const [details, setDetails] = useState({});
  const [errors, setErrors] = useState({});
  const [hidden, setHidden] = useState(new Set());
  const [adding, setAdding] = useState("");

  useEffect(() => {
    let live = true;
    for (const id of ids)
      if (!details[id]) get(`/api/sequences/${id}`)
        .then(d => { if (live) setDetails(old => ({ ...old, [id]: d })); })
        .catch(e => { if (live) setErrors(old => ({ ...old, [id]: e.message })); });
    return () => { live = false; };
  }, [ids.join(",")]);

  const add = event => {
    event.preventDefault();
    const id = Number(adding.trim());
    if (Number.isInteger(id) && id > 0 && !ids.includes(id) && ids.length < maxCompared) location.hash = compareHref([...ids, id]);
    setAdding("");
  };
  const toggle = id => setHidden(old => { const next = new Set(old); next.has(id) ? next.delete(id) : next.add(id); return next; });

  // Records in the order of the URL, each with its colour; only those loaded and switched on are drawn.
  const records = ids.map((id, i) => ({ id, colour: colours[i], detail: details[id], error: errors[id] }));
  const shown = records.filter(r => r.detail && !hidden.has(r.id));

  return html`<h1>Compare</h1>
    <div class="legend compare-legend">
      ${records.map(r => html`<span key=${r.id} class="chip">
        <button type="button" class="swatch" aria-pressed=${!hidden.has(r.id)} onClick=${() => toggle(r.id)}
          title=${hidden.has(r.id) ? "Show on every plot" : "Hide from every plot"} style=${{ "--swatch": `var(${r.colour})` }}></button>
        ${r.detail ? html`<a href="#/sequences/${r.id}">${r.detail.sequence.name || r.detail.sequence.recordId}</a>
            <span class="hint">${r.detail.sequence.kind} · ${formatBp(r.detail.sequence.lengthBp)}</span>`
          : r.error ? html`<span class="message error">${r.id}: ${r.error}</span>` : html`<span class="hint">${r.id}…</span>`}
        <a class="hint" href=${compareHref(ids.filter(x => x !== r.id))} title="Remove from the comparison">remove</a>
      </span>`)}
      ${ids.length < maxCompared && html`<form class="inline" onSubmit=${add}>
        <input type="number" min="1" value=${adding} onInput=${e => setAdding(e.target.value)} placeholder="record id" aria-label="Add a record by id" />
        <button type="submit">Add</button></form>`}
    </div>
    ${ids.length < 2 && html`<${Message}>Add a second record (up to ${maxCompared}) by its id, or use "Compare" on a record's nearest relatives.<//>`}
    ${shown.length > 0 && html`<div class="columns">
      <section class="panel"><h2>Z-curve paths</h2><${ComparePaths} records=${shown} /></section>
      <section class="panel"><h2>Signature</h2><${CompareSignatures} records=${shown} /></section>
    </div>
    <section class="panel compare-genes"><h2>Genes</h2><${CompareGenes} records=${records.filter(r => r.detail)} hidden=${hidden} /></section>`}`;
}

// One rotatable canvas for several curves: drag, the view presets, redraw on resize, and on hover the nearest point's curve label and
// describe(point index) text.
export function CurveCanvas({ curves, unit, label }) {
  const [view, setView] = useState("3d");
  const canvas = useRef(null), tip = useRef(null), screens = useRef([]), angles = useRef(views["3d"]), drag = useRef(null);
  const redraw = () => {
    if (canvas.current) screens.current = drawCurves(canvas.current, curves, { yaw: angles.current[0], pitch: angles.current[1], unit }).screens;
  };
  useEffect(() => { if (view !== "free") angles.current = views[view]; redraw(); }, [curves, view]);
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
    let best = null, bestDistance = 14 * 14;
    screens.current.forEach((screen, c) => screen.forEach((p, k) => {
      const d = (p[0] - mx) ** 2 + (p[1] - my) ** 2;
      if (d < bestDistance) { bestDistance = d; best = [c, k]; }
    }));
    if (!best) { tip.current.hidden = true; return; }
    const [c, k] = best, curve = curves[c], point = screens.current[c][k];
    tip.current.textContent = [curve.label, curve.describe?.(k)].filter(Boolean).join(" · ");
    tip.current.hidden = false;
    tip.current.style.left = `${Math.max(6, Math.min(point[0] + 12, r.width - tip.current.offsetWidth - 6))}px`;
    tip.current.style.top = `${Math.max(6, point[1] - 14)}px`;
  };
  return html`<${Segmented} label="View" value=${view} onChange=${setView} options=${[["3d", "3-D"], ["xy", "x–y"], ["xz", "x–z"], ["yz", "y–z"]]} />
    <div class="zpath"><canvas ref=${canvas} role="img" aria-label=${label}
      onPointerDown=${down} onPointerMove=${move} onPointerUp=${up} onPointerCancel=${up}
      onPointerLeave=${() => { if (!drag.current) tip.current.hidden = true; }}></canvas>
      <div class="ztip" ref=${tip} hidden></div></div>`;
}

const nameOf = r => r.detail.sequence.name || r.detail.sequence.recordId;
const aminoAcidAt = (shape, i) => `about amino acid ${Math.round(shape.at[i] / 2).toLocaleString()}`;

const removeDrift = (sums, fits, at) => sums.map((p, i) => p.map((v, a) => v - (fits[a][0] + fits[a][1] * at(i))));

function ComparePaths({ records }) {
  const anyOrigin = records.some(r => r.detail.levels.some(l => l.methodVersion.endsWith("-origin")));
  const [layout, setLayout] = useState("origin");
  const [scale, setScale] = useState("length");
  const [drift, setDrift] = useState("removed");
  const [profiles, setProfiles] = useState({});
  const [error, setError] = useState(null);

  // A chromosome with an origin-rotated level is drawn from its replication origin (2 kbp windows); anything else as deposited (500 bp).
  const wanted = records.map(r => {
    const origin = layout === "origin" && r.detail.levels.some(l => l.methodVersion.endsWith("-origin") && l.windowSize === 2000);
    return { id: r.id, origin, key: `${r.id}-${origin ? "origin" : "deposited"}` };
  });
  useEffect(() => {
    let live = true;
    setError(null);
    for (const w of wanted)
      if (!profiles[w.key]) getBytes(`/api/sequences/${w.id}/profile?${w.origin ? "windowSize=2000&origin=true" : "windowSize=500"}`)
        .then(buffer => { if (live) setProfiles(old => ({ ...old, [w.key]: parse(buffer) })); })
        .catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [wanted.map(w => w.key).join(",")]);

  const curves = useMemo(() => wanted.flatMap((w, i) => {
    const profile = profiles[w.key];
    if (!profile) return [];
    const shape = path(profile), length = records[i].detail.sequence.lengthBp;
    let points = Array.from({ length: profile.n + 1 }, (_, k) => shape.sums.map(s => s[k]));   // sums are typed arrays, one per axis
    if (drift === "removed") points = removeDrift(points, shape.fits, k => k);
    if (scale === "length") points = points.map(p => p.map(v => v * 1000 / length));
    const describe = k => `${formatBp(Math.min(k * profile.size, length))}${w.origin ? " from the origin" : ""} of ${formatBp(length)}`;
    return [{ points, from: records[i].colour, to: records[i].colour, label: nameOf(records[i]), describe }];
  }), [profiles, wanted.map(w => w.key).join(","), drift, scale]);

  if (error) return html`<${Message}>${error}<//>`;
  const deposited = wanted.filter((w, i) => layout === "origin" && !w.origin).map(w => records.find(r => r.id === w.id).detail.sequence);
  return html`<div class="toolbar">
      ${anyOrigin && html`<${Segmented} label="Layout" value=${layout} onChange=${setLayout} options=${[["origin", "From origin"], ["deposited", "As deposited"]]} />`}
      <${Segmented} label="Scale" value=${scale} onChange=${setScale} options=${[["length", "Same length"], ["bases", "Bases"]]} />
      <${Segmented} label="Drift" value=${drift} onChange=${setDrift} options=${[["removed", "Removed"], ["kept", "Kept"]]} />
    </div>
    <${CurveCanvas} curves=${curves} unit=${scale === "length" ? "bases of excess per kbp of genome" : "bases of excess"}
      label="Z-curve paths of the compared records on one set of axes; drag to rotate" />
    <p class="hint">Each record's running sums in its colour; the grey ring is each start. Same length divides each path by its genome's length,
      so a phage and a bacterium can share the axes. ${anyOrigin && layout === "origin"
        ? "Circular chromosomes start at their replication origin (2 kbp windows), which is what lines related genomes up: as deposited, each starts wherever its submitter chose."
        : "As deposited, each path starts where its record does (500 bp windows)."}
      ${deposited.length > 0 && ` Drawn as deposited, with no stored origin: ${deposited.map(s => s.name || s.recordId).join(", ")}.`}</p>`;
}

function CompareSignatures({ records }) {
  const [view, setView] = useState("karlin");
  const [values, setValues] = useState({});
  useEffect(() => {
    let live = true;
    for (const r of records) {
      const key = `${r.id}-${view}`;
      if (values[key]) continue;
      const load = view === "genes"
        ? loadGenes(r.id).then(l => (l.geneCount ? { values: l.meanWindow } : { error: "no stored genes" }))
        : get(`/api/sequences/${r.id}/signature?view=${view}`).then(s => ({ values: s.values }));
      load.catch(e => ({ error: e.message })).then(v => { if (live) setValues(old => ({ ...old, [key]: v })); });
    }
    return () => { live = false; };
  }, [records.map(r => r.id).join(","), view]);

  const rows = records.map(r => ({ ...r, result: values[`${r.id}-${view}`] }));
  const names = rows.find(r => r.result?.values)?.result.values.map(v => v.name) ?? [];
  const peak = Math.max(...rows.flatMap(r => r.result?.values?.map(v => Math.abs(v.value)) ?? []), 1e-9);
  return html`<${Segmented} label="View" value=${view} onChange=${setView}
      options=${[["karlin", "Karlin ρ*"], ["zcurve", "Z-curve coefficients"], ["genes", "Gene layer"]]} />
    ${rows.filter(r => r.result?.error).map(r => html`<p class="hint" key=${r.id}>
      <span class="dot" style=${{ "--swatch": `var(${r.colour})` }}></span>${r.detail.sequence.name || r.detail.sequence.recordId}: ${r.result.error}</p>`)}
    <div class="bars">${names.map((name, i) => html`<div class="bar gbar" key=${name}>
      <span class="mono">${name}</span>
      <div class="gtracks">${rows.filter(r => r.result?.values).map(r => {
        const v = r.result.values[i].value;
        return html`<div class="track" key=${r.id} title=${`${r.detail.sequence.name || r.id}: ${v.toFixed(3)}`}>
          <i class=${v < 0 ? "negative" : "positive"} style=${{ width: `${(Math.abs(v) / peak) * 50}%`, background: `var(${r.colour})` }}></i></div>`;
      })}</div>
    </div>`)}</div>
    <p class="hint">${view === "genes" ? "Each record's mean gene window: x, y, z over codon positions 1 and 2, their trend and curvature."
      : "The dinucleotide signature, one bar per record in its colour on a shared scale; hover a bar for its value."}</p>`;
}

function CompareGenes({ records, hidden }) {
  const [lists, setLists] = useState({});
  const [text, setText] = useState("");
  const [picked, setPicked] = useState({});   // record id -> gene index
  const [windows, setWindows] = useState({});
  const [drift, setDrift] = useState("kept");

  useEffect(() => {
    let live = true;
    for (const r of records)
      if (!lists[r.id]) loadGenes(r.id).then(l => { if (live) setLists(old => ({ ...old, [r.id]: l })); })
        .catch(e => { if (live) setLists(old => ({ ...old, [r.id]: { error: e.message, genes: [], geneCount: 0 } })); });
    return () => { live = false; };
  }, [records.map(r => r.id).join(",")]);

  const loaded = records.filter(r => lists[r.id]);
  const withGenes = loaded.filter(r => lists[r.id].geneCount > 0);
  const needle = text.trim().toLowerCase();
  const matches = g => !needle || (g.product ?? "").toLowerCase().includes(needle) || (g.familyId ?? "").toLowerCase().includes(needle);

  // Families every record with genes carries, with how many copies each has.
  const common = useMemo(() => {
    if (withGenes.length < 2) return [];
    const counts = withGenes.map(r => {
      const m = new Map();
      for (const g of lists[r.id].genes) if (g.familyId) m.set(g.familyId, (m.get(g.familyId) ?? 0) + 1);
      return m;
    });
    return [...counts[0].keys()].filter(f => counts.every(c => c.has(f)))
      .map(f => ({ familyId: f, product: lists[withGenes[0].id].genes.find(g => g.familyId === f)?.product, copies: counts.map(c => c.get(f)) }));
  }, [withGenes.map(r => r.id).join(","), Object.keys(lists).length]);
  const shownFamilies = common.filter(f => !needle || (f.product ?? "").toLowerCase().includes(needle) || f.familyId.toLowerCase().includes(needle));

  const pickFamily = familyId => setPicked(Object.fromEntries(withGenes.map(r => [r.id, lists[r.id].genes.find(g => g.familyId === familyId)?.index]).filter(([, i]) => i != null)));

  useEffect(() => {
    let live = true;
    for (const [id, index] of Object.entries(picked)) {
      const key = `${id}-${index}`;
      if (!windows[key]) get(`/api/sequences/${id}/genes/${index}`).then(g => { if (live) setWindows(old => ({ ...old, [key]: g })); }).catch(() => {});
    }
    return () => { live = false; };
  }, [JSON.stringify(picked)]);

  const selected = withGenes.filter(r => !hidden.has(r.id) && picked[r.id] != null && windows[`${r.id}-${picked[r.id]}`])
    .map(r => ({ ...r, gene: windows[`${r.id}-${picked[r.id]}`] }));
  const curves = useMemo(() => selected.map(s => {
    const shape = genePath(s.gene.windows, s.gene.windowPositions);
    const points = drift === "kept" ? shape.points : removeDrift(shape.points, shape.fits, i => shape.at[i]);
    return { points, from: s.colour, to: s.colour, marks: shape.marks, label: nameOf(s), describe: i => aminoAcidAt(shape, i) };
  }), [selected.map(s => `${s.id}-${s.gene.gene.index}`).join(","), drift]);
  const peak = Math.max(...selected.flatMap(s => s.gene.windows.flat().map(Math.abs)), 1e-9);

  if (loaded.length < records.length) return html`<${Message}>Loading the genes…<//>`;
  if (withGenes.length < 2) return html`<${Message}>At least two of the compared records need stored genes (BV-BRC genomes).<//>`;
  return html`<label class="field"><span>Product or family contains</span>
      <input value=${text} onInput=${e => setText(e.target.value)} placeholder="e.g. gyrase" /></label>
    <p class="hint">${common.length.toLocaleString()} protein families are in all ${withGenes.length} records${needle ? `; ${shownFamilies.length} match` : ""}. Pick one to
      draw each record's copy, or choose any gene per record below.</p>
    <div class="tablewrap"><table>
      <thead><tr><th>Family</th><th>Product</th>${withGenes.map(r => html`<th class="num" key=${r.id}><span class="dot" style=${{ "--swatch": `var(${r.colour})` }}></span>copies</th>`)}</tr></thead>
      <tbody>${shownFamilies.slice(0, familiesShown).map(f => html`<tr key=${f.familyId} class="row" tabIndex="0"
          onClick=${() => pickFamily(f.familyId)} onKeyDown=${e => { if (e.key === "Enter") pickFamily(f.familyId); }}>
        <td class="mono">${f.familyId}</td><td>${f.product ?? "–"}</td>${f.copies.map((c, i) => html`<td class="num" key=${i}>${c}</td>`)}
      </tr>`)}</tbody>
    </table></div>
    ${shownFamilies.length > familiesShown && html`<p class="hint">The first ${familiesShown}; narrow the filter to see others.</p>`}
    <div class="genepicks">${withGenes.map(r => {
      const list = lists[r.id].genes, options = list.filter(g => matches(g) || g.index === picked[r.id]).slice(0, genesListed);
      return html`<label class="field" key=${r.id}><span><span class="dot" style=${{ "--swatch": `var(${r.colour})` }}></span>${r.detail.sequence.name || r.detail.sequence.recordId}</span>
        <select value=${picked[r.id] ?? ""} onChange=${e => setPicked(old => ({ ...old, [r.id]: e.target.value === "" ? undefined : Number(e.target.value) }))}>
          <option value="">none</option>
          ${options.map(g => html`<option key=${g.index} value=${g.index}>${`#${g.index + 1} · ${g.product ?? g.familyId ?? "no family"} · ${g.residues} aa`}</option>`)}
        </select></label>`;
    })}</div>
    ${selected.length > 0 && html`<div class="toolbar"><${Segmented} label="Drift" value=${drift} onChange=${setDrift} options=${[["kept", "Kept"], ["removed", "Removed"]]} /></div>
      <${CurveCanvas} curves=${curves} unit="bases" label="Z-curves of the picked genes over codon positions 1 and 2; drag to rotate" />
      <p class="hint">Each picked gene's Z-curve over codon positions 1 and 2, redrawn from its stored windows; rings mark window ends.
        Genes of different lengths have different window sizes (5 windows under 500 aa). Near-identical genes draw on top of each other, the
        last record's on top: switch records off in the legend to see the ones beneath.</p>
      <div class="heatgrid">${selected.map(s => html`<div key=${s.id}>
        <p class="heatlabel"><span class="dot" style=${{ "--swatch": `var(${s.colour})` }}></span>${`${s.gene.gene.product ?? s.gene.gene.familyId ?? "Gene"} · ` +
          `${s.gene.gene.residues} aa · ${s.detail.sequence.name || s.detail.sequence.recordId}`}</p>
        <${HeatStrip} windows=${s.gene.windows} peak=${peak} /></div>`)}</div>`}`;
}
