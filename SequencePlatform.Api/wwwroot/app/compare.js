import { get, getBytes, loadGenes } from "./api.js";
import { Message, Pager, Segmented, formatBp, html, useEffect, useMemo, usePaged, useRef, useState } from "./html.js";
import { HeatStrip, genePath } from "./genes.js";
import { drawCurves, genePathOf, lineFits, parse, path, useZoom, views, windowRanges } from "./zcurve3d.js";

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
  const [picked, setPicked] = useState({});   // record id -> gene index, chosen in the genes panel and highlighted on the paths

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
      <section class="panel"><h2>Z-curve paths</h2><${ComparePaths} records=${shown} picked=${picked} /></section>
      <section class="panel"><h2>Signature</h2><${CompareSignatures} records=${shown} /></section>
    </div>
    <section class="panel compare-genes"><h2>Genes</h2><${CompareGenes} records=${records.filter(r => r.detail)} hidden=${hidden} picked=${picked} setPicked=${setPicked} /></section>`}`;
}

// One rotatable canvas for several curves: drag, the view presets, redraw on resize, and on hover the nearest point's curve label and
// describe(point index) text.
export function CurveCanvas({ curves, unit, label }) {
  const [view, setView] = useState("3d");
  const canvas = useRef(null), tip = useRef(null), screens = useRef([]), angles = useRef(views["3d"]), drag = useRef(null);
  const redraw = () => {
    if (canvas.current) screens.current = drawCurves(canvas.current, curves, { yaw: angles.current[0], pitch: angles.current[1], zoom: zoomer.zoom.current, unit }).screens;
  };
  const zoomer = useZoom(canvas, redraw, [curves]);
  useEffect(() => { if (view !== "free") angles.current = views[view]; redraw(); }, [curves, view]);
  useEffect(() => {
    if (!canvas.current) return;
    const observer = new ResizeObserver(redraw);
    observer.observe(canvas.current);
    return () => observer.disconnect();
  });
  const down = e => {
    tip.current.hidden = true;
    if (zoomer.down(e)) return;
    drag.current = { x: e.clientX, y: e.clientY, angles: angles.current }; e.currentTarget.setPointerCapture(e.pointerId);
  };
  const up = () => { drag.current = null; zoomer.up(); };
  const move = e => {
    if (zoomer.move(e)) return;
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
      <div class="ztip" ref=${tip} hidden></div></div>
    <p class="hint">Drag to rotate, scroll to zoom, shift-drag to pan, double-click to reset.</p>`;
}

const nameOf = r => r.detail.sequence.name || r.detail.sequence.recordId;
const aminoAcidAt = (shape, i) => `about amino acid ${Math.round(shape.at[i] / 2).toLocaleString()}`;

const removeDrift = (sums, fits, at) => sums.map((p, i) => p.map((v, a) => v - (fits[a][0] + fits[a][1] * at(i))));
const withoutDrift = (points, fits) => removeDrift(points, fits, k => k);

const sequencesOf = r => r.detail.annotations?.ncbi?.gene_layer?.sequences ?? [];
const hasChromosomesOf = r => r.detail.annotations?.ncbi?.chromosome_levels?.status === "stored" && r.detail.annotations.ncbi.chromosome_levels.sequences > 0;
const chromosomeOverviewBp = 20_000_000;   // as on the record page: longer chromosomes are drawn from their 32 kbp level

// Aligned layout: each record's 500 bp windows started and oriented where they best follow the first record's. Genomes deposited from
// different starts (circularly permuted phages) or on the other strand then line up, as find-similar now scores viruses wrapped around
// their end (WrapAround, results/CircularPermutation_Check.csv). Linear chromosomes (eukaryotes) aren't moved.
const alignable = r => r.detail.sequence.kind !== "eukaryote" && r.detail.levels.some(l => l.windowSize === 500 && l.methodVersion === "1");
const depositedUrl = r => `/api/sequences/${r.id}/profile?windowSize=500`;

// Per-window mean axes, x y z interleaved, optionally reversed onto the other strand (order reversed, x and y negated).
function meansOf(p, reverse = false) {
  const m = new Float64Array(3 * p.n);
  for (let i = 0; i < p.n; i++) {
    const j = reverse ? p.n - 1 - i : i, sign = reverse ? -1 : 1;
    m[3 * i] = sign * p.x[j] / p.size; m[3 * i + 1] = sign * p.y[j] / p.size; m[3 * i + 2] = p.z[j] / p.size;
  }
  return m;
}
// Every f windows averaged, so a long genome's starts are searched on a few thousand points first.
function coarse(m, f) {
  const n = Math.floor(m.length / 3 / f), c = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) for (let k = 0; k < f; k++) for (let a = 0; a < 3; a++) c[3 * i + a] += m[3 * (i * f + k) + a] / f;
  return c;
}
// Mean per-window distance between the reference and `other` started at window `shift`, wrapping around other's end.
function distanceAt(ref, other, shift) {
  const nr = ref.length / 3, no = other.length / 3, m = Math.min(nr, no);
  let d = 0;
  for (let i = 0; i < m; i++) {
    const j = 3 * ((i + shift) % no), k = 3 * i;
    d += Math.hypot(ref[k] - other[j], ref[k + 1] - other[j + 1], ref[k + 2] - other[j + 2]);
  }
  return d / Math.max(1, m);
}
function bestShift(ref, other) {
  const n = other.length / 3, f = Math.max(1, Math.ceil(Math.max(ref.length, other.length) / 3 / 1500));
  let best = { shift: 0, distance: Infinity };
  if (f === 1) {
    for (let r = 0; r < n; r++) { const d = distanceAt(ref, other, r); if (d < best.distance) best = { shift: r, distance: d }; }
    return best;
  }
  const cr = coarse(ref, f), co = coarse(other, f);
  let rough = 0, roughDistance = Infinity;
  for (let r = 0; r < co.length / 3; r++) { const d = distanceAt(cr, co, r); if (d < roughDistance) { roughDistance = d; rough = r; } }
  for (let r = rough * f - f; r <= rough * f + f; r++) {
    const shift = ((r % n) + n) % n, d = distanceAt(ref, other, shift);
    if (d < best.distance) best = { shift, distance: d };
  }
  return best;
}
// The start (in windows) and strand that make `profile` follow `reference` most closely, with the distance as deposited for comparison.
function alignTo(reference, profile) {
  const ref = meansOf(reference);
  const forward = bestShift(ref, meansOf(profile)), reverse = bestShift(ref, meansOf(profile, true));
  const pick = reverse.distance < forward.distance ? { ...reverse, reverse: true } : { ...forward, reverse: false };
  return { ...pick, deposited: distanceAt(ref, meansOf(profile), 0) };
}
// A parsed profile turned by an alignment: window i is the original window (i + shift) mod n, of the other strand when reverse.
function turnProfile(p, { shift, reverse }) {
  const n = p.n, take = (values, sign) => Int32Array.from({ length: n }, (_, i) => {
    const j = (i + shift) % n;
    return sign * values[reverse ? n - 1 - j : j];
  });
  return { ...p, x: take(p.x, reverse ? -1 : 1), y: take(p.y, reverse ? -1 : 1), z: take(p.z, 1) };
}
// Where a deposited base position lands once the windows are turned (positions past the last whole window wrap with them).
const turnedBp = (bp, { shift, reverse }, size, n) => {
  const length = n * size, mirrored = reverse ? length - bp : bp;
  return (((mirrored - shift * size) % length) + length) % length;
};
// A gene path turned the same way: genes reordered by where they start once turned, and on the other strand every step's x and y negated
// (each gene's own steps stay in coding orientation, as the path is built).
function turnGenePath(gp, genes, alignment, size, n) {
  const blocks = [...gp.spans.entries()].map(([index, [first, last]]) => {
    const g = genes[index];
    return { index, first, last, at: turnedBp(alignment.reverse ? g.endBp : g.startBp, alignment, size, n) };
  }).sort((a, b) => a.at - b.at);
  const sign = alignment.reverse ? [-1, -1, 1] : [1, 1, 1], points = [[0, 0, 0]], spans = new Map(), starts = [], order = [];
  for (const b of blocks) {
    starts.push(points.length - 1); order.push(b.index);
    for (let k = b.first; k < b.last; k++) {
      const prev = points[points.length - 1];
      points.push(prev.map((v, a) => v + sign[a] * (gp.points[k + 1][a] - gp.points[k][a])));
    }
    spans.set(b.index, [starts[starts.length - 1], points.length - 1]);
  }
  const geneAt = k => {
    let lo = 0, hi = starts.length - 1;
    while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (starts[mid] < k) lo = mid; else hi = mid - 1; }
    return order[lo];
  };
  return { points, spans, geneAt };
}

// The records' Z-curve paths on one canvas: from their stored windows, or (Source: Genes) built from their stored genes alone. The genes
// picked below (picked: record id -> gene index) are highlighted on every path. A eukaryote is one chromosome at a time, chosen per record
// by name (and moved to a picked gene's): its windows from that chromosome's own level, its gene path from the genes on it; never joined.
function ComparePaths({ records, picked }) {
  const anyOrigin = records.some(r => r.detail.levels.some(l => l.methodVersion.endsWith("-origin")));
  const withGenes = records.filter(r => r.detail.sequence.provider === "BV-BRC" || r.detail.annotations?.ncbi?.gene_layer != null);
  const withWindows = records.filter(r => r.detail.levels.some(l => l.windowSize === 500 && l.methodVersion === "1") || hasChromosomesOf(r));
  const [source, setSource] = useState(withWindows.length > 0 || withGenes.length === 0 ? "windows" : "genes");
  const [layout, setLayout] = useState(anyOrigin ? "origin" : "deposited");
  const [scale, setScale] = useState("length");
  const [drift, setDrift] = useState("removed");
  const [chosen, setChosen] = useState({});   // record id -> name of the sequence drawn, for records with several
  const [chromosomes, setChromosomes] = useState({});   // record id -> its stored chromosome levels
  const [profiles, setProfiles] = useState({});
  const [lists, setLists] = useState({});
  const [errors, setErrors] = useState({});
  const fromGenes = source === "genes";

  // A picked gene on another chromosome moves that record's path to it.
  useEffect(() => {
    const moves = Object.fromEntries(withGenes.filter(r => sequencesOf(r).length > 0 && picked[r.id] != null && lists[r.id])
      .map(r => [r.id, sequencesOf(r)[lists[r.id].genes[picked[r.id]]?.sequence]]).filter(([id, n]) => n != null && n !== chosen[id]));
    if (Object.keys(moves).length) setChosen(old => ({ ...old, ...moves }));
  }, [JSON.stringify(picked), Object.keys(lists).length]);
  useEffect(() => {
    let live = true;
    for (const r of records.filter(hasChromosomesOf))
      if (!chromosomes[r.id]) get(`/api/sequences/${r.id}/chromosomes`).then(l => { if (live) setChromosomes(old => ({ ...old, [r.id]: l.chromosomes })); })
        .catch(e => { if (live) setErrors(old => ({ ...old, [`${r.id}-chromosomes`]: `${nameOf(r)}: ${e.message}` })); });
    return () => { live = false; };
  }, [records.map(r => r.id).join(",")]);

  // A chromosome with an origin-rotated level is drawn from its replication origin (2 kbp windows, or its genes from there); anything
  // else as deposited (500 bp windows, or its genes from 0). Records without genes have no gene path, and without windows no window path.
  const wanted = (fromGenes ? withGenes : withWindows).flatMap(r => {
    const origin = layout === "origin" && r.detail.levels.some(l => l.methodVersion.endsWith("-origin") && l.windowSize === 2000);
    const perChromosome = !fromGenes && hasChromosomesOf(r);
    if (perChromosome && !chromosomes[r.id]) return [];
    const names = fromGenes ? sequencesOf(r) : perChromosome ? chromosomes[r.id].map(c => c.name) : [];
    const name = names.length === 0 ? null : names.includes(chosen[r.id]) ? chosen[r.id] : names[0];
    if (fromGenes) {
      const sequence = name == null ? null : names.indexOf(name);
      return [{ id: r.id, record: r, origin, names, name, sequence,
        url: `/api/sequences/${r.id}/genes/path?origin=${origin}${sequence == null ? "" : `&sequence=${sequence}`}` }];
    }
    if (perChromosome) {
      const c = chromosomes[r.id].find(x => x.name === name), size = c.lengthBp > chromosomeOverviewBp ? 32000 : 500;
      return [{ id: r.id, record: r, origin: false, names, name, chromosome: c,
        url: `/api/sequences/${r.id}/chromosomes/${c.index}/profile?windowSize=${size}` }];
    }
    return [{ id: r.id, record: r, origin, names, name,
      url: `/api/sequences/${r.id}/profile?${origin ? "windowSize=2000&origin=true" : "windowSize=500"}` }];
  }).map(w => ({ ...w, key: w.url }));
  useEffect(() => {
    let live = true;
    for (const w of wanted) {
      if (profiles[w.key] || errors[w.key]) continue;
      const load = fromGenes ? get(w.url).then(genePathOf) : getBytes(w.url).then(parse);
      load.then(p => { if (live) setProfiles(old => ({ ...old, [w.key]: p })); })
        .catch(e => { if (live) setErrors(old => ({ ...old, [w.key]: `${nameOf(w.record)}: ${e.message}` })); });
    }
    for (const r of withGenes)
      if (!lists[r.id]) loadGenes(r.id).then(l => { if (live) setLists(old => ({ ...old, [r.id]: l })); }).catch(() => {});
    if (layout === "aligned")
      for (const w of wanted.filter(w => alignable(w.record))) {
        const url = depositedUrl(w.record);
        if (!profiles[url] && !errors[url]) getBytes(url).then(parse).then(p => { if (live) setProfiles(old => ({ ...old, [url]: p })); })
          .catch(e => { if (live) setErrors(old => ({ ...old, [url]: `${nameOf(w.record)}: ${e.message}` })); });
      }
    return () => { live = false; };
  }, [wanted.map(w => w.key).join(","), layout]);

  // Aligned layout: each alignable record's start and strand against the first alignable record drawn.
  const reference = layout === "aligned" ? wanted.find(w => alignable(w.record)) : null;
  const alignments = useMemo(() => {
    const base = reference && profiles[depositedUrl(reference.record)];
    if (!base) return {};
    return Object.fromEntries(wanted.filter(w => w !== reference && alignable(w.record) && profiles[depositedUrl(w.record)])
      .map(w => [w.id, alignTo(base, profiles[depositedUrl(w.record)])]));
  }, [reference?.id, wanted.map(w => w.key).join(","), Object.keys(profiles).length]);

  // The picked gene of a record, when it lies on what that record's path draws.
  const geneOn = w => {
    const g = lists[w.id]?.genes[picked[w.id]];
    if (!g) return null;
    if (fromGenes) return profiles[w.key]?.spans.has(g.index) ? g : null;
    return w.chromosome && sequencesOf(w.record)[g.sequence] !== w.chromosome.name ? null : g;
  };
  const curves = useMemo(() => wanted.flatMap(w => {
    let loaded = profiles[w.key];
    if (!loaded) return [];
    const r = w.record, list = lists[r.id], gene = geneOn(w);
    // Aligned: wait for the alignment (and, for genes, the gene list) before drawing a record.
    const alignment = alignments[w.id], windows = profiles[depositedUrl(r)];
    if (layout === "aligned" && alignable(r) && w !== reference && (!alignment || (fromGenes && !list))) return [];
    if (alignment) loaded = fromGenes ? turnGenePath(loaded, list.genes, alignment, windows.size, windows.n) : turnProfile(loaded, alignment);
    // A chromosome's windows are on the same length by the chromosome's; one chromosome's gene path by how far its genes reach, the
    // only length of it stored.
    const length = w.chromosome ? w.chromosome.lengthBp : w.sequence == null ? r.detail.sequence.lengthBp
      : Math.max(1, ...(list?.genes ?? []).filter(g => g.sequence === w.sequence).map(g => g.endBp));
    let points, highlights = [], describe;
    if (fromGenes) {
      points = drift === "removed" ? withoutDrift(loaded.points, lineFits(loaded.points)) : loaded.points;
      if (gene) highlights = [loaded.spans.get(gene.index)];
      describe = k => {
        const g = list?.genes[loaded.geneAt(k)];
        return g ? `gene #${g.index + 1}${w.sequence == null ? "" : ` on ${w.name}`} · ${g.product ?? g.familyId ?? "no family"}` : "";
      };
    } else {
      const shape = path(loaded);
      points = Array.from({ length: loaded.n + 1 }, (_, k) => shape.sums.map(s => s[k]));   // sums are typed arrays, one per axis
      if (drift === "removed") points = withoutDrift(points, shape.fits.map(f => f));
      if (gene && alignment) {
        const a = turnedBp(alignment.reverse ? gene.endBp : gene.startBp, alignment, loaded.size, loaded.n);
        highlights = windowRanges(a, a + (gene.endBp - gene.startBp), 0, loaded.n * loaded.size, loaded.size, loaded.n);
      } else if (gene) highlights = windowRanges(gene.startBp, gene.endBp, w.origin ? r.detail.originBp ?? 0 : 0, length, loaded.size, loaded.n);
      const at = k => !alignment ? Math.min(k * loaded.size, length) : alignment.reverse
        ? loaded.n * loaded.size - ((k + alignment.shift) % loaded.n) * loaded.size : ((k + alignment.shift) % loaded.n) * loaded.size;
      describe = k => `${formatBp(at(k))}${w.origin ? " from the origin" : ""} of ${formatBp(length)}` +
        (alignment?.reverse ? " (other strand)" : "") + (w.chromosome ? ` (${w.chromosome.name}, ${formatBp(loaded.size)} windows)` : "");
    }
    if (scale === "length") points = points.map(p => p.map(v => v * 1000 / length));
    return [{ points, from: r.colour, to: r.colour, label: nameOf(r), describe, highlights }];
  }), [profiles, lists, wanted.map(w => w.key).join(","), drift, scale, JSON.stringify(picked), alignments, layout]);

  const moved = wanted.filter(w => alignments[w.id]).map(w => {
    const a = alignments[w.id], p = profiles[depositedUrl(w.record)];
    const where = a.shift === 0 && !a.reverse ? "kept as deposited" : `${a.reverse ? "the other strand, " : ""}started ${formatBp(a.shift * p.size)} on`;
    return `${nameOf(w.record)} ${where} (mean window distance ${a.distance.toFixed(3)}, ${a.deposited.toFixed(3)} as deposited)`;
  });
  const alignedNote = reference ? "Each path is started, and turned onto the other strand when that fits better, where it best follows " +
    `${nameOf(reference.record)}'s 500 bp windows, wrapping around its end: genomes deposited from different starts, as circularly permuted ` +
    `phages often are, then line up. ${moved.length ? moved.join("; ") + ". " : ""}` +
    (wanted.some(w => !alignable(w.record)) ? "Linear chromosomes (eukaryotes) aren't moved. " : "") : "";
  // Only circular chromosomes have an origin to start from; a eukaryote is always drawn as deposited.
  const deposited = wanted.filter(w => anyOrigin && layout === "origin" && !w.origin && w.record.detail.sequence.kind !== "eukaryote").map(w => nameOf(w.record));
  const windowsLabel = withWindows.every(hasChromosomesOf) ? "Chromosome windows" : withWindows.some(hasChromosomesOf) ? "Genome / chromosome windows" : "Genome windows";
  const missing = records.filter(r => !(fromGenes ? withGenes : withWindows).includes(r)).map(nameOf);
  const failed = [...wanted.map(w => errors[w.key]), ...records.map(r => errors[`${r.id}-chromosomes`])].filter(Boolean);
  const pickedNames = wanted.map(geneOn).filter(Boolean);
  const pickers = wanted.filter(w => w.names.length > 1);
  return html`<div class="toolbar">
      ${withGenes.length > 0 && withWindows.length > 0 && html`<${Segmented} label="Source" value=${source} onChange=${setSource} options=${[["windows", windowsLabel], ["genes", "Genes"]]} />`}
      ${(anyOrigin || records.filter(alignable).length > 1) && html`<${Segmented} label="Layout" value=${layout} onChange=${setLayout}
        options=${[...(anyOrigin ? [["origin", "From origin"]] : []), ...(records.filter(alignable).length > 1 ? [["aligned", "Aligned to first"]] : []),
          ["deposited", "As deposited"]]} />`}
      <${Segmented} label="Scale" value=${scale} onChange=${setScale} options=${[["length", "Same length"], ["bases", "Bases"]]} />
      <${Segmented} label="Drift" value=${drift} onChange=${setDrift} options=${[["removed", "Removed"], ["kept", "Kept"]]} />
      ${pickers.map(w => html`<label class="field" key=${w.id}><span><span class="dot" style=${{ "--swatch": `var(${w.record.colour})` }}></span>Sequence</span>
        <select value=${w.name} onChange=${e => setChosen(old => ({ ...old, [w.id]: e.target.value }))}>
          ${w.names.map((n, i) => html`<option key=${n} value=${n}>${`${i + 1}. ${n}${w.chromosome ? ` · ${formatBp(chromosomes[w.id][i].lengthBp)}` : ""}`}</option>`)}
        </select></label>`)}
    </div>
    ${wanted.length === 0 && records.every(r => !hasChromosomesOf(r) || chromosomes[r.id]) && html`<${Message}>${fromGenes
      ? "None of these records has stored genes." : "None of these records has stored windows: eukaryotes keep them per chromosome once backfilled."}<//>`}
    ${failed.map(f => html`<p class="hint" key=${f}>${f}</p>`)}
    <${CurveCanvas} curves=${curves} unit=${scale === "length" ? "bases of excess per kbp of genome" : "bases of excess"}
      label="Z-curve paths of the compared records on one set of axes; drag to rotate" />
    ${pickedNames.length > 0 && html`<p class="hint">${`Highlighted: ${pickedNames.map(g => `${g.product ?? g.familyId ?? "gene"} (#${g.index + 1})`).join(", ")}; ` +
      (fromGenes ? "each an exact stretch of its record's gene path." : "on a whole chromosome a gene is a window or two, so it shows as a marked spot.")}</p>`}
    <p class="hint">${(fromGenes
        ? "Built from the gene layer alone: each record's stored genes in chromosome order, each window one step of its means over codon positions 1 and 2, " +
          "reverse-strand genes flipped onto the forward strand; the DNA between genes and every third codon position are left out. Relatives " +
          "line up a little less well on it than on the genome windows (results/GenePath_Check.csv). "
        : "Each record's running sums in its colour; the grey ring is each start. ") +
      "Same length divides each path by its genome's length, so a phage and a bacterium can share the axes. " +
      (layout === "aligned" ? alignedNote
        : anyOrigin && layout === "origin"
        ? "Circular chromosomes start at their replication origin, which is what lines related genomes up: as deposited, each starts wherever its submitter chose."
        : "As deposited, each path starts where its record does.") +
      (deposited.length > 0 ? ` Drawn as deposited, with no stored origin: ${deposited.join(", ")}.` : "") +
      (pickers.length > 0 ? (fromGenes
        ? " A eukaryote's gene path is one chromosome at a time, put on the same length by how far its genes reach; it follows which strand " +
          "each gene is on, so between relatives it shows gene order more than composition."
        : " A eukaryote is drawn one chromosome at a time from that chromosome's own windows (32 kbp above 20 Mbp, else 500 bp), on the same " +
          "length by the chromosome's.") : "") +
      (missing.length > 0 ? (fromGenes ? ` No stored genes, so no gene path: ${missing.join(", ")}.`
        : ` No stored windows (a eukaryote's chromosome levels not yet backfilled), so not drawn here: ${missing.join(", ")}.`) : "")}</p>`;
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

function CompareGenes({ records, hidden, picked, setPicked }) {
  const [lists, setLists] = useState({});
  const [text, setText] = useState("");
  const [windows, setWindows] = useState({});
  const [drift, setDrift] = useState("kept");
  const [matched, setMatched] = useState({});        // per record, the match found for another record's picked gene

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
  const shownFamilies = useMemo(() => common.filter(f => !needle || (f.product ?? "").toLowerCase().includes(needle) || f.familyId.toLowerCase().includes(needle)),
    [common, needle]);
  const familyPages = usePaged(shownFamilies, { size: familiesShown });

  // The picked gene of one record matched in every other record with genes (GeneMatching: its reciprocal best hit when there is one).
  const matchFrom = source => {
    const index = picked[source.id];
    setMatched(Object.fromEntries(withGenes.filter(r => r.id !== source.id).map(r => [r.id, { pending: true }])));
    Promise.all(withGenes.filter(r => r.id !== source.id).map(r =>
      get(`/api/sequences/${source.id}/genes/${index}/matches?in=${r.id}&top=1`).then(m => [r.id, m]).catch(e => [r.id, { error: e.message }])))
      .then(results => {
        setPicked(old => ({ ...old, ...Object.fromEntries(results.filter(([, m]) => m.matches?.length).map(([id, m]) => [id, m.matches[0].gene.index])) }));
        setMatched(Object.fromEntries(results.map(([id, m]) => [id, { ...m, source }])));
      });
  };
  const matchNote = r => {
    const m = matched[r.id];
    if (!m) return null;
    if (m.pending) return "Matching…";
    if (m.error) return m.error;
    const from = `${nameOf(m.source)} #${m.gene.index + 1}`;
    if (!m.matches.length) return `No gene here is comparable with ${from} (${m.gene.residues} aa): none within ±10% of its length${m.gene.residues >= 500 ? " of 500 aa or more" : ""}.`;
    const best = m.matches[0];
    if (picked[r.id] !== best.gene.index) return null;   // picked by hand since
    return `Best match of ${from} among ${m.compared.toLocaleString()} comparable genes, score ${best.score.toFixed(3)}: ` +
      (best.reciprocal ? "a reciprocal best hit (likely the same gene)." : `not reciprocal: it matches another gene of ${nameOf(m.source)} better, so read it as a relative, not the same gene.`);
  };

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
  if (withGenes.length < 2) return html`<${Message}>At least two of the compared records need stored genes (BV-BRC genomes or NCBI-annotated eukaryotes).<//>`;
  return html`<label class="field"><span>Product or family contains</span>
      <input value=${text} onInput=${e => setText(e.target.value)} placeholder="e.g. gyrase" /></label>
    <p class="hint">${withGenes.every(r => !lists[r.id].genes.some(g => g.familyId))
      ? "These genes have no protein families (NCBI's eukaryote annotations carry none, so copies can't be matched by family): pick a gene in one record, then “Find this gene” picks its best match in the others."
      : html`${common.length.toLocaleString()} protein families are in all ${withGenes.length} records${needle ? `; ${shownFamilies.length} match` : ""}. Pick one to
      draw each record's copy, or choose any gene per record below.`}</p>
    <div class="tablewrap"><table>
      <thead><tr><th>Family</th><th>Product</th>${withGenes.map(r => html`<th class="num" key=${r.id}><span class="dot" style=${{ "--swatch": `var(${r.colour})` }}></span>copies</th>`)}</tr></thead>
      <tbody>${familyPages.rows.map(f => html`<tr key=${f.familyId} class="row" tabIndex="0"
          onClick=${() => pickFamily(f.familyId)} onKeyDown=${e => { if (e.key === "Enter") pickFamily(f.familyId); }}>
        <td class="mono">${f.familyId}</td><td>${f.product ?? "–"}</td>${f.copies.map((c, i) => html`<td class="num" key=${i}>${c}</td>`)}
      </tr>`)}</tbody>
    </table></div>
    <${Pager} paged=${familyPages} noun="families" />
    <div class="genepicks">${withGenes.map(r => {
      const list = lists[r.id].genes, options = list.filter(g => matches(g) || g.index === picked[r.id]).slice(0, genesListed);
      return html`<label class="field" key=${r.id}><span><span class="dot" style=${{ "--swatch": `var(${r.colour})` }}></span>${r.detail.sequence.name || r.detail.sequence.recordId}</span>
        <select value=${picked[r.id] ?? ""} onChange=${e => setPicked(old => ({ ...old, [r.id]: e.target.value === "" ? undefined : Number(e.target.value) }))}>
          <option value="">none</option>
          ${options.map(g => html`<option key=${g.index} value=${g.index}>${`#${g.index + 1} · ${g.product ?? g.familyId ?? "no family"} · ${g.residues} aa`}</option>`)}
        </select>
        ${picked[r.id] != null && html`<button type="button" class="quiet" onClick=${() => matchFrom(r)}>Find this gene in the other ${withGenes.length > 2 ? "records" : "record"}</button>`}
        ${matchNote(r) && html`<span class="hint">${matchNote(r)}</span>`}</label>`;
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
