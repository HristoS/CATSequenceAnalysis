import { get, getBytes, loadGenes } from "./api.js";
import { Message, Segmented, formatBp, html, useEffect, useRef, useState } from "./html.js";

// The record's Z-curve as a 3-D path, rebuilt from its stored 500 bp level: each window's mean step times its size is an exact
// count of bases, so the running sums are the classic Z-curve at window resolution, with no DNA involved.
export const views = { "3d": [-0.65, 0.42], xy: [0, 0], xz: [0, Math.PI / 2], yz: [Math.PI / 2, 0] };
const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// ZCWP (docs/database-decision.md section 4.3): "ZCWP", uint16 version, int32 window size at 8, int32 window count at 12, then x, y and
// z as little-endian float32 blocks.
export function parse(buffer) {
  const view = new DataView(buffer);
  if (view.getUint32(0, true) !== 0x5057435a) throw new Error("The stored profile isn't in ZCWP format.");
  const size = view.getInt32(8, true), n = view.getInt32(12, true);
  const axis = a => Int32Array.from({ length: n }, (_, i) => Math.round(view.getFloat32(16 + 4 * (a * n + i), true) * size));
  return { size, n, x: axis(0), y: axis(1), z: axis(2) };
}

// Running sums at window boundaries, the least-squares line of each axis (Zhang's drift), and where the running GC skew (x - y counts
// 2 (G - C)) is lowest and highest: the replication origin and terminus of a bacterial chromosome.
export function path(profile) {
  const { n } = profile;
  const sums = ["x", "y", "z"].map(a => { const s = new Float64Array(n + 1); for (let i = 0; i < n; i++) s[i + 1] = s[i] + profile[a][i]; return s; });
  const fits = sums.map(s => {
    const m = s.length, sk = (m - 1) * m / 2, skk = (m - 1) * m * (2 * m - 1) / 6;
    let sp = 0, skp = 0;
    for (let k = 0; k < m; k++) { sp += s[k]; skp += k * s[k]; }
    const b = (m * skp - sk * sp) / (m * skk - sk * sk);
    return [(sp - b * sk) / m, b];
  });
  let skew = 0, min = 0, max = 0, minAt = 0, maxAt = 0;
  for (let i = 0; i < n; i++) {
    skew += profile.x[i] - profile.y[i];
    if (skew < min) { min = skew; minAt = i + 1; }
    if (skew > max) { max = skew; maxAt = i + 1; }
  }
  return { sums, fits, minAt, maxAt };
}

function draw(canvas, profile, shape, { drift, yaw, pitch, zoom, chromosome, highlights = [] }) {
  const n = profile.n, points = new Array(n + 1);
  for (let k = 0; k <= n; k++)
    points[k] = shape.sums.map((s, a) => (drift ? s[k] : s[k] - (shape.fits[a][0] + shape.fits[a][1] * k)));
  const { ctx, screen, width } = drawCurve(canvas, points, { yaw, pitch, zoom, highlights });
  if (chromosome) {
    ring(ctx, screen[shape.minAt], 6, css("--ink")); label(ctx, screen[shape.minAt], `skew min ${formatBp(shape.minAt * profile.size)}`, width);
    ring(ctx, screen[shape.maxAt], 4, css("--ink")); label(ctx, screen[shape.maxAt], `skew max ${formatBp(shape.maxAt * profile.size)}`, width);
  }
  return screen;
}

// A 3-D path ([x, y, z] points, in bases) coloured light to dark from start to end, with the axes and a scale bar; marks are point
// indices to ring (window boundaries, for a gene). Returns the canvas context and each point's screen position for marks and tooltips.
export function drawCurve(canvas, points, { yaw, pitch, zoom, marks = [], highlights = [] }) {
  const { ctx, screens, width } = drawCurves(canvas, [{ points, from: "--seq-a", to: "--seq-b", marks, highlights }], { yaw, pitch, zoom });
  return { ctx, screen: screens[0], width };
}

// Several paths on one set of axes and one scale (the compare page): each curve is { points, from, to, marks }, with from and to the
// CSS colour variables at its start and end; zoom is useZoom's { k, x, y }. Returns each curve's screen points.
export function drawCurves(canvas, curves, { yaw, pitch, zoom = unzoomed(), unit = "bases of excess" }) {
  const width = canvas.clientWidth, height = canvas.clientHeight, dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = css("--panel"); ctx.fillRect(0, 0, width, height);

  const all = curves.flatMap(c => c.points);
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const p of all) for (let a = 0; a < 3; a++) { lo[a] = Math.min(lo[a], p[a]); hi[a] = Math.max(hi[a], p[a]); }
  const centre = all.length ? lo.map((v, a) => (v + hi[a]) / 2) : [0, 0, 0];
  let radius = 0;
  for (const p of all) radius = Math.max(radius, Math.hypot(p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]));
  radius ||= 1;
  const scale = (Math.min(width, height) / 2 - 30) / radius * zoom.k;   // one base of excess is the same length on every axis
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const turn = v => { const x1 = cy * v[0] + sy * v[2], z1 = -sy * v[0] + cy * v[2]; return [x1, cp * v[1] - sp * z1]; };
  const screens = curves.map(c => c.points.map(p => {
    const [u, v] = turn([p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]]);
    return [width / 2 + zoom.x + u * scale, height / 2 + zoom.y - v * scale];
  }));

  // Coloured by position, from the curve's start colour to its end colour.
  ctx.lineWidth = 2; ctx.lineJoin = "round"; ctx.lineCap = "round";
  curves.forEach((c, i) => {
    const screen = screens[i], n = screen.length - 1, from = rgb(css(c.from)), to = rgb(css(c.to)), chunks = Math.min(96, n);
    for (let q = 0; q < chunks; q++) {
      const a = Math.floor(q * n / chunks), b = Math.floor((q + 1) * n / chunks), t = q / Math.max(1, chunks - 1);
      ctx.strokeStyle = `rgb(${from.map((v, k) => Math.round(v + (to[k] - v) * t)).join(",")})`;
      ctx.beginPath(); ctx.moveTo(...screen[a]);
      for (let k = a + 1; k <= b; k++) ctx.lineTo(...screen[k]);
      ctx.stroke();
    }
    for (const k of c.marks ?? []) ring(ctx, screen[k], 2.5, css("--ink-2"));
    if (screen.length) ring(ctx, screen[0], 3, css("--muted"));
  });

  // Highlighted stretches ([first, last] point indices) on top of every curve: thick, in the curve's colour over an ink halo, with a
  // ring where each starts, so even a one-window stretch on a whole chromosome shows.
  curves.forEach((c, i) => {
    const screen = screens[i];
    for (const [first, last] of c.highlights ?? []) {
      const a = Math.max(0, Math.min(first, screen.length - 1)), b = Math.max(a, Math.min(last, screen.length - 1));
      for (const [width, colour] of [[8, css("--ink")], [4, css(c.to)]]) {
        ctx.lineWidth = width; ctx.strokeStyle = colour;
        ctx.beginPath(); ctx.moveTo(...screen[a]);
        for (let k = a + 1; k <= b; k++) ctx.lineTo(...screen[k]);
        if (a === b) ctx.lineTo(screen[a][0] + 0.1, screen[a][1]);
        ctx.stroke();
      }
      ring(ctx, screen[a], 7, css("--ink"));
    }
  });
  ctx.lineWidth = 2;

  const ox = width - 46, oy = height - 42;
  ctx.font = `11px ${css("--mono")}`;
  for (const [name, v] of [["x", [1, 0, 0]], ["y", [0, 1, 0]], ["z", [0, 0, 1]]]) {
    const [u, w] = turn(v);
    ctx.strokeStyle = css("--muted"); ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox + u * 24, oy - w * 24); ctx.stroke();
    ctx.fillStyle = css("--ink-2"); ctx.fillText(name, ox + u * 31 - 3, oy - w * 31 + 4);
  }
  const bar = nice(radius / 3 / zoom.k);
  ctx.strokeStyle = css("--ink-2"); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(12, 16); ctx.lineTo(12 + bar * scale, 16); ctx.moveTo(12, 12); ctx.lineTo(12, 20);
  ctx.moveTo(12 + bar * scale, 12); ctx.lineTo(12 + bar * scale, 20); ctx.stroke();
  ctx.fillStyle = css("--ink-2"); ctx.fillText(`${bar.toLocaleString()} ${unit}`, 12, 33);
  return { ctx, screens, width };
}

function ring(ctx, p, size, colour) {
  ctx.beginPath(); ctx.arc(p[0], p[1], size, 0, 2 * Math.PI);
  ctx.strokeStyle = css("--panel"); ctx.lineWidth = 5; ctx.stroke();
  ctx.strokeStyle = colour; ctx.lineWidth = 2; ctx.stroke();
}

function label(ctx, p, text, width) {
  const w = ctx.measureText(text).width, x = Math.min(p[0] + 10, width - w - 6), y = p[1] - 9;
  ctx.globalAlpha = 0.85; ctx.fillStyle = css("--panel"); ctx.fillRect(x - 3, y - 11, w + 6, 15); ctx.globalAlpha = 1;
  ctx.fillStyle = css("--ink"); ctx.fillText(text, x, y);
}

// Zoom and pan for a path canvas: the wheel zooms about the cursor, shift-drag pans, a double-click resets, and so does a change in
// any of `resetOn` (new data, drift, view). zoom.current is { k, x, y }: the scale factor and the screen offset in CSS pixels. The
// caller's pointer handlers ask down/move first and skip their rotation when these return true.
export const unzoomed = () => ({ k: 1, x: 0, y: 0 });
export function useZoom(canvas, redraw, resetOn) {
  const zoom = useRef(unzoomed()), pan = useRef(null), latest = useRef(redraw);
  latest.current = redraw;
  useEffect(() => { zoom.current = unzoomed(); }, resetOn);
  useEffect(() => {
    const element = canvas.current;
    if (!element) return;
    const wheel = e => {
      e.preventDefault();
      const r = element.getBoundingClientRect(), mx = e.clientX - r.left - r.width / 2, my = e.clientY - r.top - r.height / 2;
      const z = zoom.current, k = Math.min(2000, Math.max(0.5, z.k * Math.exp(-e.deltaY * 0.0015)));
      zoom.current = { k, x: mx - (mx - z.x) * k / z.k, y: my - (my - z.y) * k / z.k };
      latest.current();
    };
    const reset = () => { zoom.current = unzoomed(); latest.current(); };
    element.addEventListener("wheel", wheel, { passive: false });
    element.addEventListener("dblclick", reset);
    return () => { element.removeEventListener("wheel", wheel); element.removeEventListener("dblclick", reset); };
  });
  return {
    zoom,
    down: e => {
      if (!e.shiftKey) return false;
      pan.current = { x: e.clientX, y: e.clientY, from: zoom.current };
      e.currentTarget.setPointerCapture(e.pointerId);
      return true;
    },
    move: e => {
      if (!pan.current) return false;
      const p = pan.current;
      zoom.current = { ...p.from, x: p.from.x + e.clientX - p.x, y: p.from.y + e.clientY - p.y };
      latest.current();
      return true;
    },
    up: () => { pan.current = null; },
  };
}

const rgb = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
const nice = v => { const p = 10 ** Math.floor(Math.log10(v)); const m = v / p; return (m >= 5 ? 5 : m >= 2 ? 2 : 1) * p; };

// The point ranges a gene covers on a window path: windows of `size` bp, the path starting `offsetBp` into the record (its origin, or
// 0). A gene across the path's start wraps into two ranges.
export function windowRanges(startBp, endBp, offsetBp, lengthBp, size, n) {
  const at = bp => ((bp - offsetBp) % lengthBp + lengthBp) % lengthBp;
  const a = Math.min(n, Math.floor(at(startBp) / size)), b = Math.min(n, Math.ceil(at(endBp) / size));
  return at(endBp) >= at(startBp) ? [[a, Math.max(a, b)]] : [[a, n], [0, b]];
}

// The gene-built path of /genes/path: its points, each gene's stretch by gene index, and the gene a point belongs to.
export function genePathOf(response) {
  const c = response.coordinates, starts = response.firstPoint, points = [];
  for (let i = 0; i < c.length; i += 3) points.push([c[i], c[i + 1], c[i + 2]]);
  const spans = new Map(response.geneIndex.map((g, i) => [g, [starts[i], i + 1 < starts.length ? starts[i + 1] : points.length - 1]]));
  const geneAt = k => {
    let lo = 0, hi = starts.length - 1;
    while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (starts[mid] < k) lo = mid; else hi = mid - 1; }
    return response.geneIndex[lo];
  };
  return { points, spans, geneAt };
}

// Each axis's least-squares line against the point index ([intercept, slope]): the drift removed from a path.
export function lineFits(points) {
  const n = points.length, mid = (n - 1) / 2;
  let sxx = 0;
  for (let k = 0; k < n; k++) sxx += (k - mid) ** 2;
  return [0, 1, 2].map(a => {
    let mean = 0, sxy = 0;
    for (const p of points) mean += p[a] / n;
    points.forEach((p, k) => { sxy += (k - mid) * (p[a] - mean); });
    const slope = sxx ? sxy / sxx : 0;
    return [mean - slope * mid, slope];
  });
}
const withoutDrift = (points, fits) => points.map((p, k) => p.map((v, a) => v - (fits[a][0] + fits[a][1] * k)));

// The record's Z-curve path: from its stored windows, or (Source: Genes) built from its stored genes alone, in chromosome order.
// A eukaryote is drawn one sequence at a time, never joined: its windows from a chromosome's own levels (`chromosomes`, migration 0022)
// and its gene path from the genes on one sequence (`sequences`, the gene layer's names); one picker, by sequence name, serves both.
// `gene`, the gene picked in the genes panel, is highlighted on either source, and picking one moves to its sequence.
const chromosomeOverviewBp = 20_000_000;   // longer chromosomes open at 32 kbp windows, a few thousand points

export function ZCurvePath({ id, chromosome, hasWindows = true, hasChromosomes = false, hasGenes, sequences = [], gene }) {
  const [profile, setProfile] = useState(null);
  const [genePath, setGenePath] = useState(null);
  const [geneList, setGeneList] = useState(null);
  const [chromosomes, setChromosomes] = useState(null);
  const [error, setError] = useState(null);
  const windowed = hasWindows || hasChromosomes;
  const [source, setSource] = useState(windowed ? "windows" : "genes");
  const [name, setName] = useState(null);           // the sequence drawn, by name, when the record has several
  const [resolution, setResolution] = useState(null);
  const [drift, setDrift] = useState("removed");
  const [view, setView] = useState("3d");
  const canvas = useRef(null), tip = useRef(null), shape = useRef(null), screen = useRef(null), angles = useRef(views["3d"]), drag = useRef(null);
  const fromGenes = source === "genes";

  // The names the picker offers for the current source, and the one drawn.
  const names = fromGenes ? sequences : hasChromosomes ? chromosomes?.map(c => c.name) ?? [] : [];
  const current = names.length === 0 ? null : names.includes(name) ? name : names[0];
  const shownChromosome = !fromGenes && hasChromosomes ? chromosomes?.find(c => c.name === current) : null;
  const size = shownChromosome ? resolution ?? (shownChromosome.lengthBp > chromosomeOverviewBp ? 32000 : 500) : 500;

  useEffect(() => {
    setProfile(null); setGenePath(null); setGeneList(null); setChromosomes(null); setError(null);
    setSource(windowed ? "windows" : "genes"); setName(null); setResolution(null);
    if (!hasChromosomes) return;
    let live = true;
    get(`/api/sequences/${id}/chromosomes`).then(list => { if (live) setChromosomes(list.chromosomes); }).catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [id, hasWindows, hasChromosomes]);

  // The windows drawn: the record's own 500 bp level, or one chromosome's level at the chosen size.
  const windowsUrl = fromGenes ? null : hasWindows ? `/api/sequences/${id}/profile?windowSize=500`
    : shownChromosome ? `/api/sequences/${id}/chromosomes/${shownChromosome.index}/profile?windowSize=${size}` : null;
  useEffect(() => {
    if (!windowsUrl) return;
    let live = true;
    setProfile(null);
    getBytes(windowsUrl).then(buffer => {
      if (!live) return;
      const parsed = parse(buffer);
      shape.current = path(parsed);
      setProfile({ ...parsed, url: windowsUrl });
    }).catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [windowsUrl]);

  useEffect(() => { if (gene && sequences.length > 0) setName(sequences[gene.sequence]); }, [gene]);
  const geneSequence = fromGenes && sequences.length > 0 ? sequences.indexOf(current) : null;
  useEffect(() => {
    if (!fromGenes || !hasGenes || (genePath && genePath.sequence === geneSequence)) return;
    let live = true;
    setGenePath(null);
    Promise.all([get(`/api/sequences/${id}/genes/path${geneSequence == null ? "" : `?sequence=${geneSequence}`}`), loadGenes(id)])
      .then(([p, list]) => { if (live) { setGenePath({ ...genePathOf(p), count: p.geneCount, sequence: geneSequence }); setGeneList(list); } })
      .catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [id, source, geneSequence, hasGenes]);
  // Gene list for highlighting on a chromosome's windows, where the picked gene is matched by its sequence's name.
  useEffect(() => { if (hasGenes && !geneList) loadGenes(id).then(setGeneList).catch(() => {}); }, [id, hasGenes]);

  const ready = fromGenes ? genePath : profile && profile.url === windowsUrl ? profile : null;
  const geneHere = gene && (sequences.length === 0 || sequences[gene.sequence] === (shownChromosome?.name ?? current));
  const redraw = () => {
    if (!canvas.current || !ready) return;
    const options = { yaw: angles.current[0], pitch: angles.current[1], zoom: zoomer.zoom.current };
    if (fromGenes) {
      const points = drift === "kept" ? genePath.points : withoutDrift(genePath.points, lineFits(genePath.points));
      const span = gene ? genePath.spans.get(gene.index) : null;
      screen.current = drawCurve(canvas.current, points, { ...options, highlights: span ? [span] : [] }).screen;
      return;
    }
    const length = shownChromosome?.lengthBp ?? profile.n * profile.size;
    const highlights = geneHere ? windowRanges(gene.startBp, gene.endBp, 0, length, profile.size, profile.n) : [];
    screen.current = draw(canvas.current, profile, shape.current, { ...options, drift: drift === "kept", chromosome, highlights });
  };
  const zoomer = useZoom(canvas, redraw, [ready, source, drift]);
  useEffect(() => { if (view !== "free") angles.current = views[view]; redraw(); }, [ready, source, drift, view, gene]);
  useEffect(() => {
    if (!canvas.current) return;
    const observer = new ResizeObserver(redraw);
    observer.observe(canvas.current);
    const scheme = window.matchMedia("(prefers-color-scheme: dark)");
    scheme.addEventListener("change", redraw);
    return () => { observer.disconnect(); scheme.removeEventListener("change", redraw); };
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
    if (!ready) return;
    const r = canvas.current.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
    let best = -1, bestDistance = 14 * 14;
    screen.current?.forEach((p, k) => { const d = (p[0] - mx) ** 2 + (p[1] - my) ** 2; if (d < bestDistance) { bestDistance = d; best = k; } });
    if (best < 0) { tip.current.hidden = true; return; }
    if (fromGenes) {
      const g = geneList?.genes[genePath.geneAt(best)];
      tip.current.textContent = g ? `gene #${g.index + 1} at ${formatBp(g.startBp)}${sequences.length ? ` of ${sequences[g.sequence]}` : ""} ` +
        `(${g.strand === "Reverse" ? "−" : "+"}) · ${g.product ?? g.familyId ?? "no family"}` : "";
    } else {
      const w = Math.min(best, profile.n - 1);
      tip.current.textContent = `${formatBp(best * profile.size)}${shownChromosome ? ` of ${shownChromosome.name}` : ""} · window GC ` +
        `${((1 + profile.z[w] / profile.size) / 2 * 100).toFixed(1)} % · G − C ${((profile.x[w] - profile.y[w]) / 2 / profile.size * 100).toFixed(1)} %`;
    }
    tip.current.hidden = false;
    const left = screen.current[best][0] + 12;
    tip.current.style.left = `${left + tip.current.offsetWidth > r.width - 6 ? screen.current[best][0] - tip.current.offsetWidth - 12 : left}px`;
    tip.current.style.top = `${Math.max(6, screen.current[best][1] - 14)}px`;
  };

  const noWindows = "No genome windows are stored for this record yet: eukaryote genomes were ingested with their dinucleotide counts only, " +
    "and its chromosome levels haven't been backfilled.";
  if (error) return html`<${Message}>${error}<//>`;
  if (!windowed && !hasGenes)
    return html`<${Message}>${noWindows} It has no gene annotation either (NCBI publishes none for this assembly), so there's no path to draw yet.<//>`;
  if (hasChromosomes && !chromosomes) return html`<${Message}>Loading the chromosomes…<//>`;
  const picked = gene && geneHere && (!fromGenes || genePath?.spans.has(gene.index)) &&
    `Highlighted: gene #${gene.index + 1}, ${gene.product ?? gene.familyId ?? "no family"}` + (fromGenes ? "." : `, ${formatBp(gene.startBp)}–${formatBp(gene.endBp)}.`);
  const elsewhere = gene && !picked && sequences.length > 0 && `The picked gene is on ${sequences[gene.sequence]}` +
    (fromGenes || chromosomes?.some(c => c.name === sequences[gene.sequence]) ? "." : ", which has no stored windows (shorter than 100 kbp).");
  return html`<div class="toolbar">
      ${windowed && hasGenes && html`<${Segmented} label="Source" value=${source} onChange=${setSource}
        options=${[["windows", hasChromosomes ? "Chromosome windows" : "Genome windows"], ["genes", "Genes"]]} />`}
      ${names.length > 1 && html`<label class="field"><span>Sequence</span>
        <select value=${current} onChange=${e => setName(e.target.value)}>
          ${names.map((n, i) => html`<option key=${n} value=${n}>${`${i + 1}. ${n}${!fromGenes && chromosomes ? ` · ${formatBp(chromosomes[i].lengthBp)}` : ""}`}</option>`)}
        </select></label>`}
      ${shownChromosome && html`<${Segmented} label="Windows" value=${String(size)} onChange=${v => setResolution(Number(v))}
        options=${[["500", "500 bp"], ["32000", "32 kbp"]]} />`}
      <${Segmented} label="Drift" value=${drift} onChange=${setDrift} options=${[["removed", "Removed"], ["kept", "Kept"]]} />
      <${Segmented} label="View" value=${view} onChange=${setView} options=${[["3d", "3-D"], ["xy", "x–y"], ["xz", "x–z"], ["yz", "y–z"]]} />
    </div>
    ${!ready && html`<${Message}>${fromGenes ? "Building the path from the stored genes…" : "Loading the stored windows…"}<//>`}
    <div class="zpath">
      <canvas ref=${canvas} role="img" aria-label=${fromGenes ? "Z-curve path built from the stored genes; drag to rotate"
        : `Z-curve path of ${profile?.n?.toLocaleString() ?? ""} windows; drag to rotate`}
        onPointerDown=${down} onPointerMove=${move} onPointerUp=${up} onPointerCancel=${up}
        onPointerLeave=${() => { if (!drag.current) tip.current.hidden = true; }}></canvas>
      <div class="ztip" ref=${tip} hidden></div>
    </div>
    ${picked && html`<p class="hint">${picked}</p>`}
    ${elsewhere && html`<p class="hint">${elsewhere}</p>`}
    ${!windowed && html`<p class="hint">${noWindows} Only the path built from its genes can be drawn.</p>`}
    ${fromGenes
      ? html`<p class="hint">${`Built from the gene layer alone: the ${genePath?.count?.toLocaleString() ?? ""} stored genes ` +
          (sequences.length ? `on ${current ?? "this sequence"}, one chromosome at a time (never joined), in order of start` : "in chromosome order") +
          ", each window one step of its means over codon positions 1 and 2, reverse-strand genes flipped onto the forward strand; the DNA " +
          "between genes and every third codon position are left out. " +
          (sequences.length ? "On a eukaryote the path follows which strand each gene is on, so between relatives it shows gene order (synteny) " +
            "rather than composition. " : "Relatives line up a little less well on it than on the genome windows (results/GenePath_Check.csv); ") +
          "A picked gene is an exact stretch of it. Hover a point to see its gene."}</p>`
      : html`<p class="hint">x purine − pyrimidine, y amino − keto, z strong − weak (GC). Running sums of the ${profile?.n?.toLocaleString() ?? ""} stored ${formatBp(size)} windows
          ${shownChromosome ? `of ${shownChromosome.name} (${formatBp(shownChromosome.lengthBp)}), one chromosome at a time; nuclear sequences under 100 kbp aren't kept,` : ""}
          light at the start and dark at the end; drag to rotate, scroll to zoom, shift-drag to pan, double-click to reset.${chromosome ? " The running GC skew (x − y) is lowest at the replication origin and highest at the terminus." : ""}${" "}
          Removing the drift subtracts each axis's best-fit line, which is the ${shownChromosome ? "chromosome" : "genome"}'s overall composition.</p>`}`;
}
