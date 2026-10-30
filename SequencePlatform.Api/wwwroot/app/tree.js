import { get, queryString } from "./api.js";
import { Message, Segmented, formatBp, formatGc, html, useEffect, useRef, useState } from "./html.js";

const svgNamespace = "http://www.w3.org/2000/svg";
const ranks = ["phylum", "class", "order", "family", "genus", "species"];

function parse(newick) {
  let i = 0;
  function node() {
    const n = { children: [], len: 0, name: null };
    if (newick[i] === "(") { i++; n.children.push(node()); while (newick[i] === ",") { i++; n.children.push(node()); } i++; }
    let start = i;
    while (i < newick.length && !":,);".includes(newick[i])) i++;
    if (i > start) n.name = newick.slice(start, i);
    if (newick[i] === ":") { i++; start = i; while (i < newick.length && !",);".includes(newick[i])) i++; n.len = parseFloat(newick.slice(start, i)); }
    return n;
  }
  return node();
}

function element(name, attributes, parent) {
  const e = document.createElementNS(svgNamespace, name);
  for (const k in attributes) e.setAttribute(k, attributes[k]);
  if (parent) parent.appendChild(e);
  return e;
}

function leavesOf(n, out = []) {
  if (!n.children.length) out.push(n);
  n.children.forEach(c => leavesOf(c, out));
  return out;
}

const order = root => new Map(leavesOf(root).map((l, i) => [l.name, i]));

// Rotates children so each subtree sits where its leaves sit in the other tree's order; rotation never changes a tree's topology.
function rotate(n, position) {
  n.children.forEach(c => rotate(c, position));
  const mean = c => { const l = leavesOf(c); return l.reduce((s, x) => s + position.get(x.name), 0) / l.length; };
  n.children.sort((a, b) => mean(a) - mean(b));
}

function prepare(root) {
  (function walk(n, depth) {
    n.depth = depth;
    n.children.forEach(c => walk(c, depth + c.len));
    n.height = n.children.length ? 1 + Math.max(...n.children.map(c => c.height)) : 0;
  })(root, 0);
}

// A taxon is kept together when some branch separates exactly its genomes, judged on the unrooted tree.
function kept(root, leafNames, rankOf) {
  const index = new Map(leafNames.map((n, i) => [n, i])), sides = new Set();
  (function collect(node, isRoot) {
    const set = leavesOf(node).map(l => index.get(l.name)).sort((a, b) => a - b);
    if (!isRoot && set.length >= 2) {
      const members = new Set(set);
      sides.add(set.join(","));
      sides.add(leafNames.map((_, i) => i).filter(i => !members.has(i)).join(","));
    }
    node.children.forEach(c => collect(c, false));
  })(root, true);
  const groups = {};
  leafNames.forEach((name, i) => { const taxon = rankOf(name); if (taxon) (groups[taxon] ??= []).push(i); });
  // A taxon holding every leaf has no branch to be judged on, and a singleton is trivially together.
  const several = Object.entries(groups).filter(([, v]) => v.length >= 2 && v.length < leafNames.length);
  return { ok: several.filter(([, v]) => sides.has(v.join(","))).map(([k]) => k), total: several.length, sides };
}

// Normalized Robinson-Foulds distance: the share of non-trivial splits the two trees don't share.
function normalizedRf(a, b, leafNames) {
  const n = leafNames.length;
  const canonical = s => {
    const set = new Set(s.split(",").map(Number));
    return set.has(0) ? leafNames.map((_, i) => i).filter(i => !set.has(i)).join(",") : s;
  };
  const valid = s => { const c = s.split(",").length; return c >= 2 && c <= n - 2; };
  const left = new Set([...a].filter(valid).map(canonical)), right = new Set([...b].filter(valid).map(canonical));
  const shared = [...left].filter(s => right.has(s)).length;
  return left.size + right.size === 0 ? 0 : (left.size + right.size - 2 * shared) / (left.size + right.size);
}

const shortName = name => {
  const s = name.replace(/ (strain|ATCC|NCTC|str\.|isolate|subtype|genotype).*$/i, "");
  return s.length > 32 ? `${s.slice(0, 31)}…` : s;
};

function place(root, leaves, x0, direction, width, top, rowHeight, toScale) {
  const maxDepth = Math.max(...leaves.map(l => l.depth)) || 1;
  const scale = n => (toScale ? n.depth / maxDepth : (root.height - n.height) / root.height);
  leaves.forEach((l, i) => { l.y = top + i * rowHeight + rowHeight / 2; });
  (function position(n) {
    n.children.forEach(position);
    n.x = x0 + direction * scale(n) * width;
    if (n.children.length) n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2;
  })(root);
}

function drawBranches(n, parent) {
  if (!n.children.length) return;
  const ys = n.children.map(c => c.y);
  element("path", { d: `M${n.x},${Math.min(...ys)}V${Math.max(...ys)}`, class: "edge" }, parent);
  n.children.forEach(c => { element("path", { d: `M${n.x},${c.y}H${c.x}`, class: "edge" }, parent); drawBranches(c, parent); });
}

// One tree, or with `right` a tanglegram: the two trees mirrored, with a line joining each genome across them.
function render(svg, tip, meta, left, right, toScale) {
  const rowHeight = 22, top = 10, width = 1000, labelWidth = 196, gap = 8, x0 = 6;
  const treeWidth = right ? 200 : 780;
  if (right) { rotate(right, order(left)); rotate(left, order(right)); rotate(right, order(left)); }
  prepare(left);
  if (right) prepare(right);
  const leftLeaves = leavesOf(left), rightLeaves = right ? leavesOf(right) : [];
  const height = top + leftLeaves.length * rowHeight + 6;
  const leftLabel = x0 + treeWidth + gap, leftLink = leftLabel + labelWidth, rightLink = width - leftLink, rightLabel = width - leftLabel, rightTree = width - x0;

  svg.replaceChildren();
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.style.minWidth = right ? "900px" : "700px";

  const counts = {};
  Object.values(meta).forEach(m => { if (m.higher) counts[m.higher] = (counts[m.higher] || 0) + 1; });
  const tops = Object.entries(counts).sort((a, b) => b[1] - a[1]).slice(0, 3).map(([k]) => k);
  const colorOf = h => (tops.includes(h) ? `var(--p${tops.indexOf(h) + 1})` : "var(--p0)");

  place(left, leftLeaves, x0, 1, treeWidth, top, rowHeight, toScale);
  if (right) place(right, rightLeaves, rightTree, -1, treeWidth, top, rowHeight, toScale);

  const edges = element("g", {}, svg), links = element("g", {}, svg), rows = element("g", {}, svg);
  drawBranches(left, edges);
  if (right) drawBranches(right, edges);

  const linkOf = new Map(), rowsOf = new Map();
  if (right) {
    const rightY = new Map(rightLeaves.map(l => [l.name, l.y]));
    const middle = (leftLink + rightLink) / 2;
    for (const l of leftLeaves)
      linkOf.set(l.name, element("path", { d: `M${leftLink - 4},${l.y}C${middle},${l.y} ${middle},${rightY.get(l.name)} ${rightLink + 4},${rightY.get(l.name)}`, class: "link" }, links));
  }

  function row(l, side) {
    const m = meta[l.name], onLeft = side === "left", x = onLeft ? leftLabel : rightLabel;
    const g = element("g", { class: "row", tabindex: 0 }, rows);
    element("rect", { x: onLeft ? leftLabel - 6 : rightLink, y: l.y - rowHeight / 2, width: labelWidth + 6, height: rowHeight, class: "hit" }, g);
    if (toScale) element("line", { x1: l.x, x2: onLeft ? leftLabel - 8 : rightLabel + 8, y1: l.y, y2: l.y, class: "guide" }, g);
    element("circle", { cx: x, cy: l.y, r: 4, fill: colorOf(m.higher) }, g);
    element("text", { x: onLeft ? x + 9 : x - 9, y: l.y + 4, "text-anchor": onLeft ? "start" : "end" }, g).textContent = shortName(m.name);
    const same = m.genus ?? `id ${l.name}`;
    (rowsOf.get(same) || rowsOf.set(same, []).get(same)).push({ g, name: l.name });

    const show = event => {
      for (const r of rowsOf.get(same)) { r.g.classList.add("on"); linkOf.get(r.name)?.classList.add("on"); }
      tip.firstChild.textContent = m.name;
      tip.lastChild.textContent = `${m.genus ?? "no genus"} · ${m.higher ?? "–"}\nGC ${formatGc(m.gc)} · ${formatBp(m.lengthBp)}`;
      tip.hidden = false;
      const r = event?.clientX ? { left: event.clientX, top: event.clientY } : g.getBoundingClientRect();
      tip.style.left = `${Math.max(8, Math.min(r.left + 14, window.innerWidth - 290))}px`;
      tip.style.top = `${Math.min(r.top + 14, window.innerHeight - 90)}px`;
    };
    const hide = () => { tip.hidden = true; svg.querySelectorAll(".on").forEach(x => x.classList.remove("on")); };
    const open = () => { location.hash = `#/sequences/${l.name}`; };
    g.addEventListener("mousemove", show);
    g.addEventListener("mouseleave", hide);
    g.addEventListener("focus", () => show(null));
    g.addEventListener("blur", hide);
    g.addEventListener("click", open);
    g.addEventListener("keydown", event => { if (event.key === "Enter") open(); });
  }
  leftLeaves.forEach(l => row(l, "left"));
  rightLeaves.forEach(l => row(l, "right"));
  return tops.map((name, i) => ({ name, count: counts[name], swatch: `p${i + 1}` }));
}

function Score({ title, result, higherName }) {
  const number = ({ ok, total }, name) => html`<div>
    <b>${total ? `${ok.length} of ${total}` : "–"}</b>
    <small>${total ? `${name} kept together` : `no ${name} to test here`}</small>
  </div>`;
  return html`<section class="panel">
    <h2>${title}</h2>
    <div class="numbers">${number(result.genus, "genera")}${number(result.higher, higherName)}</div>
  </section>`;
}

export function TreePage({ params }) {
  const ids = params.get("ids") ?? "";
  const rank = params.get("rank") ?? "";
  const taxon = params.get("taxon") ?? "";
  const [draftRank, setDraftRank] = useState(rank || "genus");
  const [draftTaxon, setDraftTaxon] = useState(taxon);
  const [method, setMethod] = useState("upgma");
  const [compare, setCompare] = useState(false);
  const [toScale, setToScale] = useState(false);
  const [trees, setTrees] = useState(null);
  const [error, setError] = useState(null);
  const [scores, setScores] = useState(null);
  const [legend, setLegend] = useState([]);
  const svg = useRef(null), tip = useRef(null);

  useEffect(() => {
    if (!ids && !taxon) { setTrees(null); return undefined; }
    let live = true;
    setError(null);
    setTrees(null);
    const query = ids ? { ids } : { rank, taxon };
    Promise.all(["upgma", "nj"].map(m => get(`/api/trees${queryString({ ...query, method: m })}`)))
      .then(([upgma, nj]) => { if (live) setTrees({ upgma, nj }); })
      .catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [ids, rank, taxon]);

  useEffect(() => {
    if (!trees) return;
    const leaves = trees.upgma.leaves;
    const viral = leaves.every(l => l.kind === "virus");
    const meta = Object.fromEntries(leaves.map(l => [String(l.id), {
      name: l.name || l.recordId, genus: l.providerLineage.genus, higher: viral ? l.providerLineage.family : l.providerLineage.phylum,
      gc: l.gcFraction, lengthBp: l.lengthBp,
    }]));
    const names = Object.keys(meta);
    const left = parse(trees[compare ? "upgma" : method].newick), right = compare ? parse(trees.nj.newick) : null;
    setLegend(render(svg.current, tip.current, meta, left, right, toScale));

    const score = root => ({ genus: kept(root, names, n => meta[n].genus), higher: kept(root, names, n => meta[n].higher) });
    const scored = { higherName: viral ? "families" : "phyla", left: score(left), right: right && score(right) };
    scored.distance = right ? normalizedRf(scored.left.genus.sides, scored.right.genus.sides, names) : null;
    setScores(scored);
  }, [trees, compare, method, toScale]);

  const submit = event => {
    event.preventDefault();
    if (draftTaxon.trim()) location.hash = `#/tree${queryString({ rank: draftRank, taxon: draftTaxon.trim() })}`;
  };

  return html`<h1>Tree of genomes</h1>
    <p class="lede">Built from dinucleotide signatures (log ρ*, the complete second-order Z-curve). This is a composition tree: it tracks relatedness
      within families well and says little across phyla, so read it as similarity, not as phylogeny.</p>
    <form class="filters" onSubmit=${submit}>
      <label class="field"><span>Rank</span>
        <select value=${draftRank} onChange=${e => setDraftRank(e.target.value)}>${ranks.map(r => html`<option key=${r}>${r}</option>`)}</select>
      </label>
      <label class="field"><span>Taxon</span><input value=${draftTaxon} onInput=${e => setDraftTaxon(e.target.value)} placeholder="e.g. Bacillus" /></label>
      <button class="primary" type="submit">Build</button>
    </form>
    <div class="controls">
      <${Segmented} label="View" value=${compare ? "compare" : "single"} onChange=${v => setCompare(v === "compare")}
        options=${[["single", "One tree"], ["compare", "UPGMA against NJ"]]} />
      ${!compare && html`<${Segmented} label="Method" value=${method} onChange=${setMethod} options=${[["upgma", "UPGMA"], ["nj", "Neighbour-joining"]]} />`}
      <${Segmented} label="Branches" value=${toScale ? "lengths" : "topology"} onChange=${v => setToScale(v === "lengths")}
        options=${[["topology", "Topology only"], ["lengths", "To scale"]]} />
    </div>
    ${!ids && !taxon && html`<${Message}>Name a taxon, or pick genomes on the map or a record's relatives. A tree holds 3 to 500 genomes.<//>`}
    ${error && html`<${Message} error>${error}<//>`}
    ${(ids || taxon) && !trees && !error && html`<${Message}>Building the trees…<//>`}
    ${trees && !error && scores && html`<div class="scores">
      <${Score} title=${compare ? "UPGMA" : method === "nj" ? "Neighbour-joining" : "UPGMA"} result=${scores.left} higherName=${scores.higherName} />
      ${scores.right && html`<section class="panel"><h2>Disagreement</h2>
        <div class="numbers"><div><b>${Math.round(scores.distance * 100)} %</b><small>of splits differ between the two trees</small></div></div>
      </section>
      <${Score} title="Neighbour-joining" result=${scores.right} higherName=${scores.higherName} />`}
    </div>`}
    <figure class="treewrap" hidden=${!trees || !!error}>
      <div class="scroll"><svg ref=${svg} role="img" aria-label="Tree of genomes"></svg></div>
      <div class="legend">
        ${legend.map(l => html`<span key=${l.name}><i class="swatch ${l.swatch}"></i>${l.name} (${l.count})</span>`)}
        <span>Hover a genome to follow its genus; click to open it.</span>
      </div>
    </figure>
    <div class="tip" ref=${tip} hidden><b></b><span></span></div>`;
}
