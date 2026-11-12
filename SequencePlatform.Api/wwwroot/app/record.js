import { get } from "./api.js";
import { ZCurvePath } from "./zcurve3d.js";
import { GeneLayerSignature, GeneView, GenesPanel } from "./genes.js";
import { compareHref } from "./compare.js";
import { Help, Message, PageHeader, PanelBar, PlotCaption, PlotPair, Segmented, Workspace, usePlotSettings, formatBp, formatGc, helpText, html, useEffect, useState } from "./html.js";

const lineageRanks = ["superkingdom", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"];

function Lineages({ detail }) {
  const sources = [["Provider", detail.sequence.providerLineage]];
  if (detail.ncbi) sources.push(["NCBI", detail.ncbi.lineage]);
  if (detail.gtdb) sources.push([`GTDB ${detail.gtdb.release}`, detail.gtdb.lineage]);
  const shown = lineageRanks.filter(rank => sources.some(([, lineage]) => lineage[rank]));
  if (shown.length === 0) return html`<${Message}>No lineage is recorded.<//>`;
  return html`<div class="tablewrap"><table>
    <thead><tr><th>Rank</th>${sources.map(([name]) => html`<th key=${name}>${name}</th>`)}</tr></thead>
    <tbody>${shown.map(rank => html`<tr key=${rank}>
      <td>${rank}</td>${sources.map(([name, lineage]) => html`<td key=${name}>${lineage[rank] ?? "–"}</td>`)}
    </tr>`)}</tbody>
  </table></div>`;
}

function Relatives({ id, kind }) {
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  useEffect(() => {
    get(`/api/sequences/${id}/relatives?k=10`).then(setResult).catch(e => setError(e.message));
  }, [id]);

  if (error) return html`<${Message}>${error}<//>`;
  if (!result) return html`<${Message}>Finding relatives…<//>`;
  const ids = [id, ...result.items.map(r => r.sequence.id)];
  return html`<ol class="relatives">
      ${result.items.map(r => html`<li key=${r.sequence.id}>
        <a href="#/sequences/${r.sequence.id}">${r.sequence.name || r.sequence.recordId}</a>
        <span class="mono">${r.distance.toFixed(3)}</span>
      </li>`)}
    </ol>
    <p class="actions">
      <a class="button" href="#/tree?ids=${ids.join(",")}">Tree of these genomes</a>
      <a class="button" href=${compareHref(ids.slice(0, 3))}>Compare with the nearest two</a>
      <a class="button" href="#/map?group=${kind === "virus" ? "viral" : kind === "eukaryote" ? "eukaryotic" : "cellular"}&id=${id}">Show on the map</a>
    </p>
    <p class="hint">Nearest by ${result.signature}; the number is the distance, smaller is closer.</p>`;
}

function Signature({ id }) {
  const [view, setView] = useState("karlin");
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  useEffect(() => {
    if (view === "genes") return;
    setError(null);
    get(`/api/sequences/${id}/signature?view=${view}`).then(setResult).catch(e => setError(e.message));
  }, [id, view]);

  const peak = result ? Math.max(...result.values.map(v => Math.abs(v.value)), 1e-9) : 1;
  const views = html`<${Segmented} label="View" value=${view} onChange=${setView} help=${helpText.signature}
      options=${[["karlin", "Karlin ρ*"], ["zcurve", "Z-curve coefficients"], ["genes", "Gene layer"]]} />`;
  if (view === "genes") return html`${views}<${GeneLayerSignature} id=${id} />`;
  return html`${views}
    ${error && html`<${Message}>${error}<//>`}
    ${result && !error && html`<div class="bars">${result.values.map(v => html`<div class="bar" key=${v.name} title=${v.meaning}>
        <span class="mono">${v.name}</span>
        <div class="track"><i class=${v.value < 0 ? "negative" : "positive"}
          style=${{ width: `${(Math.abs(v.value) / peak) * 50}%` }}></i></div>
        <span class="mono num">${v.value.toFixed(3)}</span>
      </div>`)}</div>
      <p class="hint">The same 16 numbers in two coordinate systems: the complete second-order Z-curve is a rotation of log ρ*.</p>`}`;
}

// E2b's host guess for a phage, from the published findings (a panel only when there is one): the five nearest stored genera by 6-mer
// composition, the recorded host beside them, and how often the same first guess was right among phages with a recorded host.
function useHostGuess(id, kind) {
  const [finding, setFinding] = useState(null);
  useEffect(() => {
    setFinding(null);
    if (kind !== "virus") return;
    get(`/api/sequences/${id}/findings`).then(r => setFinding(r.items.find(f => f.experiment === "e2b-phage-host") ?? null)).catch(() => setFinding(null));
  }, [id, kind]);
  return finding;
}

function HostGuess({ finding }) {
  const p = finding.payload;
  const first = p.guesses[0];
  return html`<ol class="relatives">${p.guesses.map(g => html`<li key=${g.genus}>
      <span>${g.genus}${g.family && html` <span class="hint">${g.family}</span>`}${g.genus === p.host ? html` <span class="ok">✓ recorded host</span>` : ""}</span>
      <span class="mono">${g.distance.toFixed(4)}</span>
    </li>`)}</ol>
    ${p.host
      ? html`<p>Recorded host: <b>${p.host}</b>${p.hostFamily ? ` (${p.hostFamily})` : ""}, ranked ${p.rank.toLocaleString()} by 6-mers.</p>`
      : html`<p>No host is recorded for this ${p.class ?? "virus"}: these are guesses.</p>`}
    <p class="hint">When ${first.genus} came first for a phage with a recorded host, it was right ${p.firstGuess.right.toLocaleString()} of
      ${p.firstGuess.guessed.toLocaleString()} times${p.firstGuess.guessed === 0 ? " (never guessed first before)" : ""}.
      <a href="#/experiments/${finding.experiment}">How these guesses do overall</a>.</p>`;
}

export function RecordPage({ id }) {
  const [detail, setDetail] = useState(null);
  const [error, setError] = useState(null);
  const [gene, setGene] = useState(null);   // picked in the genes panel, highlighted on the path
  const plots = usePlotSettings("removed");   // View and Drift, shared by the path and the gene plot
  useEffect(() => {
    get(`/api/sequences/${id}`).then(setDetail).catch(e => setError(e.message));
  }, [id]);
  const hostGuess = useHostGuess(id, detail?.sequence.kind);

  if (error) return html`<${Message} error>${error}<//>`;
  if (!detail) return html`<${Message}>Loading…<//>`;
  const s = detail.sequence;
  const hasSignature = s.kind !== "plasmid";
  const hasWindows = detail.levels.some(l => l.windowSize === 500 && l.methodVersion === "1");
  const geneLayer = detail.annotations?.ncbi?.gene_layer;
  const hasGenes = s.provider === "BV-BRC" || geneLayer != null;
  const eukaryote = s.kind === "eukaryote";
  const chromosomeLevels = detail.annotations?.ncbi?.chromosome_levels;
  const hasChromosomes = chromosomeLevels?.status === "stored" && chromosomeLevels.sequences > 0;

  const lede = `${s.kind} · ${formatBp(s.lengthBp)} · GC ${formatGc(s.gcFraction)}${s.accession ? ` · ${s.accession}` : ""}${detail.owned ? " · your upload" : ""}`;
  const mapGroup = s.kind === "virus" ? "viral" : eukaryote ? "eukaryotic" : "cellular";
  const actions = html`<a class="button" href=${compareHref([s.id])}>Compare with other records</a>
    ${hasSignature && html`<a class="button" href="#/map?group=${mapGroup}&id=${s.id}">Show on the map</a>`}`;

  // The genome's path and the picked gene are the two plots of one card beside the gene list (UI-LAYOUT.md G5): picking a gene marks it on
  // the path and draws it, both in view, under one View and Drift.
  const hasPath = hasWindows || eukaryote;
  const pathSlot = hasPath && html`<${ZCurvePath} id=${s.id} chromosome=${s.kind === "bacterium" || s.kind === "archaeon"} hasWindows=${hasWindows}
    hasChromosomes=${hasChromosomes} hasGenes=${hasGenes} sequences=${geneLayer?.sequences ?? []} gene=${gene}
    viewSet=${plots.viewSet} drift=${plots.drift} />`;
  const geneSlot = gene
    ? html`<${PlotCaption} title=${`Gene #${gene.index + 1}: ${gene.product ?? gene.familyId ?? "no family"}`}>
        <button type="button" class="quiet" onClick=${() => setGene(null)}>Clear</button><//>
      <${GeneView} key=${gene.index} id=${s.id} index=${gene.index} viewSet=${plots.viewSet} drift=${plots.drift} />`
    : html`<${PlotCaption} title="Gene" /><div class="zpath empty">Click a gene in the list to draw its Z-curve here and mark it on the path.</div>`;
  const genesList = html`<section class="panel"><${PanelBar} title="Genes"
      help="The genes stored for this record (50 aa and up), each as Z-curve windows over codon positions 1 and 2." />
    <${GenesPanel} id=${s.id} selected=${gene?.index ?? null} onSelect=${setGene} /></section>`;

  return html`<${PageHeader} title=${s.name || s.recordId} lede=${lede} actions=${actions} />
    <div class="summarygrid">
      <section class="panel"><${PanelBar} title="Lineage" /><${Lineages} detail=${detail} /></section>
      <section class="panel"><${PanelBar} title="Nearest relatives" help=${"The 10 records whose whole-genome dinucleotide signature is closest to " +
          "this one's. That reflects composition, so close relatives come first, but distant ones can look alike."} />
        ${hasSignature ? html`<${Relatives} id=${s.id} kind=${s.kind} />` : html`<${Message}>Plasmids aren't compared by genome signature.<//>`}
      </section>
      ${hasSignature && html`<section class="panel"><${PanelBar} title="Dinucleotide signature" /><${Signature} id=${s.id} /></section>`}
      ${hostGuess && html`<section class="panel"><${PanelBar} title="Host guess" help=${"The five stored bacterial or archaeal genera whose genomes " +
          "are closest to this phage by 6-mer composition (both strands), closest first; the number is the distance. Composition points to a " +
          "host only loosely: on phages with a recorded host, the right genus comes first about one time in five."} />
        <${HostGuess} finding=${hostGuess} /></section>`}
    </div>
    ${hasGenes
      ? html`<${Workspace} list=${genesList} views=${html`<${PlotPair} toolbar=${plots.toolbar} left=${pathSlot || geneSlot} right=${pathSlot ? geneSlot : null} />`} />`
      : html`${pathSlot && html`<${PlotPair} title="Z-curve path" toolbar=${plots.toolbar} left=${pathSlot} />`}
        ${eukaryote && html`<section class="panel"><${PanelBar} title="Genes" />
          <${Message}>${detail.owned ? "No genes are stored for an upload: it carries no gene annotation, and the platform doesn't predict genes."
            : "No genes are stored for this genome: NCBI publishes no gene annotation for this assembly, and the platform doesn't predict genes."}<//></section>`}`}
    <section class="panel"><${PanelBar} title="Stored profile levels" help=${"The Z-curve profiles kept for this record, one per window " +
          "size: what searches compare. The DNA itself isn't kept."} />
        ${eukaryote && !detail.owned && html`<p class="hint">${chromosomeLevels?.status === "stored"
          ? `Chromosome levels: ${chromosomeLevels.sequences.toLocaleString()} nuclear sequences of at least ${formatBp(chromosomeLevels.min_length_bp)} ` +
            `(${formatBp(chromosomeLevels.sequences_bp)}), each at 500 bp and 32 kbp windows, drawn one at a time.` +
            (chromosomeLevels.left_out > 0 ? ` ${chromosomeLevels.left_out.toLocaleString()} shorter sequences (${formatBp(chromosomeLevels.left_out_bp)}) aren't kept.` : "")
          : chromosomeLevels?.status === "mismatch"
            ? "No chromosome levels: the genome file NCBI serves now isn't the one this record was made from, so nothing was stored from it."
            : "No chromosome levels yet: this genome hasn't been backfilled."}</p>`}
        ${detail.levels.length === 0
          ? html`<${Message}>${eukaryote ? "No whole-record levels: a eukaryote is one record per assembly with its nuclear dinucleotide counts " +
              "(the signature, relatives, map and trees); its windows are kept per chromosome instead, and it isn't searched as one genome." : "None are stored for this record."}<//>`
          : html`<div class="tablewrap"><table>
          <thead><tr><th class="num">Window</th><th class="num">Windows</th><th>Version</th></tr></thead>
          <tbody>${detail.levels.map(l => html`<tr key=${`${l.methodVersion}-${l.windowSize}`}>
            <td class="num">${formatBp(l.windowSize)}</td><td class="num">${l.windowCount.toLocaleString()}</td><td class="mono">${l.methodVersion}</td>
          </tr>`)}</tbody>
        </table></div>`}
        <p class="hint">Source ${detail.provenance.url ?? "upload"}, retrieved ${new Date(detail.provenance.retrievedAt).toLocaleDateString()}. The DNA itself was deleted after profiling.</p>
      </section>`;
}
