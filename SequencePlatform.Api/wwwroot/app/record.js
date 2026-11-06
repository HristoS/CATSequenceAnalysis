import { get } from "./api.js";
import { ZCurvePath } from "./zcurve3d.js";
import { GeneLayerSignature, GenesPanel } from "./genes.js";
import { compareHref } from "./compare.js";
import { Message, Segmented, formatBp, formatGc, html, useEffect, useState } from "./html.js";

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
    <p class="hint">Nearest by ${result.signature}.</p>`;
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
  const views = html`<${Segmented} label="View" value=${view} onChange=${setView}
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

export function RecordPage({ id }) {
  const [detail, setDetail] = useState(null);
  const [error, setError] = useState(null);
  useEffect(() => {
    get(`/api/sequences/${id}`).then(setDetail).catch(e => setError(e.message));
  }, [id]);

  if (error) return html`<${Message} error>${error}<//>`;
  if (!detail) return html`<${Message}>Loading…<//>`;
  const s = detail.sequence;
  const hasSignature = s.kind !== "plasmid";
  const hasWindows = detail.levels.some(l => l.windowSize === 500 && l.methodVersion === "1");

  return html`<h1>${s.name || s.recordId}</h1>
    <p class="lede">${s.kind} · ${formatBp(s.lengthBp)} · GC ${formatGc(s.gcFraction)}${s.accession ? ` · ${s.accession}` : ""}${detail.owned ? " · your upload" : ""} · <a href=${compareHref([s.id])}>Compare with other records</a></p>
    <div class="columns">
      <section class="panel"><h2>Lineage</h2><${Lineages} detail=${detail} /></section>
      <section class="panel"><h2>Nearest relatives</h2>
        ${hasSignature ? html`<${Relatives} id=${s.id} kind=${s.kind} />` : html`<${Message}>Plasmids aren't compared by genome signature.<//>`}
      </section>
      ${hasSignature && html`<section class="panel"><h2>Dinucleotide signature</h2><${Signature} id=${s.id} /></section>`}
      ${hasWindows && html`<section class="panel"><h2>Z-curve path</h2>
        <${ZCurvePath} id=${s.id} chromosome=${s.kind === "bacterium" || s.kind === "archaeon"} /></section>`}
      ${s.provider === "BV-BRC" && html`<section class="panel"><h2>Genes</h2><${GenesPanel} id=${s.id} /></section>`}
      <section class="panel"><h2>Stored profile levels</h2>
        <div class="tablewrap"><table>
          <thead><tr><th class="num">Window</th><th class="num">Windows</th><th>Version</th></tr></thead>
          <tbody>${detail.levels.map(l => html`<tr key=${`${l.methodVersion}-${l.windowSize}`}>
            <td class="num">${formatBp(l.windowSize)}</td><td class="num">${l.windowCount.toLocaleString()}</td><td class="mono">${l.methodVersion}</td>
          </tr>`)}</tbody>
        </table></div>
        <p class="hint">Source ${detail.provenance.url ?? "upload"}, retrieved ${new Date(detail.provenance.retrievedAt).toLocaleDateString()}. The DNA itself was deleted after profiling.</p>
      </section>
    </div>`;
}
