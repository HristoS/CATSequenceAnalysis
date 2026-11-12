import { get, queryString, send } from "./api.js";
import { CurveCanvas, compareHref, maxCompared } from "./compare.js";
import { HeatStrip, genePath, loadGenes } from "./genes.js";
import { FilterBox, Message, Pager, Segmented, formatBp, html, useEffect, useMemo, usePaged, useState } from "./html.js";

const pollMilliseconds = 500;
const hitsPerPage = 50;
const maxHits = 1000;   // a search keeps at most 1,000 ranked hits (MaxTop), and one page may hold them all
const lists = [["Ranked", "Ranked"], ["Identical", "Identical"], ["NotComparable", "Not comparable"]];
// Eukaryotes are left out: they're stored without profile levels or genes, so no search can reach them.
const searchKinds = [["", "any"], ["prokaryote", "prokaryotes: bacteria and archaea"], ["bacterium", "bacterium"], ["archaeon", "archaeon"],
  ["virus", "virus"], ["plasmid", "plasmid"]];

export function SearchPage({ params }) {
  const [mode, setMode] = useState(params.get("mode") ?? "similar");
  const [file, setFile] = useState(null);
  const [text, setText] = useState("");
  const [kind, setKind] = useState("");
  const [ratio, setRatio] = useState("0.5");
  const [top, setTop] = useState("");
  const [record, setRecord] = useState(params.get("record") ?? "");
  const [code, setCode] = useState("11");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState(null);

  async function submit(event) {
    event.preventDefault();
    setError(null);
    const pasted = text.trim();
    if (!file && !pasted) { setError(`Choose a FASTA file or paste a ${mode === "protein" ? "protein" : "sequence"}.`); return; }
    setBusy(true);
    try {
      const body = file ?? (pasted.startsWith(">") ? pasted : `>query\n${pasted}`);
      const query = mode === "similar" ? { kind, ratio, top }
        : mode === "protein" ? { kind, code, top }
        : { kind: record ? "" : kind, record, top };
      const accepted = await send(`/api/searches/${mode}${queryString(query)}`, { method: "POST", body, contentType: "text/x-fasta" });
      location.hash = `#/searches/${accepted.id}`;
    } catch (e) {
      setError(e.message);
      setBusy(false);
    }
  }

  return html`<h1>Search</h1>
    <p class="lede">${mode === "protein"
      ? "The protein is turned into its expected Z-curve windows and discarded: a search keeps those windows and its hash for 30 days, never the letters."
      : "The query is profiled and discarded: a search keeps both strands' profile and hash for 30 days, never the DNA."}</p>
    <form class="stack narrow" onSubmit=${submit}>
      <${Segmented} label="Question" value=${mode} onChange=${setMode}
        options=${[["similar", "Which records resemble this?"], ["localize", "Where does this fragment come from?"],
          ["protein", "Which genes encode this protein?"]]} />
      <label class="field"><span>FASTA file</span>
        <input type="file" accept=".fa,.fasta,.fna,.faa,.txt" onChange=${e => setFile(e.target.files[0] ?? null)} />
      </label>
      <label class="field"><span>${mode === "protein" ? "Or paste one protein (at least 50 aa)" : "Or paste one sequence"}</span>
        <textarea rows="5" spellCheck="false" disabled=${!!file} value=${text} onInput=${e => setText(e.target.value)}></textarea>
      </label>
      ${mode === "localize" && html`<label class="field"><span>Inside record (id, optional; required under 4 kbp)</span>
        <input type="number" min="1" value=${record} onInput=${e => setRecord(e.target.value)} />
      </label>`}
      ${(mode === "similar" || !record) && html`<label class="field"><span>Kind</span>
        <select value=${kind} onChange=${e => setKind(e.target.value)}>
          ${searchKinds.map(([k, text]) => html`<option key=${k} value=${k}>${text}</option>`)}
        </select>
      </label>`}
      ${mode === "protein" && html`<label class="field"><span>Genetic code</span>
        <select value=${code} onChange=${e => setCode(e.target.value)}>
          <option value="11">11, standard</option>
          <option value="4">4, Mycoplasma and relatives (TGA = Trp)</option>
        </select>
      </label>`}
      ${mode === "similar" && html`<label class="field"><span>Length ratio (0.25–1)</span>
        <input type="number" min="0.25" max="1" step="0.05" value=${ratio} onInput=${e => setRatio(e.target.value)} />
      </label>`}
      <label class="field"><span>Hits to keep (${mode === "localize" ? "1–100, default 20" : "1–1000, default 50"})</span>
        <input type="number" min="1" max=${mode === "localize" ? 100 : 1000} value=${top} onInput=${e => setTop(e.target.value)} />
      </label>
      <button class="primary" type="submit" disabled=${busy}>${busy ? "Sending…" : "Search"}</button>
    </form>
    ${error && html`<${Message} error>${error}<//>`}`;
}

const expected = log10 => {
  if (log10 == null) return "–";
  const value = 10 ** log10;
  return value >= 0.01 && value < 1000 ? value.toPrecision(2) : `1e${Math.round(log10)}`;
};

function GeneHits({ id, run, items }) {
  const [picked, setPicked] = useState(null);
  return html`<div class="tablewrap"><table>
    <thead><tr><th class="num">#</th><th>Record</th><th>Gene</th><th>Family</th><th class="num">Score</th><th class="num">E</th></tr></thead>
    <tbody>${items.map(h => html`<tr key=${h.position} class=${picked?.position === h.position ? "row picked" : "row"} tabIndex="0"
        onClick=${() => setPicked(h)} onKeyDown=${e => { if (e.key === "Enter") setPicked(h); }}>
      <td class="num">${h.position}</td>
      <td>${h.sequence ? html`<a href="#/sequences/${h.sequenceId}" onClick=${e => e.stopPropagation()}>${h.sequence.name || h.sequence.recordId}</a>` : `deleted record ${h.sequenceId}`}</td>
      <td class="mono">${h.geneSequence ? `${h.geneSequence}: ` : ""}${(h.geneStartBp + 1).toLocaleString()}–${h.geneEndBp.toLocaleString()} ${h.strand === "Reverse" ? "(−)" : "(+)"}</td>
      <td>${h.product ?? "–"}${h.familyId && html` <span class="hint mono">${h.familyId}</span>`}</td>
      <td class="num mono">${h.score == null ? "–" : h.score.toFixed(3)}</td>
      <td class="num mono" title=${h.pKind === "Bound" ? "An upper bound: fewer than 10 protein-gene pairs of the null sample scored this high" : ""}>${
        h.pKind === "Bound" ? `< ${expected(h.log10E)}` : expected(h.log10E)}</td>
    </tr>`)}</tbody>
  </table></div>
  <p class="hint">E: how many genes of other protein families would score this high in a search this size, from real proteins scored against
    the gene layer (build-protein-null). ${"\"< E\""} is a bound: no more than a few proteins of that sample met another family's gene this
    high. Well under 1, the family is very likely right; near 1 or above, treat the hit as unconfirmed.</p>
  ${picked ? html`<${HitView} key=${picked.position} id=${id} hit=${picked} />`
    : html`<p class="hint">Pick a hit to draw the protein's expected Z-curve over the gene's.</p>`}`;
}

// The query protein's expected windows against one hit gene's stored windows: the comparison the score summarizes.
function HitView({ id, hit }) {
  const [query, setQuery] = useState(null);
  const [gene, setGene] = useState(null);
  const [error, setError] = useState(null);
  const [drift, setDrift] = useState("kept");
  useEffect(() => {
    let live = true;
    const fail = e => { if (live) setError(e.message); };
    get(`/api/searches/${id}/query`).then(q => { if (live) setQuery(q); }).catch(fail);
    loadGenes(hit.sequenceId).then(list => {
      const g = list.genes.find(x => x.startBp === hit.geneStartBp && x.endBp === hit.geneEndBp && x.strand === hit.strand);
      if (!g) throw new Error("That gene is no longer stored.");
      return get(`/api/sequences/${hit.sequenceId}/genes/${g.index}`);
    }).then(g => { if (live) setGene(g); }).catch(fail);
    return () => { live = false; };
  }, [id, hit.position]);

  const curves = useMemo(() => {
    if (!query || !gene) return [];
    return [[query, "--p4", `query protein, ${query.lengthAa} aa`], [gene, "--p1", `${gene.gene.product ?? "hit gene"}, ${gene.gene.residues} aa`]].map(([w, colour, label]) => {
      const shape = genePath(w.windows, w.windowPositions);
      const points = drift === "kept" ? shape.points
        : shape.points.map((p, i) => p.map((v, a) => v - (shape.fits[a][0] + shape.fits[a][1] * shape.at[i])));
      return { points, from: colour, to: colour, marks: shape.marks, label, describe: i => `about amino acid ${Math.round(shape.at[i] / 2).toLocaleString()}` };
    });
  }, [query, gene, drift]);

  if (error) return html`<${Message}>${error}<//>`;
  if (!query || !gene) return html`<${Message}>Loading the windows…<//>`;
  const peak = Math.max(...query.windows.flat().concat(gene.windows.flat()).map(Math.abs), 1e-9);
  return html`<div class="hitview"><h3 class="genetitle">Hit ${hit.position}: ${gene.gene.product ?? "gene"} <span class="mono">${gene.gene.familyId ?? ""}</span></h3>
    <p class="legend"><span><span class="dot" style=${{ "--swatch": "var(--p4)" }}></span>the protein's expected curve</span>
      <span><span class="dot" style=${{ "--swatch": "var(--p1)" }}></span>the stored gene in ${hit.sequence?.name ?? `record ${hit.sequenceId}`}</span></p>
    <div class="toolbar"><${Segmented} label="Drift" value=${drift} onChange=${setDrift} options=${[["kept", "Kept"], ["removed", "Removed"]]} /></div>
    <${CurveCanvas} curves=${curves} unit="bases" label="The query protein's expected Z-curve over the hit gene's; drag to rotate" />
    <p class="hint">Both over codon positions 1 and 2. The protein's curve averages each amino acid's synonymous codons, weighted equally;
      the gene's comes from its stored windows. The score is 1 minus their mean per-window distance over the windows they share.</p>
    <div class="heatgrid">
      <div><p class="heatlabel"><span class="dot" style=${{ "--swatch": "var(--p4)" }}></span>${`query, ${query.lengthAa} aa`}</p>
        <${HeatStrip} windows=${query.windows} peak=${peak} /></div>
      <div><p class="heatlabel"><span class="dot" style=${{ "--swatch": "var(--p1)" }}></span>${`hit gene, ${gene.gene.residues} aa`}</p>
        <${HeatStrip} windows=${gene.windows} peak=${peak} /></div>
    </div></div>`;
}

// Up to three distinct records from the current page of hits, for the compare page.
function CompareTop({ run }) {
  const ids = [...new Set(run.hits.items.filter(h => h.sequence).map(h => h.sequenceId))].slice(0, maxCompared);
  return ids.length >= 2 && html`<p class="actions"><a class="button" href=${compareHref(ids)}>Compare the top ${ids.length}</a></p>`;
}

// One list's hits, all loaded, paged and filtered in the page: by record name or accession, and for proteins by product, family or sequence.
function HitList({ id, run }) {
  const text = h => [h.sequence?.name, h.sequence?.recordId, h.product, h.familyId, h.geneSequence].filter(Boolean).join(" ");
  const paged = usePaged(run.hits.items, { size: hitsPerPage, text });
  return html`<${FilterBox} paged=${paged} label="Filter these hits"
      placeholder=${run.kind === "Protein" ? "record, product, family or sequence" : "record name or accession"} />
    <${Hits} id=${id} run=${run} items=${paged.rows} />
    ${paged.filtered.length === 0 && html`<${Message}>No hit matches the filter.<//>`}
    <${Pager} paged=${paged} noun="hits" />`;
}

function Hits({ id, run, items }) {
  if (run.kind === "Protein") return html`<${GeneHits} id=${id} run=${run} items=${items} />`;
  const isSimilar = run.kind === "Similar";
  return html`<div class="tablewrap"><table>
    <thead><tr><th class="num">#</th><th>Record</th><th class="num">Length</th><th class="num">Score</th>
      ${isSimilar ? html`<th class="num">Shared</th>` : html`<th class="num">Offset</th><th>Strand</th>`}
      <th class="num">E</th></tr></thead>
    <tbody>${items.map(h => html`<tr key=${h.position}>
      <td class="num">${h.position}</td>
      <td>${h.sequence ? html`<a href="#/sequences/${h.sequenceId}">${h.sequence.name || h.sequence.recordId}</a>` : `deleted record ${h.sequenceId}`}</td>
      <td class="num">${h.sequence ? formatBp(h.sequence.lengthBp) : "–"}</td>
      <td class="num mono">${h.score == null ? "–" : h.score.toFixed(3)}</td>
      ${isSimilar
        ? html`<td class="num">${h.sharedNeighbours ?? "–"}</td>`
        : html`<td class="num">${h.offsetBp == null ? "–" : h.offsetBp.toLocaleString()}</td><td>${h.strand ?? "–"}</td>`}
      <td class="num mono" title=${h.pKind === "Bound" ? "An upper bound: the null sample has too few scores this high" : ""}>${expected(h.log10E)}</td>
    </tr>`)}</tbody>
  </table></div>`;
}

function Summary({ run }) {
  const s = run.summary;
  if (!s) return null;
  const parts = [];
  if (s.windowSize) parts.push(`${formatBp(s.windowSize)} windows`);
  if (s.bandMinBp != null) parts.push(`candidates ${formatBp(s.bandMinBp)}${s.bandMaxBp != null ? `–${formatBp(s.bandMaxBp)}` : " and up"}`);
  if (s.rankedTotal != null) parts.push(`${s.rankedTotal.toLocaleString()} ranked`);
  if (s.candidates != null) parts.push(`${s.candidates.toLocaleString()} ${run.kind === "Protein" ? "genes scored" : "scanned"}`);
  if (s.neighbourRanked) parts.push("top hits ordered by shared neighbours");
  return html`<p class="hint">${parts.join(" · ")}</p>
    ${s.weakSignal && html`<${Message}>This query is short, so the ranking reflects composition more than order. Treat it as a hint
      and check it with an aligner.<//>`}`;
}

export function SearchRunPage({ id }) {
  const [list, setList] = useState("Ranked");
  const [run, setRun] = useState(null);
  const [error, setError] = useState(null);
  const [version, setVersion] = useState(0);

  useEffect(() => {
    let live = true;
    let timer;
    async function poll() {
      try {
        const result = await get(`/api/searches/${id}${queryString({ list, limit: maxHits })}`);
        if (!live) return;
        setRun(result);
        if (result.status === "Queued" || result.status === "Running") timer = setTimeout(poll, pollMilliseconds);
      } catch (e) {
        if (live) setError(e.message);
      }
    }
    poll();
    return () => { live = false; clearTimeout(timer); };
  }, [id, list, version]);

  async function drop() {
    try {
      await send(`/api/searches/${id}`, { method: "DELETE" });
      setVersion(version + 1);
    } catch (e) {
      setError(e.message);
    }
  }

  if (error) return html`<${Message} error>${error}<//>`;
  if (!run) return html`<${Message}>Loading…<//>`;
  const title = { Similar: "Find similar", Localize: "Localize", Protein: "Protein search" }[run.kind];
  const seconds = run.completedAt ? (new Date(run.completedAt) - new Date(run.createdAt)) / 1000 : null;

  return html`<h1>${title}</h1>
    <p class="lede">${run.queryLengthAa != null ? `${run.queryLengthAa.toLocaleString()} aa` : formatBp(run.queryLengthBp)} query · ${run.status}${seconds != null ? ` in ${seconds.toFixed(2)} s` : ""}</p>
    ${run.error && html`<${Message} error>${run.error}<//>`}
    ${(run.status === "Queued" || run.status === "Running") && html`<${Message}>Waiting for the worker…<//>`}
    ${run.status === "Expired" && html`<${Message}>The results of this search were dropped. Its parameters and summary remain.<//>`}
    <${Summary} run=${run} />
    ${run.hits && html`
      ${run.kind !== "Protein" && html`<${Segmented} label="List" value=${list} options=${lists} onChange=${setList} />`}
      ${run.hits.items.length ? html`<${CompareTop} run=${run} /><${HitList} key=${list} id=${id} run=${run} />` : html`<${Message}>This list is empty.<//>`}`}
    ${run.status === "Done" && html`<p class="actions"><button type="button" onClick=${drop}>Delete results now</button></p>`}`;
}
