import { get, queryString, send } from "./api.js";
import { CurveCanvas, compareHref, maxCompared } from "./compare.js";
import { HeatStrip, genePath, loadGenes } from "./genes.js";
import { FilterBox, FormLayout, Help, HowTo, Message, PageHeader, PanelBar, Pager, Segmented, Workspace, formatBp, helpText, html, useEffect, useMemo, usePaged, useState } from "./html.js";

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

  const form = html`<form class="stack" onSubmit=${submit}>
      <${Segmented} label="Question" value=${mode} onChange=${setMode}
        help="Pick what you want to know; the fields below change with it. Hover a choice for what it expects."
        options=${[["similar", "Which records resemble this?", "A genome or a long fragment: records of similar length ranked by their Z-curve profiles."],
          ["localize", "Where does this fragment come from?", "A fragment: where in which records it fits best, on which strand."],
          ["protein", "Which genes encode this protein?", "A protein of at least 50 aa: the stored genes whose Z-curve matches what it implies."]]} />
      <label class="field"><span>FASTA file</span>
        <input type="file" accept=".fa,.fasta,.fna,.faa,.txt" onChange=${e => setFile(e.target.files[0] ?? null)} />
      </label>
      <label class="field"><span>${mode === "protein" ? "Or paste one protein (at least 50 aa)" : "Or paste one sequence"}<${Help}>${file
          ? "A file is chosen, so this box is off: clear the file to paste instead."
          : "Plain letters, or one record with its FASTA header line."}<//></span>
        <textarea rows="5" spellCheck="false" disabled=${!!file} value=${text} onInput=${e => setText(e.target.value)}
          placeholder=${mode === "protein" ? "MKV…" : "ACGT…"}></textarea>
      </label>
      ${mode === "localize" && html`<label class="field"><span>Inside record (optional)<${Help}>A record id from the catalogue. Fragments
          of 4 kbp and more are searched over every record at least as long; shorter ones only inside the record named here.<//></span>
        <input type="number" min="1" value=${record} onInput=${e => setRecord(e.target.value)} placeholder="required under 4 kbp" />
      </label>`}
      ${(mode === "similar" || !record) && html`<label class="field"><span>Kind<${Help}>Only records of this kind are searched.<//></span>
        <select value=${kind} onChange=${e => setKind(e.target.value)}>
          ${searchKinds.map(([k, text]) => html`<option key=${k} value=${k}>${text}</option>`)}
        </select>
      </label>`}
      ${mode === "protein" && html`<label class="field"><span>Genetic code<${Help}>How codons map to amino acids when the protein is turned
          into the codon composition it implies. Code 4 reads TGA as tryptophan, as Mycoplasma and its relatives do.<//></span>
        <select value=${code} onChange=${e => setCode(e.target.value)}>
          <option value="11">11, standard</option>
          <option value="4">4, Mycoplasma and relatives (TGA = Trp)</option>
        </select>
      </label>`}
      ${mode === "similar" && html`<label class="field"><span>Length ratio<${Help}>Only records whose length is between ratio × L and L ÷ ratio
          are compared, L being the query's length: 0.5 means from half to twice the query. From 0.25 to 1; lower widens the range.<//></span>
        <input type="number" min="0.25" max="1" step="0.05" value=${ratio} onInput=${e => setRatio(e.target.value)} />
      </label>`}
      <label class="field"><span>Hits to keep<${Help}>${mode === "localize" ? "1 to 100; 20 when left empty." : "1 to 1,000; 50 when left empty."}<//></span>
        <input type="number" min="1" max=${mode === "localize" ? 100 : 1000} value=${top} onInput=${e => setTop(e.target.value)}
          placeholder=${mode === "localize" ? "20" : "50"} />
      </label>
      <button class="primary" type="submit" disabled=${busy}>${busy ? "Sending…" : "Search"}</button>
    </form>
    ${error && html`<${Message} error>${error}<//>`}`;

  return html`<${PageHeader} title="Search" lede=${questions[mode]} />
    <${FormLayout} form=${form} guide=${html`<${SearchGuide} mode=${mode} />`} />`;
}

// The guide beside the search form: what the chosen question returns, what it accepts, and what is kept.
function SearchGuide({ mode }) {
  const kept = mode === "protein"
    ? "The protein is turned into its expected Z-curve windows and discarded: the search keeps those windows and a hash for 30 days, never the letters."
    : "The query is profiled and discarded: the search keeps both strands' profile and a hash for 30 days, never the DNA.";
  return html`<section class="panel"><h2>What you get</h2>${guides[mode].get}</section>
    <section class="panel"><h2>What it accepts</h2><ul>${guides[mode].accepts.map(a => html`<li key=${a}>${a}</li>`)}</ul></section>
    <section class="panel"><h2>What is kept</h2><p>${kept}</p>
      <p>The search waits in a queue; its page fills in by itself, and you can delete the results there at any time.</p></section>`;
}

const guides = {
  similar: {
    get: html`<p>Records ranked by how alike their Z-curve profiles are, each with a score and an E value saying how likely so good a hit
      is by chance. Two more lists hold records identical to the query and records that couldn't be compared.</p>`,
    accepts: ["One genome or long fragment, as a FASTA file or pasted.", "Bacteria, archaea, viruses and plasmids are searched; eukaryotes aren't.",
      "Short queries rank by composition more than by order: treat their hits as hints."],
  },
  localize: {
    get: html`<p>Where the fragment fits best: which records, at which offset and on which strand, best first.</p>`,
    accepts: ["A fragment of 4 kbp or more is searched over every record at least as long.",
      "A shorter fragment needs a record to search inside: give its id."],
  },
  protein: {
    get: html`<p>The stored genes whose Z-curve matches what the protein's codons would draw, with each gene's protein family where known, a
      score and an E value.</p>`,
    accepts: ["One protein of at least 50 amino acids.", "From 500 aa it is compared with the other long genes; shorter ones with genes within ±10% of its length.",
      "Genes are stored for BV-BRC genomes and NCBI-annotated eukaryotes."],
  },
};

const questions = {
  similar: "Paste or upload a genome, or a long piece of one, to rank the catalogue's records of similar length by how alike their Z-curve profiles are.",
  localize: "Paste or upload a fragment to find where it fits best: in which records, at which offset and on which strand.",
  protein: "Paste or upload a protein to find the stored genes that could encode it, scored by the Z-curve its codons would draw.",
};

const scoreHelp = "How alike the two profiles are, 1 being identical. Compare scores within one search: their scale depends on the window size used.";
const eHelp = "How many hits this good a search this size would find by chance among unrelated records. Well under 1: very likely real; " +
  "near 1 or above: treat as unconfirmed. “<” marks an upper bound.";
const geneEHelp = "How many genes of other protein families would score this high in a search this size, from real proteins scored against " +
  "the gene layer. Well under 1: the family is very likely right; near 1 or above: unconfirmed. “< E” is a bound: no more than a few " +
  "proteins of that sample met another family's gene this high.";

const sharedHelp = "For the top hits: how many of the query's 20 best hits are also among this hit's own 20 nearest records. The top 100 " +
  "are ordered by it, which favours hits that sit among the query's other relatives.";
const listHelp = "Ranked: the scored hits. Identical: records with exactly the query's sequence, on either strand (by hash, not by score). " +
  "Not comparable: records in the length range with no stored profile at the window size the query needs.";

const expected = log10 => {
  if (log10 == null) return "–";
  const value = 10 ** log10;
  return value >= 0.01 && value < 1000 ? value.toPrecision(2) : `1e${Math.round(log10)}`;
};

// The protein hits table; the hit list holds which hit is picked and draws it beside the table.
function GeneHits({ items, picked, setPicked }) {
  return html`<div class="tablewrap"><table>
    <thead><tr><th class="num">#</th><th>Record</th><th>Gene</th><th>Family</th>
      <th class="num">Score<${Help}>1 minus the mean distance, per window, between the protein's expected Z-curve and the gene's stored one,
        over the windows they share. 1 is a perfect match.<//></th>
      <th class="num">E<${Help}>${geneEHelp}<//></th></tr></thead>
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
  </table></div>`;
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
  return html`<div class="hitview"><h3 class="genetitle">${gene.gene.product ?? "gene"} <span class="mono">${gene.gene.familyId ?? ""}</span></h3>
    <p class="legend"><span><span class="dot" style=${{ "--swatch": "var(--p4)" }}></span>the protein's expected curve</span>
      <span><span class="dot" style=${{ "--swatch": "var(--p1)" }}></span>the stored gene in ${hit.sequence?.name ?? `record ${hit.sequenceId}`}</span></p>
    <${CurveCanvas} curves=${curves} unit="bases" label="The query protein's expected Z-curve over the hit gene's; drag to rotate"
      controls=${html`<${Segmented} label="Drift" value=${drift} onChange=${setDrift} help=${helpText.geneDrift} options=${[["kept", "Kept"], ["removed", "Removed"]]} />`} />
    <${HowTo}><p>Both curves are over codon positions 1 and 2. The protein's averages each amino acid's synonymous codons, weighted equally;
      the gene's comes from its stored windows. The score is 1 minus their mean per-window distance over the windows they share.</p>
      <p>The heat strips show the 9 numbers per window: blue positive, orange negative, darker further from zero.</p><//>
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
  return ids.length >= 2 && html`<a class="button" href=${compareHref(ids)}>Compare the top ${ids.length}</a>`;
}

// One list's hits, all loaded, paged and filtered in the page: by record name or accession, and for proteins by product, family or sequence.
// Protein hits open into a plot, so they're a workspace (docs/ui-layout.md G5): the list on one side, the picked hit drawn on the other.
function HitList({ id, run }) {
  const text = h => [h.sequence?.name, h.sequence?.recordId, h.product, h.familyId, h.geneSequence].filter(Boolean).join(" ");
  const paged = usePaged(run.hits.items, { size: hitsPerPage, text });
  const [picked, setPicked] = useState(null);
  const list = html`<${FilterBox} paged=${paged} label="Filter these hits"
      placeholder=${run.kind === "Protein" ? "record, product, family or sequence" : "record name or accession"} />
    <${Hits} id=${id} run=${run} items=${paged.rows} picked=${picked} setPicked=${setPicked} />
    ${paged.filtered.length === 0 && html`<${Message}>No hit matches the filter.<//>`}
    <${Pager} paged=${paged} noun="hits" />`;
  if (run.kind !== "Protein") return list;

  const onPage = picked && paged.rows.some(h => h.position === picked.position);
  const views = html`<section class="panel"><${PanelBar} title=${picked ? `Hit ${picked.position}` : "Picked hit"}>
      ${picked && html`<button type="button" class="quiet" onClick=${() => setPicked(null)}>Clear</button>`}<//>
    ${picked ? html`${!onPage && html`<p class="showing">This hit isn't on the list's current page.
        <button type="button" onClick=${() => paged.locate(h => h.position === picked.position)}>Show it in the list</button></p>`}
      <${HitView} key=${picked.position} id=${id} hit=${picked} />`
    : html`<p class="placeholder">Click a hit to draw the protein's expected Z-curve over the gene's.</p>`}</section>`;
  return html`<${Workspace} split="even" list=${html`<section class="panel">${list}</section>`} views=${views} />`;
}

function Hits({ run, items, picked, setPicked }) {
  if (run.kind === "Protein") return html`<${GeneHits} items=${items} picked=${picked} setPicked=${setPicked} />`;
  const isSimilar = run.kind === "Similar";
  return html`<div class="tablewrap"><table>
    <thead><tr><th class="num">#</th><th>Record</th><th class="num">Length</th><th class="num">Score<${Help}>${scoreHelp}<//></th>
      ${isSimilar ? html`<th class="num">Shared<${Help}>${sharedHelp}<//></th>`
        : html`<th class="num">Offset<${Help}>Where in the record the fragment fits best, in bases from the record's start.<//></th>
          <th>Strand<${Help}>Which strand of the record the fragment matched.<//></th>`}
      <th class="num">E<${Help}>${eHelp}<//></th></tr></thead>
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

  const lede = `${run.queryLengthAa != null ? `${run.queryLengthAa.toLocaleString()} aa` : formatBp(run.queryLengthBp)} query · ${run.status}` +
    (seconds != null ? ` in ${seconds.toFixed(2)} s` : "");
  const actions = html`${run.hits?.items.length > 0 && html`<${CompareTop} run=${run} />`}
    <a class="button" href=${`#/search?mode=${{ Similar: "similar", Localize: "localize", Protein: "protein" }[run.kind]}`}>New search</a>
    ${run.status === "Done" && html`<button type="button" onClick=${drop} title="Results are dropped after 30 days anyway">Delete results now</button>`}`;

  return html`<${PageHeader} title=${title} lede=${lede} actions=${actions} />
    ${run.error && html`<${Message} error>${run.error}<//>`}
    ${(run.status === "Queued" || run.status === "Running") && html`<${Message}>Waiting for the worker…<//>`}
    ${run.status === "Expired" && html`<${Message}>The results of this search were dropped. Its parameters and summary remain.<//>`}
    <${Summary} run=${run} />
    ${run.hits && html`
      ${run.kind !== "Protein" && html`<${Segmented} label="List" value=${list} options=${lists} onChange=${setList} help=${listHelp} />`}
      ${run.hits.items.length ? html`<${HitList} key=${list} id=${id} run=${run} />` : html`<${Message}>This list is empty.<//>`}`}`;
}
