import { get, queryString, send } from "./api.js";
import { Message, Segmented, formatBp, html, useEffect, useState } from "./html.js";

const pollMilliseconds = 500;
const hitsPerPage = 50;
const lists = [["Ranked", "Ranked"], ["Identical", "Identical"], ["NotComparable", "Not comparable"]];

export function SearchPage({ params }) {
  const [mode, setMode] = useState(params.get("mode") ?? "similar");
  const [file, setFile] = useState(null);
  const [text, setText] = useState("");
  const [kind, setKind] = useState("");
  const [ratio, setRatio] = useState("0.5");
  const [top, setTop] = useState("");
  const [record, setRecord] = useState(params.get("record") ?? "");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState(null);

  async function submit(event) {
    event.preventDefault();
    setError(null);
    const pasted = text.trim();
    if (!file && !pasted) { setError("Choose a FASTA file or paste a sequence."); return; }
    setBusy(true);
    try {
      const body = file ?? (pasted.startsWith(">") ? pasted : `>query\n${pasted}`);
      const query = mode === "similar" ? { kind, ratio, top } : { kind: record ? "" : kind, record, top };
      const accepted = await send(`/api/searches/${mode}${queryString(query)}`, { method: "POST", body, contentType: "text/x-fasta" });
      location.hash = `#/searches/${accepted.id}`;
    } catch (e) {
      setError(e.message);
      setBusy(false);
    }
  }

  return html`<h1>Search</h1>
    <p class="lede">The query is profiled and discarded: a search keeps both strands' profile and hash for 30 days, never the DNA.</p>
    <form class="stack narrow" onSubmit=${submit}>
      <${Segmented} label="Question" value=${mode} onChange=${setMode}
        options=${[["similar", "Which records resemble this?"], ["localize", "Where does this fragment come from?"]]} />
      <label class="field"><span>FASTA file</span>
        <input type="file" accept=".fa,.fasta,.fna,.txt" onChange=${e => setFile(e.target.files[0] ?? null)} />
      </label>
      <label class="field"><span>Or paste one sequence</span>
        <textarea rows="5" spellCheck="false" disabled=${!!file} value=${text} onInput=${e => setText(e.target.value)}></textarea>
      </label>
      ${mode === "localize" && html`<label class="field"><span>Inside record (id, optional; required under 4 kbp)</span>
        <input type="number" min="1" value=${record} onInput=${e => setRecord(e.target.value)} />
      </label>`}
      ${(mode === "similar" || !record) && html`<label class="field"><span>Kind</span>
        <select value=${kind} onChange=${e => setKind(e.target.value)}>
          ${["", "bacterium", "archaeon", "virus", "plasmid"].map(k => html`<option key=${k} value=${k}>${k || "any"}</option>`)}
        </select>
      </label>`}
      ${mode === "similar" && html`<label class="field"><span>Length ratio (0.25–1)</span>
        <input type="number" min="0.25" max="1" step="0.05" value=${ratio} onInput=${e => setRatio(e.target.value)} />
      </label>`}
      <label class="field"><span>Hits to keep (${mode === "similar" ? "1–1000, default 50" : "1–100, default 20"})</span>
        <input type="number" min="1" max=${mode === "similar" ? 1000 : 100} value=${top} onInput=${e => setTop(e.target.value)} />
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

function Hits({ run }) {
  const isSimilar = run.kind === "Similar";
  return html`<div class="tablewrap"><table>
    <thead><tr><th class="num">#</th><th>Record</th><th class="num">Length</th><th class="num">Score</th>
      ${isSimilar ? html`<th class="num">Shared</th>` : html`<th class="num">Offset</th><th>Strand</th>`}
      <th class="num">E</th></tr></thead>
    <tbody>${run.hits.items.map(h => html`<tr key=${h.position}>
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
  if (s.candidates != null) parts.push(`${s.candidates.toLocaleString()} scanned`);
  if (s.neighbourRanked) parts.push("top hits ordered by shared neighbours");
  return html`<p class="hint">${parts.join(" · ")}</p>
    ${s.weakSignal && html`<${Message}>This query is short, so the ranking reflects composition more than order. Treat it as a hint and check it with an aligner.<//>`}`;
}

export function SearchRunPage({ id }) {
  const [list, setList] = useState("Ranked");
  const [after, setAfter] = useState(0);
  const [run, setRun] = useState(null);
  const [error, setError] = useState(null);
  const [version, setVersion] = useState(0);

  useEffect(() => {
    let live = true;
    let timer;
    async function poll() {
      try {
        const result = await get(`/api/searches/${id}${queryString({ list, after: after || "", limit: hitsPerPage })}`);
        if (!live) return;
        setRun(result);
        if (result.status === "Queued" || result.status === "Running") timer = setTimeout(poll, pollMilliseconds);
      } catch (e) {
        if (live) setError(e.message);
      }
    }
    poll();
    return () => { live = false; clearTimeout(timer); };
  }, [id, list, after, version]);

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
  const title = run.kind === "Similar" ? "Find similar" : "Localize";
  const seconds = run.completedAt ? (new Date(run.completedAt) - new Date(run.createdAt)) / 1000 : null;

  return html`<h1>${title}</h1>
    <p class="lede">${formatBp(run.queryLengthBp)} query · ${run.status}${seconds != null ? ` in ${seconds.toFixed(2)} s` : ""}</p>
    ${run.error && html`<${Message} error>${run.error}<//>`}
    ${(run.status === "Queued" || run.status === "Running") && html`<${Message}>Waiting for the worker…<//>`}
    ${run.status === "Expired" && html`<${Message}>The results of this search were dropped. Its parameters and summary remain.<//>`}
    <${Summary} run=${run} />
    ${run.hits && html`
      <${Segmented} label="List" value=${list} options=${lists} onChange=${v => { setList(v); setAfter(0); }} />
      ${run.hits.items.length ? html`<${Hits} run=${run} />` : html`<${Message}>This list is empty.<//>`}
      <div class="pager">
        <button type="button" disabled=${after === 0} onClick=${() => setAfter(Math.max(0, after - hitsPerPage))}>Previous</button>
        <button type="button" disabled=${run.hits.nextAfter == null} onClick=${() => setAfter(run.hits.nextAfter)}>Next</button>
      </div>`}
    ${run.status === "Done" && html`<p class="actions"><button type="button" onClick=${drop}>Delete results now</button></p>`}`;
}
