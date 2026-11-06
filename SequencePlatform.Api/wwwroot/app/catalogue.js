import { get, queryString } from "./api.js";
import { Message, formatBp, formatGc, html, useEffect, useState } from "./html.js";

const pageSize = 50;
const ranks = ["superkingdom", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"];
const empty = { name: "", kind: "", provider: "", rank: "genus", taxon: "", minLength: "", maxLength: "" };
const groupLabels = { prokaryote: "prokaryotes: bacteria and archaea" };
const counted = (label, count) => `${label} (${count.toLocaleString()})`;

// Kind and provider choices come from what the catalogue holds (/api/sequences/facets), with their record counts.
function choices(facets) {
  if (!facets) return { kinds: [["", "any"]], providers: [["", "any"]] };
  return {
    kinds: [["", counted("any", facets.total)],
      ...facets.kindGroups.map(g => [g.value, counted(groupLabels[g.value] ?? g.value, g.count)]),
      ...facets.kinds.map(k => [k.value, counted(k.value, k.count)])],
    providers: [["", "any"], ...facets.providers.map(p => [p.value, counted(p.value, p.count)])],
  };
}

export function CataloguePage() {
  const [draft, setDraft] = useState(empty);
  const [filter, setFilter] = useState(empty);
  const [starts, setStarts] = useState([0]);   // each visited page's "after", so Previous can go back
  const [page, setPage] = useState(null);
  const [error, setError] = useState(null);
  const [facets, setFacets] = useState(null);
  const [facetsError, setFacetsError] = useState(null);
  const after = starts[starts.length - 1];

  useEffect(() => {
    let live = true;
    get("/api/sequences/facets")
      .then(result => { if (live) setFacets(result); })
      .catch(e => { if (live) setFacetsError(e.message); });
    return () => { live = false; };
  }, []);
  const { kinds, providers } = choices(facets);

  useEffect(() => {
    let live = true;
    setError(null);
    const { rank, taxon, ...rest } = filter;
    get(`/api/sequences${queryString({ ...rest, rank: taxon ? rank : "", taxon, after: after || "", limit: pageSize })}`)
      .then(result => { if (live) setPage(result); })
      .catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [filter, after]);

  const change = name => event => setDraft({ ...draft, [name]: event.target.value });
  const apply = event => { event.preventDefault(); setStarts([0]); setFilter(draft); };

  return html`<h1>Catalogue</h1>
    <form class="filters" onSubmit=${apply}>
      <label class="field"><span>Name contains</span><input value=${draft.name} onInput=${change("name")} placeholder="e.g. Escherichia" /></label>
      <label class="field"><span>Kind</span>
        <select value=${draft.kind} onChange=${change("kind")}>${kinds.map(([v, text]) => html`<option key=${v} value=${v}>${text}</option>`)}</select>
      </label>
      <label class="field"><span>Source</span>
        <select value=${draft.provider} onChange=${change("provider")}>${providers.map(([v, text]) => html`<option key=${v} value=${v}>${text}</option>`)}</select>
      </label>
      <label class="field"><span>Rank</span>
        <select value=${draft.rank} onChange=${change("rank")}>${ranks.map(r => html`<option key=${r}>${r}</option>`)}</select>
      </label>
      <label class="field"><span>Taxon</span><input value=${draft.taxon} onInput=${change("taxon")} placeholder="e.g. Bacillus" /></label>
      <label class="field"><span>Min length (bp)</span><input type="number" min="0" value=${draft.minLength} onInput=${change("minLength")} /></label>
      <label class="field"><span>Max length (bp)</span><input type="number" min="0" value=${draft.maxLength} onInput=${change("maxLength")} /></label>
      <button class="primary" type="submit">Filter</button>
    </form>
    ${facetsError && html`<${Message}>The filter choices couldn't be loaded (${facetsError}); "any" still works.<//>`}
    ${error && html`<${Message} error>${error}<//>`}
    ${page && !error && html`<div class="tablewrap"><table>
        <thead><tr><th>Name</th><th>Kind</th><th>Genus</th><th class="num">Length</th><th class="num">GC</th><th>Accession</th></tr></thead>
        <tbody>${page.items.map(r => html`<tr key=${r.id}>
          <td><a href="#/sequences/${r.id}">${r.name || r.recordId}</a></td>
          <td>${r.kind}</td>
          <td>${r.providerLineage.genus ?? "–"}</td>
          <td class="num">${formatBp(r.lengthBp)}</td>
          <td class="num">${formatGc(r.gcFraction)}</td>
          <td class="mono">${r.accession ?? "upload"}</td>
        </tr>`)}</tbody>
      </table></div>
      ${page.items.length === 0 && html`<${Message}>Nothing matches these filters.<//>`}
      <div class="pager">
        <button type="button" disabled=${starts.length === 1} onClick=${() => setStarts(starts.slice(0, -1))}>Previous</button>
        <span>Page ${starts.length}</span>
        <button type="button" disabled=${page.nextAfter == null} onClick=${() => setStarts([...starts, page.nextAfter])}>Next</button>
      </div>`}`;
}
