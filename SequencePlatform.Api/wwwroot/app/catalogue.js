import { get, queryString } from "./api.js";
import { compareHref, maxCompared } from "./compare.js";
import { Help, Message, PageHeader, RangeField, formatBp, formatGc, html, recordKinds, useEffect, useState } from "./html.js";

const pageSize = 50;
const ranks = ["superkingdom", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"];
const empty = { name: "", kinds: [], provider: "", rank: "genus", taxon: "", minLength: "", maxLength: "" };
const counted = (label, count) => `${label} (${count.toLocaleString()})`;

// Kind and provider choices come from what the catalogue holds (/api/sequences/facets), with their record counts: the kinds a record can
// be (the upload page offers the same), then any other kind stored.
function choices(facets) {
  const count = new Map((facets?.kinds ?? []).map(k => [k.value, k.count]));
  const extra = [...count.keys()].filter(k => !recordKinds.some(([v]) => v === k)).map(k => [k, k]);
  return {
    kinds: [...recordKinds, ...extra].map(([value, label]) => ({ value, label, count: facets ? count.get(value) ?? 0 : null })),
    providers: [["", "any"], ...(facets?.providers ?? []).map(p => [p.value, counted(p.value, p.count)])],
  };
}

export function CataloguePage() {
  const [draft, setDraft] = useState(empty);
  const [filter, setFilter] = useState(empty);
  // Each visited page's "after" (the server pages by key, so only pages already reached can be jumped to) and which of them is shown.
  const [starts, setStarts] = useState([0]);
  const [at, setAt] = useState(0);
  const [page, setPage] = useState(null);
  const [error, setError] = useState(null);
  const [facets, setFacets] = useState(null);
  const [facetsError, setFacetsError] = useState(null);
  const [chosen, setChosen] = useState([]);   // records picked for the compare page, kept across pages and filters
  const after = starts[at];

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
    const { rank, taxon, kinds: picked, ...rest } = filter;
    get(`/api/sequences${queryString({ ...rest, kind: picked.join(","), rank: taxon ? rank : "", taxon, after: after || "", limit: pageSize })}`)
      .then(result => { if (live) setPage(result); })
      .catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [filter, after]);

  const change = name => event => setDraft({ ...draft, [name]: event.target.value });
  const apply = event => { event.preventDefault(); setStarts([0]); setAt(0); setFilter(draft); };
  const reset = () => { setDraft(empty); setStarts([0]); setAt(0); setFilter(empty); };
  const next = () => {
    if (at + 1 < starts.length && starts[at + 1] === page.nextAfter) setAt(at + 1);
    else { setStarts([...starts.slice(0, at + 1), page.nextAfter]); setAt(at + 1); }
  };
  // Ticking a kind applies at once, with the other filters as last applied.
  const toggleKind = value => {
    const picked = draft.kinds.includes(value) ? draft.kinds.filter(k => k !== value) : [...draft.kinds, value];
    setDraft({ ...draft, kinds: picked });
    setStarts([0]); setAt(0); setFilter({ ...filter, kinds: picked });
  };
  const filtered = Object.keys(empty).some(k => k !== "rank" && JSON.stringify(filter[k]) !== JSON.stringify(empty[k]));

  return html`<${PageHeader} title="Catalogue"
      lede="Every record the platform holds. Open one for its lineage, relatives and Z-curve, or tick two or three to compare them." />
    <fieldset class="kindpicks">
      <legend>Kinds<${Help}>Tick the kinds to list; with none ticked, every kind is listed. The numbers are how many records of each kind
        you can see.<//></legend>
      ${kinds.map(k => html`<label key=${k.value} class="kindpick">
        <input type="checkbox" checked=${draft.kinds.includes(k.value)} onChange=${() => toggleKind(k.value)} />
        ${k.label}${k.count != null && html`<span class="hint">${k.count.toLocaleString()}</span>`}</label>`)}
    </fieldset>
    <form class="filterbar" onSubmit=${apply}>
      <label class="field"><span>Name contains</span><input value=${draft.name} onInput=${change("name")} placeholder="e.g. Escherichia" /></label>
      <label class="field"><span>Source<${Help}>Where the record came from: a public database, or an upload.<//></span>
        <select value=${draft.provider} onChange=${change("provider")}>${providers.map(([v, text]) => html`<option key=${v} value=${v}>${text}</option>`)}</select>
      </label>
      <div class="pair">
        <label class="field"><span>Rank<${Help}>The level of the lineage the taxon is matched at. Only used when a taxon is given.<//></span>
          <select value=${draft.rank} onChange=${change("rank")}>${ranks.map(r => html`<option key=${r}>${r}</option>`)}</select>
        </label>
        <label class="field"><span>Taxon<${Help}>The exact name at the chosen rank, as the provider's lineage spells it: e.g. genus Bacillus, or family Enterobacteriaceae.<//></span>
          <input value=${draft.taxon} onInput=${change("taxon")} placeholder="e.g. Bacillus" /></label>
      </div>
      <${RangeField} label="Length (bp)" min=${draft.minLength} max=${draft.maxLength}
        onMin=${v => setDraft(d => ({ ...d, minLength: v }))} onMax=${v => setDraft(d => ({ ...d, maxLength: v }))} />
      <div class="actions formbuttons">
        <button class="primary" type="submit">Filter</button>
        ${filtered && html`<button type="button" onClick=${reset}>Clear filters</button>`}
      </div>
    </form>
    ${facetsError && html`<${Message}>The filter choices couldn't be loaded (${facetsError}); "any" still works.<//>`}
    ${error && html`<${Message} error>${error}<//>`}
    ${page && !error && html`<div class="tablewrap"><table>
        <thead><tr><th><span class="visually-hidden">Compare</span></th><th>Name</th><th>Kind</th><th>Genus</th><th class="num">Length</th><th class="num">GC</th><th>Accession</th></tr></thead>
        <tbody>${page.items.map(r => {
          const on = chosen.some(c => c.id === r.id);
          return html`<tr key=${r.id}>
          <td><input type="checkbox" checked=${on} disabled=${!on && chosen.length >= maxCompared} aria-label=${`Compare ${r.name || r.recordId}`}
            title=${!on && chosen.length >= maxCompared ? `Up to ${maxCompared} records` : "Choose for the compare page"}
            onChange=${() => setChosen(on ? chosen.filter(c => c.id !== r.id) : [...chosen, { id: r.id, name: r.name || r.recordId }])} /></td>
          <td><a href="#/sequences/${r.id}">${r.name || r.recordId}</a></td>
          <td>${r.kind}</td>
          <td>${r.providerLineage.genus ?? "–"}</td>
          <td class="num">${formatBp(r.lengthBp)}</td>
          <td class="num">${formatGc(r.gcFraction)}</td>
          <td class="mono">${r.accession ?? "upload"}</td>
        </tr>`;
        })}</tbody>
      </table></div>
      ${page.items.length === 0 && html`<${Message}>Nothing matches these filters.<//>`}
      <div class="pager">
        <button type="button" disabled=${at === 0} onClick=${() => setAt(0)} title="First page" aria-label="First page">«</button>
        <button type="button" disabled=${at === 0} onClick=${() => setAt(at - 1)}>‹ Previous</button>
        <span class="pagelabel">Page
          <select class="pagechoice" aria-label="Page" value=${at} onChange=${e => setAt(Number(e.target.value))}>
            ${starts.map((_, i) => html`<option key=${i} value=${i}>${i + 1}</option>`)}
          </select></span>
        <button type="button" disabled=${page.nextAfter == null} onClick=${next}>Next ›</button>
        <span class="hint">${page.nextAfter == null ? "the last page" : `${pageSize} records a page; the list holds the pages visited so far`}</span>
      </div>`}
    ${chosen.length > 0 && html`<div class="chosen tray" role="region" aria-label="Records chosen to compare">
      <span>${chosen.length} of ${maxCompared} chosen: ${chosen.map(c => c.name).join(", ")}</span>
      ${chosen.length >= 2 ? html`<a class="button primary" href=${compareHref(chosen.map(c => c.id))}>Compare</a>` : html`<span class="hint">choose one more to compare</span>`}
      <button type="button" onClick=${() => setChosen([])}>Clear</button></div>`}`;
}
