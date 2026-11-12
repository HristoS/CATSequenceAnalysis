import { queryString, send } from "./api.js";
import { FormLayout, Help, Message, PageHeader, formatBp, html, recordKinds, useState } from "./html.js";

const maxBytes = 32 * 1024 * 1024;

export function UploadPage() {
  const [file, setFile] = useState(null);
  const [kind, setKind] = useState("virus");
  const [name, setName] = useState("");
  const [taxid, setTaxid] = useState("");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState(null);
  const [records, setRecords] = useState(null);

  async function submit(event) {
    event.preventDefault();
    setError(null);
    setRecords(null);
    if (file.size > maxBytes) { setError("An upload is limited to 32 MiB."); return; }
    setBusy(true);
    try {
      // The file itself is the body, which the browser streams: the API profiles it as it arrives and keeps no DNA.
      const result = await send(`/api/sequences${queryString({ kind, name, taxid })}`, { method: "POST", body: file, contentType: "text/x-fasta" });
      setRecords(result.records);
    } catch (e) {
      setError(e.message);
    } finally {
      setBusy(false);
    }
  }

  const form = html`<form class="stack" onSubmit=${submit}>
      <label class="field"><span>FASTA file (up to 32 MiB)</span>
        <input type="file" accept=".fa,.fasta,.fna,.txt" required onChange=${e => setFile(e.target.files[0] ?? null)} />
      </label>
      <label class="field"><span>Kind<${Help}>What the sequences are. It decides which searches, maps and comparisons they join.<//></span>
        <select value=${kind} onChange=${e => setKind(e.target.value)}>
          ${recordKinds.map(([value, label]) => html`<option key=${value} value=${value}>${label}</option>`)}
        </select>
      </label>
      <label class="field"><span>Name (optional)</span><input value=${name} onInput=${e => setName(e.target.value)} /></label>
      <label class="field"><span>NCBI taxid (optional)<${Help}>The NCBI Taxonomy id of the organism, if you know it: the record page then
          shows NCBI's lineage for it.<//></span><input type="number" min="1" value=${taxid} onInput=${e => setTaxid(e.target.value)} /></label>
      <button class="primary" type="submit" disabled=${busy || !file}>${busy ? "Uploading…" : "Upload"}</button>
    </form>
    ${error && html`<${Message} error>${error}<//>`}
    ${records && html`<section class="panel"><h2>Stored ${records.length} ${records.length === 1 ? "record" : "records"}</h2>
      <ul class="plain">${records.map(r => html`<li key=${r.id}>
        <a href="#/sequences/${r.id}">${r.recordId}</a> · ${formatBp(r.lengthBp)} · ${r.windowCount.toLocaleString()} windows
      </li>`)}</ul>
    </section>`}`;

  const guide = html`<section class="panel"><h2>What happens</h2>
      <p>Each record of the FASTA file becomes one record of yours, profiled as the file streams in. They're private: only you see them in the
        catalogue, and only your searches can find them.</p></section>
    <section class="panel"><h2>What is kept</h2>
      <p>Only each record's Z-curve profile, its hashes and the details you give here. The DNA itself is discarded as it is read.</p></section>
    <section class="panel"><h2>Limits</h2><ul>
      <li>One file of up to 32 MiB: a bacterial genome, a set of viruses, or a few chromosomes of a small eukaryote.</li>
      <li>Every record of the file gets the kind chosen here.</li></ul></section>`;

  return html`<${PageHeader} title="Upload" lede="Add your own sequences: a FASTA file becomes private records you can search, compare and open like any other." />
    <${FormLayout} form=${form} guide=${guide} />`;
}
