import { queryString, send } from "./api.js";
import { Message, formatBp, html, useState } from "./html.js";

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

  return html`<h1>Upload</h1>
    <p class="lede">A FASTA file becomes your own private records. Each is profiled as it streams in; only the profile, hashes and metadata are kept, never the DNA.</p>
    <form class="stack narrow" onSubmit=${submit}>
      <label class="field"><span>FASTA file (up to 32 MiB)</span>
        <input type="file" accept=".fa,.fasta,.fna,.txt" required onChange=${e => setFile(e.target.files[0] ?? null)} />
      </label>
      <label class="field"><span>Kind</span>
        <select value=${kind} onChange=${e => setKind(e.target.value)}>
          ${["virus", "bacterium", "archaeon", "plasmid"].map(k => html`<option key=${k}>${k}</option>`)}
        </select>
      </label>
      <label class="field"><span>Name (optional)</span><input value=${name} onInput=${e => setName(e.target.value)} /></label>
      <label class="field"><span>NCBI taxid (optional)</span><input type="number" min="1" value=${taxid} onInput=${e => setTaxid(e.target.value)} /></label>
      <button class="primary" type="submit" disabled=${busy || !file}>${busy ? "Uploading…" : "Upload"}</button>
    </form>
    ${error && html`<${Message} error>${error}<//>`}
    ${records && html`<section class="panel"><h2>Stored ${records.length} ${records.length === 1 ? "record" : "records"}</h2>
      <ul class="plain">${records.map(r => html`<li key=${r.id}>
        <a href="#/sequences/${r.id}">${r.recordId}</a> · ${formatBp(r.lengthBp)} · ${r.windowCount.toLocaleString()} windows
      </li>`)}</ul>
    </section>`}`;
}
