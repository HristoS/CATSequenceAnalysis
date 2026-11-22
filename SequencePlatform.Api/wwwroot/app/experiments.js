import { download, get, send } from "./api.js";
import { ChartPanel, ResultTables } from "./chart.js";
import { FilterBox, Help, HowTo, Message, PageHeader, Pager, PanelBar, Segmented, html, useEffect, useMemo, usePaged, useState } from "./html.js";

const kingdoms = [["Metazoa", "Animals", "p1"], ["Fungi", "Fungi", "p2"], ["Viridiplantae", "Plants", "p3"]];
const swatchOf = kingdom => kingdoms.find(([k]) => k === kingdom)?.[2] ?? "p0";
const verdicts = { pass: "✓ Pass", fail: "✗ Fail", context: "Measurement" };
const verdictHelp = "Pass and Fail are judged against a rule written down before the experiment was run, so the result can't be bent to fit " +
  "it afterwards. A Measurement has no rule: it describes how the method behaves, it doesn't decide anything.";
const minPerClass = 3;

// How many genomes the experiment analysed, against how many findings it lists: the two are easy to confuse. More listed than analysed
// means some got a finding with no answer to check it against (e2b's phages without a recorded host). Null when it lists none and says
// nothing about coverage (the article experiments published from their CSVs).
function coverage(run) {
  const analysed = run.summary.analysed, listed = run.findings;
  if (analysed == null) return listed > 0 ? `${listed.toLocaleString()} listed` : null;
  if (listed > analysed)
    return `${listed.toLocaleString()} listed: ${analysed.toLocaleString()} checked against a known answer, ${(listed - analysed).toLocaleString()} more without one`;
  return `analysed ${analysed.toLocaleString()} genomes · ${listed.toLocaleString()} listed`;
}

// The registry's article groups, in the order the list shows them, under names that say what they're about (the article codes mean
// nothing to a reader); runs without one are the experiments that publish their own findings.
const sections = [["A-genome", "Comparing whole genomes"], ["A-gene", "Genes and proteins"], ["B", "Searching the catalogue"],
  ["P9", "What composition says about biology"], [null, "Findings across the catalogue"]];

// What an experiment tests (summary.method, from the registry): the platform's own Z-curve profiles, of genomes and of genes, are marked
// so they stand out from tests of classic signatures the platform also stores. [badge, class, title of the panel on what the test
// shows about the method]. The runs that publish themselves carry no method and count as other.
const methods = {
  zcurve: ["Our Z-curve method · genomes", "zcurve", "What it shows about the Z-curve method"],
  "zcurve-gene": ["Our Z-curve method · genes", "zcurve", "What it shows about the Z-curve method"],
  signature: ["Classic composition signature", "signature", "Which method this tests"],
  cat: ["Earlier method (CAT)", "cat", "Which method this tests"],
  platform: ["Platform infrastructure", "platform", "What this measures"],
};
const methodOf = run => run.summary.method ?? run.summary.explain?.method;
const isZCurve = run => methods[methodOf(run)]?.[1] === "zcurve";
function MethodBadge({ method }) {
  const [label, kind] = methods[method] ?? [];
  return label ? html`<span class=${`methodbadge ${kind}`}>${label}</span>` : null;
}
const paragraphs = text => (Array.isArray(text) ? text : [text]).map((p, i) => html`<p key=${i}>${p}</p>`);
const where = run => (run.runnable ? "Rerunnable from the database" : "Computed offline");

// Shares and scores (|x| < 1) keep three decimals, as the P9 pages always showed them; counts and timings read as written.
const formatKey = value => (Math.abs(value) < 1 && value !== 0 ? value.toFixed(Math.abs(value) < 0.01 ? 5 : 3)
  : value.toLocaleString(undefined, { maximumFractionDigits: 2 }));
const formatSize = bytes => (bytes >= 1024 ? `${(bytes / 1024).toFixed(1)} KB` : `${bytes} B`);

function VerdictBadge({ verdict }) {
  return html`<span class=${`verdict ${verdict}`}>${verdicts[verdict] ?? verdict}</span>`;
}

export function ExperimentsPage() {
  const [runs, setRuns] = useState(null);
  const [error, setError] = useState(null);
  const [shown, setShown] = useState("all");
  useEffect(() => { get("/api/experiments").then(r => setRuns(r.items)).catch(e => setError(e.message)); }, []);

  if (error) return html`<${Message} error>${error}<//>`;
  if (!runs) return html`<${Message}>Loading…<//>`;
  // Oldest first within a section: publish-experiments publishes in the registry's order.
  const ordered = [...runs].sort((a, b) => new Date(a.publishedAt) - new Date(b.publishedAt))
    .filter(run => shown === "all" || (shown === "zcurve") === isZCurve(run));
  const zcurveCount = runs.filter(isZCurve).length;
  return html`<${PageHeader} title="Experiments"
      lede="Every experiment behind the method: what it asks, how it was tested, and what came out." />
    <p class="hint">Each was run offline, on a frozen copy of the catalogue or on the bundled genomes. ${verdictHelp}</p>
    <div class="toolbar"><${Segmented} label="Method" value=${shown} onChange=${setShown}
      help="Our Z-curve method: the platform's own profiles of genomes (windows along the genome) and of genes (windows from the start codon). Other: classic composition signatures, the earlier CAT method and the platform's infrastructure."
      options=${[["all", `All (${runs.length})`], ["zcurve", `Our Z-curve method (${zcurveCount})`], ["other", `Other (${runs.length - zcurveCount})`]]} /></div>
    ${runs.length === 0 && html`<${Message}>No experiment has been published yet.<//>`}
    ${sections.map(([article, heading]) => {
      const inSection = ordered.filter(run => (run.article ?? null) === article);
      return inSection.length > 0 && html`<section class="stack" key=${heading}>
        <h2>${heading} <span class="hint">${inSection.length}</span></h2>
        <div class="scores">${inSection.map(run => html`<a class="panel card" key=${run.name} href="#/experiments/${encodeURIComponent(run.name)}">
          <h3>${run.title}</h3>
          <p class="badges"><${VerdictBadge} verdict=${run.verdict} /><${MethodBadge} method=${methodOf(run)} /></p>
          <p class="hint">${run.question}</p>
          <p class="hint">${[coverage(run), run.files.length > 0 && `${run.files.length} file${run.files.length > 1 ? "s" : ""}`, where(run)].filter(Boolean).join(" · ")}</p>
        </a>`)}</div>
      </section>`;
    })}`;
}

function CopyCommand({ command }) {
  const [copied, setCopied] = useState(false);
  const copy = () => navigator.clipboard?.writeText(command).then(() => { setCopied(true); setTimeout(() => setCopied(false), 1500); });
  return html`<div class="actions"><code class="command">${command}</code>
    <button type="button" onClick=${copy}>${copied ? "Copied" : "Copy"}</button></div>`;
}

function DownloadButton({ file }) {
  const [state, setState] = useState(null);
  const start = () => {
    setState("…");
    download(file.href, file.name).then(() => setState(null)).catch(e => setState(e.message));
  };
  return html`<button type="button" onClick=${start} disabled=${state === "…"} title=${state && state !== "…" ? state : undefined}>
    Download ${file.name} <span class="hint">${state && state !== "…" ? "failed" : formatSize(file.byteSize)}</span></button>`;
}

const back = html`<a class="button" href="#/experiments">All experiments</a>`;
const provenance = run => html`<p class="hint">Computed on <code>${run.snapshot}</code> by <code>${run.script}</code>${run.commit && html` at <code>${run.commit}</code>`}${run.machine && ` on ${run.machine}`}${run.summary.section && `; article section ${run.summary.section}`}. Published ${new Date(run.publishedAt).toLocaleString()}.</p>`;

// An experiment the registry doesn't explain yet: its verdict, numbers and note as the article states them.
function RunHeader({ run }) {
  const s = run.summary;
  const covered = coverage(run);
  return html`<${PageHeader} title=${run.title} lede=${run.question} actions=${back} />
    <section class="panel stack">
      <p><${VerdictBadge} verdict=${run.verdict} /> ${s.rule && html`<span class="hint">${s.rule}</span>`}</p>
      ${s.keyNumbers && html`<${KeyNumbers} numbers=${s.keyNumbers} />`}
      ${s.note && html`<p>${s.note}</p>`}
      ${covered && html`<p><b>${covered}</b>${s.analysedNote && html`<span class="hint">: ${s.analysedNote}</span>`}</p>`}
      ${provenance(run)}
    </section>
    <${ResultsPanel} run=${run} />`;
}

const KeyNumbers = ({ numbers }) => html`<div class="numbers">${numbers.map(k => html`<div key=${k.label}><b>${formatKey(k.value)}</b><small>${k.label}</small></div>`)}</div>`;

// An explained experiment (summary.explain, from the registry), for a biologist who doesn't know the method: the answer and why, what it
// asks, how it was tested and what it means, what it shows about the method, the results drawn from its own tables, the terms it uses.
function ExplainedRun({ run }) {
  const s = run.summary, e = s.explain;
  const [, kind, methodTitle] = methods[methodOf(run)] ?? [];
  return html`<${PageHeader} title=${run.title} lede=${run.question} actions=${back} />
    <section class="panel stack">
      <${PanelBar} title="The answer" help=${verdictHelp} />
      <p class="badges"><${VerdictBadge} verdict=${run.verdict} /><${MethodBadge} method=${methodOf(run)} /></p>
      ${paragraphs(e.verdict)}
      ${s.keyNumbers && html`<${KeyNumbers} numbers=${s.keyNumbers} />`}
    </section>
    <div class="summarygrid">
      <section class="panel stack"><h2>What it asks</h2>${paragraphs(e.asks)}</section>
      <section class="panel stack"><h2>How it was tested</h2>${paragraphs(e.how)}</section>
      <section class="panel stack"><h2>What it means</h2>${paragraphs(e.meaning)}</section>
    </div>
    ${e.methodNote && html`<section class=${`panel stack methodnote ${kind ?? ""}`}><h2>${methodTitle ?? "Method"}</h2>${paragraphs(e.methodNote)}</section>`}
    ${(e.charts ?? []).map(spec => html`<${ChartPanel} key=${spec.title} spec=${spec} files=${run.files} />`)}
    ${e.terms && html`<section class="panel stack"><h2>Terms used here</h2>
      <dl class="terms">${e.terms.map(([term, meaning]) => html`<div key=${term}><dt>${term}</dt><dd>${meaning}</dd></div>`)}</dl></section>`}
    <${ResultsPanel} run=${run} details=${true} />`;
}

// The run's tables, readable in the page and downloadable, with how to reproduce them; an explained run keeps its provenance and the
// article's own note here, under the tables, since they're for whoever checks the work.
function ResultsPanel({ run, details = false }) {
  const s = run.summary;
  if (run.files.length === 0 && !run.command) return null;
  return html`<section class="panel stack">
    <${PanelBar} title="Result tables" help=${"The tables behind the numbers above, exactly as the experiment wrote them. " +
      (run.runnable ? "This experiment needs only the database, so it can be rerun on it." : "This experiment needs DNA or outside data, so it is rerun offline only.")} />
    <${ResultTables} files=${run.files} />
    ${run.files.length > 0 && html`<div class="actions">${run.files.map(f => html`<${DownloadButton} key=${f.name} file=${f} />`)}</div>`}
    ${run.command && html`<p class="hint">Reproduce from the repository root (${where(run).toLowerCase()}):</p><${CopyCommand} command=${run.command} />`}
    ${details && html`<${HowTo} title="Technical details">
      ${s.rule && html`<p>Rule as written: ${s.rule}</p>`}
      ${s.note && html`<p>Note from the run: ${s.note}</p>`}
      ${coverage(run) && html`<p>${coverage(run)}${s.analysedNote && `: ${s.analysedNote}`}</p>`}
      ${provenance(run)}
    <//>`}
  </section>`;
}

const jobActive = job => job && (job.status === "Queued" || job.status === "Running");
const pollMilliseconds = 3000;
const minutes = (from, to) => Math.max(1, Math.round((new Date(to) - new Date(from)) / 60000));

function jobText(job) {
  if (!job) return null;
  const by = `by ${job.requestedBy}`;
  switch (job.status) {
    case "Queued": return `Queued ${by} at ${new Date(job.createdAt).toLocaleString()}; it starts when no search is waiting.`;
    case "Running": return `Running since ${new Date(job.startedAt).toLocaleTimeString()} (started ${by}).`;
    case "Done": return `Last rerun ${by} finished ${new Date(job.completedAt).toLocaleString()}, in about ${minutes(job.startedAt, job.completedAt)} min.`;
    default: return `Last rerun ${by} failed: ${job.error ?? "no reason given"}`;
  }
}

// A database experiment rerun on today's catalogue (migration 0024): the latest rerun's numbers beside the article's, which stay as
// published. Only an admin sees the button; the experiment worker runs the job and the panel polls until it ends.
function RerunPanel({ article }) {
  const name = article.name;
  const [reruns, setReruns] = useState(null);
  const [server, setServer] = useState(undefined);
  const [error, setError] = useState(null);
  const [starting, setStarting] = useState(false);
  const latest = reruns?.jobs[0];

  const loadServer = () => get(`/api/experiments/${encodeURIComponent(name)}?origin=server&limit=1`)
    .then(page => setServer(page.experiment))
    .catch(e => (e.status === 404 ? setServer(null) : setError(e.message)));
  useEffect(() => {
    get(`/api/experiments/${encodeURIComponent(name)}/runs`).then(setReruns).catch(e => setError(e.message));
    loadServer();
  }, [name]);
  useEffect(() => {
    if (!jobActive(latest)) return undefined;
    const timer = setTimeout(() => get(latest.href).then(job => {
      setReruns(r => ({ ...r, jobs: [job, ...r.jobs.slice(1)] }));
      if (job.status === "Done") loadServer();
    }).catch(e => setError(e.message)), pollMilliseconds);
    return () => clearTimeout(timer);
  }, [latest]);

  if (!reruns || server === undefined) return null;
  if (!reruns.rerunnable && !server) return null;
  const start = () => {
    setStarting(true);
    setError(null);
    send(`/api/experiments/${encodeURIComponent(name)}/runs`, { method: "POST" })
      .then(job => setReruns(r => ({ ...r, jobs: [job, ...r.jobs] })))
      .catch(e => setError(e.message))
      .finally(() => setStarting(false));
  };
  return html`<section class="panel stack">
    <${PanelBar} title="Rerun on this database"
      help=${"The same experiment run again on the catalogue as it is today, by the platform's own code. The article's numbers above stay as " +
        "published; the rerun's sit beside them. Only admins can start one, and it waits while any search is queued or running."}>
      ${reruns.canRun && html`<button type="button" class="primary" onClick=${start} disabled=${starting || jobActive(latest)}>
        ${jobActive(latest) ? "Rerun in progress…" : "Rerun now"}</button>`}
    <//>
    ${error && html`<${Message} error>${error}<//>`}
    ${latest && html`<p class=${latest.status === "Failed" ? "message error" : "hint"}>${jobText(latest)}</p>`}
    ${server ? html`<${RerunComparison} article=${article} server=${server} />`
      : !jobActive(latest) && html`<${Message}>Not rerun on this database yet.<//>`}
  </section>`;
}

// One row per key number either run reports, by label: a rerun reports the article's labels where it computes the same number.
function RerunComparison({ article, server }) {
  const articleNumbers = article.summary.keyNumbers ?? [];
  const serverNumbers = server.summary.keyNumbers ?? [];
  const labels = [...new Set([...articleNumbers, ...serverNumbers].map(k => k.label))];
  const value = (numbers, label) => numbers.find(k => k.label === label)?.value;
  const cell = v => (v == null ? "–" : formatKey(v));
  return html`<div class="tablewrap"><table>
      <thead><tr><th></th><th class="num">Article run</th><th class="num">This database</th></tr></thead>
      <tbody>
        <tr><td>Verdict</td><td class="num"><${VerdictBadge} verdict=${article.verdict} /></td><td class="num"><${VerdictBadge} verdict=${server.verdict} /></td></tr>
        ${labels.map(label => html`<tr key=${label}><td>${label}</td>
          <td class="num mono">${cell(value(articleNumbers, label))}</td><td class="num mono">${cell(value(serverNumbers, label))}</td></tr>`)}
      </tbody>
    </table></div>
    ${server.summary.rule && html`<p class="hint">${server.summary.rule}</p>`}
    ${coverage(server) && html`<p><b>${coverage(server)}</b>${server.summary.analysedNote && html`<span class="hint">: ${server.summary.analysedNote}</span>`}</p>`}
    <p class="hint">Computed on <code>${server.snapshot}</code> by <code>${server.script}</code>${server.commit && html` at <code>${server.commit}</code>`}${server.machine && ` on ${server.machine}`}. Published ${new Date(server.publishedAt).toLocaleString()}.</p>
    ${server.files.length > 0 && html`<div class="actions">${server.files.map(f => html`<${DownloadButton} key=${f.name} file=${f} />`)}</div>`}`;
}

// Each experiment's own view, and how many findings the page needs up front (a view that loads its own asks for one).
const views = { "e4-cpg": [() => CpgView, 1000], "e1-outliers": [() => OutlierView, 1000], "e2b-phage-host": [() => PhageHostView, 1] };

export function ExperimentPage({ name }) {
  const [page, setPage] = useState(null);
  const [error, setError] = useState(null);
  const [view, limit] = views[name] ?? [() => null, 1000];
  useEffect(() => { get(`/api/experiments/${encodeURIComponent(name)}?order=-score&limit=${limit}`).then(setPage).catch(e => setError(e.message)); }, [name]);

  if (error) return html`<${Message} error>${error}<//>`;
  if (!page) return html`<${Message}>Loading…<//>`;
  const View = view();
  return html`${page.experiment.summary.explain ? html`<${ExplainedRun} run=${page.experiment} />` : html`<${RunHeader} run=${page.experiment} />`}
    <${RerunPanel} article=${page.experiment} />
    ${View ? html`<${View} page=${page} />` : page.total > 0 && html`<${Message}>${page.total.toLocaleString()} findings; this experiment has no page of its own yet.<//>`}`;
}

// The findings, already in the order wanted, paged and filtered in the page by genome name, accession or any text in the finding
// (family, class, ...); columns are [header, cell(finding), numeric?].
const findingText = f => [f.sequence.name, f.sequence.recordId, ...Object.values(f.payload ?? {}).filter(v => typeof v === "string")].join(" ");
function FindingsTable({ findings, columns, pageSize = 50, text = findingText, placeholder = "genome, accession or family" }) {
  const paged = usePaged(findings, { size: pageSize, text });
  const shown = paged.rows;
  return html`<${FilterBox} paged=${paged} label="Filter these findings" placeholder=${placeholder} />
    <div class="tablewrap"><table>
      <thead><tr>${columns.map(([header, , numeric]) => html`<th key=${header} class=${numeric ? "num" : undefined}>${header}</th>`)}</tr></thead>
      <tbody>${shown.map(f => html`<tr key=${`${f.kind}-${f.sequence.id}`}>
        ${columns.map(([header, cell, numeric]) => html`<td key=${header} class=${numeric ? "num mono" : undefined}>${cell(f)}</td>`)}
      </tr>`)}</tbody>
    </table></div>
    ${paged.filtered.length === 0 && html`<${Message}>${findings.length ? "No finding matches the filter." : "None."}<//>`}
    <${Pager} paged=${paged} noun="findings" />`;
}

const genomeLink = f => html`<a href="#/sequences/${f.sequence.id}">${f.sequence.name || f.sequence.recordId}</a>`;
const keyNames = { gtdbFamilyDiffers: "GTDB family differs", checkmPoor: "CheckM: contamination ≥ 5 % or completeness < 90 %", ncbiFamilyDiffers: "NCBI family changed since" };
const groupNames = { cellular: "Bacteria and archaea", viral: "Viruses" };

function OutlierView({ page }) {
  const [group, setGroup] = useState("cellular");
  const s = page.experiment.summary;
  const findings = page.items.filter(f => f.payload.group === group);
  const columns = [
    ["Genome", genomeLink],
    ["Family", f => f.payload.family],
    ["Same-family neighbours", f => `${(f.payload.share * 10).toFixed(0)} of 10 (family median ${(f.payload.familyMedianShare * 10).toFixed(0)})`],
    ["Neighbours' family", f => f.payload.majorityFamily ? `${f.payload.majorityFamily} ×${f.payload.majorityCount}` : "–"],
    [group === "cellular" ? "GTDB family" : "NCBI family now", f => (group === "cellular" ? f.payload.gtdbFamily : f.payload.ncbiFamily) ?? "–"],
    ...(group === "cellular" ? [["CheckM", f => f.payload.checkmCompleteness == null ? "–" : `${f.payload.checkmCompleteness} / ${f.payload.checkmContamination}`, true]] : []),
    ["Score", f => f.score.toFixed(2), true],
  ];
  return html`<section class="panel stack">
      <h2>Does the score find mislabels better than GC and length?</h2>
      <div class="tablewrap"><table>
        <thead><tr><th>Group</th><th>Answer key<${Help}>The outside evidence a flagged genome is checked against.<//></th><th>Space</th>
          <th class="num">Genomes</th><th class="num">Base rate<${Help}>The share of all genomes the answer key flags: what picking at random would find.<//></th>
          <th class="num">AUC (95 %)<${Help}>How well the score ranks flagged genomes above the others: 0.5 is chance, 1 perfect.<//></th>
          <th class="num">Top 1 % rate<${Help}>The share flagged by the answer key among the 1 % highest scores; compare it with the base rate.<//></th></tr></thead>
        <tbody>${s.results.map(r => html`<tr key=${`${r.group}-${r.key}-${r.space}`}>
          <td>${groupNames[r.group] ?? r.group}</td><td>${keyNames[r.key] ?? r.key}</td><td>${r.space === "signature" ? "signature" : "GC + length"}</td>
          <td class="num">${r.withKey.toLocaleString()}</td><td class="num mono">${(r.baseRate * 100).toFixed(1)} %</td>
          <td class="num mono">${r.auc.toFixed(3)} (${r.aucLow.toFixed(3)}–${r.aucHigh.toFixed(3)})</td><td class="num mono">${(r.top1PctRate * 100).toFixed(1)} %</td>
        </tr>`)}</tbody>
      </table></div>
      ${s.reading && html`<p>${s.reading}</p>`}
      ${s.afterTheRun && html`<p class="hint">Found after the run: of ${s.afterTheRun.cellularTopWithGtdb} flagged bacteria and archaea with a GTDB family, GTDB agrees with the provider for ${s.afterTheRun.gtdbAgreesWithProvider} and with the neighbours for ${s.afterTheRun.gtdbAgreesWithNeighbours}.</p>`}
    </section>
    <section class="panel stack">
      <h2>The 200 highest scores per group</h2>
      <p class="hint">Genomes whose 10 nearest by signature mostly belong to another family, although their family's genomes usually keep together. Read as "composition unusual for its family", not as a mislabel.</p>
      <${Segmented} label="Group" value=${group} onChange=${setGroup} options=${[["cellular", "Bacteria and archaea"], ["viral", "Viruses"]]} />
      <${FindingsTable} findings=${findings} columns=${columns} />
    </section>`;
}

function CpgView({ page }) {
  const [measure, setMeasure] = useState("rhoCg");
  return html`<section class="panel stack">
      <h2>${measure === "rhoCg" ? "ρ*(CG)" : "ρ*(TA)"} per genome, by class</h2>
      <${Segmented} label="Dinucleotide" value=${measure} onChange=${setMeasure} options=${[["rhoCg", "CpG"], ["rhoTa", "TpA"]]} />
      <${ClassStrip} findings=${page.items} measure=${measure} />
    </section>
    <section class="panel"><h2>Summary by group</h2><${GroupTable} groups=${page.experiment.summary.groups} measure=${measure} /></section>
    <section class="panel stack">
      <${PanelBar} title=${`Every genome analysed (${page.items.length.toLocaleString()})`}
        help="Every eukaryote the run analysed, including those in classes with fewer than 3 genomes, which the strip above leaves out." />
      <${FindingsTable} findings=${page.items} columns=${cpgColumns} placeholder="genome, accession, kingdom or class" />
    </section>`;
}

const cpgColumns = [
  ["Genome", genomeLink],
  ["Kingdom", f => kingdoms.find(([k]) => k === f.payload.kingdom)?.[1] ?? f.payload.kingdom ?? "–"],
  ["Class", f => f.payload.class ?? "–"],
  ["ρ*(CG)", f => f.payload.rhoCg.toFixed(3), true],
  ["ρ*(TA)", f => f.payload.rhoTa.toFixed(3), true],
];

const median = values => { const s = [...values].sort((a, b) => a - b); return s[Math.floor(s.length / 2)]; };

// One row per class with enough genomes, ordered by median; a dot per genome. 1 is "as composition predicts"; below it the pair is avoided.
function ClassStrip({ findings, measure }) {
  const [tip, setTip] = useState(null);
  const rows = useMemo(() => {
    const byClass = new Map();
    for (const f of findings) {
      const key = f.payload.class;
      if (!key) continue;
      if (!byClass.has(key)) byClass.set(key, { name: key, kingdom: f.payload.kingdom, items: [] });
      byClass.get(key).items.push(f);
    }
    return [...byClass.values()].filter(r => r.items.length >= minPerClass)
      .map(r => ({ ...r, median: median(r.items.map(f => f.payload[measure])) }))
      .sort((a, b) => a.median - b.median);
  }, [findings, measure]);

  const width = 760, left = 190, right = 16, top = 26, rowHeight = 24;
  const height = top + rows.length * rowHeight + 10;
  const max = Math.max(1.25, ...findings.map(f => f.payload[measure])) * 1.02;
  const x = v => left + (v / max) * (width - left - right);
  const ticks = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75].filter(t => t <= max);
  const shown = kingdoms.filter(([k]) => rows.some(r => r.kingdom === k));

  return html`<figure class="strip"><div class="scroll">
      <svg viewBox="0 0 ${width} ${height}" role="img" aria-label=${`ρ* per genome for ${rows.length} classes`}>
        ${ticks.map(t => html`<g key=${t}>
          <line class=${t === 1 ? "reference" : "grid"} x1=${x(t)} x2=${x(t)} y1=${top - 6} y2=${height - 6} />
          <text class="axis" x=${x(t)} y=${top - 12} text-anchor="middle">${t}</text>
        </g>`)}
        ${rows.map((r, row) => {
          const y = top + row * rowHeight + rowHeight / 2;
          return html`<g key=${r.name}>
            <text class="label" x=${left - 10} y=${y + 4} text-anchor="end">${r.name} <tspan class="count">${r.items.length}</tspan></text>
            <line class="median" x1=${x(r.median)} x2=${x(r.median)} y1=${y - 9} y2=${y + 9} />
            ${r.items.map((f, i) => {
              const cy = y + (((i * 37) % 7) - 3) * 1.6;
              const show = event => setTip({ f, x: event.clientX, y: event.clientY });
              return html`<a key=${f.sequence.id} href="#/sequences/${f.sequence.id}" aria-label=${`${f.sequence.name}: ${f.payload[measure].toFixed(3)}`}
                  onMouseEnter=${show} onMouseMove=${show} onMouseLeave=${() => setTip(null)} onFocus=${e => setTip({ f, x: e.target.getBoundingClientRect().x, y: e.target.getBoundingClientRect().y })} onBlur=${() => setTip(null)}>
                <circle class="hit" cx=${x(f.payload[measure])} cy=${cy} r="8" />
                <circle class=${`dot ${swatchOf(f.payload.kingdom)}`} cx=${x(f.payload[measure])} cy=${cy} r="4" />
              </a>`;
            })}
          </g>`;
        })}
      </svg></div>
      <figcaption class="legend">
        ${shown.map(([k, label, swatch]) => html`<span key=${k}><i class=${`swatch ${swatch}`}></i>${label}</span>`)}
        <span><i class="tick"></i>median</span>
        <span>1 = as composition predicts; below 1 the pair is avoided</span>
      </figcaption>
      ${tip && html`<div class="tip" style=${{ left: `${tip.x + 14}px`, top: `${tip.y + 14}px` }}>
        <b>${tip.f.sequence.name}</b>
        <span>${tip.f.payload.class} · ${tip.f.payload.kingdom}${"\n"}ρ*(CG) ${tip.f.payload.rhoCg.toFixed(3)}   ρ*(TA) ${tip.f.payload.rhoTa.toFixed(3)}</span>
      </div>`}
    </figure>`;
}

function GroupTable({ groups, measure }) {
  const key = measure === "rhoCg" ? "rhoCg" : "rhoTa";
  const sets = groups.filter(g => g.set !== "class");
  return html`<div class="tablewrap"><table>
    <thead><tr><th>Group</th><th class="num">Genomes</th><th class="num">25th</th><th class="num">Median</th><th class="num">75th</th></tr></thead>
    <tbody>${sets.map(g => html`<tr key=${g.set}>
      <td>${g.set}</td><td class="num">${g.genomes.toLocaleString()}</td>
      ${g[key].map((v, i) => html`<td class="num mono" key=${i}>${v.toFixed(3)}</td>`)}
    </tr>`)}</tbody>
  </table></div>`;
}

// E2b: a phage's host genus from 6-mer composition. The accuracy against the round-1 method and the plain baselines, then every phage's
// guesses, loaded a page at a time because there are thousands.
const predictorNames = { a6: "6-mer composition", signature: "log ρ* (round 1)", gc: "GC alone", mostCommon: "always the commonest host" };
const guessText = f => [f.sequence.name, f.sequence.recordId, f.payload.host, f.payload.class, ...f.payload.guesses.map(g => g.genus)].join(" ");
const firstGuessRate = f => `${f.payload.firstGuess.right.toLocaleString()} of ${f.payload.firstGuess.guessed.toLocaleString()}`;

function PhageHostView({ page }) {
  const s = page.experiment.summary;
  const [scope, setScope] = useState("all");
  const [kind, setKind] = useState("recorded");
  const counts = Object.fromEntries((s.classes ?? []).filter(c => c.qualifies).map(c => [c.name, c.unrecorded]));
  return html`<section class="panel stack">
      <${PanelBar} title="How often is the recorded host found?"
        help="For each phage with a recorded host, the stored genera are ranked by their genome nearest to the phage; a hit is the recorded genus first (top 1) or among the first five (top 5)." />
      <${Segmented} label="Phages counted" value=${scope} onChange=${setScope}
        help=${`"Averaged over host genera" gives each of the ${s.macroGenera} host genera with at least 5 phages the same weight, so the many Escherichia phages don't decide it.`}
        options=${[["all", `All ${s.analysed.toLocaleString()} phages`], ["macro", "Averaged over host genera"]]} />
      <div class="accuracypair">
        <${AccuracyBars} title="Top 1" predictors=${s.predictors} value=${p => (scope === "all" ? p.top1 : p.macroTop1)} reference=${scope === "all" ? s.shuffled : null} />
        <${AccuracyBars} title="Top 5" predictors=${s.predictors} value=${p => (scope === "all" ? p.top5 : p.macroTop5)} />
      </div>
      <${HowTo}>
        <p>${s.reading}</p>
        <p>6-mers against log ρ*: ${differenceText(s.differences, "signature", scope)}; against always guessing the commonest host: ${differenceText(s.differences, "mostCommon", scope)}
          (difference in top-1, 95 % interval). ${scope === "all" ? `With host labels shuffled across phages, 6-mers reach at most ${(s.shuffled * 100).toFixed(1)} % (the grey line).` : ""}</p>
        <p>${s.classRule} Qualifying: ${Object.entries(counts).map(([name, n]) => `${name} (${n.toLocaleString()} without a host)`).join(", ")}.</p>
      <//>
    </section>
    <section class="panel stack">
      <${PanelBar} title="Every phage's five nearest host genera"
        help="Ordered by the distance to the nearest genus, closest first. A guess is composition, not evidence of infection: check how often the same first guess was right." />
      <${Segmented} label="Phages" value=${kind} onChange=${setKind}
        options=${[["recorded", "With a recorded host"], ["unrecorded", "Without a recorded host"]]} />
      <${GuessList} key=${kind} name=${page.experiment.name} kind=${kind} />
    </section>`;
}

function differenceText(differences, baseline, scope) {
  const d = differences.find(x => x.baseline === baseline);
  const [value, low, high] = scope === "all" ? [d.overall, d.overallLow, d.overallHigh] : [d.macro, d.macroLow, d.macroHigh];
  const points = v => `${v >= 0 ? "+" : ""}${(v * 100).toFixed(1)}`;
  return `${points(value)} points (${points(low)} to ${points(high)})`;
}

// One bar per predictor from 0, the 6-mer bar in the series colour and the baselines grey; each value printed beside its bar.
function AccuracyBars({ title, predictors, value, reference = null }) {
  const max = Math.max(...predictors.map(value), reference ?? 0, 0.05);
  const scale = v => `${(v / max) * 100}%`;
  return html`<figure class="accuracy">
    <figcaption>${title}</figcaption>
    ${predictors.map(p => html`<div class="bar" key=${p.predictor} title=${`${predictorNames[p.predictor] ?? p.predictor}: ${(value(p) * 100).toFixed(1)} %`}>
      <span>${predictorNames[p.predictor] ?? p.predictor}</span>
      <div class="track from0">
        <i class=${p.predictor === "a6" ? "lead" : "base"} style=${{ width: scale(value(p)) }}></i>
        ${reference != null && html`<b class="reference" style=${{ left: scale(reference) }}></b>`}
      </div>
      <span class="mono num">${(value(p) * 100).toFixed(1)} %</span>
    </div>`)}
  </figure>`;
}

// All findings of one kind, fetched 1,000 at a time and shown as they arrive.
function GuessList({ name, kind }) {
  const [findings, setFindings] = useState([]);
  const [total, setTotal] = useState(null);
  const [error, setError] = useState(null);
  useEffect(() => {
    let cancelled = false;
    (async () => {
      const loaded = [];
      for (let offset = 0; ; offset += 1000) {
        const page = await get(`/api/experiments/${encodeURIComponent(name)}?kind=${kind}&order=score&offset=${offset}&limit=1000`);
        if (cancelled) return;
        loaded.push(...page.items);
        setTotal(page.total);
        setFindings([...loaded]);
        if (loaded.length >= page.total || page.items.length === 0) break;
      }
    })().catch(e => setError(e.message));
    return () => { cancelled = true; };
  }, [name, kind]);

  if (error) return html`<${Message} error>${error}<//>`;
  if (total == null) return html`<${Message}>Loading…<//>`;
  const first = f => f.payload.guesses[0];
  const columns = kind === "recorded"
    ? [
        ["Phage", genomeLink],
        ["Recorded host", f => f.payload.host],
        ["First guess", f => html`${first(f).genus}${first(f).genus === f.payload.host ? html` <span class="ok">✓</span>` : ""}`],
        ["Recorded host's rank", f => f.payload.rank.toLocaleString(), true],
        ["First guess right for", f => firstGuessRate(f), true],
        ["Distance", f => first(f).distance.toFixed(4), true],
      ]
    : [
        ["Phage", genomeLink],
        ["Class", f => f.payload.class],
        ["First guess", f => `${first(f).genus} (${first(f).family ?? "–"})`],
        ["First guess right for", f => firstGuessRate(f), true],
        ["Next guesses", f => f.payload.guesses.slice(1).map(g => g.genus).join(", ")],
        ["Distance", f => first(f).distance.toFixed(4), true],
      ];
  return html`${findings.length < total && html`<p class="hint">Loaded ${findings.length.toLocaleString()} of ${total.toLocaleString()}…</p>`}
    <${FindingsTable} findings=${findings} columns=${columns} text=${guessText} placeholder="phage, host or guessed genus" />`;
}
