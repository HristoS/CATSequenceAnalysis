import { get } from "./api.js";
import { FilterBox, Help, Message, PageHeader, Pager, Segmented, html, useEffect, useMemo, usePaged, useState } from "./html.js";

const kingdoms = [["Metazoa", "Animals", "p1"], ["Fungi", "Fungi", "p2"], ["Viridiplantae", "Plants", "p3"]];
const swatchOf = kingdom => kingdoms.find(([k]) => k === kingdom)?.[2] ?? "p0";
const verdicts = { pass: "✓ Pass", fail: "✗ Fail", context: "Context" };
const minPerClass = 3;

// How many genomes the experiment analysed, against how many findings it lists: the two are easy to confuse.
const coverage = run => (run.summary.analysed != null
  ? `analysed ${run.summary.analysed.toLocaleString()} genomes · ${run.findings.toLocaleString()} listed`
  : `${run.findings.toLocaleString()} listed`);

function VerdictBadge({ verdict }) {
  return html`<span class=${`verdict ${verdict}`}>${verdicts[verdict] ?? verdict}</span>`;
}

export function ExperimentsPage() {
  const [runs, setRuns] = useState(null);
  const [error, setError] = useState(null);
  useEffect(() => { get("/api/experiments").then(r => setRuns(r.items)).catch(e => setError(e.message)); }, []);

  if (error) return html`<${Message} error>${error}<//>`;
  if (!runs) return html`<${Message}>Loading…<//>`;
  return html`<${PageHeader} title="Experiments"
      lede=${"What the catalogue shows beyond \"related genomes look alike\". Each experiment runs offline on a frozen copy of the catalogue, " +
        "against a rule fixed before the run."} />
    ${runs.length === 0 && html`<${Message}>No experiment has been published yet.<//>`}
    <div class="scores">${runs.map(run => html`<a class="panel card" key=${run.name} href="#/experiments/${run.name}">
      <h2>${run.title}</h2>
      <p><${VerdictBadge} verdict=${run.verdict} /></p>
      <p class="hint">${run.question}</p>
      <p class="hint">${coverage(run)} · ${new Date(run.publishedAt).toLocaleDateString()}</p>
    </a>`)}</div>`;
}

function RunHeader({ run }) {
  const s = run.summary;
  return html`<${PageHeader} title=${run.title} lede=${run.question} actions=${html`<a class="button" href="#/experiments">All experiments</a>`} />
    <section class="panel stack">
      <p><${VerdictBadge} verdict=${run.verdict} /> ${s.rule && html`<span class="hint">${s.rule}</span>`}</p>
      ${s.keyNumbers && html`<div class="numbers">${s.keyNumbers.map(k => html`<div key=${k.label}><b>${k.value.toFixed(3)}</b><small>${k.label}</small></div>`)}</div>`}
      <p><b>${coverage(run)}</b>${s.analysedNote && html`<span class="hint">: ${s.analysedNote}</span>`}</p>
      <p class="hint">Computed on <code>${run.snapshot}</code> by <code>${run.script}</code>, published ${new Date(run.publishedAt).toLocaleString()}.</p>
    </section>`;
}

export function ExperimentPage({ name }) {
  const [page, setPage] = useState(null);
  const [error, setError] = useState(null);
  useEffect(() => { get(`/api/experiments/${encodeURIComponent(name)}?order=-score&limit=1000`).then(setPage).catch(e => setError(e.message)); }, [name]);

  if (error) return html`<${Message} error>${error}<//>`;
  if (!page) return html`<${Message}>Loading…<//>`;
  const views = { "e4-cpg": CpgView, "e1-outliers": OutlierView };
  const View = views[name];
  return html`<${RunHeader} run=${page.experiment} />
    ${View ? html`<${View} page=${page} />` : html`<${Message}>${page.total.toLocaleString()} findings; this experiment has no page of its own yet.<//>`}`;
}

// The findings, already in the order wanted, paged and filtered in the page by genome name, accession or any text in the finding
// (family, class, ...); columns are [header, cell(finding), numeric?].
const findingText = f => [f.sequence.name, f.sequence.recordId, ...Object.values(f.payload ?? {}).filter(v => typeof v === "string")].join(" ");
function FindingsTable({ findings, columns, pageSize = 50 }) {
  const paged = usePaged(findings, { size: pageSize, text: findingText });
  const shown = paged.rows;
  return html`<${FilterBox} paged=${paged} label="Filter these findings" placeholder="genome, accession or family" />
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
    <section class="panel"><h2>Summary by group</h2><${GroupTable} groups=${page.experiment.summary.groups} measure=${measure} /></section>`;
}

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
