import { getBytes } from "./api.js";
import { FilterBox, HowTo, Message, Pager, PanelBar, Segmented, html, useEffect, useMemo, usePaged, useState } from "./html.js";

// An experiment's result files, read and drawn in the page: every one but a dot-plot table is under 140 KB, so each is fetched whole,
// once, and shared by its charts and its table. A chart is described in the registry (results/experiments.json, explain.charts), not
// coded per experiment.
const texts = new Map();
function loadText(href) {
  if (!texts.has(href)) texts.set(href, getBytes(href).then(bytes => new TextDecoder().decode(bytes)).catch(e => { texts.delete(href); throw e; }));
  return texts.get(href);
}

function useTable(file) {
  const [table, setTable] = useState(null);
  const [error, setError] = useState(null);
  useEffect(() => {
    if (!file) return;
    setTable(null);
    setError(null);
    loadText(file.href).then(text => setTable(parseCsv(text))).catch(e => setError(e.message));
  }, [file?.href]);
  return [table, error];
}

// The CSVs as the experiments wrote them: comma-separated, fields padded with spaces and sometimes quoted. '#' lines are notes. After the
// main table, blocks separated by blank lines are either further tables (three or more columns under a header) or "measure, value" lines
// (a rule, a verdict, a total), kept as pairs. A chart reads tables[spec.table ?? 0].
function parseCsv(text) {
  const notes = [], blocks = [];
  let block = [];
  for (const line of text.replace(/^﻿/, "").split(/\r?\n/)) {
    if (line.startsWith("#")) notes.push(line.replace(/^#\s*/, ""));
    else if (!line.trim()) { if (block.length > 0) blocks.push(block); block = []; }
    else block.push(splitLine(line));
  }
  if (block.length > 0) blocks.push(block);
  // A row wider than its header has an unquoted comma in its last field (a note): those cells go back together.
  const fit = (r, width) => (width > 0 && r.length > width ? [...r.slice(0, width - 1), r.slice(width - 1).join(", ")] : r);
  const table = ([header = [], ...lines]) => {
    const rows = lines.map(r => fit(r, header.length));
    return { header, rows, records: rows.map(r => Object.fromEntries(header.map((h, i) => [h, r[i] ?? ""]))) };
  };
  const [first = [], ...rest] = blocks;
  const tables = [table(first)], pairs = [];
  for (const b of rest) {
    if (b.length > 1 && b[0].length > 2 && b.every(r => r.length === b[0].length)) tables.push(table(b));
    else pairs.push(...b);
  }
  return { ...tables[0], tables, pairs, notes };
}

function splitLine(line) {
  const fields = [];
  let field = "", quoted = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (quoted) {
      if (c === '"' && line[i + 1] === '"') { field += '"'; i++; }
      else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"' && field.trim() === "") { quoted = true; field = ""; }
    else if (c === ",") { fields.push(field.trim()); field = ""; }
    else field += c;
  }
  fields.push(field.trim());
  return fields;
}

// "85.0%" reads as 0.85 and "23/30" as 0.767, so a column written as a percentage or a count of trials draws on the same scale as one
// written as a share.
const toNumber = text => {
  if (text == null || text === "") return null;
  const fraction = /^(\d+)\s*\/\s*(\d+)$/.exec(text);
  if (fraction) return Number(fraction[2]) > 0 ? Number(fraction[1]) / Number(fraction[2]) : null;
  const percent = text.endsWith("%");
  const value = Number(percent ? text.slice(0, -1) : text);
  return Number.isFinite(value) ? (percent ? value / 100 : value) : null;
};

// Three significant digits for measured quantities (times, rates), whatever their size.
const significant = v => (v === 0 ? "0" : Math.abs(v) >= 1000 ? Math.round(v).toLocaleString() : Number(v.toPrecision(3)).toString());
const formats = {
  percent: v => `${(v * 100).toFixed(1)} %`,
  signed: v => `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(3)}`,
  decimal: v => v.toFixed(3),
  number: significant,
  bp: v => (v >= 1e9 ? `${significant(v / 1e9)} Gbp` : v >= 1e6 ? `${significant(v / 1e6)} Mbp` : v >= 1e3 ? `${significant(v / 1e3)} kbp` : `${v} bp`),
};
const formatOf = (spec, key = "format") => formats[spec[key]] ?? formats.decimal;

// Series take the categorical colours in a fixed order (styles.css --p1 …), never by rank, so a series keeps its colour on every chart.
const swatches = ["p1", "p2", "p3", "p4"];

const matches = (rule, value) => (Array.isArray(rule) ? rule.map(String).includes(value) : String(rule) === value);
const chosen = (records, { where = {}, except = {} }) => records.filter(r =>
  Object.entries(where).every(([column, rule]) => matches(rule, r[column])) && !Object.entries(except).some(([column, rule]) => matches(rule, r[column])));
const labelOf = (spec, column, value) => spec.labels?.[column]?.[value] ?? value;
// "{Rank} · {Stratum}" names a category from several columns; a bare name is one column.
const categoryOf = (spec, record) => (spec.category.includes("{")
  ? spec.category.replace(/\{([^}]+)\}/g, (_, column) => labelOf(spec, column, record[column]))
  : labelOf(spec, spec.category, record[spec.category]));

// Ticks at 1, 2 or 5 times a power of ten: the step that gives closest to `count` ticks across the domain.
function niceTicks(lo, hi, count = 5) {
  const power = 10 ** Math.floor(Math.log10((hi - lo) / count));
  const step = [1, 2, 5, 10].map(m => m * power).reduce((best, s) => (Math.abs((hi - lo) / s - count) < Math.abs((hi - lo) / best - count) ? s : best));
  const ticks = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) ticks.push(Math.round(t / step) * step);
  return ticks;
}

// On a log axis, the powers of ten across the domain, with 2 and 5 between them while it spans no more than three decades.
function logTicks(lo, hi) {
  const ticks = [];
  for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++)
    for (const m of Math.log10(hi / lo) <= 3 ? [1, 2, 5] : [1]) { const t = m * 10 ** e; if (t >= lo * 0.999 && t <= hi * 1.001) ticks.push(t); }
  return ticks;
}

function useTip() {
  const [tip, setTip] = useState(null);
  const show = (title, lines) => event => setTip({ title, lines, x: event.clientX, y: event.clientY });
  const focus = (title, lines) => event => { const r = event.currentTarget.getBoundingClientRect(); setTip({ title, lines, x: r.right, y: r.top }); };
  const hide = () => setTip(null);
  const view = tip && html`<div class="tip" style=${{ left: `${tip.x + 14}px`, top: `${tip.y + 14}px` }}><b>${tip.title}</b><span>${tip.lines.join("\n")}</span></div>`;
  return { show, focus, hide, view };
}

const recordsOf = (table, spec) => table.tables[spec.table ?? 0]?.records ?? [];

export function ChartPanel({ spec, files }) {
  const file = files.find(f => f.name === spec.file);
  const [table, error] = useTable(file);
  return html`<section class="panel stack">
    <${PanelBar} title=${spec.title} help=${spec.help} />
    ${!file ? html`<${Message}>${spec.file} isn't published with this run.<//>`
      : error ? html`<${Message} error>${error}<//>`
      : !table ? html`<${Message}>Loading…<//>`
      : spec.kind === "scatter" ? html`<${Scatter} spec=${spec} records=${chosen(recordsOf(table, spec), spec)} />`
      : spec.kind === "line" ? html`<${Lines} spec=${spec} records=${chosen(recordsOf(table, spec), spec)} />`
      : html`<${Bars} spec=${spec} records=${chosen(recordsOf(table, spec), spec)} />`}
    ${spec.read && html`<${HowTo}><p>${spec.read}</p><//>`}
  </section>`;
}

// The bars' series: listed columns (spec.series), or one value column split by another column (spec.seriesBy + spec.value).
function barGroups(spec, records) {
  const series = spec.series
    ? spec.series.map(s => ({ ...s, pick: () => true }))
    : [...new Set(records.map(r => r[spec.seriesBy]))].map(v => ({ column: spec.value, low: spec.low, high: spec.high, label: labelOf(spec, spec.seriesBy, v), pick: r => r[spec.seriesBy] === v }));
  const groups = new Map();
  // aggregate "share-true": one bar per category, the share of its rows whose column reads True (a table of one row per trial).
  if (spec.aggregate === "share-true") {
    for (const record of records) {
      const category = categoryOf(spec, record);
      if (!groups.has(category)) groups.set(category, series.map(() => ({ hits: 0, rows: 0 })));
      groups.get(category).forEach((count, i) => { count.rows++; if (record[series[i].column] === "True") count.hits++; });
    }
    return { series, groups: [...groups.entries()].map(([category, counts]) => ({ category, values: counts.map(c => ({ value: c.hits / c.rows, low: null, high: null })) })) };
  }
  for (const record of records) {
    const category = categoryOf(spec, record);
    if (!groups.has(category)) groups.set(category, series.map(() => null));
    const values = groups.get(category);
    series.forEach((s, i) => {
      const value = s.pick(record) ? toNumber(record[s.column]) : null;
      if (value != null) values[i] = { value, low: s.low ? toNumber(record[s.low]) : null, high: s.high ? toNumber(record[s.high]) : null };
    });
  }
  return { series, groups: [...groups.entries()].map(([category, values]) => ({ category, values })) };
}

// Horizontal bars from a baseline (0, or 0.5 for an AUC), squared there and rounded at the data end; an interval as a whisker; values
// printed at the bar ends while there are few enough to read; a reference line where the rule or chance sits.
function Bars({ spec, records }) {
  const tip = useTip();
  const { series, groups } = useMemo(() => barGroups(spec, records), [spec, records]);
  if (groups.length === 0) return html`<${Message}>No rows match this chart.<//>`;
  const format = formatOf(spec), baseline = spec.baseline ?? 0;
  const bars = groups.flatMap(g => g.values.filter(Boolean));
  const extent = [baseline, spec.reference?.value, ...bars.flatMap(b => [b.value, b.low, b.high])].filter(v => v != null);
  const [lo, hi] = spec.domain ?? [Math.min(...extent), Math.max(...extent)];
  const span = hi - lo || 1;
  const domainLo = spec.domain ? lo : lo - span * 0.04, domainHi = spec.domain ? hi : hi + span * 0.04;
  const width = 820, left = 250, right = 72, top = 34, barHeight = 14, barGap = 2, groupGap = 12;
  const groupHeight = series.length * barHeight + (series.length - 1) * barGap;
  const height = top + groups.length * groupHeight + (groups.length - 1) * groupGap + 10;
  const x = v => left + ((v - domainLo) / (domainHi - domainLo)) * (width - left - right);
  const labelled = bars.length <= 24;

  return html`<figure class="chart">
    <svg viewBox="0 0 ${width} ${height}" role="img" aria-label=${spec.title}>
      ${niceTicks(domainLo, domainHi).map(t => html`<g key=${t}>
        <line class="grid" x1=${x(t)} x2=${x(t)} y1=${top - 6} y2=${height - 6} />
        <text class="axis" x=${x(t)} y=${top - 12} text-anchor="middle">${format(t)}</text></g>`)}
      ${spec.reference && html`<line class="reference" x1=${x(spec.reference.value)} x2=${x(spec.reference.value)} y1=${top - 6} y2=${height - 4} />`}
      ${groups.map((g, row) => {
        const y0 = top + row * (groupHeight + groupGap);
        return html`<g key=${g.category}>
          <text class="label" x=${left - 10} y=${y0 + groupHeight / 2 + 4} text-anchor="end">${g.category}</text>
          ${g.values.map((b, i) => {
            if (!b) return null;
            const y = y0 + i * (barHeight + barGap), lines = [`${series[i].label}: ${format(b.value)}`,
              ...(b.low != null && b.high != null ? [`${spec.intervalLabel ?? "95 % interval"} ${format(b.low)} to ${format(b.high)}`] : [])];
            const end = Math.max(x(b.value), x(baseline), b.high != null ? x(b.high) : 0);
            return html`<g key=${i} class="row" tabIndex="0" aria-label=${`${g.category}, ${lines.join(", ")}`}
                onMouseEnter=${tip.show(g.category, lines)} onMouseMove=${tip.show(g.category, lines)} onMouseLeave=${tip.hide}
                onFocus=${tip.focus(g.category, lines)} onBlur=${tip.hide}>
              <rect class="hit" x=${left} y=${y - barGap / 2} width=${width - left} height=${barHeight + barGap} />
              <path class=${`mark ${swatches[i % swatches.length]}`} d=${barPath(x(baseline), x(b.value), y, barHeight)} />
              ${b.low != null && b.high != null && html`<path class="ci" d=${`M${x(b.low)},${y + barHeight / 2}H${x(b.high)}M${x(b.low)},${y + 3}V${y + barHeight - 3}M${x(b.high)},${y + 3}V${y + barHeight - 3}`} />`}
              ${labelled && html`<text class="value" x=${end + 6} y=${y + barHeight / 2 + 4}>${format(b.value)}</text>`}
            </g>`;
          })}
        </g>`;
      })}
      <line class="zero" x1=${x(baseline)} x2=${x(baseline)} y1=${top - 6} y2=${height - 4} />
    </svg>
    <figcaption class="legend">
      ${series.length > 1 && series.map((s, i) => html`<span key=${s.label}><i class=${`swatch ${swatches[i % swatches.length]}`}></i>${s.label}</span>`)}
      ${baseline !== 0 && html`<span><i class="tick"></i>${spec.baselineLabel ?? `bars start at ${format(baseline)}`}</span>`}
      ${spec.reference && html`<span><i class="dash"></i>${spec.reference.label}</span>`}
    </figcaption>
    ${tip.view}
  </figure>`;
}

// Square at the baseline, rounded (4 px) at the data end, whichever side of the baseline the value falls.
function barPath(x0, x1, y, h, radius = 4) {
  const direction = x1 >= x0 ? 1 : -1, r = Math.min(radius, Math.abs(x1 - x0), h / 2);
  return `M${x0},${y}H${x1 - direction * r}Q${x1},${y} ${x1},${y + r}V${y + h - r}Q${x1},${y + h} ${x1 - direction * r},${y + h}H${x0}Z`;
}

// Lines over a numeric x (a length, a mutation rate), one per series: listed columns, or one value column split by another. x may be
// read through xValues ("Mutation_0.10" → 0.10), and either axis may be logarithmic; y may be scaled (bases/s → M bases/s).
function Lines({ spec, records }) {
  const tip = useTip();
  const xOf = r => (spec.xValues ? spec.xValues[r[spec.x]] ?? null : toNumber(r[spec.x]));
  const yOf = (r, column) => { const v = toNumber(r[column]); return v == null ? null : v * (spec.scale ?? 1); };
  const series = (spec.series
    ? spec.series.map(s => ({ label: s.label, points: records.map(r => ({ x: xOf(r), y: yOf(r, s.column), r })) }))
    : [...new Set(records.map(r => r[spec.seriesBy]))].map(v => ({ label: labelOf(spec, spec.seriesBy, v),
        points: records.filter(r => r[spec.seriesBy] === v).map(r => ({ x: xOf(r), y: yOf(r, spec.value), r })) })))
    .map(s => ({ ...s, points: s.points.filter(p => p.x != null && p.y != null && (!spec.logX || p.x > 0) && (!spec.logY || p.y > 0)).sort((a, b) => a.x - b.x) }))
    .filter(s => s.points.length > 0);
  if (series.length === 0) return html`<${Message}>No rows match this chart.<//>`;
  const format = formatOf(spec), xFormat = formatOf(spec, "xFormat");
  const xs = series.flatMap(s => s.points.map(p => p.x)), ys = [...series.flatMap(s => s.points.map(p => p.y)), ...(spec.reference ? [spec.reference.value] : [])];
  const [xlo, xhi] = [Math.min(...xs), Math.max(...xs)];
  let [ylo, yhi] = spec.domain ?? [Math.min(...ys), Math.max(...ys)];
  if (!spec.domain && !spec.logY) { const pad = (yhi - ylo || Math.abs(yhi) || 1) * 0.06; ylo -= pad; yhi += pad; }
  if (!spec.domain && spec.logY) { ylo /= 1.3; yhi *= 1.3; }
  const width = 760, height = 420, left = 64, right = 20, top = 16, bottom = 50;
  const along = (v, lo, hi, log) => (log ? Math.log10(v / lo) / Math.log10(hi / lo) : (v - lo) / (hi - lo || 1));
  const sx = v => left + along(v, xlo, xhi, spec.logX) * (width - left - right);
  const sy = v => height - bottom - along(v, ylo, yhi, spec.logY) * (height - top - bottom);
  const xTicks = spec.xTicks ?? (spec.logX ? logTicks(xlo, xhi) : niceTicks(xlo, xhi));
  const yTicks = spec.logY ? logTicks(ylo, yhi) : niceTicks(ylo, yhi);

  return html`<figure class="chart line">
    <svg viewBox="0 0 ${width} ${height}" role="img" aria-label=${spec.title}>
      ${yTicks.map(t => html`<g key=${`y${t}`}><line class="grid" x1=${left} x2=${width - right} y1=${sy(t)} y2=${sy(t)} />
        <text class="axis" x=${left - 8} y=${sy(t) + 4} text-anchor="end">${format(t)}</text></g>`)}
      ${xTicks.map(t => html`<text key=${`x${t}`} class="axis" x=${sx(t)} y=${height - bottom + 16} text-anchor="middle">${spec.xTickLabels?.[t] ?? xFormat(t)}</text>`)}
      <line class="zero" x1=${left} x2=${width - right} y1=${height - bottom} y2=${height - bottom} />
      ${spec.reference && html`<line class="reference" x1=${left} x2=${width - right} y1=${sy(spec.reference.value)} y2=${sy(spec.reference.value)} />`}
      <text class="axistitle" x=${(left + width - right) / 2} y=${height - 10} text-anchor="middle">${spec.xLabel ?? spec.x}</text>
      <text class="axistitle" transform=${`translate(14 ${(top + height - bottom) / 2}) rotate(-90)`} text-anchor="middle">${spec.yLabel ?? spec.value ?? ""}</text>
      ${series.map((s, i) => html`<g key=${s.label}>
        <path class=${`trace ${swatches[i % swatches.length]}`} d=${s.points.map((p, j) => `${j ? "L" : "M"}${sx(p.x)},${sy(p.y)}`).join("")} />
        ${s.points.map((p, j) => {
          const title = `${s.label}`, lines = [`${spec.xLabel ?? spec.x}: ${xFormat(p.x)}`, `${spec.yLabel ?? spec.value}: ${format(p.y)}`];
          return html`<g key=${j} class="row" tabIndex="0" aria-label=${`${title}, ${lines.join(", ")}`}
              onMouseEnter=${tip.show(title, lines)} onMouseMove=${tip.show(title, lines)} onMouseLeave=${tip.hide} onFocus=${tip.focus(title, lines)} onBlur=${tip.hide}>
            <circle class="hit" cx=${sx(p.x)} cy=${sy(p.y)} r="9" />
            <circle class=${`mark dot ${swatches[i % swatches.length]}`} cx=${sx(p.x)} cy=${sy(p.y)} r="4" />
          </g>`;
        })}
      </g>`)}
    </svg>
    <figcaption class="legend">
      ${series.length > 1 && series.map((s, i) => html`<span key=${s.label}><i class=${`swatch ${swatches[i % swatches.length]}`}></i>${s.label}</span>`)}
      ${spec.reference && html`<span><i class="dash"></i>${spec.reference.label}</span>`}
      ${(spec.logX || spec.logY) && html`<span>${spec.logX && spec.logY ? "both axes" : spec.logX ? "horizontal axis" : "vertical axis"} on a log scale</span>`}
    </figcaption>
    ${tip.view}
  </figure>`;
}

// One dot per row on two axes of the same quantity, with the diagonal where they are equal; split colours the dots by which side of a
// value a column falls.
function Scatter({ spec, records }) {
  const tip = useTip();
  const points = records.map(r => ({ record: r, x: toNumber(r[spec.x]), y: toNumber(r[spec.y]) })).filter(p => p.x != null && p.y != null);
  if (points.length === 0) return html`<${Message}>No rows match this chart.<//>`;
  const format = formatOf(spec);
  const values = points.flatMap(p => [p.x, p.y]);
  const [lo, hi] = spec.domain ?? [Math.min(...values), Math.max(...values)];
  const width = 640, height = 480, left = 56, right = 16, top = 16, bottom = 50;
  const sx = v => left + ((v - lo) / (hi - lo)) * (width - left - right), sy = v => height - bottom - ((v - lo) / (hi - lo)) * (height - top - bottom);
  const ticks = niceTicks(lo, hi);
  const side = p => (spec.split ? (toNumber(p.record[spec.split.column]) > spec.split.at ? "above" : "below") : "above");
  const swatch = { above: "p1", below: "p0" };
  const counts = { above: points.filter(p => side(p) === "above").length, below: points.filter(p => side(p) === "below").length };

  return html`<figure class="chart scatter">
    <svg viewBox="0 0 ${width} ${height}" role="img" aria-label=${spec.title}>
      ${ticks.map(t => html`<g key=${t}>
        <line class="grid" x1=${sx(t)} x2=${sx(t)} y1=${top} y2=${height - bottom} />
        <line class="grid" x1=${left} x2=${width - right} y1=${sy(t)} y2=${sy(t)} />
        <text class="axis" x=${sx(t)} y=${height - bottom + 16} text-anchor="middle">${format(t)}</text>
        <text class="axis" x=${left - 8} y=${sy(t) + 4} text-anchor="end">${format(t)}</text></g>`)}
      ${spec.diagonal && html`<line class="reference" x1=${sx(lo)} y1=${sy(lo)} x2=${sx(hi)} y2=${sy(hi)} />`}
      <text class="axistitle" x=${(left + width - right) / 2} y=${height - 10} text-anchor="middle">${spec.xLabel ?? spec.x}</text>
      <text class="axistitle" transform=${`translate(14 ${(top + height - bottom) / 2}) rotate(-90)`} text-anchor="middle">${spec.yLabel ?? spec.y}</text>
      ${points.map((p, i) => {
        const title = p.record[spec.name] ?? `Row ${i + 1}`;
        const lines = [`${spec.yLabel ?? spec.y}: ${format(p.y)}`, `${spec.xLabel ?? spec.x}: ${format(p.x)}`, ...(spec.tip ?? []).map(c => `${c}: ${p.record[c]}`)];
        return html`<g key=${i} class="row" tabIndex="0" aria-label=${`${title}, ${lines.join(", ")}`}
            onMouseEnter=${tip.show(title, lines)} onMouseMove=${tip.show(title, lines)} onMouseLeave=${tip.hide} onFocus=${tip.focus(title, lines)} onBlur=${tip.hide}>
          <circle class="hit" cx=${sx(p.x)} cy=${sy(p.y)} r="8" />
          <circle class=${`mark dot ${swatch[side(p)]}`} cx=${sx(p.x)} cy=${sy(p.y)} r="4.5" />
        </g>`;
      })}
    </svg>
    <figcaption class="legend">
      ${spec.split && ["above", "below"].map(s => html`<span key=${s}><i class=${`swatch ${swatch[s]}`}></i>${spec.split[s]} (${counts[s]})</span>`)}
      ${spec.diagonal && html`<span><i class="dash"></i>${spec.diagonalLabel ?? "equal on both axes"}</span>`}
    </figcaption>
    ${tip.view}
  </figure>`;
}

// Every result file readable in the page, a file at a time, paged and filtered like the other tables.
export function ResultTables({ files }) {
  const [name, setName] = useState(files[0]?.name);
  if (files.length === 0) return null;
  const file = files.find(f => f.name === name) ?? files[0];
  const choice = files.length > 4
    ? html`<label class="field"><span>File</span><select value=${file.name} onChange=${e => setName(e.target.value)}>
        ${files.map(f => html`<option key=${f.name} value=${f.name}>${f.name}</option>`)}</select></label>`
    : files.length > 1 && html`<${Segmented} label="File" value=${file.name} onChange=${setName} options=${files.map(f => [f.name, f.name])} />`;
  return html`${choice}<${CsvTable} key=${file.href} file=${file} />`;
}

const pageRows = 25;
function CsvTable({ file }) {
  const [table, error] = useTable(file);
  const rows = table?.rows ?? none;
  const paged = usePaged(rows, { size: pageRows, text: r => r.join(" ") });
  if (error) return html`<${Message} error>${error}<//>`;
  if (!table) return html`<${Message}>Loading ${file.name}…<//>`;
  return html`${rows.length > pageRows && html`<div class="toolbar"><${FilterBox} paged=${paged} label="Filter rows" placeholder="any text in a row" /></div>`}
    <${Grid} header=${table.header} rows=${paged.rows} all=${rows} />
    <${Pager} paged=${paged} noun="rows" />
    ${table.tables.slice(1).map((t, i) => html`<${Grid} key=${i} header=${t.header} rows=${t.rows} all=${t.rows} />`)}
    ${table.pairs.length > 0 && html`<div class="tablewrap"><table class="csv pairs"><tbody>
      ${table.pairs.map((p, i) => html`<tr key=${i}><th scope="row">${p[0]}</th><td>${p.slice(1).join(", ")}</td></tr>`)}</tbody></table></div>`}
    ${table.notes.length > 0 && html`<p class="hint">${table.notes.join(" ")}</p>`}`;
}

// A column is numeric, and right-aligned, when every filled cell in it reads as a number.
function Grid({ header, rows, all }) {
  const numeric = header.map((_, i) => all.some(r => r[i]) && all.every(r => !r[i] || toNumber(r[i]) != null));
  return html`<div class="tablewrap"><table class="csv">
    <thead><tr>${header.map((h, i) => html`<th key=${i} class=${numeric[i] ? "num" : undefined}>${h}</th>`)}</tr></thead>
    <tbody>${rows.map((r, j) => html`<tr key=${j}>${header.map((_, i) => html`<td key=${i} class=${numeric[i] ? "num mono" : undefined}>${r[i] ?? ""}</td>`)}</tr>`)}</tbody>
  </table></div>`;
}
const none = [];
