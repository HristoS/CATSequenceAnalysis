import { getBytes } from "./api.js";
import { Message, Segmented, formatBp, html, useEffect, useRef, useState } from "./html.js";

// The record's Z-curve as a 3-D path, rebuilt from its stored 500 bp level: each window's mean step times its size is an exact
// count of bases, so the running sums are the classic Z-curve at window resolution, with no DNA involved.
const views = { "3d": [-0.65, 0.42], xy: [0, 0], xz: [0, Math.PI / 2], yz: [Math.PI / 2, 0] };
const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// ZCWP (docs/database-decision.md section 4.3): "ZCWP", uint16 version, int32 window size at 8, int32 window count at 12, then x, y and
// z as little-endian float32 blocks.
function parse(buffer) {
  const view = new DataView(buffer);
  if (view.getUint32(0, true) !== 0x5057435a) throw new Error("The stored profile isn't in ZCWP format.");
  const size = view.getInt32(8, true), n = view.getInt32(12, true);
  const axis = a => Int32Array.from({ length: n }, (_, i) => Math.round(view.getFloat32(16 + 4 * (a * n + i), true) * size));
  return { size, n, x: axis(0), y: axis(1), z: axis(2) };
}

// Running sums at window boundaries, the least-squares line of each axis (Zhang's drift), and where the running GC skew (x - y counts
// 2 (G - C)) is lowest and highest: the replication origin and terminus of a bacterial chromosome.
function path(profile) {
  const { n } = profile;
  const sums = ["x", "y", "z"].map(a => { const s = new Float64Array(n + 1); for (let i = 0; i < n; i++) s[i + 1] = s[i] + profile[a][i]; return s; });
  const fits = sums.map(s => {
    const m = s.length, sk = (m - 1) * m / 2, skk = (m - 1) * m * (2 * m - 1) / 6;
    let sp = 0, skp = 0;
    for (let k = 0; k < m; k++) { sp += s[k]; skp += k * s[k]; }
    const b = (m * skp - sk * sp) / (m * skk - sk * sk);
    return [(sp - b * sk) / m, b];
  });
  let skew = 0, min = 0, max = 0, minAt = 0, maxAt = 0;
  for (let i = 0; i < n; i++) {
    skew += profile.x[i] - profile.y[i];
    if (skew < min) { min = skew; minAt = i + 1; }
    if (skew > max) { max = skew; maxAt = i + 1; }
  }
  return { sums, fits, minAt, maxAt };
}

function draw(canvas, profile, shape, { drift, yaw, pitch, chromosome }) {
  const width = canvas.clientWidth, height = canvas.clientHeight, dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = css("--panel"); ctx.fillRect(0, 0, width, height);

  const n = profile.n, points = new Array(n + 1);
  for (let k = 0; k <= n; k++)
    points[k] = shape.sums.map((s, a) => (drift ? s[k] : s[k] - (shape.fits[a][0] + shape.fits[a][1] * k)));
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const p of points) for (let a = 0; a < 3; a++) { lo[a] = Math.min(lo[a], p[a]); hi[a] = Math.max(hi[a], p[a]); }
  const centre = lo.map((v, a) => (v + hi[a]) / 2);
  let radius = 0;
  for (const p of points) radius = Math.max(radius, Math.hypot(p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]));
  radius ||= 1;
  const scale = (Math.min(width, height) / 2 - 30) / radius;   // one base of excess is the same length on every axis
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const turn = v => { const x1 = cy * v[0] + sy * v[2], z1 = -sy * v[0] + cy * v[2]; return [x1, cp * v[1] - sp * z1]; };
  const screen = points.map(p => { const [u, v] = turn([p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]]); return [width / 2 + u * scale, height / 2 - v * scale]; });

  // Coloured by position, light at the start and dark at the end.
  const from = rgb(css("--seq-a")), to = rgb(css("--seq-b")), chunks = Math.min(96, n);
  ctx.lineWidth = 2; ctx.lineJoin = "round"; ctx.lineCap = "round";
  for (let q = 0; q < chunks; q++) {
    const a = Math.floor(q * n / chunks), b = Math.floor((q + 1) * n / chunks), t = q / Math.max(1, chunks - 1);
    ctx.strokeStyle = `rgb(${from.map((c, i) => Math.round(c + (to[i] - c) * t)).join(",")})`;
    ctx.beginPath(); ctx.moveTo(...screen[a]);
    for (let k = a + 1; k <= b; k++) ctx.lineTo(...screen[k]);
    ctx.stroke();
  }
  ring(ctx, screen[0], 3, css("--muted"));
  if (chromosome) {
    ring(ctx, screen[shape.minAt], 6, css("--ink")); label(ctx, screen[shape.minAt], `skew min ${formatBp(shape.minAt * profile.size)}`, width);
    ring(ctx, screen[shape.maxAt], 4, css("--ink")); label(ctx, screen[shape.maxAt], `skew max ${formatBp(shape.maxAt * profile.size)}`, width);
  }

  const ox = width - 46, oy = height - 42;
  ctx.font = `11px ${css("--mono")}`;
  for (const [name, v] of [["x", [1, 0, 0]], ["y", [0, 1, 0]], ["z", [0, 0, 1]]]) {
    const [u, w] = turn(v);
    ctx.strokeStyle = css("--muted"); ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox + u * 24, oy - w * 24); ctx.stroke();
    ctx.fillStyle = css("--ink-2"); ctx.fillText(name, ox + u * 31 - 3, oy - w * 31 + 4);
  }
  const bar = nice(radius / 3);
  ctx.strokeStyle = css("--ink-2"); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(12, 16); ctx.lineTo(12 + bar * scale, 16); ctx.moveTo(12, 12); ctx.lineTo(12, 20);
  ctx.moveTo(12 + bar * scale, 12); ctx.lineTo(12 + bar * scale, 20); ctx.stroke();
  ctx.fillStyle = css("--ink-2"); ctx.fillText(`${bar.toLocaleString()} bases of excess`, 12, 33);
  return screen;
}

function ring(ctx, p, size, colour) {
  ctx.beginPath(); ctx.arc(p[0], p[1], size, 0, 2 * Math.PI);
  ctx.strokeStyle = css("--panel"); ctx.lineWidth = 5; ctx.stroke();
  ctx.strokeStyle = colour; ctx.lineWidth = 2; ctx.stroke();
}

function label(ctx, p, text, width) {
  const w = ctx.measureText(text).width, x = Math.min(p[0] + 10, width - w - 6), y = p[1] - 9;
  ctx.globalAlpha = 0.85; ctx.fillStyle = css("--panel"); ctx.fillRect(x - 3, y - 11, w + 6, 15); ctx.globalAlpha = 1;
  ctx.fillStyle = css("--ink"); ctx.fillText(text, x, y);
}

const rgb = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
const nice = v => { const p = 10 ** Math.floor(Math.log10(v)); const m = v / p; return (m >= 5 ? 5 : m >= 2 ? 2 : 1) * p; };

export function ZCurvePath({ id, chromosome }) {
  const [profile, setProfile] = useState(null);
  const [error, setError] = useState(null);
  const [drift, setDrift] = useState("removed");
  const [view, setView] = useState("3d");
  const canvas = useRef(null), tip = useRef(null), shape = useRef(null), screen = useRef(null), angles = useRef(views["3d"]), drag = useRef(null);

  useEffect(() => {
    let live = true;
    setProfile(null); setError(null);
    getBytes(`/api/sequences/${id}/profile?windowSize=500`).then(buffer => {
      if (!live) return;
      const parsed = parse(buffer);
      shape.current = path(parsed);
      setProfile(parsed);
    }).catch(e => { if (live) setError(e.message); });
    return () => { live = false; };
  }, [id]);

  const redraw = () => {
    if (!profile || !canvas.current) return;
    screen.current = draw(canvas.current, profile, shape.current,
      { drift: drift === "kept", yaw: angles.current[0], pitch: angles.current[1], chromosome });
  };
  useEffect(() => { if (view !== "free") angles.current = views[view]; redraw(); }, [profile, drift, view]);
  useEffect(() => {
    if (!canvas.current) return;
    const observer = new ResizeObserver(redraw);
    observer.observe(canvas.current);
    const scheme = window.matchMedia("(prefers-color-scheme: dark)");
    scheme.addEventListener("change", redraw);
    return () => { observer.disconnect(); scheme.removeEventListener("change", redraw); };
  });

  const down = e => { drag.current = { x: e.clientX, y: e.clientY, angles: angles.current }; e.currentTarget.setPointerCapture(e.pointerId); tip.current.hidden = true; };
  const up = () => { drag.current = null; };
  const move = e => {
    if (drag.current) {
      const d = drag.current;
      angles.current = [d.angles[0] + (e.clientX - d.x) * 0.01, Math.max(-1.55, Math.min(1.55, d.angles[1] + (e.clientY - d.y) * 0.01))];
      if (view !== "free") setView("free"); else redraw();
      return;
    }
    const r = canvas.current.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
    let best = -1, bestDistance = 14 * 14;
    screen.current?.forEach((p, k) => { const d = (p[0] - mx) ** 2 + (p[1] - my) ** 2; if (d < bestDistance) { bestDistance = d; best = k; } });
    if (best < 0) { tip.current.hidden = true; return; }
    const w = Math.min(best, profile.n - 1);
    tip.current.textContent = `${formatBp(best * profile.size)} · window GC ${((1 + profile.z[w] / profile.size) / 2 * 100).toFixed(1)} % · G − C ${((profile.x[w] - profile.y[w]) / 2 / profile.size * 100).toFixed(1)} %`;
    tip.current.hidden = false;
    const left = screen.current[best][0] + 12;
    tip.current.style.left = `${left + tip.current.offsetWidth > r.width - 6 ? screen.current[best][0] - tip.current.offsetWidth - 12 : left}px`;
    tip.current.style.top = `${Math.max(6, screen.current[best][1] - 14)}px`;
  };

  if (error) return html`<${Message}>${error}<//>`;
  if (!profile) return html`<${Message}>Loading the stored windows…<//>`;
  return html`<div class="toolbar">
      <${Segmented} label="Drift" value=${drift} onChange=${setDrift} options=${[["removed", "Removed"], ["kept", "Kept"]]} />
      <${Segmented} label="View" value=${view} onChange=${setView} options=${[["3d", "3-D"], ["xy", "x–y"], ["xz", "x–z"], ["yz", "y–z"]]} />
    </div>
    <div class="zpath">
      <canvas ref=${canvas} role="img" aria-label=${`Z-curve path of ${profile.n.toLocaleString()} windows; drag to rotate`}
        onPointerDown=${down} onPointerMove=${move} onPointerUp=${up} onPointerCancel=${up}
        onPointerLeave=${() => { if (!drag.current) tip.current.hidden = true; }}></canvas>
      <div class="ztip" ref=${tip} hidden></div>
    </div>
    <p class="hint">x purine − pyrimidine, y amino − keto, z strong − weak (GC). Running sums of the ${profile.n.toLocaleString()} stored ${formatBp(profile.size)} windows,
      light at the start and dark at the end; drag to rotate.${chromosome ? " The running GC skew (x − y) is lowest at the replication origin and highest at the terminus." : ""}${" "}
      Removing the drift subtracts each axis's best-fit line, which is the genome's overall composition.</p>`;
}
