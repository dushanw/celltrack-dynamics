/* celltrack-dynamics browser viewer. Pure static: loads data/<tag>/{meta,tracks}.json, frames/*.jpg,
   masks/*.png produced by `python -m celltrack.view export-web`. No build step, no dependencies. */
"use strict";

const $ = (id) => document.getElementById(id);
const canvas = $("canvas"), ctx = canvas.getContext("2d");

const S = {                       // viewer state
  meta: null, tracks: [], events: [], byFrame: [],   // byFrame[t] = [{track, i}] points visible at frame t
  frames: [], masks: [],          // Image objects (lazy)
  t: 0, playing: false, timer: null,
  view: { scale: 1, tx: 0, ty: 0 },
  drag: null, hover: null, selected: null,          // selected = lineage id
  dpr: window.devicePixelRatio || 1,
};
const UI = {
  image: $("l-image"), masks: $("l-masks"), maskOpacity: $("mask-opacity"), trails: $("l-trails"),
  markers: $("l-markers"), events: $("l-events"), ids: $("l-ids"), lowconf: $("l-lowconf"), conf: $("conf"),
  tail: $("tail"), time: $("time"), play: $("play"), fps: $("fps"), lag: $("lag"), dataset: $("dataset"),
};
const trailColorMode = () => document.querySelector('input[name="trail-color"]:checked').value;
const rgb = (c, a = 1) => `rgba(${c[0]},${c[1]},${c[2]},${a})`;

// ----------------------------------------------------------------------------- data
async function loadIndex() {
  const [idx, plate] = await Promise.all([
    fetch("data/index.json").then((r) => r.json()),
    fetch("data/plate.json").then((r) => (r.ok ? r.json() : null)).catch(() => null),
  ]);
  S.index = idx; S.plate = plate;
  S.tagByWell = new Map(idx.filter((e) => e.well).map((e) => [e.well, e.tag]));
  const opt = (e) => `<option value="${e.tag}">${e.label}</option>`;
  const wells = idx.filter((e) => e.well).sort((a, b) => a.well.localeCompare(b.well, undefined, { numeric: true })), other = idx.filter((e) => !e.well);
  UI.dataset.innerHTML = (wells.length ? `<optgroup label="Plate wells (6–10 h excerpts)">${wells.map(opt).join("")}</optgroup>` : "") +
    (other.length ? `<optgroup label="Other excerpts">${other.map(opt).join("")}</optgroup>` : "");
  const q = new URLSearchParams(location.search), want = q.get("d") || (q.get("w") && S.tagByWell.get(q.get("w").toUpperCase()));
  if (want && idx.some((e) => e.tag === want)) UI.dataset.value = want;
  else if (S.tagByWell.has("E1")) UI.dataset.value = S.tagByWell.get("E1");   // default: a tier-A control well
  if (plate) initPlate();
  await loadDataset(UI.dataset.value);
}

async function loadDataset(tag) {
  stop();
  $("loading").hidden = false;
  const base = `data/${tag}/`;
  const [meta, tr] = await Promise.all([fetch(base + "meta.json").then((r) => r.json()), fetch(base + "tracks.json").then((r) => r.json())]);
  S.meta = meta; S.tracks = tr.tracks; S.events = tr.events; S.selected = null; S.hover = null;
  S.trackById = new Map(S.tracks.map((k) => [k.id, k]));
  S.byFrame = Array.from({ length: meta.n_frames }, () => []);
  for (const k of S.tracks) k.pts.forEach((p, i) => S.byFrame[p[0]].push({ track: k, i }));
  S.frames = new Array(meta.n_frames).fill(null); S.masks = new Array(meta.n_frames).fill(null);
  history.replaceState(null, "", meta.well ? `?w=${meta.well}` : `?d=${tag}`);
  $("title").textContent = `${meta.source}  ·  frames ${meta.start}–${meta.stop - 1}  ·  crop ${meta.crop[2]} px at (y ${meta.crop[0]}, x ${meta.crop[1]})  ·  ${meta.n_tracks} tracks`;
  UI.time.max = meta.n_frames - 1; S.t = 0; UI.time.value = 0;
  renderLegend(); renderStats(); renderMatrix(); renderDwell(); renderSelection(); renderWell(); markCurrentWell();
  resetView();
  // preload frames progressively (first the current one)
  frameImage(0); await awaitImg(S.frames[0]);
  $("loading").hidden = true;
  draw();
  (async () => { for (let i = 1; i < meta.n_frames; i++) { if (S.meta !== meta) return; frameImage(i); await awaitImg(S.frames[i]); } })();
}

/* image cache: arr[i] is an Image; returns it when decoded, else null (and redraws when it arrives) */
function getImg(arr, i, src) {
  if (!arr[i]) {
    const im = new Image();
    im.onload = () => { if (i === S.t) draw(); };
    im.src = src; arr[i] = im;
  }
  const im = arr[i];
  return im.complete && im.naturalWidth ? im : null;
}
const frameImage = (i) => getImg(S.frames, i, `data/${S.meta.tag}/frames/${String(i).padStart(3, "0")}.jpg`);
const maskImage = (i) => getImg(S.masks, i, `data/${S.meta.tag}/masks/${String(i).padStart(3, "0")}.png`);
function awaitImg(im) { return im.complete && im.naturalWidth ? Promise.resolve(im) : new Promise((res) => { const f = im.onload; im.onload = (e) => { f && f(e); res(im); }; }); }

// ----------------------------------------------------------------------------- view transforms
function resize() {
  const r = canvas.getBoundingClientRect();
  canvas.width = Math.round(r.width * S.dpr); canvas.height = Math.round(r.height * S.dpr);
  draw();
}
function resetView() {
  const r = canvas.getBoundingClientRect(), n = S.meta ? S.meta.size : 768;
  S.view.scale = Math.min(r.width, r.height - 40) / n * 0.98;
  S.view.tx = (r.width - n * S.view.scale) / 2; S.view.ty = (r.height - 40 - n * S.view.scale) / 2;
  draw();
}
const toImage = (px, py) => [(px - S.view.tx) / S.view.scale, (py - S.view.ty) / S.view.scale];

// ----------------------------------------------------------------------------- drawing
function draw() {
  if (!S.meta) return;
  const { scale, tx, ty } = S.view, n = S.meta.size, t = S.t;
  ctx.setTransform(S.dpr, 0, 0, S.dpr, 0, 0);
  ctx.fillStyle = "#000"; ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.save(); ctx.translate(tx, ty); ctx.scale(scale, scale);
  ctx.imageSmoothingEnabled = scale < 2;
  ctx.beginPath(); ctx.rect(0, 0, n, n); ctx.clip();
  const fr = frameImage(t);
  if (UI.image.checked && fr) ctx.drawImage(fr, 0, 0); else { ctx.fillStyle = "#222"; ctx.fillRect(0, 0, n, n); }
  if (UI.masks.checked) {
    const m = maskImage(t);
    if (m) { ctx.globalAlpha = +UI.maskOpacity.value; ctx.drawImage(m, 0, 0); ctx.globalAlpha = 1; }
  }
  const lw = Math.max(0.8, 1.6 / scale), confThr = +UI.conf.value, dimLow = UI.lowconf.checked;
  const sel = S.selected, hov = S.hover ? S.hover.track : null, colors = S.meta.state_colors, states = S.meta.states;
  const tail = +UI.tail.value, mode = trailColorMode();
  // trails: for every track with a point at frame <= t within tail, draw its segments in (t - tail, t]
  if (UI.trails.checked) {
    ctx.lineCap = "round"; ctx.lineJoin = "round";
    for (const k of S.tracks) {
      const pts = k.pts; if (pts[0][0] > t || pts[pts.length - 1][0] < t - tail) continue;
      const faded = (sel !== null && k.lineage !== sel);
      ctx.lineWidth = (hov === k) ? lw * 2.2 : lw;
      for (let i = 1; i < pts.length; i++) {
        const a = pts[i - 1], b = pts[i]; if (b[0] > t || b[0] <= t - tail) continue;
        const age = (t - b[0]) / tail, conf = b[5];
        let alpha = (1 - 0.7 * age) * (faded ? 0.12 : 1); if (dimLow && conf !== null && conf < confThr) alpha *= 0.25;
        const c = mode === "id" ? k.color : colors[states[b[3]]];
        ctx.strokeStyle = rgb(c, alpha); ctx.beginPath(); ctx.moveTo(a[1], a[2]); ctx.lineTo(b[1], b[2]); ctx.stroke();
        if (dimLow && conf !== null && conf < confThr) { ctx.setLineDash([3 / scale, 3 / scale]); ctx.stroke(); ctx.setLineDash([]); }
      }
      // parent link (division): dotted line from mother's last point to daughter's first
      if (k.parent !== null && pts[0][0] <= t && pts[0][0] > t - tail) {
        const mo = S.trackById.get(k.parent); if (mo) { const e = mo.pts[mo.pts.length - 1];
          ctx.strokeStyle = rgb(colors.dividing, faded ? 0.15 : 0.9); ctx.setLineDash([2 / scale, 2 / scale]);
          ctx.beginPath(); ctx.moveTo(e[1], e[2]); ctx.lineTo(pts[0][1], pts[0][2]); ctx.stroke(); ctx.setLineDash([]); }
      }
    }
  }
  // markers at current frame
  const cur = S.byFrame[t];
  if (UI.markers.checked || UI.ids.checked) {
    const r = 6 / Math.sqrt(scale);
    ctx.font = `${Math.max(9, 11 / scale)}px sans-serif`; ctx.textBaseline = "bottom";
    for (const { track: k, i } of cur) {
      const p = k.pts[i], faded = (sel !== null && k.lineage !== sel);
      if (UI.markers.checked) {
        ctx.strokeStyle = rgb(colors[states[p[3]]], faded ? 0.2 : 1); ctx.lineWidth = (hov === k ? 2.5 : 1.4) / scale;
        ctx.beginPath(); ctx.arc(p[1], p[2], r, 0, 2 * Math.PI); ctx.stroke();
      }
      if (UI.ids.checked && !faded) { ctx.fillStyle = "rgba(255,255,255,.85)"; ctx.fillText(k.id, p[1] + r, p[2] - r); }
    }
  }
  if (UI.events.checked) {
    ctx.lineWidth = 2 / scale;
    for (const e of S.events) {
      if (Math.abs(e[0] - t) > 1) continue;
      ctx.strokeStyle = rgb(colors[e[3]], e[0] === t ? 1 : 0.4); ctx.beginPath(); ctx.arc(e[1], e[2], 16, 0, 2 * Math.PI); ctx.stroke();
    }
  }
  ctx.restore();
  // scale bar / time
  ctx.fillStyle = "rgba(255,255,255,.9)"; ctx.font = "12px sans-serif"; ctx.textBaseline = "top";
  const bar = 100 * scale; ctx.fillRect(12, 12, bar, 2); ctx.fillText(S.meta.pixel_size_um ? `${(100 * S.meta.pixel_size_um).toFixed(0)} µm` : "100 px", 12, 18);
  const hrs = ((S.meta.start + t) * S.meta.frame_interval_min / 60);
  $("time-label").textContent = `frame ${S.meta.start + t}  ·  ${hrs.toFixed(2)} h  ·  ${cur.length} cells`;
}

// ----------------------------------------------------------------------------- side panels
function renderLegend() {
  const c = S.meta.state_colors;
  $("legend").innerHTML = S.meta.states.map((s) => `<div><span class="swatch" style="background:${rgb(c[s])}"></span>${s}` +
    (s === "stationary" || s === "slow" || s === "fast" ? ` <span class="muted small">${speedRange(s)}</span>` : "") + `</div>`).join("") +
    `<div class="muted small" style="margin-top:6px">motion state = step length between consecutive frames; thresholds fitted per movie (3-component Gaussian mixture).</div>`;
}
function speedRange(s) {
  const T = S.meta.speed_thresholds_px_per_frame; if (!T) return "";
  const th = Array.isArray(T) ? T : [T.threshold_stationary_slow, T.threshold_slow_fast];
  if (th.some((v) => typeof v !== "number")) return "";
  const u = S.meta.pixel_size_um ? (v) => `${(v * S.meta.pixel_size_um / S.meta.frame_interval_min).toFixed(1)} µm/min` : (v) => `${v.toFixed(0)} px/frame`;
  return s === "stationary" ? `< ${u(th[0])}` : s === "slow" ? `${u(th[0])} – ${u(th[1])}` : `> ${u(th[1])}`;
}
function kv(obj) { return Object.entries(obj).map(([k, v]) => `<span>${k}</span><span>${v}</span>`).join(""); }
function renderStats() {
  const m = S.meta, f = m.state_fraction;
  $("stats").innerHTML = kv({
    "frames": `${m.n_frames} (${(m.n_frames * m.frame_interval_min / 60).toFixed(1)} h, ${m.frame_interval_min} min/frame)`,
    "cell-frames": m.n_cells, "tracks": m.n_tracks, "lineages": m.n_lineages,
    "divisions": m.n_divisions, "deaths": m.n_deaths,
    "stationary / slow / fast": `${pct(f.stationary)} / ${pct(f.slow)} / ${pct(f.fast)}`,
  });
}
const pct = (x) => `${(100 * (x || 0)).toFixed(0)}%`;
function renderMatrix() {
  const M = S.meta.transitions[UI.lag.value]; if (!M) { $("matrix").textContent = "n/a"; return; }
  const cell = (p) => { const v = Math.min(1, p / 0.7); return `background: rgba(76,153,255,${(0.08 + 0.85 * v).toFixed(2)})`; };
  let h = `<table class="matrix"><tr><th></th>${M.cols.map((c) => `<th class="rot">${c}</th>`).join("")}<th class="rot">n</th></tr>`;
  M.rows.forEach((r, i) => { h += `<tr><th class="row">${r}</th>${M.probs[i].map((p) => `<td style="${cell(p)}" title="${r} → ${M.cols[M.probs[i].indexOf(p)]}: ${p}">${p >= 0.005 ? p.toFixed(2) : "·"}</td>`).join("")}<td class="muted">${M.n ? M.n[i].toLocaleString() : ""}</td></tr>`; });
  $("matrix").innerHTML = h + "</table>";
}
function renderDwell() {
  const d = S.meta.dwell_minutes; if (!d) { $("dwell").textContent = "n/a"; return; }
  // summary.json layout: {"50%": {state: minutes}, "75%": {...}, "mean": {...}, "count": {...}}
  const med = d["50%"] || d.median || {}, q75 = d["75%"] || {}, mean = d.mean || {};
  const rows = {}; for (const s of ["stationary", "slow", "fast"]) if (s in med) rows[s] = `${med[s]} <span class="muted">(75%: ${q75[s] ?? "–"}, mean ${mean[s] ?? "–"})</span>`;
  $("dwell").innerHTML = kv(rows);
}
function renderSelection() {
  const sec = $("selection"); if (S.selected === null) { sec.hidden = true; return; }
  const ks = S.tracks.filter((k) => k.lineage === S.selected).sort((a, b) => a.pts[0][0] - b.pts[0][0]);
  const n = ks.reduce((s, k) => s + k.pts.length, 0), f0 = ks[0].pts[0][0], f1 = Math.max(...ks.map((k) => k.pts[k.pts.length - 1][0]));
  const speeds = ks.flatMap((k) => k.pts.map((p) => p[4]).filter((v) => v !== null));
  const mean = speeds.length ? (speeds.reduce((a, b) => a + b, 0) / speeds.length) : 0;
  $("sel").innerHTML = kv({
    "lineage": S.selected, "tracks": ks.map((k) => `${k.id}${k.parent !== null ? ` (← ${k.parent})` : ""}`).join(", "),
    "frames": `${S.meta.start + f0} – ${S.meta.start + f1}`, "cell-frames": n, "mean step": `${mean.toFixed(1)} px/frame`,
    "divisions": ks.filter((k) => k.end === "divided").length, "ends": ks.map((k) => k.end || "outside excerpt").join(", "),
  });
  sec.hidden = false;
}

// ----------------------------------------------------------------------------- plate map + plate-wide results
const PLATE_ROWS = "ABCDEFGH", PLATE_COLS = [1, 2, 3, 4, 5, 6];
const WELL_METRICS = ["cells_per_mm2", "hc_median_step_um_min", "msd_10min_um2", "D_um2_per_min", "msd_alpha", "persistence_time_min",
  "median_straightness", "frac_stationary", "frac_fast", "cluster_fraction", "division_rate_per_1000cells_h", "count_fold_change",
  "median_area_um2", "median_eccentricity", "bidir_fraction_full"];
const wellOf = (meta) => meta.well || ((meta.source || "").match(/well_([A-H]\d+)/) || [])[1] || null;
const plateWell = (w) => S.plate && S.plate.wells.find((x) => x.well === w);
function controlOf(w) { return w.group === "drug" ? plateWell(w.row === "A" ? "E1" : "F1") : plateWell(`${w.row}1`); }
function metricVal(w, metric, win) { const m = w.m[metric]; if (!m) return undefined; return m[win] !== undefined ? m[win] : m.all; }
function fmtVal(v, dec) { return v === undefined || v === null || Number.isNaN(v) ? "–" : (+v).toFixed(dec); }
function condLabel(w) { return `${w.condition}${w.group === "coating" ? ` · ${w.density / 1000}k` : ""}`; }
/* 3-stop colour map (dark blue → accent → amber), u in [0,1] */
function cmap(u) {
  const st = [[42, 58, 85], [76, 153, 255], [255, 211, 106]]; u = Math.max(0, Math.min(1, u)) * 2;
  const i = Math.min(1, Math.floor(u)), f = u - i, a = st[i], b = st[i + 1];
  return `rgb(${a.map((x, k) => Math.round(x + (b[k] - x) * f)).join(",")})`;
}

function initPlate() {
  const P = S.plate, pm = $("plate-metric"), pw = $("plate-window");
  pm.innerHTML = Object.entries(P.metrics).filter(([k]) => WELL_METRICS.includes(k)).map(([k, v]) => `<option value="${k}">${v[0]}</option>`).join("");
  pm.value = "hc_median_step_um_min";
  pw.innerHTML = P.windows.map((w) => `<option value="${w}">${w === "all" ? "whole movie" : w.replace("h", " h")}</option>`).join("");
  pw.value = "6-12h";
  pm.onchange = pw.onchange = () => { drawPlate(); renderWell(); };
  // static geometry: 6 columns × 8 rows, labels on the top/left
  const svg = $("plate"), X0 = 24, Y0 = 16, DX = 39, DY = 30, R = 12.5;
  let h = PLATE_COLS.map((c, j) => `<text x="${X0 + j * DX + DX / 2}" y="11" text-anchor="middle">${c}</text>`).join("");
  [...PLATE_ROWS].forEach((r, i) => { h += `<text x="10" y="${Y0 + i * DY + DY / 2 + 3.5}" text-anchor="middle">${r}</text>`; });
  for (const w of P.wells) {
    const cx = X0 + (w.col - 1) * DX + DX / 2, cy = Y0 + PLATE_ROWS.indexOf(w.row) * DY + DY / 2;
    h += `<g data-well="${w.well}"><circle class="cur" cx="${cx}" cy="${cy}" r="${R + 3.5}" fill="none" stroke="none" stroke-width="2"/>` +
      `<circle class="well tier-${w.tier || "x"}" data-well="${w.well}" cx="${cx}" cy="${cy}" r="${R}" fill="#333"/>` +
      `<text class="lbl" x="${cx}" y="${cy + 3}">${w.well}</text></g>`;
  }
  svg.innerHTML = h;
  svg.onclick = (e) => { const g = e.target.closest("[data-well]"); if (!g) return; const tag = S.tagByWell.get(g.dataset.well); if (tag) { UI.dataset.value = tag; loadDataset(tag); } };
  const tip = $("plate-tip"), sec = $("plate-section");
  svg.onmousemove = (e) => {
    const g = e.target.closest("[data-well]"); if (!g) { tip.hidden = true; return; }
    const w = plateWell(g.dataset.well), m = P.metrics[pm.value], v = metricVal(w, pm.value, pw.value);
    tip.innerHTML = `<b>${w.well}</b> · ${condLabel(w)} · tier ${w.tier || "?"}<br><span class="muted">${m[0]}:</span> ${fmtVal(v, m[2])} ${m[1]}` +
      (S.tagByWell.has(w.well) ? "" : `<br><span class="muted">no excerpt exported</span>`);
    const r = sec.getBoundingClientRect(); tip.style.left = `${Math.min(e.clientX - r.left + 12, r.width - 170)}px`; tip.style.top = `${e.clientY - r.top + 12}px`; tip.hidden = false;
  };
  svg.onmouseleave = () => { tip.hidden = true; };
  drawPlate();
}
function drawPlate() {
  const P = S.plate; if (!P) return;
  const metric = $("plate-metric").value, win = $("plate-window").value, m = P.metrics[metric];
  const vals = P.wells.map((w) => metricVal(w, metric, win)).filter((v) => typeof v === "number" && !Number.isNaN(v));
  const lo = Math.min(...vals), hi = Math.max(...vals), span = hi - lo || 1;
  for (const w of P.wells) {
    const v = metricVal(w, metric, win), el = $("plate").querySelector(`circle.well[data-well="${w.well}"]`);
    el.setAttribute("fill", typeof v === "number" && !Number.isNaN(v) ? cmap((v - lo) / span) : "#333");
  }
  const onlyAll = vals.length && P.wells.every((w) => !w.m[metric] || w.m[metric][win] === undefined);
  $("plate-scale").innerHTML = `<span>${fmtVal(lo, m[2])}</span><span class="bar"></span><span>${fmtVal(hi, m[2])}</span>` +
    `<span>${m[1]}${onlyAll && win !== "all" ? " (whole movie)" : ""}</span>`;
}
function markCurrentWell() {
  if (!S.plate) return;
  const cur = S.meta && wellOf(S.meta);
  for (const g of $("plate").querySelectorAll("g[data-well]")) g.querySelector(".cur").setAttribute("stroke", g.dataset.well === cur ? "#fff" : "none");
}
function renderWell() {
  const sec = $("well-section"), w = S.plate && S.meta && plateWell(wellOf(S.meta));
  if (!w) { sec.hidden = true; return; }
  const P = S.plate, ctrl = controlOf(w), isCtrl = ctrl && ctrl.well === w.well, win = $("plate-window").value;
  $("well-title").textContent = `· well ${w.well}`;
  $("well-info").innerHTML = kv({
    "condition": w.condition, "seeding": `${w.density.toLocaleString()} cells/well`,
    "group": w.group === "drug" ? "drug panel (no coating)" : "coating × density", "replicate": `${w.replicate} of 2`,
  });
  const tierCls = { A: "tier-a", B: "tier-b", C: "tier-c" }[w.tier] || "muted";
  $("well-tier").innerHTML = `<span class="${tierCls}">QC tier ${w.tier || "?"}</span> <span class="muted">— ${w.tier_reason || ""}</span>`;
  $("well-ctrl").textContent = ctrl ? `${ctrl.well}, ${condLabel(ctrl)}` : "n/a";
  const winLbl = win === "all" ? "whole movie" : win.replace("h", " h"), allCol = win !== "all";
  let h = `<tr><th>metric</th><th>${winLbl}</th>${allCol ? "<th>whole movie</th>" : ""}<th>${isCtrl ? "" : "vs control"}</th></tr>`;
  for (const k of WELL_METRICS) {
    const m = P.metrics[k]; if (!m || !w.m[k]) continue;
    const vw = w.m[k][win], va = w.m[k].all, vc = ctrl ? metricVal(ctrl, k, win) : undefined, v = vw !== undefined ? vw : va;
    const ratio = !isCtrl && typeof v === "number" && typeof vc === "number" && vc !== 0 ? v / vc : undefined;
    const rc = ratio === undefined ? "na" : ratio > 1.15 ? "up" : ratio < 0.85 ? "down" : "";
    h += `<tr><td title="${k}">${m[0]}${m[1] ? ` <span class="muted">(${m[1]})</span>` : ""}</td><td${vw === undefined ? ' class="na"' : ""}>${fmtVal(vw, m[2])}</td>${allCol ? `<td>${fmtVal(va, m[2])}</td>` : ""}` +
      `<td class="${rc}">${isCtrl ? '<span class="muted">control</span>' : ratio === undefined ? "–" : `${ratio.toFixed(2)}×`}</td></tr>`;
  }
  $("well-metrics").innerHTML = h;
  sec.hidden = false;
}

// ----------------------------------------------------------------------------- interaction
function setT(t) { S.t = Math.max(0, Math.min(S.meta.n_frames - 1, t)); UI.time.value = S.t; draw(); }
function play() { if (!S.meta) return; S.playing = true; UI.play.innerHTML = "&#10074;&#10074;"; schedule(); }
function schedule() { clearTimeout(S.timer); if (!S.playing) return; S.timer = setTimeout(() => { setT(S.t + 1 >= S.meta.n_frames ? 0 : S.t + 1); schedule(); }, 1000 / +UI.fps.value); }
function stop() { S.playing = false; clearTimeout(S.timer); UI.play.innerHTML = "&#9654;"; }

UI.play.onclick = () => (S.playing ? stop() : play());
UI.time.oninput = () => { stop(); setT(+UI.time.value); };
UI.fps.onchange = () => schedule();
$("reset").onclick = resetView;
UI.dataset.onchange = () => loadDataset(UI.dataset.value);
UI.lag.onchange = renderMatrix;
UI.tail.oninput = () => { $("tail-v").textContent = UI.tail.value; draw(); };
UI.conf.oninput = () => { $("conf-v").textContent = (+UI.conf.value).toFixed(2); draw(); };
for (const el of [UI.image, UI.masks, UI.trails, UI.markers, UI.events, UI.ids, UI.lowconf, UI.maskOpacity, ...document.querySelectorAll('input[name="trail-color"]')]) el.addEventListener("input", draw);
$("clear-sel").onclick = () => { S.selected = null; renderSelection(); draw(); };

canvas.addEventListener("wheel", (e) => {
  e.preventDefault(); const r = canvas.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top;
  const f = Math.exp(-e.deltaY * 0.0015), [ix, iy] = toImage(px, py);
  S.view.scale = Math.max(0.1, Math.min(40, S.view.scale * f));
  S.view.tx = px - ix * S.view.scale; S.view.ty = py - iy * S.view.scale; draw();
}, { passive: false });
canvas.addEventListener("mousedown", (e) => { S.drag = { x: e.clientX, y: e.clientY, tx: S.view.tx, ty: S.view.ty, moved: false }; canvas.classList.add("dragging"); });
window.addEventListener("mouseup", (e) => {
  if (S.drag && !S.drag.moved && S.hover) { S.selected = (S.selected === S.hover.track.lineage) ? null : S.hover.track.lineage; renderSelection(); draw(); }
  S.drag = null; canvas.classList.remove("dragging");
});
canvas.addEventListener("mousemove", (e) => {
  const r = canvas.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top;
  if (S.drag) { const dx = e.clientX - S.drag.x, dy = e.clientY - S.drag.y; if (Math.hypot(dx, dy) > 3) S.drag.moved = true;
    S.view.tx = S.drag.tx + dx; S.view.ty = S.drag.ty + dy; draw(); return; }
  if (!S.meta) return;
  const [ix, iy] = toImage(px, py); let best = null, bd = 14 / Math.sqrt(S.view.scale) + 4;
  for (const h of S.byFrame[S.t]) { const p = h.track.pts[h.i], d = Math.hypot(p[1] - ix, p[2] - iy); if (d < bd) { bd = d; best = h; } }
  const changed = (best ? best.track : null) !== (S.hover ? S.hover.track : null); S.hover = best;
  const tip = $("tooltip");
  if (best) {
    const k = best.track, p = k.pts[best.i], m = S.meta;
    tip.innerHTML = `<b>track ${k.id}</b> · lineage ${k.lineage}${k.parent !== null ? ` · daughter of ${k.parent}` : ""}<br>` +
      `state <span style="color:${rgb(m.state_colors[m.states[p[3]]])}">${m.states[p[3]]}</span>` +
      (p[4] !== null ? ` · step ${p[4]} px` : "") + (p[5] !== null ? ` · link conf ${p[5]}` : "") + `<br>` +
      `<span class="muted">frames ${m.start + k.pts[0][0]}–${m.start + k.pts[k.pts.length - 1][0]} (${k.pts.length})` +
      ` · start: ${k.start || "—"} · end: ${k.end || "continues outside excerpt"}</span>`;
    tip.style.left = `${px + 14}px`; tip.style.top = `${py + 14}px`; tip.hidden = false;
  } else tip.hidden = true;
  if (changed) draw();
});
canvas.addEventListener("mouseleave", () => { S.hover = null; $("tooltip").hidden = true; draw(); });
window.addEventListener("keydown", (e) => {
  if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;
  if (e.code === "Space") { e.preventDefault(); S.playing ? stop() : play(); }
  else if (e.key === "ArrowRight") { stop(); setT(S.t + 1); } else if (e.key === "ArrowLeft") { stop(); setT(S.t - 1); }
  else if (e.key === "r" || e.key === "R") resetView();
  else if (e.key === "m") { UI.masks.checked = !UI.masks.checked; draw(); }
});
new ResizeObserver(resize).observe($("stage"));
loadIndex().catch((err) => { $("loading").textContent = `failed to load data: ${err}`; console.error(err); });
