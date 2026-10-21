"""Phase 4b: state-transition matrices, dwell times, event rates and the HTML report.

Rows   = state of a cell at frame t      : stationary, slow, fast, dividing, (entered, appeared for lag 1)
Columns= state of its descendant at t+L  : stationary, slow, fast, dividing, left, died, lost
A division within the lag contributes one transition per daughter.

Usage:
    python -m celltrack.transition --name trackastra_ilp
"""

from __future__ import annotations

import argparse
import base64
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .common import Units, ensure_dir, load_config, setup_logging

log = setup_logging("transition")

ROW_STATES = ["stationary", "slow", "fast", "dividing"]
ENTRY_STATES = ["entered", "appeared"]
COL_STATES = ["stationary", "slow", "fast", "dividing", "left", "died", "merged", "lost"]
TERMINAL = {"left", "died", "merged", "lost"}
# variant where 'in_cluster' replaces the motion state of cells inside clusters
ROW_STATES_CL = ["stationary", "slow", "fast", "in_cluster", "dividing"]
COL_STATES_CL = ["stationary", "slow", "fast", "in_cluster", "dividing", "left", "died", "merged", "lost"]


def cluster_as_state(d: pd.DataFrame) -> pd.DataFrame:
    """Copy of the state table in which cells in a cluster are in state 'in_cluster'
    (unless they are dividing/dying at that frame)."""
    e = d.copy()
    m = e["in_cluster"].astype(bool) & e["state"].isin(["stationary", "slow", "fast"])
    e.loc[m, "state"] = "in_cluster"
    return e


# --------------------------------------------------------------------------- transition enumeration
def enumerate_transitions(d: pd.DataFrame, lag: int, conf_thr: float, use_bidir: bool,
                          row_states=ROW_STATES) -> pd.DataFrame:
    """One row per (source node, descendant-or-terminal at +lag)."""
    state = d.set_index("node_id")["state"].replace({"dying": "died"}).to_dict()
    frame = d.set_index("node_id")["frame"].to_dict()
    end_type = d.set_index("node_id")["end_type"].to_dict()
    conf = d.set_index("node_id")["edge_conf_in"].to_dict()
    bidir = d.set_index("node_id")["bidir"].to_dict() if ("bidir" in d and use_bidir) else {}
    children: dict[int, list[int]] = {}
    for c, p in zip(d.node_id.to_numpy(), d.parent_node.to_numpy()):
        if p != -1:
            children.setdefault(int(p), []).append(int(c))
    T1 = int(d.frame.max())
    rows = []
    src_cols = d[["node_id", "state", "frame", "in_cluster", "lineage_id"]].to_numpy()
    for nid, s, f, clus, lin in src_cols:
        if s not in row_states:
            continue
        # walk `lag` steps; stack of (node, steps_left, min_conf_on_path, all_bidir)
        stack = [(int(nid), lag, 1.0, True)]
        while stack:
            n, k, mc, bd = stack.pop()
            if k == 0:
                rows.append((nid, s, f, clus, lin, state[n], mc, bd))
                continue
            kids = children.get(n, [])
            if not kids:
                et = end_type.get(n, "")
                if et in TERMINAL:
                    rows.append((nid, s, f, clus, lin, et, mc, bd))
                # censored / still running at the end of the movie: no transition
                continue
            for c in kids:
                w = conf.get(c, np.nan)
                w = 1.0 if not np.isfinite(w) else w
                b = bool(bidir.get(c, True)) if bidir else True
                gap = frame[c] - frame[n]
                stack.append((c, max(0, k - gap), min(mc, w), bd and b))
    tr = pd.DataFrame(rows, columns=["src", "from", "frame", "in_cluster", "lineage_id", "to", "min_conf", "bidir"])
    tr["high_conf"] = (tr.min_conf >= conf_thr) & (tr.bidir if use_bidir else True)
    return tr


def entry_transitions(d: pd.DataFrame) -> pd.DataFrame:
    e = d[d.start_type.isin(ENTRY_STATES)]
    return pd.DataFrame(dict(src=e.node_id, **{"from": e.start_type}, frame=e.frame, in_cluster=e.in_cluster,
                             lineage_id=e.lineage_id, to=e.state.replace({"dying": "died"}), min_conf=1.0, bidir=True, high_conf=True))


def count_matrix(tr: pd.DataFrame, rows=None, cols=COL_STATES) -> pd.DataFrame:
    rows = rows or (ROW_STATES + ENTRY_STATES)
    m = pd.crosstab(tr["from"], tr["to"]).reindex(index=rows, columns=cols, fill_value=0)
    return m.loc[m.sum(1) > 0].copy()


def prob_matrix(cm: pd.DataFrame) -> pd.DataFrame:
    return cm.div(cm.sum(1).replace(0, np.nan), axis=0)


def bootstrap_ci(tr: pd.DataFrame, n_boot: int, rng, rows, cols=COL_STATES):
    """Resample lineages with replacement; return (lo, hi) matrices of the row-normalised probabilities."""
    lin = tr["lineage_id"].to_numpy()
    uniq, inv = np.unique(lin, return_inverse=True)
    fi = pd.Categorical(tr["from"], categories=rows).codes
    ti = pd.Categorical(tr["to"], categories=cols).codes
    ok = (fi >= 0) & (ti >= 0)
    fi, ti, inv = fi[ok], ti[ok], inv[ok]
    flat = fi * len(cols) + ti
    mats = np.empty((n_boot, len(rows), len(cols)))
    for b in range(n_boot):
        w_lin = np.bincount(rng.integers(0, len(uniq), len(uniq)), minlength=len(uniq))
        w = w_lin[inv]
        cm = np.bincount(flat, weights=w, minlength=len(rows) * len(cols)).reshape(len(rows), len(cols))
        mats[b] = cm / np.clip(cm.sum(1, keepdims=True), 1, None)
    lo = pd.DataFrame(np.nanpercentile(mats, 2.5, 0), index=rows, columns=cols)
    hi = pd.DataFrame(np.nanpercentile(mats, 97.5, 0), index=rows, columns=cols)
    return lo, hi


# --------------------------------------------------------------------------- dwell times and rates
def dwell_times(d: pd.DataFrame, units: Units) -> pd.DataFrame:
    rows = []
    for tid, g in d.sort_values(["track_id", "frame"]).groupby("track_id"):
        st = g["motion"].to_numpy()
        fr = g["frame"].to_numpy()
        i = 0
        while i < len(st):
            j = i
            while j + 1 < len(st) and st[j + 1] == st[i]:
                j += 1
            censored = (i == 0) or (j == len(st) - 1)
            rows.append((tid, st[i], fr[j] - fr[i] + 1, censored))
            i = j + 1
    dw = pd.DataFrame(rows, columns=["track_id", "state", "frames", "censored"])
    dw["minutes"] = units.frames_to_minutes(dw["frames"])
    return dw


def event_rates(d: pd.DataFrame, units: Units) -> dict:
    cell_hours = float(units.frames_to_hours(len(d)))
    ev = dict(
        cell_hours=cell_hours,
        divisions=int((d.n_children == 2).sum()),
        divisions_isolated_mother=int((d.get("division_quality", pd.Series("", index=d.index)) == "isolated").sum()),
        divisions_in_cluster_lowconf=int((d.get("division_quality", pd.Series("", index=d.index)) == "in_cluster").sum()),
        deaths=int((d.end_type == "died").sum()),
        left=int((d.end_type == "left").sum()),
        merged=int((d.end_type == "merged").sum()),
        lost=int((d.end_type == "lost").sum()),
        entered=int((d.start_type == "entered").sum()),
        appeared=int((d.start_type == "appeared").sum()),
    )
    for k in ["divisions", "divisions_isolated_mother", "divisions_in_cluster_lowconf", "deaths", "left", "merged", "lost", "entered", "appeared"]:
        ev[f"{k}_per_cell_hour"] = ev[k] / cell_hours if cell_hours else np.nan
    n_first = int((d.frame == d.frame.min()).sum())
    n_last = int((d.frame == d.frame.max()).sum())
    hours = float(units.frames_to_hours(d.frame.max() - d.frame.min()))
    ev.update(n_cells_first_frame=n_first, n_cells_last_frame=n_last, duration_hours=hours,
              population_doubling_time_hours=float(hours * np.log(2) / np.log(n_last / n_first)) if n_last > n_first else np.nan,
              division_doubling_time_hours=float(np.log(2) / ev["divisions_per_cell_hour"]) if ev["divisions_per_cell_hour"] else np.nan)
    return ev


# --------------------------------------------------------------------------- plotting
def heatmap(P: pd.DataFrame, C: pd.DataFrame, title: str, path: Path, lo=None, hi=None):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1.3 * len(P.columns) + 2, 0.8 * len(P.index) + 1.8))
    im = ax.imshow(P.to_numpy(dtype=float), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(P.columns)), P.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(P.index)), [f"{r}  (n={int(C.loc[r].sum())})" for r in P.index])
    for i, r in enumerate(P.index):
        for j, c in enumerate(P.columns):
            p = P.loc[r, c]
            if np.isfinite(p):
                txt = f"{p:.2f}" if p >= 0.005 else f"{p:.1e}" if p > 0 else "0"
                if lo is not None and p >= 0.005:
                    txt += f"\n[{lo.loc[r, c]:.2f},{hi.loc[r, c]:.2f}]"
                ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="white" if p < 0.6 else "black")
    ax.set_xlabel("state at t + lag"); ax.set_ylabel("state at t"); ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.03)
    plt.tight_layout(); plt.savefig(path, dpi=90); plt.close(fig)


# --------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--name", default="hoct")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    units = Units.from_config(cfg)
    scfg = cfg["states"]
    out = ensure_dir(Path(cfg["paths"]["analysis_dir"]) / f"transition_{args.name}")
    d = pd.read_parquet(Path(cfg["paths"]["analysis_dir"]) / f"states_{args.name}.parquet")
    conf_thr = cfg["tracking"]["high_confidence_prob"]
    use_bidir = "bidir" in d.columns and d["bidir"].notna().any()
    rng = np.random.default_rng(0)
    summary = dict(name=args.name, units=dict(frame_interval_min=units.frame_interval_min, pixel_size_um=units.pixel_size_um),
                   high_conf_rule=f"all edges on path prob >= {conf_thr}" + (" and bidirectionally consistent" if use_bidir else ""),
                   matrices={}, figures=[])

    entries = entry_transitions(d)
    win_h = scfg["time_window_hours"]
    d["window"] = (d["hours"] // win_h).astype(int)
    d_cl = cluster_as_state(d)
    entries_cl = entry_transitions(d_cl)
    for lag in scfg["lags"]:
        tr = enumerate_transitions(d, lag, conf_thr, use_bidir)
        tr = tr.merge(d[["node_id", "window"]].rename(columns={"node_id": "src"}), on="src", how="left")
        tr_cl = enumerate_transitions(d_cl, lag, conf_thr, use_bidir, row_states=ROW_STATES_CL)
        if lag == 1:
            entries["window"] = (d.set_index("node_id").loc[entries.src, "window"].to_numpy())
            tr = pd.concat([tr, entries], ignore_index=True)
            tr_cl = pd.concat([tr_cl, entries_cl], ignore_index=True)
        variants = {
            "all": tr,
            "highconf": tr[tr.high_conf],
            "isolated": tr[~tr.in_cluster.astype(bool)],
            "in_cluster": tr[tr.in_cluster.astype(bool)],
            "cluster_state": tr_cl,  # 'in_cluster' as its own state
        }
        for w, g in tr.groupby("window"):
            variants[f"window{w}_{w*win_h:.0f}-{(w+1)*win_h:.0f}h"] = g
        for vname, g in variants.items():
            if len(g) < 20:
                continue
            if vname == "cluster_state":
                C = count_matrix(g, rows=ROW_STATES_CL + ENTRY_STATES, cols=COL_STATES_CL)
            else:
                C = count_matrix(g)
            P = prob_matrix(C)
            tag = f"lag{lag}_{vname}"
            C.to_csv(out / f"counts_{tag}.csv")
            P.to_csv(out / f"probs_{tag}.csv")
            lo = hi = None
            if vname in ("all", "highconf", "cluster_state"):
                lo, hi = bootstrap_ci(g, scfg["bootstrap_n"], rng, list(C.index), cols=list(C.columns))
                lo.to_csv(out / f"probs_{tag}_ci_lo.csv"); hi.to_csv(out / f"probs_{tag}_ci_hi.csv")
            minutes = lag * units.frame_interval_min
            fig_p = out / f"heatmap_{tag}.png"
            heatmap(P, C, f"{args.name}: lag {lag} frames ({minutes:g} min), {vname}, n={len(g)}", fig_p, lo, hi)
            summary["matrices"][tag] = dict(n_transitions=int(len(g)), rows=list(C.index), cols=list(C.columns),
                                            counts=C.to_numpy().tolist(), probs=P.round(4).to_numpy().tolist())
            if vname in ("all", "cluster_state", "highconf", "isolated", "in_cluster") or vname.startswith("window"):
                summary["figures"].append(fig_p.name)

    dw = dwell_times(d, units)
    dw.to_parquet(out / "dwell_times.parquet", index=False)
    dws = dw[~dw.censored].groupby("state")["minutes"].describe()[["count", "mean", "50%", "75%"]]
    dws.to_csv(out / "dwell_time_summary.csv")
    summary["dwell_minutes_uncensored"] = dws.round(2).to_dict()
    summary["event_rates"] = event_rates(d, units)
    summary["state_fraction"] = d["state"].value_counts(normalize=True).round(4).to_dict()
    st_info = json.load(open(Path(cfg["paths"]["analysis_dir"]) / f"states_{args.name}.json"))
    summary["speed_thresholds_px_per_frame"] = st_info["thresholds"]["thresholds"]
    summary["divisions_linker"] = st_info.get("divisions_linker", {})
    with open(out / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    write_report(cfg, args.name, summary, out, d, units)
    log.info("Wrote %s", out)


def _validation_html(val: dict) -> str:
    """Compact tables for the validation results."""
    parts = []
    bk = val.get("bookkeeping", {})
    if bk:
        keys = ["n_first", "n_last", "divisions_total", "interior_appearances_total", "interior_ends_total",
                "border_appearances_total", "border_ends_total", "interior_appearance_rate_per_cell_frame",
                "interior_end_rate_per_cell_frame", "median_track_length_frames", "mean_track_length_frames",
                "outlier_frames_count_dev_gt_8pct"]
        parts.append("<h3>Population bookkeeping (raw linker graph, before refinement)</h3><ul>" +
                     "".join(f"<li>{k}: {bk[k]:.4g}</li>" if isinstance(bk.get(k), float) else f"<li>{k}: {bk.get(k)}</li>" for k in keys if k in bk) + "</ul>")
    tf = val.get("track_free", {})
    if tf:
        parts.append("<h3>Tracking-free step-length estimates (px/frame)</h3><ul>" +
                     "".join(f"<li>{k}: {v:.3g}</li>" for k, v in tf.items() if isinstance(v, (int, float))) + "</ul>"
                     "<p>The pair-displacement excess is an upper bound on the true step length (neighbours within a cluster also "
                     "contribute); linked steps below it indicate that the longest jumps are the ones most often broken.</p>")
    fs = val.get("frame_skip", {})
    if fs:
        rows = [f"<tr><td>{a}-{b}</td><td>{p.get('agreement_overall', float('nan')):.3f}</td>"
                + "".join(f"<td>{p.get('agreement_by_density', {}).get(k, float('nan')):.2f}</td>" for k in ["dense", "medium", "sparse"]) + "</tr>"
                for (a, b), p in zip(fs["windows"], fs["per_window"]) if p]
        parts.append(f"<h3>Frame-skip stability (same linker at half frame rate reproduces 2-step links): overall {fs['agreement_overall']:.3f}</h3>"
                     "<table><tr><th>frames</th><th>agreement</th><th>dense</th><th>medium</th><th>sparse</th></tr>" + "".join(rows) + "</table>")
    trf = val.get("time_reversal_full")
    tr = val.get("time_reversal", {})
    if trf:
        parts.append(f"<h3>Whole-movie time-reversal consistency: {trf['bidirectional_fraction']:.3f} of {trf['n_edges']} links</h3>"
                     "<table><tr><th>by density</th>" + "".join(f"<td>{k}: {v:.2f}</td>" for k, v in trf["bidir_by_density"].items()) + "</tr>"
                     "<tr><th>by time</th>" + "".join(f"<td>{k}: {v:.2f}</td>" for k, v in trf["bidir_by_6h_window"].items()) + "</tr></table>"
                     "<p>This flag ('bidir') defines the high-confidence link subset used in the 'highconf' matrices.</p>")
    elif tr:
        parts.append(f"<h3>Time-reversal consistency (3 windows): {tr['bidirectional_fraction']:.3f}</h3>")
    mg = val.get("mini_gt", {})
    if mg:
        head = "<tr><th>crop</th><th>run</th><th>edge F1</th><th>precision</th><th>recall</th><th>CTC TRA</th><th>CTC LNK</th><th>pred/consensus divisions</th></tr>"
        rows = []
        for crop, r in mg.items():
            for nm, s in r["scores"].items():
                rows.append(f"<tr><td>{crop} (n={r['n_nodes']}, consensus edges={r['n_consensus_edges']})</td><td>{nm}</td>"
                            f"<td>{s['edge_f1']:.3f}</td><td>{s['edge_precision']:.3f}</td><td>{s['edge_recall']:.3f}</td>"
                            f"<td>{s.get('ctc_TRA', float('nan')):.3f}</td><td>{s.get('ctc_LNK', float('nan')):.3f}</td><td>{s['n_pred_div']}/{s['n_gt_div']}</td></tr>")
        parts.append("<h3>Consensus mini ground truth (majority of Trackastra greedy / ILP / reversed, LapTrack, HOCT on 512 px crops, 20 frames)</h3>"
                     "<p>Biased towards the voters (three are Trackastra variants) and only covers links the linkers agree on: use as a ranking, "
                     "not an absolute accuracy. 'full:' rows are the full-movie results restricted to the crop.</p>"
                     "<table>" + head + "".join(rows) + "</table>")
    return "\n".join(parts)


def write_report(cfg, name, summary, out: Path, d: pd.DataFrame, units: Units):
    def img(p: Path):
        if not p.exists():
            return ""
        b = base64.b64encode(p.read_bytes()).decode()
        return f'<img src="data:image/png;base64,{b}" style="max-width:100%">'

    val_dir = Path(cfg["paths"]["analysis_dir"]) / "validation"
    val_json = val_dir / f"validation_{name}.json"
    val = json.load(open(val_json)) if val_json.exists() else {}
    mot_dir = Path(cfg["paths"]["analysis_dir"]) / f"motility_{name}"
    mot_json = mot_dir / "summary.json"
    mot = json.load(open(mot_json)) if mot_json.exists() else {}
    ev = summary["event_rates"]
    thr = summary["speed_thresholds_px_per_frame"]
    unit_note = ("Pixel size unknown: lengths in px, speeds in px/frame (1 frame = "
                 f"{units.frame_interval_min} min).") if not units.has_pixel_size else f"Pixel size {units.pixel_size_um} um."
    rows = []
    rows.append(f"<h1>Cell tracking report: {name}</h1><p>{unit_note}</p>")
    rows.append("<h2>Population and event rates</h2><ul>")
    for k in ["n_cells_first_frame", "n_cells_last_frame", "duration_hours", "population_doubling_time_hours",
              "divisions", "divisions_isolated_mother", "divisions_in_cluster_lowconf", "divisions_per_cell_hour", "division_doubling_time_hours",
              "deaths", "deaths_per_cell_hour", "left", "entered", "merged", "lost", "appeared", "cell_hours"]:
        v = ev.get(k)
        rows.append(f"<li>{k}: {v:.4g}</li>" if isinstance(v, float) else f"<li>{k}: {v}</li>")
    rows.append("</ul>")
    rows.append(f"<p>Motion-state thresholds (px/frame): stationary &lt; {thr['threshold_stationary_slow']:.1f} &lt; slow &lt; "
                f"{thr['threshold_slow_fast']:.1f} &lt; fast  (method: {thr['method']}). State fractions: {summary['state_fraction']}</p>")
    rows.append(f"<p>High-confidence subset rule: {summary['high_conf_rule']}.</p>")
    rows.append("<h2>State-transition matrices</h2><p>Row-normalised probabilities; brackets are 95% bootstrap CIs over lineages. "
                "Rows 'entered'/'appeared' (lag 1 only) show the first state after a cell enters the field or appears unexplained. "
                "Columns 'left'/'died'/'merged'/'lost' are absorbing track ends within the lag. Variants: 'all'; 'cluster_state' "
                "(cells touching >= 2 neighbours are in state 'in_cluster' instead of a motion state); 'highconf' (only high-confidence links); "
                "'isolated'/'in_cluster' (source cell outside / inside a cluster); and per time window.</p>")
    for f in summary["figures"]:
        rows.append(img(out / f))
    rows.append("<h2>Dwell times (uncensored runs, minutes)</h2>")
    rows.append(pd.DataFrame(summary["dwell_minutes_uncensored"]).to_html(float_format=lambda x: f"{x:.1f}"))
    if mot:
        rows.append("<h2>Motility</h2><pre>" + json.dumps(mot, indent=1, default=float)[:6000] + "</pre>")
        for f in ["msd.png", "speed_by_time.png", "turning_angles.png", "velocity_autocorrelation.png"]:
            rows.append(img(mot_dir / f))
    if val:
        rows.append("<h2>Validation (no manual ground truth)</h2>")
        rows.append(_validation_html(val))
        rows.append("<details><summary>raw validation JSON</summary><pre>" + json.dumps(val, indent=1, default=float) + "</pre></details>")
        for f in ["bookkeeping.png", "step_distribution_vs_trackfree.png", "reversal_full_density.png", "reversal_full_time.png"]:
            rows.append(img(val_dir / f))
        for pat in ["frame_skip_density_f*.png", "frame_skip_step_f*.png", "reversal_density_f*.png"]:
            for f in sorted(val_dir.glob(pat), key=lambda q: int(q.stem.rsplit("_f", 1)[1])):
                rows.append(img(f))
    cmp_p = Path(cfg["paths"]["analysis_dir"]) / "linker_comparison_f0-113.csv"
    if cmp_p.exists():
        rows.append("<h2>Linker comparison (frames 0-112, after identical refinement)</h2>"
                    "<p>Fewer tracks / longer median track / fewer unexplained appearances = less fragmentation of the same detections.</p>"
                    + pd.read_csv(cmp_p).to_html(index=False))
    movies = sorted((Path(cfg["paths"]["qc_dir"]) / "movies").glob(f"{name}_*.mp4"))
    if movies:
        rows.append("<h2>Overlay movies</h2><ul>" + "".join(f"<li>{m.relative_to(Path(cfg['paths']['qc_dir']).parent)}</li>" for m in movies) + "</ul>"
                    "<p>Trails = last 8 frames; marker colour = state (blue stationary, yellow slow, red fast, green dividing, magenta dying).</p>")
    qc = Path(cfg["paths"]["qc_dir"]) / "segmentation_qc.png"
    rows.append("<h2>Segmentation QC</h2>" + img(qc))
    dl = summary.get("divisions_linker", {})
    rows.append("<h2>Caveats</h2><ul>"
                "<li>Divisions whose mother touched another cell ('in_cluster', edge confidence 0.3) are low-confidence: visual spot checks suggest roughly half are "
                "cluster splits. They are excluded from the 'highconf' matrices; the 'isolated' mothers are the reliable set. Expected true count lies between the two "
                "(population bookkeeping: divisions - deaths = growth - entered + left).</li>"
                "<li>Single-frame links inside dense clusters are the least reliable; prefer lag-4/lag-12, the 'isolated' and the high-confidence matrices.</li>"
                "<li>Linking is 1-to-1 (every linker 'division' is detached and re-examined). Divisions are then proposed only when an orphan "
                "object appears next to a rounded, isolated, normally sized mother that has been tracked for several frames and the two daughters "
                "are of similar, sub-mother size, lie close together and conserve the mother's area. This is deliberately conservative: "
                f"{dl.get('divisions', {}).get('accepted', '?')} divisions accepted out of {dl.get('divisions', {}).get('candidates', '?')} orphan candidates "
                f"({dl.get('linker_divisions_detached', '?')} linker 'divisions' detached, {dl.get('merge_split_repaired', '?')} merge/split artefacts repaired by gap edges).</li>"
                "<li>'merged' = the track ends because the segmentation fused the cell with a touching neighbour; 'lost' = unexplained interior end (segmentation miss or linking failure).</li>"
                "<li>'died' is a heuristic: an isolated cell that stopped moving and shows at least two more apoptotic signatures (shrinkage, brightening, rounding, nearby debris). "
                "Corpses that remain visible are counted 'stationary' until no longer detected. Expect under-detection.</li>"
                "<li>Pixel size is unknown: all lengths are in pixels. Supply acquisition.pixel_size_um in configs/default.yaml to get micrometres.</li></ul>")
    (out / "report.html").write_text("<html><head><meta charset='utf-8'><style>body{font-family:sans-serif;max-width:1200px;margin:auto}"
                                     "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:3px 6px}</style></head><body>"
                                     + "\n".join(rows) + "</body></html>")


if __name__ == "__main__":
    main()
