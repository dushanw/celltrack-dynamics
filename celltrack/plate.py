"""Plate-level aggregation of per-well results.

Reads the synced per-well outputs (outputs/wells/<well>/...) and the experiment plate map
(xlsx from the data share), and writes tidy tables + plate figures to outputs/wells/plate/:

    plate_map.csv                 well -> group / condition / density / replicate / control well
    well_summary.csv              one row per well x time window, calibrated to um and minutes
    transition_global_long.csv    per-well transition matrices (stationary/slow/fast) at lag 1/4/12
                                  using ONE pooled speed-threshold pair for the whole plate
    speed_vs_density.csv          median step vs local density (nearest-neighbour distance bins
                                  and cells/mm2 bins) per well
    frame_series.csv              per-well per-frame cell count, median area/intensity (drift QC)
    plate_heatmap_<metric>.png    8 x 6 plate maps
    condition_<metric>.png        replicate-level dot plots per condition (coating x density; drugs)
    plate_thresholds.json         pooled thresholds + calibration used

The per-well pipeline fits speed-state thresholds per well, which makes "fast fraction" and
transition matrices incomparable across wells; here they are recomputed with pooled thresholds.
Run:  python -m celltrack.plate --wells outputs/wells --platemap <xlsx> --pixel-size 0.5
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .common import ensure_dir, setup_logging
from .states import fit_speed_thresholds

log = setup_logging("plate")

WINDOWS = {"0-6h": (0, 6), "6-12h": (6, 12), "12-18h": (12, 18), "18-24h": (18, 24), "all": (0, 1e9)}
STATES = ["stationary", "slow", "fast"]
LAGS = (1, 4, 12)
NODE_COLS = ["frame", "track_id", "step_px", "edge_conf_in", "in_cluster", "touching", "area", "eccentricity",
             "solidity", "circularity", "intensity_mean", "nn_dist", "near_border", "division_confirmed",
             "end_type", "state", "bidir", "hours", "n_children"]


# --------------------------------------------------------------------------- plate map
def parse_platemap(xlsx: Path) -> pd.DataFrame:
    """Parse the lab's plate-map sheet: a grid with row letters in one column, condition names in
    the six columns to its right, and a seeding-density label in the column to its left that applies
    to that row and the following one. Rows whose conditions are drugs are the treatment group."""
    import openpyxl

    ws = openpyxl.load_workbook(xlsx, data_only=True).worksheets[0]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    recs, density, notes = [], None, []
    coat_words = ("coating",)
    for r in rows:
        cells = ["" if c is None else str(c).strip() for c in r]
        letter_idx = next((i for i, c in enumerate(cells) if re.fullmatch(r"[A-H]", c)), None)
        if letter_idx is None:
            notes += [c for c in cells if len(c) > 12]
            continue
        label = cells[letter_idx - 1] if letter_idx > 0 else ""
        m = re.search(r"(\d+)\s*cell", label.replace(",", ""))
        if m:
            density = int(m.group(1))
        row_letter = cells[letter_idx]
        conds = cells[letter_idx + 1: letter_idx + 7]
        is_coating = all(any(w in c.lower() for w in coat_words) for c in conds if c)
        for j, cond in enumerate(conds, start=1):
            if not cond:
                continue
            recs.append(dict(well=f"well_{row_letter}{j}", row=row_letter, col=j, density=density,
                             group="coating" if is_coating else "drug",
                             condition=cond if cond.lower().startswith("no ") else cond.replace(" coating", "").strip()))
    pm = pd.DataFrame(recs).sort_values(["row", "col"]).reset_index(drop=True)
    # replicate index: wells sharing (group, condition, density) ordered by row
    pm["replicate"] = pm.groupby(["group", "condition", "density"]).cumcount() + 1
    # matched control: "No coating" wells at the same seeding density (E1/F1 for the drug rows)
    ctrl = pm[(pm.group == "coating") & (pm.condition.str.lower().str.startswith("no "))]
    ctrl_by_density = {d: ",".join(sorted(g.well)) for d, g in ctrl.groupby("density")}
    pm["control_wells"] = pm.density.map(ctrl_by_density).astype(object)
    pm.loc[pm.well.isin(ctrl.well), "control_wells"] = ""
    pm.attrs["notes"] = notes
    return pm.sort_values(["row", "col"]).reset_index(drop=True)


# --------------------------------------------------------------------------- per-well loading
def load_nodes(well_dir: Path) -> pd.DataFrame:
    import pyarrow.parquet as pq

    p = well_dir / "analysis" / "states_hoct.parquet"
    have = set(pq.ParquetFile(p).schema.names)
    d = pd.read_parquet(p, columns=[c for c in NODE_COLS + ["start_type"] if c in have])
    if "bidir" in d:  # object column with True/False/None -> float with NaN
        d["bidir"] = pd.to_numeric(d.bidir.map({True: 1.0, False: 0.0}), errors="coerce")
    return d


def _json(p: Path, default=None):
    try:
        return json.load(open(p))
    except FileNotFoundError:
        return default


def _lin_trend(y: np.ndarray) -> float:
    """Relative linear change over the series (end - start) / mean, from a least-squares fit."""
    y = np.asarray(y, float)
    if len(y) < 4 or not np.isfinite(y).all() or y.mean() == 0:
        return float("nan")
    x = np.arange(len(y))
    slope = np.polyfit(x, y, 1)[0]
    return float(slope * (len(y) - 1) / y.mean())


# --------------------------------------------------------------------------- metrics
def well_metrics(well: str, nodes: pd.DataFrame, well_dir: Path, t1: float, t2: float, px: float, dt: float,
                 fov_mm2: float, hc_prob: float) -> tuple[list[dict], list[dict], list[dict], pd.DataFrame]:
    """Per-window summary rows, long-format global transition rows, speed-vs-density rows, frame series."""
    um_per_min = px / dt
    nodes = nodes.copy()
    nodes["gstate"] = pd.cut(nodes.step_px, [-np.inf, t1, t2, np.inf], labels=STATES)
    # state at t+lag for the same track (for transition matrices)
    key = nodes.set_index(["track_id", "frame"])["gstate"]
    for lag in LAGS:
        idx = pd.MultiIndex.from_arrays([nodes.track_id, nodes.frame + lag])
        nodes[f"gstate_lag{lag}"] = key.reindex(idx).to_numpy()

    summary = _json(well_dir / "tracks" / "hoct_summary.json", {})
    val = _json(well_dir / "analysis" / "validation" / "validation_hoct.json", {})
    mot = _json(well_dir / "analysis" / "motility_hoct" / "summary.json", {})
    sjson = _json(well_dir / "analysis" / "states_hoct.json", {})
    thr = sjson.get("thresholds", {}).get("thresholds", {})
    per_track = pd.read_parquet(well_dir / "analysis" / "motility_hoct" / "per_track.parquet")
    book = pd.read_csv(well_dir / "analysis" / "validation" / "bookkeeping_per_frame.csv")

    frames_per_h = 60 / dt
    counts = nodes.groupby("frame").size()
    frame_series = pd.DataFrame({
        "well": well, "frame": counts.index, "hours": counts.index * dt / 60, "n_cells": counts.values,
        "median_area_um2": nodes.groupby("frame").area.median().values * px ** 2,
        "median_intensity": nodes.groupby("frame").intensity_mean.median().values,
        "median_step_um_min": nodes.groupby("frame").step_px.median().values * um_per_min,
    })

    rows, trans_rows, dens_rows = [], [], []
    for wname, (h0, h1) in WINDOWS.items():
        w = nodes[(nodes.hours >= h0) & (nodes.hours < h1)]
        if len(w) == 0:
            continue
        nfr = w.frame.nunique()
        cells_per_frame = len(w) / nfr
        cell_hours = len(w) * dt / 60
        steps = w.step_px.dropna()
        hc = w[(w.edge_conf_in >= hc_prob)].step_px.dropna()
        r = dict(well=well, window=wname, n_frames=nfr, cells_per_frame=cells_per_frame,
                 cells_per_mm2=cells_per_frame / fov_mm2,
                 cluster_fraction=float(w.in_cluster.mean()), touching_fraction=float(w.touching.mean()),
                 median_nn_dist_um=float(w.nn_dist.median() * px),
                 # movement
                 median_step_um_min=float(steps.median() * um_per_min),
                 mean_step_um_min=float(steps.mean() * um_per_min),
                 hc_median_step_um_min=float(hc.median() * um_per_min) if len(hc) else np.nan,
                 hc_step_fraction=float(len(hc) / max(1, len(steps))),
                 p90_step_um_min=float(steps.quantile(0.9) * um_per_min),
                 # global-threshold occupancy
                 frac_stationary=float((w.gstate == "stationary").mean()),
                 frac_slow=float((w.gstate == "slow").mean()),
                 frac_fast=float((w.gstate == "fast").mean()),
                 # per-well-threshold occupancy (as produced by the per-well pipeline)
                 frac_fast_wellthr=float((w.state == "fast").mean()),
                 # events
                 divisions=int(w.division_confirmed.astype(bool).sum()),
                 deaths=int((w.end_type == "died").sum()),
                 lost_ends=int((w.end_type == "lost").sum()),
                 appeared_interior=int((w.start_type == "appeared").sum()),
                 # morphology
                 median_area_um2=float(w.area.median() * px ** 2), median_eccentricity=float(w.eccentricity.median()),
                 median_solidity=float(w.solidity.median()), median_circularity=float(w.circularity.median()),
                 median_intensity=float(w.intensity_mean.median()),
                 # QC
                 bidir_fraction=float(w.bidir.mean()) if "bidir" in w and w.bidir.notna().any() else np.nan,
                 mean_edge_conf=float(w.edge_conf_in.mean()),
                 border_fraction=float(w.near_border.mean()),
                 count_trend_rel=_lin_trend(w.groupby("frame").size().values),
                 intensity_trend_rel=_lin_trend(w.groupby("frame").intensity_mean.median().values),
                 area_trend_rel=_lin_trend(w.groupby("frame").area.median().values),
                 )
        r["division_rate_per_1000cells_h"] = r["divisions"] / max(1e-9, cell_hours) * 1000
        r["death_rate_per_1000cells_h"] = r["deaths"] / max(1e-9, cell_hours) * 1000
        r["lost_rate_per_cell_frame"] = r["lost_ends"] / max(1, len(w))
        # MSD / D / alpha from the per-well motility tables (ensemble MSD, px^2 -> um^2)
        msd_file = well_dir / "analysis" / "motility_hoct" / ("msd_all.csv" if wname == "all" else f"msd_window_{wname}.csv")
        if msd_file.exists():
            m = pd.read_csv(msd_file)
            m = m[m.n_pairs > 50]
            for lag, nm in ((1, "msd_2p5min_um2"), (4, "msd_10min_um2"), (12, "msd_30min_um2")):
                v = m.loc[m.lag_frames == lag, "msd_px2"]
                r[nm] = float(v.iloc[0] * px ** 2) if len(v) else np.nan
            fit = m[(m.lag_frames >= 1) & (m.lag_frames <= 10)]
            if len(fit) >= 3:
                a, _ = np.polyfit(np.log(fit.lag_frames), np.log(fit.msd_px2), 1)
                r["msd_alpha"] = float(a)
                slope = np.polyfit(fit.lag_frames * dt, fit.msd_px2 * px ** 2, 1)[0]  # um^2 / min
                r["D_um2_per_min"] = float(slope / 4)
        pw = (mot.get("persistence") or {}).get("all" if wname == "all" else f"window_{wname}") or {}
        r["persistence_time_min"] = pw.get("persistence_minutes", np.nan)
        # per-track persistence / straightness for tracks starting in the window
        pt = per_track[(per_track.start_hours >= h0) & (per_track.start_hours < h1)]
        if len(pt):
            r["median_straightness"] = float(pt.straightness.median())
            r["median_dir_persistence"] = float(pt.directional_persistence.median())
            r["n_tracks_ge5"] = int(len(pt))
        # local-density dependence of speed (nn-distance bins) within the window
        if wname == "all":
            bins_um = [0, 10, 20, 30, 50, 80, 1e9]
            nb = pd.cut(w.nn_dist * px, bins_um, labels=["0-10", "10-20", "20-30", "30-50", "50-80", ">80"])
            for b, g in w.groupby(nb, observed=True):
                s = g.step_px.dropna()
                if len(s) >= 200:
                    dens_rows.append(dict(well=well, kind="nn_dist_um", bin=str(b), n_steps=len(s),
                                          median_step_um_min=float(s.median() * um_per_min),
                                          frac_fast=float((g.gstate == "fast").mean())))
            fc = w.groupby("frame").size()
            fd = (fc / fov_mm2)
            dbins = [0, 100, 200, 300, 400, 600, 800, 1200, 1e9]
            lab = pd.cut(fd, dbins, labels=[f"{dbins[i]}-{dbins[i+1]}" if dbins[i+1] < 1e8 else f">{dbins[i]}" for i in range(len(dbins) - 1)])
            w2 = w.assign(dbin=w.frame.map(lab))
            for b, g in w2.groupby("dbin", observed=True):
                s = g.step_px.dropna()
                if len(s) >= 200:
                    dens_rows.append(dict(well=well, kind="cells_per_mm2", bin=str(b), n_steps=len(s),
                                          median_step_um_min=float(s.median() * um_per_min),
                                          frac_fast=float((g.gstate == "fast").mean())))
        # transition matrices with pooled thresholds
        for lag in LAGS:
            a, b = w.gstate, w[f"gstate_lag{lag}"]
            ok = a.notna() & pd.notna(b)
            ct = pd.crosstab(a[ok], pd.Categorical(b[ok], categories=STATES)).reindex(index=STATES, columns=STATES, fill_value=0)
            P = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0)
            for s_from in STATES:
                for s_to in STATES:
                    trans_rows.append(dict(well=well, window=wname, lag_frames=lag, lag_min=lag * dt, from_state=s_from,
                                           to_state=s_to, n=int(ct.loc[s_from, s_to]), p=float(P.loc[s_from, s_to])))
        rows.append(r)

    # whole-movie only extras
    allrow = next(r for r in rows if r["window"] == "all")
    tf = val.get("track_free", {})
    allrow.update(
        trackfree_pair_excess_step_um_min=tf.get("pair_excess_median_step_px", np.nan) * um_per_min,
        trackfree_track_step_um_min=tf.get("track_median_step_px", np.nan) * um_per_min,
        bidir_fraction_full=val.get("time_reversal_full", {}).get("bidirectional_fraction", np.nan),
        linked_fraction=summary.get("linked_fraction", np.nan),
        n_tracks=summary.get("track_starts", np.nan),
        wellthr_t1_um_min=thr.get("threshold_stationary_slow", np.nan) * um_per_min,
        wellthr_t2_um_min=thr.get("threshold_slow_fast", np.nan) * um_per_min,
        n_first=int(book.n.iloc[0]), n_last=int(book.n.iloc[-1]),
        outlier_frames=len(val.get("bookkeeping", {}).get("outlier_frames_count_dev_gt_8pct", [])),
    )
    # doubling time from log-linear fit of the cell count (includes entering/leaving -> apparent)
    bk = book[book.n > 0]
    if len(bk) > 10:
        slope = np.polyfit(bk.hours, np.log(bk.n), 1)[0]
        allrow["apparent_doubling_time_h"] = float(np.log(2) / slope) if slope > 0 else np.inf
        allrow["count_fold_change"] = float(bk.n.iloc[-1] / bk.n.iloc[0])
    return rows, trans_rows, dens_rows, frame_series


# --------------------------------------------------------------------------- figures
def _plate_grid(summary_all: pd.DataFrame, pm: pd.DataFrame, metric: str) -> np.ndarray:
    g = np.full((8, 6), np.nan)
    d = summary_all.set_index("well")[metric]
    for _, r in pm.iterrows():
        if r.well in d.index:
            g["ABCDEFGH".index(r.row), r.col - 1] = d[r.well]
    return g


def plate_heatmaps(summary_all: pd.DataFrame, pm: pd.DataFrame, out: Path, metrics: dict[str, str]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for metric, label in metrics.items():
        if metric not in summary_all:
            continue
        g = _plate_grid(summary_all, pm, metric)
        fig, ax = plt.subplots(figsize=(7.2, 7.8))
        im = ax.imshow(g, cmap="viridis", aspect="equal")
        ax.set_xticks(range(6)); ax.set_xticklabels(range(1, 7))
        ax.set_yticks(range(8)); ax.set_yticklabels(list("ABCDEFGH"))
        for i in range(8):
            for j in range(6):
                if np.isfinite(g[i, j]):
                    v = g[i, j]
                    txt = f"{v:.2f}" if abs(v) < 10 else f"{v:.0f}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                            color="white" if (v - np.nanmin(g)) / (np.nanmax(g) - np.nanmin(g) + 1e-9) < 0.6 else "black")
        # condition labels under columns and density/group labels at rows
        cols_coat = pm[pm.group == "coating"].drop_duplicates("col").set_index("col").condition
        cols_drug = pm[pm.group == "drug"].drop_duplicates("col").set_index("col").condition
        ax.set_xticks(range(6))
        ax.set_xticklabels([f"{j+1}\n{cols_coat.get(j+1, '')}\n({cols_drug.get(j+1, '')})" for j in range(6)], fontsize=7)
        row_lab = {r: f"{r}  {pm[pm.row == r].group.iloc[0]} {pm[pm.row == r].density.iloc[0]//1000}k" for r in "ABCDEFGH" if (pm.row == r).any()}
        ax.set_yticklabels([row_lab.get(r, r) for r in "ABCDEFGH"], fontsize=8)
        ax.set_title(f"{label}\n(whole movie, 0-24 h; rows C-H coating x density, rows A-B drugs at 10k)", fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.7, label=label)
        fig.tight_layout()
        fig.savefig(out / f"plate_heatmap_{metric}.png", dpi=130)
        plt.close(fig)


def condition_plots(summary_all: pd.DataFrame, pm: pd.DataFrame, out: Path, metrics: dict[str, str]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = summary_all if "group" in summary_all else summary_all.merge(pm, on="well")
    coat = d[d.group == "coating"]
    drug = d[d.group == "drug"]
    ctrl10k = coat[(coat.density == 10000) & coat.condition.str.lower().str.startswith("no")]
    coat_order = list(pm[pm.group == "coating"].sort_values("col").condition.unique())
    drug_order = ["Control (E1,F1)"] + list(pm[pm.group == "drug"].sort_values("col").condition.unique())
    for metric, label in metrics.items():
        if metric not in d:
            continue
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[6, 6.5]))
        ax = axes[0]
        for k, dens in enumerate(sorted(coat.density.unique())):
            sub = coat[coat.density == dens]
            means = sub.groupby("condition")[metric].mean().reindex(coat_order)
            xs = np.arange(len(coat_order)) + (k - 1) * 0.22
            ax.plot(xs, means.values, "-", alpha=0.5, color=f"C{k}")
            for i, c in enumerate(coat_order):
                ys = sub[sub.condition == c][metric].values
                ax.scatter(np.full(len(ys), xs[i]), ys, color=f"C{k}", s=28, label=f"{dens//1000}k cells/well" if i == 0 else None, zorder=3)
        ax.set_xticks(range(len(coat_order))); ax.set_xticklabels(coat_order, rotation=20, ha="right", fontsize=8)
        ax.set_ylabel(label); ax.set_title("Coating x seeding density (points = replicate wells)", fontsize=10)
        ax.legend(fontsize=8, frameon=False); ax.grid(axis="y", alpha=0.3)
        ax = axes[1]
        groups = [("Control (E1,F1)", ctrl10k[metric].values)] + [(c, drug[drug.condition == c][metric].values) for c in drug_order[1:]]
        for i, (c, ys) in enumerate(groups):
            ax.scatter(np.full(len(ys), i), ys, color="C3" if i else "C0", s=28, zorder=3)
            if len(ys):
                ax.hlines(np.mean(ys), i - 0.25, i + 0.25, color="C3" if i else "C0", alpha=0.7)
        ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[0] for g in groups], rotation=20, ha="right", fontsize=8)
        ax.set_title("Drug treatments at 10k cells/well, no coating (vs matched control)", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        fig.suptitle(f"{label} - whole movie (0-24 h)", fontsize=11)
        fig.tight_layout()
        fig.savefig(out / f"condition_{metric}.png", dpi=130)
        plt.close(fig)


PLOT_METRICS = {
    "cells_per_mm2": "Cell density (cells / mm2, mean over movie)",
    "hc_median_step_um_min": "Median step speed, high-confidence links (um/min)",
    "msd_10min_um2": "Ensemble MSD at 10 min lag (um2)",
    "D_um2_per_min": "Diffusion coefficient D from MSD slope (um2/min)",
    "msd_alpha": "MSD exponent alpha (lags 2.5-25 min)",
    "frac_fast": "Fraction of cell-frames in 'fast' state (pooled thresholds)",
    "frac_stationary": "Fraction of cell-frames 'stationary' (pooled thresholds)",
    "cluster_fraction": "Fraction of cell-frames in clusters",
    "division_rate_per_1000cells_h": "Confirmed divisions per 1000 cells per hour",
    "median_area_um2": "Median cell area (um2)",
    "median_straightness": "Median track straightness (net / path)",
    "bidir_fraction_full": "Tracking QC: forward/reverse link agreement",
    "count_fold_change": "Cell count fold change (last / first frame)",
    "trackfree_pair_excess_step_um_min": "Tracking-free step estimate, pair-excess (um/min)",
}


# --------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wells", default="outputs/wells", help="directory with well_*/ folders")
    ap.add_argument("--platemap", required=True, help="plate-map xlsx")
    ap.add_argument("--out", default=None, help="default: <wells>/plate")
    ap.add_argument("--pixel-size", type=float, default=0.5, help="um per pixel")
    ap.add_argument("--frame-interval", type=float, default=2.5, help="minutes per frame")
    ap.add_argument("--hc-prob", type=float, default=0.5, help="edge probability for the high-confidence subset")
    ap.add_argument("--sample-per-well", type=int, default=30000, help="steps per well for the pooled threshold fit")
    ap.add_argument("--only", nargs="*", default=None, help="restrict to these wells (debug)")
    args = ap.parse_args(argv)

    wells_dir = Path(args.wells)
    out = ensure_dir(args.out or wells_dir / "plate")
    px, dt = args.pixel_size, args.frame_interval
    pm = parse_platemap(Path(args.platemap))
    pm.to_csv(out / "plate_map.csv", index=False)
    log.info("plate map: %d wells, groups %s", len(pm), pm.group.value_counts().to_dict())

    well_dirs = sorted(p for p in wells_dir.glob("well_*") if (p / "analysis" / "states_hoct.parquet").exists())
    if args.only:
        well_dirs = [p for p in well_dirs if p.name in args.only]
    meta = _json(well_dirs[0] / "ingest_meta.json", {})
    H, W = meta.get("shape_yx", [3072, 3072])
    fov_mm2 = (H * px / 1000) * (W * px / 1000)
    log.info("%d wells; FOV %.3f mm2; %.2f um/px, %.1f min/frame", len(well_dirs), fov_mm2, px, dt)

    # pass 1: pooled speed thresholds (equal weight per well)
    rng = np.random.default_rng(0)
    pooled = []
    for wd in well_dirs:
        s = pd.read_parquet(wd / "analysis" / "states_hoct.parquet", columns=["step_px"]).step_px.dropna().to_numpy()
        pooled.append(rng.choice(s, min(args.sample_per_well, len(s)), replace=False))
    t1, t2, info = fit_speed_thresholds(np.concatenate(pooled))
    info.update(pixel_size_um=px, frame_interval_min=dt, t1_um_per_min=t1 * px / dt, t2_um_per_min=t2 * px / dt,
                fov_mm2=fov_mm2, wells=[p.name for p in well_dirs])
    json.dump(info, open(out / "plate_thresholds.json", "w"), indent=2)
    log.info("pooled thresholds: stationary<%.1f px (%.2f um/min) < slow < %.1f px (%.2f um/min) < fast  [%s]",
             t1, t1 * px / dt, t2, t2 * px / dt, info["method"])

    # pass 2: per-well metrics
    rows, trans, dens, series = [], [], [], []
    for wd in well_dirs:
        nodes = load_nodes(wd)
        r, t, d, fs = well_metrics(wd.name, nodes, wd, t1, t2, px, dt, fov_mm2, args.hc_prob)
        rows += r; trans += t; dens += d; series.append(fs)
        a = r[-1]
        log.info("%s: %.0f cells/frame, step %.2f um/min (hc %.2f), fast %.2f, div %.1f/1000c/h, bidir %.2f",
                 wd.name, a["cells_per_frame"], a["median_step_um_min"], a["hc_median_step_um_min"], a["frac_fast"],
                 a["division_rate_per_1000cells_h"], a.get("bidir_fraction_full", np.nan))
    summary = pd.DataFrame(rows).merge(pm, on="well", how="left")
    summary.to_csv(out / "well_summary.csv", index=False)
    pd.DataFrame(trans).merge(pm[["well", "group", "condition", "density", "replicate"]], on="well").to_csv(out / "transition_global_long.csv", index=False)
    pd.DataFrame(dens).merge(pm[["well", "group", "condition", "density"]], on="well").to_csv(out / "speed_vs_density.csv", index=False)
    pd.concat(series).to_csv(out / "frame_series.csv", index=False)

    summary_all = summary[summary.window == "all"]
    plate_heatmaps(summary_all, pm, out, PLOT_METRICS)
    condition_plots(summary_all, pm, out, PLOT_METRICS)
    log.info("Wrote %s (%d rows), figures for %d metrics", out / "well_summary.csv", len(summary), len(PLOT_METRICS))


if __name__ == "__main__":
    main()
