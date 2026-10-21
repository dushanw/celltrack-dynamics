"""Phase 4a: per-frame cell states and events from the lineage graph + detection features.

Motion states (from the step into the frame, px/frame):
    stationary | slow | fast          thresholds fitted from the step-speed distribution
Modifier:
    in_cluster                        >= k touching neighbours for >= m consecutive frames
Event states (override motion state on the event frame):
    dividing                          division parent frame (+/- window on parent/daughters flagged)
    dying                             last frame of a track that ends in the field with a death signature
Track boundary events (stored per node as start_type / end_type):
    start: initial | entered | divided_from | appeared
    end:   censored | left | divided | died | lost

Usage:
    python -m celltrack.states --name trackastra_ilp
Writes outputs/analysis/states_<name>.parquet and thresholds JSON.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .common import Units, ensure_dir, load_config, setup_logging
from .graphio import load_tracks

log = setup_logging("states")

MOTION = ["stationary", "slow", "fast"]


def fit_speed_thresholds(speed: np.ndarray) -> tuple[float, float, dict]:
    """Two thresholds from a 3-component Gaussian mixture on log(speed+1).
    Falls back to terciles if the mixture is degenerate."""
    v = speed[np.isfinite(speed) & (speed >= 0)]
    x = np.log1p(v)[:, None]
    info = {"method": "gmm3_log1p", "n": int(len(v))}
    try:
        from sklearn.mixture import GaussianMixture

        gm = GaussianMixture(3, random_state=0, n_init=3).fit(x)
        order = np.argsort(gm.means_.ravel())
        mu = gm.means_.ravel()[order]
        sd = np.sqrt(gm.covariances_.ravel()[order])
        w = gm.weights_[order]
        info.update(means_px=list(np.expm1(mu)), weights=list(w))
        grid = np.linspace(x.min(), x.max(), 2000)

        def comp(i):
            return w[i] * np.exp(-0.5 * ((grid - mu[i]) / sd[i]) ** 2) / sd[i]

        # thresholds = posterior crossover between neighbouring components
        t1 = grid[np.argmax(comp(1) > comp(0))]
        t2 = grid[np.argmax(comp(2) > comp(1))]
        t1, t2 = float(np.expm1(t1)), float(np.expm1(t2))
        if not (0 < t1 < t2) or w.min() < 0.03:
            raise ValueError("degenerate mixture")
    except Exception as e:  # noqa: BLE001
        t1, t2 = map(float, np.quantile(v, [1 / 3, 2 / 3]))
        info.update(method="terciles", reason=str(e)[:100])
    info.update(threshold_stationary_slow=t1, threshold_slow_fast=t2)
    return t1, t2, info


def consecutive_runs(flags: np.ndarray) -> np.ndarray:
    """Length of the run of True ending at each position."""
    out = np.zeros(len(flags), dtype=int)
    run = 0
    for i, f in enumerate(flags):
        run = run + 1 if f else 0
        out[i] = run
    return out


def assign_states(tracks: pd.DataFrame, det: pd.DataFrame, debris: pd.DataFrame | None, cfg, units: Units):
    scfg = cfg["states"]
    margin = cfg["segmentation"]["border_margin_px"]
    T0, T1 = int(tracks.frame.min()), int(tracks.frame.max())
    from .refine import mark_merges, refine_graph

    n_div_raw = int((tracks.n_children == 2).sum())
    tracks, div_stats = refine_graph(tracks, det, cfg)
    div_stats["linker_raw_divisions"] = n_div_raw
    log.info("Graph refinement: %s", div_stats)
    merged_flag = mark_merges(tracks, det, cfg)
    tracks["merged_end"] = merged_flag.to_numpy()
    d = tracks.merge(det, on=["frame", "label", "x", "y"], how="left", suffixes=("", "_det"))
    d = d.sort_values(["track_id", "frame"]).reset_index(drop=True)

    # ---------------------------------------------------------------- step speed
    pos = d.set_index("node_id")[["x", "y", "frame"]]
    has_parent = d.parent_node != -1
    pp = pos.reindex(d.loc[has_parent, "parent_node"].to_numpy())
    step = np.full(len(d), np.nan)
    dt = np.full(len(d), np.nan)
    step[has_parent.to_numpy()] = np.hypot(d.loc[has_parent, "x"].to_numpy() - pp["x"].to_numpy(),
                                           d.loc[has_parent, "y"].to_numpy() - pp["y"].to_numpy())
    dt[has_parent.to_numpy()] = d.loc[has_parent, "frame"].to_numpy() - pp["frame"].to_numpy()
    d["step_px"] = step
    d["speed_px_per_frame"] = step / dt
    # first node of a track: borrow the next step so it gets a motion state too
    nxt = d.groupby("track_id")["speed_px_per_frame"].shift(-1)
    d["speed_filled"] = d["speed_px_per_frame"].fillna(nxt)
    d["speed_um_per_min"] = units.speed_px_per_frame_to_um_per_min(d["speed_px_per_frame"])

    thr = scfg.get("speed_thresholds_px_per_frame")
    if thr:
        t1, t2 = thr
        thr_info = {"method": "config", "threshold_stationary_slow": t1, "threshold_slow_fast": t2}
    else:
        t1, t2, thr_info = fit_speed_thresholds(d["speed_px_per_frame"].to_numpy())
    log.info("Speed thresholds (px/frame): stationary<%.1f<slow<%.1f<fast  [%s]", t1, t2, thr_info["method"])
    d["motion"] = pd.cut(d["speed_filled"], [-np.inf, t1, t2, np.inf], labels=MOTION).astype(object)
    d.loc[d["speed_filled"].isna(), "motion"] = "stationary"  # isolated single-frame detections

    # ---------------------------------------------------------------- cluster modifier
    d["touching"] = d["n_touching"].fillna(0) >= scfg["cluster_min_neighbours"]
    run_len = d.groupby("track_id")["touching"].transform(
        lambda s: pd.Series(consecutive_runs(s.to_numpy()), index=s.index))
    d["in_cluster"] = run_len >= scfg["cluster_min_frames"]

    # ---------------------------------------------------------------- track starts / ends
    g = d.groupby("track_id")
    first_idx = g.head(1).index
    last_idx = g.tail(1).index
    d["start_type"] = ""
    d["end_type"] = ""
    f = d.loc[first_idx]
    st = np.where(f.parent_node != -1, "divided_from",
         np.where(f.frame == T0, "initial",
         np.where(f.dist_to_border < margin, "entered", "appeared")))
    d.loc[first_idx, "start_type"] = st
    l = d.loc[last_idx]
    et = np.where(l.n_children == 2, "divided",
         np.where(l.frame == T1, "censored",
         np.where(l.dist_to_border < margin, "left",
         np.where(l.merged_end, "merged", "interior_end"))))
    d.loc[last_idx, "end_type"] = et

    # ---------------------------------------------------------------- death signature for interior ends
    look = scfg["death_lookback_frames"]
    deb_by_frame = {f: gg[["x", "y"]].to_numpy() for f, gg in debris.groupby("frame")} if debris is not None and len(debris) else {}
    from scipy.spatial import cKDTree

    deb_tree = {f: cKDTree(v) for f, v in deb_by_frame.items()}
    death_rows = []
    for i in last_idx[d.loc[last_idx, "end_type"] == "interior_end"]:
        tid = d.at[i, "track_id"]
        tr = d[(d.track_id == tid)].tail(look + 1)
        if len(tr) < 3:
            score = dict(shrink=False, bright=False, still=False, round=False, debris=False)
        else:
            a = tr["area"].to_numpy()
            inten = tr["intensity_mean"].to_numpy()
            sp = tr["speed_px_per_frame"].to_numpy()
            circ = tr["circularity"].to_numpy()
            ref = slice(0, max(1, len(tr) // 2))
            score = dict(
                shrink=bool(a[-1] < 0.7 * np.nanmax(a[ref])),
                bright=bool(inten[-1] > 1.15 * np.nanmean(inten[ref])),
                still=bool(np.nanmean(sp[-2:]) < t1) if np.isfinite(np.nanmean(sp[-2:])) else False,
                round=bool(circ[-1] > np.nanmean(circ[ref]) + 0.1),
                debris=False,
            )
        fr, x, y = int(d.at[i, "frame"]), d.at[i, "x"], d.at[i, "y"]
        for ff in range(fr + 1, fr + 4):
            if ff in deb_tree and deb_tree[ff].query([x, y])[0] < 30:
                score["debris"] = True
                break
        # isolation: a cell that vanishes while touching / next to another cell has most
        # likely been merged into it by the segmentation, not died
        isolated = bool(d.at[i, "n_touching"] == 0 and d.at[i, "nn_dist"] > scfg.get("death_min_nn_dist_px", 45))
        n = sum(score.values())
        death_rows.append((i, n, score, isolated))
    # death = isolated cell that stopped moving AND shows >= 2 further apoptotic signatures
    # (shrinkage, brightening/condensation, rounding, debris nearby). Unexplained interior
    # ends are far more often linking failures than deaths, hence the conservative rule.
    min_score = scfg.get("death_min_score", 3)
    d["death_score"] = np.nan
    for i, n, score, isolated in death_rows:
        d.at[i, "death_score"] = n
        d.at[i, "end_type"] = "died" if (n >= min_score and score["still"] and isolated) else "lost"

    # ---------------------------------------------------------------- frames around divisions
    w = scfg["division_window_frames"]
    d["near_division"] = False
    for tid, gg in d.groupby("track_id"):
        if gg["end_type"].iloc[-1] == "divided":
            d.loc[gg.index[-w:], "near_division"] = True
        if gg["start_type"].iloc[0] == "divided_from":
            d.loc[gg.index[:w], "near_division"] = True

    # ---------------------------------------------------------------- final per-frame state
    d["state"] = d["motion"]
    d.loc[d.n_children == 2, "state"] = "dividing"
    d.loc[d.end_type == "died", "state"] = "dying"
    d["hours"] = units.frames_to_hours(d["frame"])
    return d, dict(thresholds=thr_info, t1=t1, t2=t2, divisions=div_stats)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--name", default="hoct")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    units = Units.from_config(cfg)
    out = ensure_dir(cfg["paths"]["analysis_dir"])
    det = pd.read_parquet(cfg["paths"]["detections"])
    deb_p = Path(cfg["paths"]["debris"])
    debris = pd.read_parquet(deb_p) if deb_p.exists() else None
    tracks, _ = load_tracks(cfg["paths"]["tracks_dir"], args.name)
    # attach bidirectional-consistency flags if the validation produced them
    conf_p = Path(cfg["paths"]["analysis_dir"]) / "validation" / f"edge_confidence_{args.name}.parquet"
    d, info = assign_states(tracks, det, debris, cfg, units)
    if conf_p.exists():
        ec = pd.read_parquet(conf_p).rename(columns={"v": "node_id"})[["node_id", "bidir", "weight_rev"]]
        d = d.merge(ec, on="node_id", how="left")
    d.to_parquet(out / f"states_{args.name}.parquet", index=False)
    summary = dict(
        thresholds=info,
        n_nodes=int(len(d)),
        state_counts=d["state"].value_counts().to_dict(),
        start_types=d.loc[d.start_type != "", "start_type"].value_counts().to_dict(),
        end_types=d.loc[d.end_type != "", "end_type"].value_counts().to_dict(),
        divisions_confirmed=int(d["division_confirmed"].sum()),
        divisions_linker=info["divisions"],
        in_cluster_fraction=float(d["in_cluster"].mean()),
    )
    with open(out / f"states_{args.name}.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    log.info("%s", json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
