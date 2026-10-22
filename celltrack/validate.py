"""Phase 3: validation without manual ground truth.

Components (each writes JSON + figures under outputs/analysis/validation):
    skip        frame-skip stability: track every 2nd frame vs every frame
    reverse     time-reversal consistency: forward vs backward links of the same linker (3 windows)
    reverse_full whole-movie reverse consistency, if a `<name>_rev` run exists (flags every edge)
    bookkeeping population bookkeeping and unexplained appearances/disappearances
    trackfree   tracking-free step-length estimates (pair-displacement excess, optical flow)
    minigt      consensus mini ground truth on crops, scored with traccuracy
    all         everything above

Usage:
    python -m celltrack.validate all --name trackastra_ilp [--start 0 --stop 576]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from .common import Units, ensure_dir, load_config, setup_logging
from .graphio import load_tracks, node_id, save_tracks, split_node_id, tables_to_graph

log = setup_logging("validate")


def _edges_df(tracks: pd.DataFrame) -> pd.DataFrame:
    e = tracks[tracks.parent_node != -1][["parent_node", "node_id", "edge_conf_in", "frame"]].copy()
    e = e.rename(columns={"parent_node": "u", "node_id": "v", "edge_conf_in": "weight"})
    return e


def _density_bins(det: pd.DataFrame):
    """Three local-density regimes from nearest-neighbour distance (px)."""
    bins = [0, 35, 60, np.inf]
    labels = ["dense (<35px)", "medium (35-60px)", "sparse (>60px)"]
    return pd.cut(det["nn_dist"], bins=bins, labels=labels)


# --------------------------------------------------------------------------- re-linking helper
def _relink(cfg, linker: str, imgs, segs, det_sub, frames, out: Path, tag: str, reverse=False, crop=None):
    """Link a (possibly subsampled / reversed) substack with the given linker and return
    a graph with canonical node ids and edges oriented forward in real time."""
    from .track import run_hoct, run_trackastra

    if linker == "hoct":
        order = np.argsort(frames)[::-1] if reverse else np.arange(len(frames))
        H = run_hoct(np.asarray(imgs)[order], np.asarray(segs)[order], det_sub, int(frames.min()), cfg,
                     out / f"relink_{tag}_work", "cpu", crop=crop, frames=np.asarray(frames)[order])
        return H, None
    mode = "greedy_nodiv" if reverse else "greedy"
    return run_trackastra(np.asarray(imgs), np.asarray(segs), det_sub, frames, cfg, mode=mode, reverse=reverse)


def _linker_of(name: str, cfg) -> str:
    return "hoct" if name.startswith("hoct") else "trackastra"


# --------------------------------------------------------------------------- frame-skip stability
def frame_skip(cfg, tracks: pd.DataFrame, det: pd.DataFrame, start: int, stop: int, out: Path, crop=None, linker="trackastra"):
    """Track frames start, start+2, ... with the same linker and compare with 2-step
    compositions of the full-rate tracks. Agreement = fraction of 2-step links reproduced."""
    from .track import load_inputs

    imgs, segs, det_sub, start, stop = load_inputs(cfg, start, stop, crop)
    frames = np.arange(start, stop)[::2]
    H2, _ = _relink(cfg, linker, np.asarray(imgs)[::2], np.asarray(segs)[::2], det_sub, frames, out, f"skip{start}", crop=crop)
    skip_edges = {(u, v) for u, v in H2.edges}

    G = tables_to_graph(tracks)
    in_region = set(node_id(f, l) for f, l in zip(det_sub.frame, det_sub.label))
    two_step = []
    for u in G.nodes:
        tu = G.nodes[u]["time"]
        if tu < start or tu >= stop or (tu - start) % 2 or u not in in_region:
            continue
        for m in G.successors(u):
            for v in G.successors(m):
                if v in in_region:
                    two_step.append((u, v, G.out_degree(u) == 2 or G.out_degree(m) == 2))
    if not two_step:
        return {}
    ts = pd.DataFrame(two_step, columns=["u", "v", "via_division"])
    ts["reproduced"] = [(u, v) in skip_edges for u, v in zip(ts.u, ts.v)]
    # step length and density at u
    d_idx = det.set_index(["frame", "label"])
    fu = ts.u // 1_000_000
    lu = ts.u % 1_000_000
    fv = ts.v // 1_000_000
    lv = ts.v % 1_000_000
    pu = d_idx.loc[list(zip(fu, lu)), ["x", "y", "nn_dist"]].to_numpy()
    pv = d_idx.loc[list(zip(fv, lv)), ["x", "y"]].to_numpy()
    ts["step2"] = np.hypot(*(pu[:, :2] - pv).T)
    ts["nn_dist"] = pu[:, 2]
    ts["density"] = pd.cut(ts["nn_dist"], [0, 35, 60, np.inf], labels=["dense", "medium", "sparse"])
    ts["step_bin"] = pd.cut(ts["step2"], [0, 20, 40, 60, 100, np.inf])
    res = dict(
        n_two_step_links=int(len(ts)),
        agreement_overall=float(ts.reproduced.mean()),
        agreement_excluding_divisions=float(ts[~ts.via_division].reproduced.mean()),
        agreement_by_density={str(k): float(v) for k, v in ts.groupby("density", observed=True).reproduced.mean().items()},
        agreement_by_step={str(k): float(v) for k, v in ts.groupby("step_bin", observed=True).reproduced.mean().items()},
        n_by_density={str(k): int(v) for k, v in ts.groupby("density", observed=True).size().items()},
    )
    ts.assign(density=ts.density.astype(str), step_bin=ts.step_bin.astype(str)).to_parquet(
        out / f"frame_skip_links_f{start}.parquet", index=False)
    _bar(res["agreement_by_density"], f"frame-skip agreement by local density (frames {start}-{stop})", out / f"frame_skip_density_f{start}.png")
    _bar(res["agreement_by_step"], f"frame-skip agreement by 2-frame step (px, frames {start}-{stop})", out / f"frame_skip_step_f{start}.png")
    return res


# --------------------------------------------------------------------------- time reversal
def time_reversal(cfg, tracks: pd.DataFrame, det: pd.DataFrame, start: int, stop: int, out: Path, crop=None, linker="trackastra"):
    """Link the time-reversed substack with the same linker; a forward edge is
    'bidirectional' when the backward run reproduces it."""
    from .track import load_inputs

    imgs, segs, det_sub, start, stop = load_inputs(cfg, start, stop, crop)
    frames = np.arange(start, stop)
    Hr, cand_r = _relink(cfg, linker, imgs, segs, det_sub, frames, out, f"rev{start}", reverse=True, crop=crop)
    rev_edges = {(u, v) for u, v in Hr.edges}
    e = _edges_df(tracks)
    e = e[(e.frame > start) & (e.frame < stop)]
    # restrict to edges whose target lies in the (cropped) re-linked region
    in_region = set(node_id(f, l) for f, l in zip(det_sub.frame, det_sub.label))
    e = e[e.v.isin(in_region) & e.u.isin(in_region)].copy()
    e["bidir"] = [(u, v) in rev_edges for u, v in zip(e.u, e.v)]
    # reverse-model probability for the forward edge (as a second confidence; hoct: 0/1)
    if cand_r is not None:
        cr = cand_r.set_index(["u", "v"])["weight"]
        e["weight_rev"] = [cr.get((u, v), 0.0) for u, v in zip(e.u, e.v)]
    else:
        e["weight_rev"] = e["bidir"].astype(float)
    d_idx = det.set_index(["frame", "label"])
    fu, lu = e.u // 1_000_000, e.u % 1_000_000
    e["nn_dist"] = d_idx.loc[list(zip(fu, lu)), "nn_dist"].to_numpy()
    e["density"] = pd.cut(e["nn_dist"], [0, 35, 60, np.inf], labels=["dense", "medium", "sparse"])
    e["conf_bin"] = pd.cut(e["weight"], [0, 0.5, 0.7, 0.9, 1.0001])
    res = dict(
        n_edges=int(len(e)),
        bidirectional_fraction=float(e.bidir.mean()),
        bidir_by_density={str(k): float(v) for k, v in e.groupby("density", observed=True).bidir.mean().items()},
        bidir_by_forward_conf={str(k): float(v) for k, v in e.groupby("conf_bin", observed=True).bidir.mean().items()},
    )
    e[["u", "v", "weight", "weight_rev", "bidir"]].to_parquet(out / f"edge_confidence_f{start}.parquet", index=False)
    _bar(res["bidir_by_density"], f"time-reversal agreement by local density (frames {start}-{stop})", out / f"reversal_density_f{start}.png")
    _bar(res["bidir_by_forward_conf"], f"time-reversal agreement by forward edge probability (frames {start}-{stop})", out / f"reversal_conf_f{start}.png")
    return res


def reverse_full(cfg, name: str, tracks: pd.DataFrame, det: pd.DataFrame, units: Units, out: Path) -> dict | None:
    """Whole-movie reverse consistency: compare the forward tracks with a complete
    time-reversed run of the same linker (`<name>_rev`), if it exists. Produces the
    bidir flag for *every* edge (used by the high-confidence matrices) and the
    agreement by local density and by time."""
    rev_p = Path(cfg["paths"]["tracks_dir"]) / f"{name}_rev_tracks.parquet"
    if not rev_p.exists():
        return None
    rev = pd.read_parquet(rev_p)
    rev_edges = set(zip(rev.parent_node[rev.parent_node != -1], rev.node_id[rev.parent_node != -1]))
    e = _edges_df(tracks)
    e["bidir"] = [(u, v) in rev_edges for u, v in zip(e.u, e.v)]
    e["weight_rev"] = e["bidir"].astype(float)
    d_idx = det.set_index(["frame", "label"])
    fu, lu = e.u // 1_000_000, e.u % 1_000_000
    e["nn_dist"] = d_idx.loc[list(zip(fu, lu)), "nn_dist"].to_numpy()
    e["density"] = pd.cut(e["nn_dist"], [0, 35, 60, np.inf], labels=["dense", "medium", "sparse"])
    e["hour_bin"] = (units.frames_to_hours(e["frame"]) // 6 * 6).astype(int)
    e[["u", "v", "weight", "weight_rev", "bidir"]].to_parquet(out / f"edge_confidence_{name}.parquet", index=False)
    res = dict(
        n_edges=int(len(e)),
        n_reverse_edges=int(len(rev_edges)),
        bidirectional_fraction=float(e.bidir.mean()),
        bidir_by_density={str(k): float(v) for k, v in e.groupby("density", observed=True).bidir.mean().items()},
        bidir_by_6h_window={f"{k}-{k+6}h": float(v) for k, v in e.groupby("hour_bin").bidir.mean().items()},
    )
    _bar(res["bidir_by_density"], "whole-movie reverse agreement by local density", out / "reversal_full_density.png")
    _bar(res["bidir_by_6h_window"], "whole-movie reverse agreement by time", out / "reversal_full_time.png")
    return res


def merge_edge_confidence(out: Path, name: str):
    """Concatenate the per-window reverse-consistency tables into one file that
    states.py can join on (edge_confidence_<name>.parquet)."""
    parts = sorted(out.glob("edge_confidence_f*.parquet"))
    if parts:
        pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True).drop_duplicates(["u", "v"]).to_parquet(
            out / f"edge_confidence_{name}.parquet", index=False)


# --------------------------------------------------------------------------- bookkeeping
def bookkeeping(cfg, tracks: pd.DataFrame, det: pd.DataFrame, units: Units, out: Path):
    margin = cfg["segmentation"]["border_margin_px"]
    d = tracks.merge(det[["frame", "label", "dist_to_border"]], on=["frame", "label"], how="left")
    T0, T1 = d.frame.min(), d.frame.max()
    per = d.groupby("frame").size().rename("n").to_frame()
    starts = d[(d.parent_node == -1) & (d.frame > T0)]
    ends = d[(d.n_children == 0) & (d.frame < T1)]
    divs = d[d.n_children == 2]
    per["appear_interior"] = starts[starts.dist_to_border >= margin].groupby("frame").size()
    per["appear_border"] = starts[starts.dist_to_border < margin].groupby("frame").size()
    per["end_interior"] = ends[ends.dist_to_border >= margin].groupby("frame").size()
    per["end_border"] = ends[ends.dist_to_border < margin].groupby("frame").size()
    per["divisions"] = divs.groupby("frame").size()
    per = per.fillna(0).astype(int)
    per["dn"] = per["n"].diff().fillna(0).astype(int)
    per["hours"] = units.frames_to_hours(per.index)
    cell_frames = per["n"].sum()
    # frames whose detection count deviates strongly from the local trend (defocus / illumination glitches)
    roll = per["n"].rolling(9, center=True, min_periods=3).median()
    dev = (per["n"] - roll) / roll
    outlier_frames = [int(f) for f in per.index[dev.abs() > 0.08]]
    res = dict(
        outlier_frames_count_dev_gt_8pct=outlier_frames,
        frames=[int(T0), int(T1)],
        cell_frames=int(cell_frames),
        n_first=int(per["n"].iloc[0]),
        n_last=int(per["n"].iloc[-1]),
        divisions_total=int(per["divisions"].sum()),
        interior_appearances_total=int(per["appear_interior"].sum()),
        interior_ends_total=int(per["end_interior"].sum()),
        border_appearances_total=int(per["appear_border"].sum()),
        border_ends_total=int(per["end_border"].sum()),
        # fragmentation: unexplained interior appearances per cell per frame (should be ~0 for perfect tracking)
        interior_appearance_rate_per_cell_frame=float(per["appear_interior"].sum() / cell_frames),
        interior_end_rate_per_cell_frame=float(per["end_interior"].sum() / cell_frames),
        count_noise_std_dn=float(per["dn"].std()),
        mean_track_length_frames=float(tracks.groupby("track_id").size().mean()),
        median_track_length_frames=float(tracks.groupby("track_id").size().median()),
    )
    per.to_csv(out / "bookkeeping_per_frame.csv")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    ax[0].plot(per.hours, per.n)
    ax[0].set_xlabel("time (h)"); ax[0].set_ylabel("cells detected")
    w = 12
    for c in ["divisions", "appear_interior", "end_interior", "appear_border", "end_border"]:
        ax[1].plot(per.hours, per[c].rolling(w, center=True).mean(), label=c)
    ax[1].set_xlabel("time (h)"); ax[1].set_ylabel(f"events / frame ({w}-frame mean)"); ax[1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(out / "bookkeeping.png", dpi=80); plt.close(fig)
    return res


# --------------------------------------------------------------------------- tracking-free estimates
def pair_displacement_excess(det: pd.DataFrame, frames, R=120.0, shift=700.0, bins=None):
    """Histogram of all t -> t+1 pair distances within R minus the same histogram
    for spatially shifted (uncorrelated) t+1 positions. The excess is the true
    step-length distribution, independent of any linking."""
    from scipy.spatial import cKDTree

    bins = bins if bins is not None else np.arange(0, R + 2, 2.0)
    h_true = np.zeros(len(bins) - 1)
    h_bg = np.zeros(len(bins) - 1)
    by = {f: g[["x", "y"]].to_numpy() for f, g in det[det.frame.isin(set(frames) | {f + 1 for f in frames})].groupby("frame")}
    W = H = 3072
    for f in frames:
        if f not in by or f + 1 not in by:
            continue
        A, B = by[f], by[f + 1]
        tA = cKDTree(A)
        dd = tA.sparse_distance_matrix(cKDTree(B), R, output_type="ndarray")["v"]
        h_true += np.histogram(dd, bins)[0]
        # background: four large shifts with wrap-around (keeps density gradient approximately)
        for dx, dy in [(shift, 0), (-shift, 0), (0, shift), (0, -shift)]:
            Bs = np.column_stack([(B[:, 0] + dx) % W, (B[:, 1] + dy) % H])
            dd = tA.sparse_distance_matrix(cKDTree(Bs), R, output_type="ndarray")["v"]
            h_bg += np.histogram(dd, bins)[0] / 4.0
    excess = np.clip(h_true - h_bg, 0, None)
    centers = 0.5 * (bins[1:] + bins[:-1])
    # the excess beyond ~60 px is dominated by background noise; restrict the summary
    p = excess / max(excess.sum(), 1)
    cdf = np.cumsum(p)
    med = float(np.interp(0.5, cdf, centers))
    mean = float((p * centers).sum())
    return dict(centers=centers, excess=excess, h_true=h_true, h_bg=h_bg, median_step=med, mean_step=mean)


def optical_flow_speed(cfg, frames, crop=(1024, 1024, 1024)):
    """Mean optical-flow magnitude inside cell masks (iLK, skimage) on a crop."""
    import zarr
    from skimage.registration import optical_flow_ilk

    flat = zarr.open_array(cfg["paths"]["flat_zarr"], mode="r")
    masks = zarr.open_array(cfg["paths"]["masks_zarr"], mode="r")
    y0, x0, s = crop
    speeds = []
    for f in frames:
        a = flat[f, y0:y0 + s, x0:x0 + s].astype(np.float32) / 255
        b = flat[f + 1, y0:y0 + s, x0:x0 + s].astype(np.float32) / 255
        m = masks[f, y0:y0 + s, x0:x0 + s] > 0
        v, u = optical_flow_ilk(a, b, radius=15, num_warp=5, prefilter=True)
        mag = np.hypot(u, v)[m]
        if mag.size:
            speeds.append(float(np.median(mag)))
    return dict(frames=list(map(int, frames)), median_flow_speed_px=float(np.median(speeds)) if speeds else np.nan,
                per_frame=speeds)


def track_free(cfg, tracks: pd.DataFrame, det: pd.DataFrame, out: Path):
    T0, T1 = int(det.frame.min()), int(det.frame.max())
    sample = np.unique(np.linspace(T0, T1 - 1, min(40, T1 - T0)).astype(int))
    pde = pair_displacement_excess(det, sample)
    e = _edges_df(tracks)
    d_idx = det.set_index(["frame", "label"])[["x", "y"]]
    pu = d_idx.loc[list(zip(e.u // 1_000_000, e.u % 1_000_000))].to_numpy()
    pv = d_idx.loc[list(zip(e.v // 1_000_000, e.v % 1_000_000))].to_numpy()
    step = np.hypot(*(pu - pv).T)
    step = step[np.isin(e.frame.to_numpy() - 1, sample)]
    hi = e.weight.to_numpy()[np.isin(e.frame.to_numpy() - 1, sample)] >= cfg["tracking"]["high_confidence_prob"]
    of = optical_flow_speed(cfg, sample[:: max(1, len(sample) // 8)][:8])
    res = dict(
        pair_excess_median_step_px=pde["median_step"],
        pair_excess_mean_step_px=pde["mean_step"],
        track_median_step_px=float(np.median(step)),
        track_mean_step_px=float(np.mean(step)),
        track_highconf_median_step_px=float(np.median(step[hi])) if hi.any() else np.nan,
        optical_flow_median_speed_px=of["median_flow_speed_px"],
        note="optical flow underestimates for displacements > ~1 cell radius; pair-excess is the reference",
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4))
    c = pde["centers"]
    ax.plot(c, pde["excess"] / max(pde["excess"].sum(), 1) / np.diff(c)[0], label="pair-displacement excess (tracking-free)")
    ax.hist(step, bins=np.arange(0, 122, 2), density=True, alpha=0.4, label="linked steps (all)")
    if hi.any():
        ax.hist(step[hi], bins=np.arange(0, 122, 2), density=True, histtype="step", label="linked steps (high conf)")
    ax.set_xlabel("step length (px / frame)"); ax.set_ylabel("density"); ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(out / "step_distribution_vs_trackfree.png", dpi=80); plt.close(fig)
    return res


# --------------------------------------------------------------------------- consensus mini ground truth
def mini_gt(cfg, names: list[str], det: pd.DataFrame, out: Path, crops, n_frames=20):
    """Build a consensus 'ground truth' on small crops from several linkers plus the
    time-reversed Trackastra run, then score every linker against it with traccuracy.

    The consensus keeps an edge when the majority of runs agree. It is biased towards
    the linkers used, so treat the resulting numbers as a *ranking* and as an
    optimistic error-rate estimate, not as an absolute accuracy."""
    from .track import load_inputs, run_hoct, run_laptrack, run_trackastra

    out_gt = ensure_dir(out / "minigt")
    results = {}
    for ci, (f0, y0, x0, s) in enumerate(crops):
        f1 = f0 + n_frames
        imgs, segs, det_sub, f0, f1 = load_inputs(cfg, f0, f1, (y0, x0, s))
        frames = np.arange(f0, f1)
        runs = {}
        runs["ta_greedy"], _ = run_trackastra(imgs, segs, det_sub, frames, cfg, mode="greedy")
        runs["ta_ilp"], _ = run_trackastra(imgs, segs, det_sub, frames, cfg, mode="ilp_nodiv", chunk=n_frames)
        runs["ta_rev"], _ = run_trackastra(imgs, segs, det_sub, frames, cfg, mode="greedy_nodiv", reverse=True)
        runs["laptrack"] = run_laptrack(det_sub, cfg)
        try:
            runs["hoct"] = run_hoct(imgs, segs, det_sub, f0, cfg, out_gt / f"crop{ci}_hoct_work", "cpu", crop=(y0, x0, s))
        except Exception as e:  # HOCT is optional for the consensus
            log.warning("HOCT failed on crop %d: %s", ci, str(e)[:200])
        for nm in names:  # full-movie results restricted to the crop
            try:
                tr, _ = load_tracks(cfg["paths"]["tracks_dir"], nm)
                G = tables_to_graph(tr)
                keep = set(node_id(f, l) for f, l in zip(det_sub.frame, det_sub.label))
                sub = G.subgraph([n for n in keep if n in G]).copy()
                if sub.number_of_nodes() > 0.5 * len(det_sub):  # run covers this crop/window
                    runs[f"full:{nm}"] = sub
            except FileNotFoundError:
                pass
        votes: dict[tuple, int] = {}
        voters = [v for v in ["ta_greedy", "ta_ilp", "ta_rev", "laptrack", "hoct"] if v in runs]
        for nm in voters:
            for u, v in runs[nm].edges:
                votes[(u, v)] = votes.get((u, v), 0) + 1
        need = len(voters) // 2 + 1  # strict majority
        gt = nx.DiGraph()
        gt.add_nodes_from(runs["ta_greedy"].nodes(data=True))
        for (u, v), c in sorted(votes.items(), key=lambda kv: -kv[1]):
            if c >= need and gt.in_degree(v) == 0 and gt.out_degree(u) < 2:
                gt.add_edge(u, v)
        save_tracks(gt, out_gt, f"crop{ci}_consensus", write_geff=False)
        for nm, G in runs.items():
            save_tracks(G, out_gt, f"crop{ci}_{nm.replace(':', '_')}", write_geff=False)
        scores = {nm: score_against(gt, G) for nm, G in runs.items()}
        results[f"crop{ci}_f{f0}_y{y0}_x{x0}"] = dict(
            n_nodes=gt.number_of_nodes(), n_consensus_edges=gt.number_of_edges(),
            consensus_fraction=gt.number_of_edges() / max(1, runs["ta_greedy"].number_of_edges()),
            scores=scores,
        )
        _disagreement_montage(cfg, runs.get(f"full:{names[0]}", runs["ta_ilp"]), gt, det_sub, (y0, x0, s), out_gt / f"crop{ci}_disagreements.png")
    return results


def score_against(gt: nx.DiGraph, pred: nx.DiGraph) -> dict:
    """Edge precision/recall/F1 and division F1 vs a reference graph on identical nodes,
    plus traccuracy CTC-style metrics when available."""
    ge, pe = set(gt.edges), set(pred.edges)
    tp = len(ge & pe)
    prec = tp / max(1, len(pe))
    rec = tp / max(1, len(ge))
    f1 = 2 * prec * rec / max(1e-9, prec + rec)
    gd = {n for n in gt.nodes if gt.out_degree(n) == 2}
    pdv = {n for n in pred.nodes if pred.out_degree(n) == 2}
    dtp = len(gd & pdv)
    dprec = dtp / max(1, len(pdv))
    drec = dtp / max(1, len(gd))
    res = dict(edge_precision=prec, edge_recall=rec, edge_f1=f1, n_pred_edges=len(pe),
               division_precision=dprec, division_recall=drec, n_gt_div=len(gd), n_pred_div=len(pdv))
    try:
        res.update(_traccuracy(gt, pred))
    except Exception as e:  # keep the simple metrics if traccuracy's API differs
        res["traccuracy_error"] = str(e)[:200]
    return res


def _traccuracy(gt: nx.DiGraph, pred: nx.DiGraph) -> dict:
    from traccuracy import TrackingGraph
    from traccuracy.matchers import PointMatcher
    from traccuracy.metrics import CTCMetrics, DivisionMetrics

    def tg(G):
        H = nx.DiGraph()
        for n, d in G.nodes(data=True):
            H.add_node(n, t=int(d["time"]), y=float(d["coords"][0]), x=float(d["coords"][1]))
        H.add_edges_from(G.edges)
        return TrackingGraph(H, frame_key="t", location_keys=("y", "x"))

    matched = PointMatcher(threshold=1.0).compute_mapping(tg(gt), tg(pred))
    out = {}
    r = CTCMetrics().compute(matched)
    out.update({f"ctc_{k}": float(v) for k, v in r.results.items() if isinstance(v, (int, float))})
    r = DivisionMetrics(max_frame_buffer=1).compute(matched)
    for k, v in r.results.items():
        if isinstance(v, dict):
            for kk, vv in v.items():
                if isinstance(vv, (int, float)):
                    out[f"div_{k}_{kk}"] = float(vv)
    return out


def _disagreement_montage(cfg, pred: nx.DiGraph, gt: nx.DiGraph, det: pd.DataFrame, crop, path: Path, n=8):
    """Side-by-side frames for links where pred and consensus disagree, for visual review."""
    import matplotlib
    import zarr

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    flat = zarr.open_array(cfg["paths"]["flat_zarr"], mode="r")
    y0, x0, s = crop
    diff = [(u, v) for u, v in pred.edges if not gt.has_edge(u, v)][:n]
    if not diff:
        return
    fig, axes = plt.subplots(len(diff), 2, figsize=(8, 4 * len(diff)))
    axes = np.atleast_2d(axes)
    for i, (u, v) in enumerate(diff):
        for j, nid in enumerate((u, v)):
            f, l = split_node_id(nid)
            y, x = pred.nodes[nid]["coords"]
            r = 60
            ya, yb, xa, xb = int(max(y - r, 0)), int(y + r), int(max(x - r, 0)), int(x + r)
            axes[i, j].imshow(flat[f, ya:yb, xa:xb], cmap="gray", vmin=0, vmax=255)
            axes[i, j].set_title(f"t={f} label={l}" + ("" if j == 0 else "  (pred link; consensus: " +
                                 (f"->{split_node_id(next(gt.successors(u)))[1]}" if gt.out_degree(u) else "end") + ")"), fontsize=8)
            axes[i, j].plot(x - xa, y - ya, "r+")
            axes[i, j].axis("off")
    plt.tight_layout(); plt.savefig(path, dpi=70); plt.close(fig)


# --------------------------------------------------------------------------- helpers
def _bar(d: dict, title: str, path: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.bar(range(len(d)), list(d.values()))
    ax.set_xticks(range(len(d)), list(d.keys()), rotation=20, fontsize=8)
    ax.set_ylim(0, 1); ax.set_title(title, fontsize=10)
    plt.tight_layout(); plt.savefig(path, dpi=80); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["skip", "reverse", "reverse_full", "bookkeeping", "trackfree", "minigt", "all"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--name", default="hoct", help="tracks name under outputs/tracks")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=10**9)
    ap.add_argument("--crop", type=int, nargs=3, default=None, metavar=("Y0", "X0", "SIZE"),
                    help="crop for the skip/reverse re-tracking runs (keeps them fast)")
    ap.add_argument("--compare", nargs="*", default=[], help="other track names to score in minigt")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    units = Units.from_config(cfg)
    out = ensure_dir(Path(cfg["paths"]["analysis_dir"]) / "validation")
    det = pd.read_parquet(cfg["paths"]["detections"])
    tracks, _ = load_tracks(cfg["paths"]["tracks_dir"], args.name)
    T0, T1 = int(tracks.frame.min()), int(tracks.frame.max()) + 1
    start, stop = max(args.start, T0), min(args.stop, T1)
    results = {"tracks": args.name, "frames": [start, stop]}
    linker = _linker_of(args.name, cfg)
    res_path = out / f"validation_{args.name}.json"
    if res_path.exists():
        results.update(json.load(open(res_path)))

    if args.what in ("bookkeeping", "all"):
        results["bookkeeping"] = bookkeeping(cfg, tracks, det, units, out)
        log.info("bookkeeping: %s", results["bookkeeping"])
    if args.what in ("trackfree", "all"):
        results["track_free"] = track_free(cfg, tracks, det, out)
        log.info("track-free: %s", results["track_free"])
    if args.what in ("skip", "all"):
        # re-tracking is expensive; use a spatial crop and a 60-frame window at three times
        wins = _windows(start, stop, n=3, length=60)
        crop = args.crop or (768, 768, 1536)
        parts = [frame_skip(cfg, tracks, det, a, b, out, crop, linker) for a, b in wins]
        results["frame_skip"] = {"windows": [list(map(int, w)) for w in wins], "per_window": parts,
                                 "agreement_overall": float(np.mean([p["agreement_overall"] for p in parts if p]))}
        log.info("frame skip: %s", results["frame_skip"]["agreement_overall"])
    if args.what in ("reverse", "all"):
        wins = _windows(start, stop, n=3, length=60)
        crop = args.crop or (768, 768, 1536)
        parts = [time_reversal(cfg, tracks, det, a, b, out, crop, linker) for a, b in wins]
        results["time_reversal"] = {"windows": [list(map(int, w)) for w in wins], "per_window": parts,
                                    "bidirectional_fraction": float(np.mean([p["bidirectional_fraction"] for p in parts]))}
        merge_edge_confidence(out, args.name)
    if args.what in ("reverse", "reverse_full", "all"):
        rf = reverse_full(cfg, args.name, tracks, det, units, out)
        if rf is not None:  # whole-movie flags supersede the window-based ones
            results["time_reversal_full"] = rf
            log.info("whole-movie time reversal: %s", rf)
        log.info("time reversal: %s", results["time_reversal"]["bidirectional_fraction"])
    if args.what in ("minigt", "all"):
        wins = _windows(start, stop, n=3, length=20)
        crops = [c for c in (_pick_crop(det, a, b, cfg["acquisition"]["shape_yx"]) for a, b in wins) if c is not None]
        results["mini_gt"] = mini_gt(cfg, [args.name, *args.compare], det, out, crops)
        log.info("mini-GT: %s", json.dumps(results["mini_gt"], indent=1, default=float)[:2000])

    with open(res_path, "w") as fh:
        json.dump(results, fh, indent=2, default=float)
    log.info("Wrote %s", res_path)


def _pick_crop(det: pd.DataFrame, a: int, b: int, shape_yx, size: int = 512, min_per_frame: float = 3.0):
    """Crop (f0, y0, x0, size) for the mini-GT window [a, b). The fixed centre crop is kept
    when it holds enough cells (so dense wells stay comparable across runs); for sparse
    wells the densest `size`² region is used instead, and a window without any cells is
    skipped (returns None)."""
    d = det[(det.frame >= a) & (det.frame < b)]
    need = min_per_frame * (b - a)
    y0 = x0 = 1280
    inside = d[(d.y >= y0) & (d.y < y0 + size) & (d.x >= x0) & (d.x < x0 + size)]
    if len(inside) >= need:
        return (a, y0, x0, size)
    if len(d) == 0:
        log.warning("mini-GT: no detections in frames %d-%d, skipping crop", a, b)
        return None
    step = size // 2  # coarse search for the densest crop: count cells per half-crop bin
    H, W = shape_yx
    iy = (d.y.values // step).astype(int)
    ix = (d.x.values // step).astype(int)
    counts = np.zeros((iy.max() + 2, ix.max() + 2), dtype=int)
    np.add.at(counts, (iy, ix), 1)
    n_y, n_x = max(1, (H - size) // step + 1), max(1, (W - size) // step + 1)  # candidate starts
    counts = np.pad(counts, ((0, max(0, n_y + 1 - counts.shape[0])), (0, max(0, n_x + 1 - counts.shape[1]))))
    full = counts[:n_y, :n_x] + counts[1:n_y + 1, :n_x] + counts[:n_y, 1:n_x + 1] + counts[1:n_y + 1, 1:n_x + 1]
    i, j = np.unravel_index(int(np.argmax(full)), full.shape)
    y0, x0 = int(i * step), int(j * step)
    log.info("mini-GT: centre crop too sparse in frames %d-%d (%d cells), using densest crop at y=%d x=%d (%d cells)",
             a, b, len(inside), y0, x0, int(full[i, j]))
    return (a, y0, x0, size)


def _windows(start, stop, n, length):
    length = min(length, stop - start)
    if stop - start <= length:
        return [(start, stop)]
    starts = np.linspace(start, stop - length, n).astype(int)
    return [(int(s), int(s + length)) for s in starts]


if __name__ == "__main__":
    main()
