"""Graph refinement: repair segmentation merge/split artifacts and propose divisions.

Why: touching cells are sometimes segmented as one object for a few frames. In the
lineage graph this is a *merge* (a track ends next to an enlarged object) followed by a
*split* (a new track starts next to it, or the linker reports a 'division'). Real
divisions are rare (~0.1 % of cell-frames), so such artifacts dominate the division
count unless repaired. Linkers are therefore run 1-to-1 where possible and divisions are
decided here, with explicit morphological and temporal criteria.

Terminology: an *orphan* is a node with no parent that is not in the first frame.
For a linker that reports divisions, the weaker child of every 2-child parent is first
detached and treated as an orphan, so all linkers go through the same logic.

resolve_merge_splits(orphans):
    orphan D at frame t; S = enlarged object at t-1 within reach of D (the merged
    object); walk S's tracklet back to its start S0; if a track X ended within
    `gap_frames` before S0, near S0, and S0 is large enough to contain X plus another
    cell, link X -> D (gap-closed through the merged object).
propose_divisions(orphans):
    orphan D at frame t; mother M at t-1 within reach with exactly one child C at t;
    accept M -> (C, D) as a division when M's tracklet is old enough, both daughters
    are small relative to M and to the frame median, similar to each other and persist.
mark_merges:
    interior track ends adjacent to an enlarged object in the next frame are labelled
    'merged' (identity lost by under-segmentation) instead of 'lost'.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .graphio import graph_to_tables, tables_to_graph


def _reach(area: float, floor: float = 40.0) -> float:
    return max(floor, 1.6 * np.sqrt(max(area, 1.0)))


class _Ctx:
    """Fast lookups shared by the refinement steps."""

    def __init__(self, G, det: pd.DataFrame):
        self.G = G
        d = det.set_index(["frame", "label"])
        self.x = d["x"].to_dict()
        self.y = d["y"].to_dict()
        self.area = d["area"].to_dict()
        self.touch = d["n_touching"].to_dict() if "n_touching" in d else {}
        self.circ = d["circularity"].to_dict() if "circularity" in d else {}
        self.med_area = det.groupby("frame")["area"].median().to_dict()
        self.T0, self.T1 = int(det.frame.min()), int(det.frame.max())
        self.by_frame = {}
        for f, g in det.groupby("frame"):
            ids = (g["frame"].to_numpy() * 1_000_000 + g["label"].to_numpy()).astype(np.int64)
            self.by_frame[f] = (cKDTree(g[["x", "y"]].to_numpy()), ids)

    def info(self, n):
        k = (n // 1_000_000, n % 1_000_000)
        return k[0], self.x[k], self.y[k], self.area[k]

    def area_of(self, n):
        return self.area[(n // 1_000_000, n % 1_000_000)]

    def near(self, f, x, y, r):
        if f not in self.by_frame:
            return []
        tree, ids = self.by_frame[f]
        return [int(ids[i]) for i in tree.query_ball_point((x, y), r)]

    def tracklet_len_back(self, n, max_len=10**6):
        """Number of nodes in the linear chain ending at n (inclusive), capped."""
        G, k = self.G, 1
        while G.in_degree(n) == 1 and k <= max_len:
            p = next(G.predecessors(n))
            if G.out_degree(p) != 1:
                break
            n, k = p, k + 1
        return k, n

    def tracklet_len_fwd(self, n, max_len=10**6):
        G, k = self.G, 1
        while G.out_degree(n) == 1 and k <= max_len:
            c = next(G.successors(n))
            if G.in_degree(c) != 1:
                break
            n, k = c, k + 1
        return k, n


def detach_weak_children(G) -> int:
    """Turn linker divisions into orphans (weaker child edge removed)."""
    n = 0
    for p in [n for n in G.nodes if G.out_degree(n) == 2]:
        kids = sorted(G.successors(p), key=lambda c: G.edges[p, c].get("weight", 1.0))
        G.remove_edge(p, kids[0])
        n += 1
    return n


def resolve_merge_splits(G, ctx: _Ctx, cfg) -> int:
    rc = cfg["states"].get("merge_split", {})
    gap_frames = rc.get("gap_frames", 3)
    area_ratio = rc.get("min_merged_area_ratio", 1.25)
    max_merged = rc.get("max_merged_frames", 12)
    repaired = 0
    orphans = [n for n in G.nodes if G.in_degree(n) == 0 and G.nodes[n]["time"] > ctx.T0]
    for D in sorted(orphans, key=lambda n: G.nodes[n]["time"]):
        fD, xD, yD, aD = ctx.info(D)
        # merged object S at t-1: enlarged, within reach, and not itself ending here
        best_S, best_d = None, np.inf
        for S in ctx.near(fD - 1, xD, yD, _reach(aD) + 20):
            _, xS, yS, aS = ctx.info(S)
            if aS < area_ratio * aD or G.out_degree(S) == 0:
                continue
            dd = np.hypot(xS - xD, yS - yD)
            if dd < best_d:
                best_S, best_d = S, dd
        if best_S is None:
            continue
        # walk back along the merged object's tracklet to the merge point: either the
        # tracklet start, or the frame where its area jumped by >= area_ratio
        cur, k = best_S, 0
        merge_node = None
        while k <= max_merged:
            prev = next(G.predecessors(cur), None)
            if prev is None:
                merge_node = cur
                break
            if ctx.area_of(cur) >= area_ratio * ctx.area_of(prev):
                merge_node = cur
                break
            cur, k = prev, k + 1
        if merge_node is None:
            continue
        fS0, xS0, yS0, aS0 = ctx.info(merge_node)
        y_prev = next(G.predecessors(merge_node), None)
        best_X, best_c = None, np.inf
        for f in range(fS0 - gap_frames, fS0):
            for X in ctx.near(f, xS0, yS0, _reach(aS0)):
                if G.out_degree(X) != 0 or X == y_prev or X == D:
                    continue
                _, xX, yX, aX = ctx.info(X)
                if aS0 < area_ratio * aX:
                    continue
                c = np.hypot(xX - xS0, yX - yS0) + 5.0 * (fS0 - f - 1)
                if c < best_c:
                    best_X, best_c = X, c
        if best_X is None:
            continue
        G.add_edge(best_X, D, weight=0.5, via_merge=int(fD - ctx.info(best_X)[0]))
        repaired += 1
    return repaired


def propose_divisions(G, ctx: _Ctx, cfg) -> dict:
    dc = cfg["states"].get("division", {})
    min_parent = dc.get("min_parent_frames", 4)
    max_frac = dc.get("max_daughter_frac", 0.8)
    rel_med = dc.get("max_daughter_rel_median", 0.85)
    min_kid = dc.get("min_daughter_frames", 2)
    max_step = dc.get("max_daughter_step_px", 45)
    min_circ = dc.get("min_mother_circularity", 0.75)
    min_rel_mother = dc.get("min_mother_rel_median", 0.9)
    sep_factor = dc.get("max_daughter_separation_factor", 2.5)  # x sqrt(mean daughter area)
    max_rel_mother = dc.get("max_mother_rel_median", 2.2)
    min_rel_daughter = dc.get("min_daughter_rel_median", 0.45)
    area_tol = dc.get("area_conservation_tol", 0.35)
    stats = dict(candidates=0, accepted=0, reasons={})
    orphans = [n for n in G.nodes if G.in_degree(n) == 0 and G.nodes[n]["time"] > ctx.T0]
    for D in sorted(orphans, key=lambda n: G.nodes[n]["time"]):
        fD, xD, yD, aD = ctx.info(D)
        # mothers at t-1 with exactly one child at t
        cands = []
        for M in ctx.near(fD - 1, xD, yD, max_step):
            if G.out_degree(M) != 1:
                continue
            C = next(G.successors(M))
            if G.nodes[C]["time"] != fD:
                continue
            cands.append((np.hypot(ctx.x[(fD - 1, M % 1_000_000)] - xD, ctx.y[(fD - 1, M % 1_000_000)] - yD), M, C))
        if not cands:
            continue
        stats["candidates"] += 1
        cands.sort()
        accepted = False
        for _, M, C in cands[:2]:
            _, xM, yM, aM = ctx.info(M)
            _, xC, yC, aC = ctx.info(C)
            reasons = []
            kM = (fD - 1, M % 1_000_000)
            n_back, M0 = ctx.tracklet_len_back(M)
            if n_back < min_parent:
                reasons.append("short_parent")
            # mitotic mother: isolated (segmentation is reliable), rounded, not small
            prevM = next(G.predecessors(M), None)
            touching = ctx.touch.get(kM, 0) > 0 or (prevM is not None and ctx.touch.get((fD - 2, prevM % 1_000_000), 0) > 0)
            if touching:
                reasons.append("mother_touching")
            if ctx.circ.get(kM, 1.0) < min_circ:
                reasons.append("mother_not_round")
            med_prev = ctx.med_area.get(fD - 1, 0)
            if aM < min_rel_mother * med_prev:
                reasons.append("mother_small")
            if aM > max_rel_mother * med_prev:
                reasons.append("mother_too_large")  # merged cluster, not one cell
            if np.hypot(xC - xD, yC - yD) > sep_factor * np.sqrt(0.5 * (aC + aD)):
                reasons.append("daughters_apart")
            if 0.5 * (aD + aC) < min_rel_daughter * ctx.med_area.get(fD, 0):
                reasons.append("daughters_tiny")  # apoptotic fragmentation / debris
            if abs(aD + aC - aM) > area_tol * aM:
                reasons.append("area_not_conserved")
            if not (aD <= max_frac * aM and aC <= max_frac * aM):
                reasons.append("daughter_too_large")
            if not (0.5 <= aD / max(aC, 1) <= 2.0):
                reasons.append("asymmetric")
            if 0.5 * (aD + aC) > rel_med * ctx.med_area.get(fD, np.inf):
                reasons.append("daughters_normal_size")
            lenD, endD = ctx.tracklet_len_fwd(D)
            lenC, endC = ctx.tracklet_len_fwd(C)
            if not ((lenD >= min_kid or G.nodes[endD]["time"] >= ctx.T1) and (lenC >= min_kid or G.nodes[endC]["time"] >= ctx.T1)):
                reasons.append("short_daughter")
            # a touching mother is the only *soft* criterion: every morphological test passed,
            # but the segmentation of touching cells is less reliable -> accept with lower confidence
            soft_only = reasons == ["mother_touching"]
            if reasons and not soft_only:
                for r in reasons:
                    stats["reasons"][r] = stats["reasons"].get(r, 0) + 1
                continue
            G.add_edge(M, D, weight=0.3 if soft_only else 0.6, division=True,
                       division_quality="in_cluster" if soft_only else "isolated")
            accepted = "in_cluster" if soft_only else "isolated"
            break
        if accepted:
            stats["accepted"] += 1
            stats[f"accepted_{accepted}"] = stats.get(f"accepted_{accepted}", 0) + 1
    return stats


def refine_graph(tracks: pd.DataFrame, det: pd.DataFrame, cfg) -> tuple[pd.DataFrame, dict]:
    """Full refinement: detach linker divisions -> repair merge/splits -> propose divisions."""
    G = tables_to_graph(tracks)
    ctx = _Ctx(G, det)
    stats = {"linker_divisions_detached": detach_weak_children(G)}
    stats["merge_split_repaired"] = resolve_merge_splits(G, ctx, cfg)
    stats["divisions"] = propose_divisions(G, ctx, cfg)
    new_tracks, _ = graph_to_tables(G)
    via = {(u, v): dd["via_merge"] for u, v, dd in G.edges(data=True) if "via_merge" in dd}
    new_tracks["gap_via_merge"] = [via.get((p, n), 0) for p, n in zip(new_tracks.parent_node, new_tracks.node_id)]
    new_tracks["division_confirmed"] = new_tracks.n_children == 2
    # quality of the division a node is the *mother* of ('isolated' | 'in_cluster' | '')
    dq = {u: dd["division_quality"] for u, v, dd in G.edges(data=True) if "division_quality" in dd}
    new_tracks["division_quality"] = [dq.get(n, "") if k == 2 else "" for n, k in zip(new_tracks.node_id, new_tracks.n_children)]
    return new_tracks, stats


def mark_merges(tracks: pd.DataFrame, det: pd.DataFrame, cfg) -> pd.Series:
    """Boolean per node: a track end whose cell was most likely absorbed into an enlarged
    neighbouring object in the next frame (under-segmentation)."""
    rc = cfg["states"].get("merge_split", {})
    area_ratio = rc.get("min_merged_area_ratio", 1.25)
    T1 = int(tracks.frame.max())
    d = det.set_index(["frame", "label"])["area"].to_dict()
    med_area = det.groupby("frame")["area"].median().to_dict()
    by_frame = {f: (cKDTree(g[["x", "y"]].to_numpy()), g["area"].to_numpy()) for f, g in det.groupby("frame")}
    ends = tracks[(tracks.n_children == 0) & (tracks.frame < T1)]
    flags = pd.Series(False, index=tracks.index)
    for i, r in ends.iterrows():
        f = int(r.frame) + 1
        if f not in by_frame:
            continue
        aX = d[(int(r.frame), int(r.label))]
        tree, areas = by_frame[f]
        idx = tree.query_ball_point((r.x, r.y), _reach(aX))
        if any(areas[j] >= area_ratio * aX and areas[j] >= 1.3 * med_area.get(f, np.inf) for j in idx):
            flags.at[i] = True
    return flags
