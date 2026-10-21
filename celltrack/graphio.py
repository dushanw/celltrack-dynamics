"""Shared lineage-graph representation and conversions.

Canonical in-memory form: a networkx.DiGraph whose nodes are detections with
attributes ``time`` (frame), ``label`` (mask label in that frame), ``coords``
(y, x) and whose edges go parent -> child with a ``weight`` in [0, 1] (the
linker's association confidence; 1.0 for linkers that do not provide one).

Canonical on-disk form: ``tracks.parquet`` (one row per detection) and
``lineages.parquet`` (one row per tracklet), plus an optional GEFF store.

A *tracklet* (track_id) is a maximal linear chain; it starts at an appearance
or right after a division and ends at a disappearance or a division. A
*lineage* (lineage_id) is a connected component of the graph.
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


def node_id(frame: int, label: int) -> int:
    """Deterministic node id so graphs from different linkers are comparable."""
    return int(frame) * 1_000_000 + int(label)


def split_node_id(nid: int) -> tuple[int, int]:
    return int(nid) // 1_000_000, int(nid) % 1_000_000


def graph_from_detections(det: pd.DataFrame) -> nx.DiGraph:
    """Edge-less graph containing every detection as a node."""
    G = nx.DiGraph()
    for f, l, y, x in zip(det["frame"].to_numpy(), det["label"].to_numpy(), det["y"].to_numpy(), det["x"].to_numpy()):
        G.add_node(node_id(f, l), time=int(f), label=int(l), coords=(float(y), float(x)))
    return G


def normalise_trackastra_graph(G_sol: nx.DiGraph, det: pd.DataFrame) -> nx.DiGraph:
    """Trackastra uses its own integer node ids; remap to node_id(frame, label)
    and add detections that received no edges as isolated nodes."""
    H = graph_from_detections(det)
    for u, v, d in G_sol.edges(data=True):
        nu, nv = G_sol.nodes[u], G_sol.nodes[v]
        a, b = node_id(nu["time"], nu["label"]), node_id(nv["time"], nv["label"])
        if a in H and b in H:
            H.add_edge(a, b, weight=float(d.get("weight", 1.0)))
    return H


def graph_to_tables(G: nx.DiGraph) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assign tracklet ids and lineage ids; return (tracks, lineages) tables."""
    # lineage = weakly connected component
    lineage_of = {}
    for i, comp in enumerate(nx.weakly_connected_components(G)):
        for n in comp:
            lineage_of[n] = i
    track_of: dict[int, int] = {}
    parent_track: dict[int, int] = {}
    next_tid = 0
    # tracklet starts: nodes with in-degree 0, or whose parent has out-degree >= 2
    starts = [n for n in G.nodes if G.in_degree(n) == 0 or G.out_degree(next(G.predecessors(n))) >= 2]
    for s in sorted(starts, key=lambda n: (G.nodes[n]["time"], n)):
        tid = next_tid
        next_tid += 1
        preds = list(G.predecessors(s))
        parent_track[tid] = track_of.get(preds[0], -1) if preds else -1
        n = s
        while True:
            track_of[n] = tid
            succ = list(G.successors(n))
            if len(succ) == 1 and G.in_degree(succ[0]) == 1:
                n = succ[0]
            else:
                break
    # (parents of division children may be assigned after children in rare orderings; fix up)
    for tid, s in enumerate(sorted(starts, key=lambda n: (G.nodes[n]["time"], n))):
        preds = list(G.predecessors(s))
        if preds:
            parent_track[tid] = track_of[preds[0]]

    rows = []
    for n, d in G.nodes(data=True):
        preds = list(G.predecessors(n))
        w_in = G.edges[preds[0], n].get("weight", 1.0) if preds else np.nan
        rows.append(
            dict(
                node_id=n,
                frame=int(d["time"]),
                label=int(d["label"]),
                y=float(d["coords"][0]),
                x=float(d["coords"][1]),
                track_id=track_of[n],
                lineage_id=lineage_of[n],
                parent_node=preds[0] if preds else -1,
                edge_conf_in=w_in,
                n_children=G.out_degree(n),
            )
        )
    cols = ["node_id", "frame", "label", "y", "x", "track_id", "lineage_id", "parent_node", "edge_conf_in", "n_children"]
    tracks = pd.DataFrame(rows, columns=cols).sort_values(["track_id", "frame"]).reset_index(drop=True)
    g = tracks.groupby("track_id")
    lineages = pd.DataFrame(
        dict(
            track_id=g.size().index,
            parent_track_id=[parent_track[t] for t in g.size().index],
            lineage_id=g["lineage_id"].first().values,
            start_frame=g["frame"].min().values,
            end_frame=g["frame"].max().values,
            length=g.size().values,
            n_children=g["n_children"].last().values,  # out-degree of the last node
            mean_edge_conf=g["edge_conf_in"].mean().values,
        )
    )
    return tracks, lineages


def tables_to_graph(tracks: pd.DataFrame) -> nx.DiGraph:
    G = nx.DiGraph()
    for r in tracks.itertuples(index=False):
        G.add_node(int(r.node_id), time=int(r.frame), label=int(r.label), coords=(float(r.y), float(r.x)))
    for r in tracks.itertuples(index=False):
        if r.parent_node != -1:
            G.add_edge(int(r.parent_node), int(r.node_id), weight=float(r.edge_conf_in))
    return G


def save_tracks(G: nx.DiGraph, out_dir: Path, name: str, write_geff: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tracks, lineages = graph_to_tables(G)
    tracks.to_parquet(out_dir / f"{name}_tracks.parquet", index=False)
    lineages.to_parquet(out_dir / f"{name}_lineages.parquet", index=False)
    if write_geff:
        try:
            import geff

            H = nx.DiGraph()
            for n, d in G.nodes(data=True):
                H.add_node(n, t=int(d["time"]), y=float(d["coords"][0]), x=float(d["coords"][1]), seg_id=int(d["label"]))
            for u, v, d in G.edges(data=True):
                H.add_edge(u, v, weight=float(d.get("weight", 1.0)))
            geff.write(H, out_dir / f"{name}.geff", axis_names=["t", "y", "x"], axis_types=["time", "space", "space"], overwrite=True)
        except Exception as e:  # GEFF is a convenience export; never fail the run on it
            print(f"[graphio] GEFF export failed: {e}")
    return tracks, lineages


def load_tracks(out_dir: Path, name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    out_dir = Path(out_dir)
    return (
        pd.read_parquet(out_dir / f"{name}_tracks.parquet"),
        pd.read_parquet(out_dir / f"{name}_lineages.parquet"),
    )


def graph_from_geff(path: Path, det: pd.DataFrame) -> nx.DiGraph:
    """Load a GEFF written by another tool (e.g. HOCT) and remap onto our node ids."""
    import geff

    G_in, _meta = geff.read(str(path))
    H = graph_from_detections(det)
    # find frame and label attribute names
    sample = next(iter(G_in.nodes(data=True)))[1]
    t_key = next(k for k in ("t", "time", "frame") if k in sample)
    l_key = next((k for k in ("seg_id", "label", "segm_id") if k in sample), None)
    if l_key is None:
        # fall back to nearest detection by position
        from scipy.spatial import cKDTree

        by_frame = {f: (g[["y", "x"]].to_numpy(), g["label"].to_numpy()) for f, g in det.groupby("frame")}
        trees = {f: cKDTree(xy) for f, (xy, _) in by_frame.items()}

        def remap(d):
            f = int(d[t_key])
            if f not in trees:
                return None
            dist, i = trees[f].query([d["y"], d["x"]])
            return node_id(f, by_frame[f][1][i]) if dist < 20 else None

    else:

        def remap(d):
            return node_id(int(d[t_key]), int(d[l_key]))

    ids = {n: remap(d) for n, d in G_in.nodes(data=True)}
    for u, v, d in G_in.edges(data=True):
        a, b = ids[u], ids[v]
        if a is not None and b is not None and a in H and b in H and a != b and H.in_degree(b) == 0:
            w = d.get("weight", d.get("prob", d.get("score", 1.0)))
            H.add_edge(a, b, weight=float(w))
    return H
