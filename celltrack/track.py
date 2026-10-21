"""Phase 2: link detections across frames.

Linkers:
    trackastra  transformer association model (primary); greedy or chunked ILP
    laptrack    distance-based LAP with splitting and gap closing (baseline)
    hoct        Higher-Order Cell Tracking transformer via its CLI (comparison)

Usage:
    python -m celltrack.track trackastra --mode greedy [--start 0 --stop 576] [--crop Y0 X0 SIZE] [--name NAME]
    python -m celltrack.track trackastra --mode ilp --chunk 24
    python -m celltrack.track laptrack
    python -m celltrack.track hoct

Every run writes outputs/tracks/<name>_tracks.parquet, <name>_lineages.parquet and <name>.geff.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import dask.array as da
import networkx as nx
import numpy as np
import pandas as pd
import zarr

from .common import ensure_dir, load_config, pick_device, setup_logging
from .graphio import graph_from_detections, node_id, normalise_trackastra_graph, save_tracks

log = setup_logging("track")


# --------------------------------------------------------------------------- data access
def load_inputs(cfg, start, stop, crop=None):
    P = cfg["paths"]
    flat = zarr.open_array(P["flat_zarr"], mode="r")
    masks = zarr.open_array(P["masks_zarr"], mode="r")
    done_path = Path(P["masks_zarr"]).parent / "masks_done.json"
    done = set(json.load(open(done_path))) if done_path.exists() else set(range(masks.shape[0]))
    stop = min(stop, masks.shape[0])
    # only use the leading contiguous block of segmented frames
    t = start
    while t < stop and t in done:
        t += 1
    stop = t
    if stop <= start:
        raise RuntimeError(f"No segmented frames available from frame {start}")
    det = pd.read_parquet(P["detections"])
    det = det[(det.frame >= start) & (det.frame < stop)]
    if crop is not None:
        y0, x0, s = crop
        imgs = da.from_zarr(flat)[start:stop, y0:y0 + s, x0:x0 + s]
        segs = da.from_zarr(masks)[start:stop, y0:y0 + s, x0:x0 + s]
        det = det[(det.y >= y0) & (det.y < y0 + s) & (det.x >= x0) & (det.x < x0 + s)]
    else:
        imgs = da.from_zarr(flat)[start:stop]
        segs = da.from_zarr(masks)[start:stop]
    log.info("Frames %d..%d, %d detections%s", start, stop - 1, len(det), f", crop {crop}" if crop else "")
    return imgs, segs, det, start, stop


# --------------------------------------------------------------------------- trackastra
_TRACKASTRA_MODEL = None


def get_trackastra(cfg):
    global _TRACKASTRA_MODEL
    if _TRACKASTRA_MODEL is None:
        from trackastra.model import Trackastra

        tcfg = cfg["tracking"]["trackastra"]
        device = pick_device(tcfg["device"])
        log.info("Loading Trackastra %s on %s", tcfg["model"], device)
        _TRACKASTRA_MODEL = Trackastra.from_pretrained(tcfg["model"], device=device)
    return _TRACKASTRA_MODEL


def run_trackastra(imgs, segs, det, frames, cfg, mode="greedy", chunk=0, reverse=False):
    """Link a (sub)stack with Trackastra.

    Args:
        imgs, segs: arrays (T', Y, X) in *local* order.
        det: detections table (global frame numbers, labels matching segs).
        frames: global frame number of each local index (len T'). Allows strided
            or offset substacks; edges are returned in global frame coordinates.
        reverse: run the model on the time-reversed stack and flip the edges back.
    Returns:
        (solution graph in canonical ids, candidate-edge table with probabilities)
    """
    from trackastra.tracking import build_graph, track_greedy

    model = get_trackastra(cfg)
    tcfg = cfg["tracking"]["trackastra"]
    frames = np.asarray(frames)
    imgs_np = np.asarray(imgs)
    segs_np = np.asarray(segs)
    if reverse:
        imgs_np, segs_np, frames = imgs_np[::-1].copy(), segs_np[::-1].copy(), frames[::-1]
    t0 = time.time()
    pred = model._predict(imgs_np, segs_np, batch_size=None)  # library default per device
    log.info("Trackastra prediction took %.1f min (%d nodes, %d candidate edges)",
             (time.time() - t0) / 60, len(pred["nodes"]), len(pred["weights"]))
    # delta_t=2 adds candidate edges across one missing frame (gap closing in the ILP)
    cand = build_graph(nodes=pred["nodes"], weights=pred["weights"], max_neighbors=10,
                       delta_t=tcfg.get("delta_t", 2) if mode.startswith("ilp") else 1)

    # Modes: ilp / greedy allow 2 children (linker divisions); ilp_nodiv / greedy_nodiv are
    # strict 1-to-1 matchings, divisions are then proposed post hoc (states.confirm_divisions).
    if mode == "ilp":
        sol = _ilp_chunked(cand, chunk or tcfg.get("ilp_chunk", 24), allow_divisions=True)
    elif mode == "ilp_nodiv":
        sol = _ilp_chunked(cand, chunk or tcfg.get("ilp_chunk", 24), allow_divisions=False)
    elif mode == "greedy_nodiv":
        sol = track_greedy(cand, allow_divisions=False)
    else:
        sol = track_greedy(cand, allow_divisions=True)

    def canon(u, v):
        nu, nv = cand.nodes[u], cand.nodes[v]
        a = node_id(frames[int(nu["time"])], nu["label"])
        b = node_id(frames[int(nv["time"])], nv["label"])
        return (b, a) if reverse else (a, b)

    H = graph_from_detections(det)
    for u, v, d in sol.edges(data=True):
        a, b = canon(u, v)
        if a in H and b in H:
            H.add_edge(a, b, weight=float(d.get("weight", 1.0)))
    cand_rows = [(*canon(u, v), float(d["weight"])) for u, v, d in cand.edges(data=True)]
    cand_df = pd.DataFrame(cand_rows, columns=["u", "v", "weight"])
    return H, cand_df


def _ilp_chunked(cand: nx.DiGraph, chunk: int, allow_divisions: bool = True) -> nx.DiGraph:
    """Solve the ILP in overlapping temporal chunks and stitch by source frame."""
    from trackastra.tracking.ilp import track_ilp

    times = nx.get_node_attributes(cand, "time")
    t_min, t_max = min(times.values()), max(times.values())
    T = t_max - t_min + 1
    if T <= chunk + 4:
        return track_ilp(cand, ilp_config="gt", allow_divisions=allow_divisions)
    overlap = 4
    sol = nx.DiGraph()
    sol.add_nodes_from(cand.nodes(data=True))
    s = t_min
    keep_from = t_min
    while s <= t_max:
        e = min(s + chunk + overlap, t_max + 1)
        keep_to = e - overlap // 2 if e <= t_max else e  # frames whose outgoing edges we keep
        nodes = [n for n, t in times.items() if s <= t < e]
        sub = cand.subgraph(nodes).copy()
        log.info("ILP chunk frames %d..%d (%d nodes, %d edges), keeping sources %d..%d",
                 s, e - 1, sub.number_of_nodes(), sub.number_of_edges(), keep_from, keep_to - 1)
        t0 = time.time()
        part = track_ilp(sub, ilp_config="gt", allow_divisions=allow_divisions)
        log.info("  solved in %.1f s", time.time() - t0)
        for u, v, d in part.edges(data=True):
            if keep_from <= times[u] < keep_to:
                if sol.in_degree(v) == 0 and sol.out_degree(u) < (2 if allow_divisions else 1):
                    sol.add_edge(u, v, **d)
        keep_from = keep_to
        s = keep_to
        if e > t_max:
            break
    return sol


# --------------------------------------------------------------------------- laptrack
def run_laptrack(det, cfg):
    from laptrack import LapTrack

    lcfg = cfg["tracking"]["laptrack"]
    lt = LapTrack(
        cutoff=lcfg["max_distance"] ** 2,
        splitting_cutoff=lcfg["splitting_max_distance"] ** 2,
        merging_cutoff=False,
        gap_closing_cutoff=lcfg["max_distance"] ** 2,
        gap_closing_max_frame_count=lcfg["gap_closing_max_frame_count"],
    )
    d = det[["frame", "label", "x", "y"]].reset_index(drop=True)
    track_df, split_df, _merge_df = lt.predict_dataframe(d, coordinate_cols=["x", "y"], frame_col="frame", only_coordinate_cols=False)
    track_df = track_df.reset_index()
    H = graph_from_detections(det)
    track_df = track_df.sort_values(["track_id", "frame"])
    for _tid, g in track_df.groupby("track_id"):
        ids = [node_id(f, l) for f, l in zip(g["frame"], g["label"])]
        for a, b in zip(ids[:-1], ids[1:]):
            H.add_edge(a, b, weight=1.0)
    if len(split_df):
        last = track_df.groupby("track_id").last()
        first = track_df.groupby("track_id").first()
        for r in split_df.itertuples(index=False):
            p, c = last.loc[r.parent_track_id], first.loc[r.child_track_id]
            H.add_edge(node_id(p.frame, p.label), node_id(c.frame, c.label), weight=1.0)
    return H


# --------------------------------------------------------------------------- hoct
def run_hoct(imgs, segs, det, start, cfg, workdir: Path, device: str, crop=None, frames=None):
    """HOCT works from files; write the (sub)stack as TIFFs and call the CLI.
    Its GEFF carries crop-local positions (no mask labels), so nodes are matched
    back to detections by position.

    `frames` (optional) gives the original frame index of each plane of `imgs`; this
    lets callers pass subsampled or time-reversed stacks (validation) and still get
    edges in canonical node ids. Default: start, start+1, ..."""
    import tifffile

    from .graphio import graph_from_geff

    imgs, segs = np.asarray(imgs), np.asarray(segs)
    det_orig = det
    frames = np.arange(start, start + len(imgs)) if frames is None else np.asarray(frames)
    if not np.array_equal(frames, np.arange(start, start + len(imgs))):
        # relabel detections to a contiguous, forward substack 0..n-1 in the given plane order
        fmap = {int(f): i for i, f in enumerate(frames)}
        det = det[det.frame.isin(fmap)].copy()
        det["frame"] = det["frame"].map(fmap)
        start = 0
    workdir = ensure_dir(workdir)
    img_p, seg_p, out_p = workdir / "imgs.tif", workdir / "segs.tif", workdir / "hoct.geff"
    # photometric/metadata are explicit so that stacks with <= 4 frames are not written as RGB(A) pages
    tifffile.imwrite(img_p, imgs, bigtiff=True, photometric="minisblack", metadata={"axes": "TYX"})
    tifffile.imwrite(seg_p, segs, bigtiff=True, photometric="minisblack", metadata={"axes": "TYX"})
    cmd = ["hoct", "track", str(img_p), str(seg_p), "-o", str(out_p), "-ow", "-d", device]
    log.info("Running: %s", " ".join(cmd))
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    (workdir / "hoct.log").write_text(res.stdout + "\n" + res.stderr)
    if res.returncode != 0:
        raise RuntimeError(f"hoct failed (see {workdir/'hoct.log'}):\n{res.stderr[-2000:]}")
    log.info("HOCT finished in %.1f min", (time.time() - t0) / 60)
    # HOCT node times are relative to the substack; shift detections accordingly
    det_shift = det.copy()
    det_shift["frame"] = det_shift["frame"] - start
    if crop is not None:
        det_shift["y"] = det_shift["y"] - crop[0]
        det_shift["x"] = det_shift["x"] - crop[1]
    G = graph_from_geff(out_p, det_shift)
    for p in (img_p, seg_p):  # the temporary TIFF copies are large; keep only the GEFF + log
        p.unlink(missing_ok=True)
    # map substack plane index -> original frame
    plane_to_frame = {i: int(f) for i, f in enumerate(frames)} if start == 0 and len(frames) == len(imgs) else None
    H = graph_from_detections(det_orig[det_orig.frame.isin(set(int(f) for f in frames))])
    for u, v, d in G.edges(data=True):
        fu, lu = u // 1_000_000, u % 1_000_000
        fv, lv = v // 1_000_000, v % 1_000_000
        if plane_to_frame is not None:
            fu, fv = plane_to_frame[fu], plane_to_frame[fv]
        else:
            fu, fv = fu + start, fv + start
        a, b = node_id(fu, lu), node_id(fv, lv)
        if fu > fv:  # time-reversed stack: orient edges forward in real time
            a, b = b, a
        H.add_edge(a, b, **d)
    return H


# --------------------------------------------------------------------------- generic chunking
def stitch_chunks(solve, det: pd.DataFrame, start: int, stop: int, chunk: int, overlap: int,
                  cache_dir: Path | None = None, tag: str = "", keep_by_target: bool = False) -> nx.DiGraph:
    """Run `solve(a, b) -> nx.DiGraph (canonical ids)` on overlapping windows [a, b) and
    stitch by *source frame*: window k contributes the outgoing edges of frames in its
    keep-range; keep-ranges partition [start, stop). Per-window edge lists are cached as
    parquet so interrupted or incremental runs resume."""
    H = graph_from_detections(det)
    frame_of = {n: d["time"] for n, d in H.nodes(data=True)}
    a, keep_from = start, start
    while a < stop:
        b = min(a + chunk + overlap, stop)
        keep_to = stop if b >= stop else b - overlap // 2
        cache = cache_dir / f"{tag}_chunk_{a}_{b}.parquet" if cache_dir else None
        if cache is not None and cache.exists():
            edges = pd.read_parquet(cache)
            log.info("chunk %d..%d loaded from cache (%d edges)", a, b - 1, len(edges))
        else:
            t0 = time.time()
            G = solve(a, b)
            edges = pd.DataFrame([(u, v, float(d.get("weight", 1.0))) for u, v, d in G.edges(data=True)],
                                 columns=["u", "v", "weight"])
            if cache is not None:
                cache.parent.mkdir(parents=True, exist_ok=True)
                edges.to_parquet(cache, index=False)
            log.info("chunk %d..%d solved in %.1f min (%d edges)", a, b - 1, (time.time() - t0) / 60, len(edges))
        for u, v, w in edges.itertuples(index=False):
            key = v if keep_by_target else u
            if u in H and v in H and keep_from <= frame_of[key] < keep_to:
                if H.in_degree(v) == 0 and H.out_degree(u) < 2:
                    H.add_edge(int(u), int(v), weight=float(w))
        keep_from = keep_to
        a = keep_to
        if b >= stop:
            break
    return H


# --------------------------------------------------------------------------- main
def summarize(H: nx.DiGraph) -> dict:
    n = H.number_of_nodes()
    e = H.number_of_edges()
    div = sum(1 for x in H.nodes if H.out_degree(x) == 2)
    starts = sum(1 for x in H.nodes if H.in_degree(x) == 0)
    ends = sum(1 for x in H.nodes if H.out_degree(x) == 0)
    return dict(nodes=n, edges=e, divisions=div, track_starts=starts, track_ends=ends,
                linked_fraction=round(e / max(1, n), 3))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("linker", choices=["trackastra", "laptrack", "hoct"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=10**9)
    ap.add_argument("--crop", type=int, nargs=3, default=None, metavar=("Y0", "X0", "SIZE"))
    ap.add_argument("--name", default=None)
    ap.add_argument("--mode", default=None, help="trackastra: greedy | greedy_nodiv | ilp | ilp_nodiv")
    ap.add_argument("--chunk", type=int, default=0, help="ILP temporal chunk size (frames)")
    ap.add_argument("--reverse", action="store_true", help="trackastra: link the time-reversed movie")
    ap.add_argument("--device", default=None)
    ap.add_argument("--model", default=None, help="trackastra pretrained model name (general_2d | ctc)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    if args.device:
        cfg["tracking"]["trackastra"]["device"] = args.device
    if args.model:
        cfg["tracking"]["trackastra"]["model"] = args.model
    out_dir = ensure_dir(cfg["paths"]["tracks_dir"])
    imgs, segs, det, start, stop = load_inputs(cfg, args.start, args.stop, args.crop)
    mode = args.mode or cfg["tracking"]["trackastra"]["mode"]
    name = args.name or (f"{args.linker}_{mode}" if args.linker == "trackastra" else args.linker)
    if args.reverse:
        name += "_rev"
    if args.crop or args.start or stop < cfg["acquisition"]["n_frames"]:
        name += f"_f{start}-{stop}" + (f"_c{args.crop[0]}-{args.crop[1]}-{args.crop[2]}" if args.crop else "")

    t0 = time.time()
    cand_df = None
    cache_dir = out_dir / "chunks"
    if args.linker == "trackastra":
        if args.chunk and not args.reverse:
            # prediction + linking per temporal chunk (bounded memory, resumable)
            def solve(a, b):
                G, _ = run_trackastra(imgs[a - start:b - start], segs[a - start:b - start],
                                      det[(det.frame >= a) & (det.frame < b)], np.arange(a, b), cfg,
                                      mode=mode, chunk=b - a)
                return G

            H = stitch_chunks(solve, det, start, stop, args.chunk, 4, cache_dir,
                              tag=f"ta_{cfg['tracking']['trackastra']['model']}_{mode}" + (f"_c{args.crop[0]}-{args.crop[1]}-{args.crop[2]}" if args.crop else ""))
        else:
            H, cand_df = run_trackastra(imgs, segs, det, np.arange(start, stop), cfg, mode=mode,
                                        chunk=args.chunk, reverse=args.reverse)
    elif args.linker == "laptrack":
        H = run_laptrack(det, cfg)
    else:
        # HOCT's JIT model contains float64 tensors, which MPS cannot load -> CPU unless CUDA
        dev = args.device or ("cuda" if pick_device("cuda") == "cuda" else "cpu")
        rev = args.reverse  # time-reversed linking (for the reverse-consistency confidence of every edge)

        def solve(a, b):
            sl = slice(a - start, b - start)
            fr = np.arange(a, b)
            if rev:
                return run_hoct(np.asarray(imgs[sl])[::-1], np.asarray(segs[sl])[::-1],
                                det[(det.frame >= a) & (det.frame < b)], a, cfg,
                                out_dir / f"{name}_work", dev, crop=args.crop, frames=fr[::-1])
            return run_hoct(imgs[sl], segs[sl], det[(det.frame >= a) & (det.frame < b)], a, cfg,
                            out_dir / f"{name}_work", dev, crop=args.crop)

        if args.chunk:
            H = stitch_chunks(solve, det, start, stop, args.chunk, 4, cache_dir,
                              tag=("hoct_rev" if rev else "hoct") + (f"_c{args.crop[0]}-{args.crop[1]}-{args.crop[2]}" if args.crop else ""),
                              keep_by_target=rev)
        else:
            H = solve(start, stop)
    tracks, lineages = save_tracks(H, out_dir, name)
    if cand_df is not None:
        cand_df.to_parquet(out_dir / f"{name}_candidates.parquet", index=False)
    summ = summarize(H)
    summ.update(dict(name=name, linker=args.linker, mode=mode, start=start, stop=stop, crop=args.crop,
                     minutes=round((time.time() - t0) / 60, 2)))
    with open(out_dir / f"{name}_summary.json", "w") as fh:
        json.dump(summ, fh, indent=2)
    log.info("Done %s: %s", name, summ)


if __name__ == "__main__":
    main()
