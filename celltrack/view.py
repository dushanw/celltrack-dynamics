"""Phase 5b: visualisation.

    python -m celltrack.view napari --name trackastra_ilp [--crop Y0 X0 SIZE] [--start 0 --stop 100]
        opens napari with image, labels, tracks (coloured by state) and lineage graph

    python -m celltrack.view movie --name trackastra_ilp --crop 1200 1200 768 --start 0 --stop 120
        writes an MP4 overlay (tracks + state colours + event markers) under outputs/qc/movies

    python -m celltrack.view export-demo --name hoct --crop 1200 1200 768 --start 0 --stop 96
        writes a small self-contained zip (image, masks, states, matrices) for demo/celltrack_demo.py

    python -m celltrack.view export-web --name hoct --crop 1200 1200 768 --start 480 --stop 576
        writes JPEG frames + mask overlays + tracks.json into demo/web/data/<tag>/ for the browser viewer
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

from .common import ensure_dir, load_config, setup_logging

log = setup_logging("view")

STATE_COLORS = {
    "stationary": (0.3, 0.6, 1.0),
    "slow": (1.0, 0.85, 0.2),
    "fast": (1.0, 0.3, 0.2),
    "dividing": (0.2, 1.0, 0.4),
    "dying": (0.9, 0.2, 0.9),
}


def load_states(cfg, name, start, stop, crop, fraction: float = 1.0, seed: int = 0):
    """Load the per-node state table for a window/crop. `fraction` < 1 keeps a random
    subset of *lineages* (so mothers and daughters stay together), reproducible via `seed`."""
    p = Path(cfg["paths"]["analysis_dir"]) / f"states_{name}.parquet"
    if p.exists():
        d = pd.read_parquet(p)
    else:  # fall back to raw tracks without states
        from .graphio import load_tracks

        d, _ = load_tracks(cfg["paths"]["tracks_dir"], name)
        d["state"] = "slow"
    if fraction < 1.0:
        key = "lineage_id" if "lineage_id" in d else "track_id"
        ids = d[key].unique()
        keep = np.random.default_rng(seed).choice(ids, size=max(1, int(round(fraction * len(ids)))), replace=False)
        d = d[d[key].isin(keep)]
        log.info("Showing %.1f%% of %s: %d of %d (seed %d)", 100 * fraction, key.replace("_id", "s"), len(keep), len(ids), seed)
    d = d[(d.frame >= start) & (d.frame < stop)]
    if crop:
        y0, x0, s = crop
        d = d[(d.y >= y0) & (d.y < y0 + s) & (d.x >= x0) & (d.x < x0 + s)].copy()
        d["y"] -= y0
        d["x"] -= x0
    d["frame_local"] = d.frame - start
    return d


def napari_tracks_arrays(d: pd.DataFrame):
    """napari Tracks layer wants [track_id, t, y, x] and a graph dict child->[parents]."""
    data = d[["track_id", "frame_local", "y", "x"]].to_numpy(dtype=float)
    graph = {}
    if "parent_node" in d:
        tid_of = d.set_index("node_id")["track_id"].to_dict()
        for tid, g in d.groupby("track_id"):
            p = g.sort_values("frame_local").iloc[0]["parent_node"]
            if p != -1 and p in tid_of and tid_of[p] != tid:
                graph[int(tid)] = [int(tid_of[p])]
    codes = {s: i for i, s in enumerate(STATE_COLORS)}  # stationary=0 .. dying=4
    props = {"state": d["state"].astype(str).to_numpy(),
             "state_code": d["state"].map(codes).fillna(1).to_numpy(dtype=float),
             "speed": d.get("speed_px_per_frame", pd.Series(np.zeros(len(d)))).fillna(0).to_numpy()}
    return data, graph, props


def build_viewer(img, lab, d: pd.DataFrame, title: str, show: bool = True):
    """Create a napari viewer with the standard layer set:
    image, masks (hidden), tracks coloured by id, tracks coloured by state (hidden), event markers.
    `img`/`lab` are (T, Y, X) arrays (numpy, zarr or dask); `d` is the output of `load_states`
    (positions already crop-local, `frame_local` starting at 0)."""
    import napari

    data, graph, props = napari_tracks_arrays(d)
    v = napari.Viewer(title=title, show=show)
    v.add_image(img, name="phase contrast (flattened)", contrast_limits=(0, 255))
    v.add_labels(lab, name="masks", opacity=0.25, visible=False)
    # one colour per track id: identity swaps show up as a trail changing colour
    v.add_tracks(data, graph=graph, properties=props, name="tracks (colour = track id)", tail_length=10, tail_width=3)
    v.add_tracks(data, graph=graph, properties=props, name="tracks (colour = state)", tail_length=10, tail_width=3,
                 color_by="state_code", colormap="turbo", visible=False)
    ev = d[d.state.isin(["dividing", "dying"])]
    if len(ev):
        v.add_points(ev[["frame_local", "y", "x"]].to_numpy(), size=30, face_color="transparent",
                     border_color=[STATE_COLORS[s] for s in ev.state], name="events (green=division, magenta=death)")
    return v


def _load_window(cfg, args):
    """Image + masks for the requested window/crop (lazy dask arrays for the full field)."""
    flat = zarr.open_array(cfg["paths"]["flat_zarr"], mode="r")
    masks = zarr.open_array(cfg["paths"]["masks_zarr"], mode="r")
    start, stop = args.start, min(args.stop, flat.shape[0])
    if args.crop:
        y0, x0, s = args.crop
        img = flat[start:stop, y0:y0 + s, x0:x0 + s]
        lab = masks[start:stop, y0:y0 + s, x0:x0 + s]
    else:
        import dask.array as da

        img = da.from_zarr(flat)[start:stop]
        lab = da.from_zarr(masks)[start:stop]
    return img, lab, start, stop


def cmd_napari(args, cfg):
    import napari

    img, lab, start, stop = _load_window(cfg, args)
    d = load_states(cfg, args.name, start, stop, args.crop, args.fraction, args.seed)
    sub = f", {args.fraction:.0%} of lineages, seed {args.seed}" if args.fraction < 1.0 else ""
    build_viewer(img, lab, d, f"celltrack {args.name}  (frames {start}-{stop - 1}, {len(d.track_id.unique())} tracks{sub})")
    napari.run()


DEMO_COLUMNS = ["node_id", "frame", "track_id", "lineage_id", "parent_node", "y", "x", "area", "circularity",
                "n_touching", "in_cluster", "speed_px_per_frame", "state", "start_type", "end_type",
                "division_quality", "edge_conf_in", "bidir", "hours"]


def cmd_export_demo(args, cfg):
    """Write a self-contained demo bundle (zip) for `demo/celltrack_demo.py`: cropped image + masks
    as .npy, the per-cell state table for the crop, the transition-matrix figures/CSVs and a
    meta.json. Everything a collaborator needs to browse the tracking without the raw data."""
    import json
    import shutil
    import tempfile
    import zipfile

    crop = args.crop or (1200, 1200, 768)
    args.crop = crop
    img, lab, start, stop = _load_window(cfg, args)
    d = load_states(cfg, args.name, start, stop, crop, args.fraction, args.seed)
    cols = [c for c in DEMO_COLUMNS if c in d.columns] + ["frame_local"]
    out = Path(args.out) if args.out else Path(cfg["paths"]["analysis_dir"]) / f"demo_{args.name}_f{start}-{stop}_y{crop[0]}_x{crop[1]}_s{crop[2]}.zip"
    out.parent.mkdir(parents=True, exist_ok=True)
    tdir = Path(tempfile.mkdtemp())
    np.save(tdir / "image.npy", np.asarray(img, dtype=np.uint8))
    np.save(tdir / "masks.npy", np.asarray(lab))
    d[cols].reset_index(drop=True).to_parquet(tdir / "states.parquet", index=False)
    tr = Path(cfg["paths"]["analysis_dir"]) / f"transition_{args.name}"
    (tdir / "transition").mkdir()
    for p in list(tr.glob("heatmap_lag*_all.png")) + list(tr.glob("heatmap_lag*_highconf.png")) + \
             list(tr.glob("probs_lag*_all.csv")) + list(tr.glob("dwell_times*.csv")) + [tr / "summary.json"]:
        if p.exists():
            shutil.copy(p, tdir / "transition" / p.name)
    meta = dict(name=args.name, start=start, stop=stop, crop=list(crop), shape=list(np.asarray(img).shape),
                frame_interval_min=cfg["acquisition"]["frame_interval_min"], pixel_size_um=cfg["acquisition"]["pixel_size_um"],
                n_tracks=int(d.track_id.nunique()), n_cells=len(d), source=Path(cfg["paths"]["raw_tif"]).name,
                states=list(STATE_COLORS))
    (tdir / "meta.json").write_text(json.dumps(meta, indent=1))
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(tdir.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(tdir))
    shutil.rmtree(tdir)
    log.info("Wrote %s (%.1f MB): %d frames, %d tracks", out, out.stat().st_size / 1e6, stop - start, meta["n_tracks"])


def cmd_movie(args, cfg):
    import imageio.v2 as imageio
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    flat = zarr.open_array(cfg["paths"]["flat_zarr"], mode="r")
    start, stop = args.start, min(args.stop, flat.shape[0])
    crop = args.crop or (1024, 1024, 1024)
    y0, x0, s = crop
    d = load_states(cfg, args.name, start, stop, crop, args.fraction, args.seed)
    out = ensure_dir(Path(cfg["paths"]["qc_dir"]) / "movies")
    sub = f"_p{int(round(100 * args.fraction))}s{args.seed}" if args.fraction < 1.0 else ""
    path = out / f"{args.name}_f{start}-{stop}_y{y0}_x{x0}_s{s}{sub}.mp4"
    from matplotlib.collections import LineCollection

    tail = 8
    big = s > 1024  # full-field overview: smaller markers, no track ids
    size_in = 8 if not big else 12
    writer = imageio.get_writer(path, fps=args.fps, codec="libx264", quality=7, macro_block_size=8)
    d = d.sort_values(["track_id", "frame"])
    by_frame = {f: g for f, g in d.groupby("frame")}
    fig = plt.figure(figsize=(size_in, size_in), dpi=100)
    canvas = FigureCanvasAgg(fig)
    for t in range(start, stop):
        fig.clf()
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(flat[t, y0:y0 + s, x0:x0 + s], cmap="gray", vmin=0, vmax=255)
        cur = by_frame.get(t, d.iloc[0:0])
        # trails: consecutive positions of each track within the last `tail` frames
        win = pd.concat([by_frame[f] for f in range(max(start, t - tail + 1), t + 1) if f in by_frame]) if t >= start else d.iloc[0:0]
        if len(win):
            win = win.sort_values(["track_id", "frame"])
            same = win.track_id.to_numpy()[1:] == win.track_id.to_numpy()[:-1]
            xy = win[["x", "y"]].to_numpy()
            segs = np.stack([xy[:-1][same], xy[1:][same]], axis=1)
            ax.add_collection(LineCollection(segs, colors="yellow", linewidths=0.8 if not big else 0.5, alpha=0.7))
        for st, col in STATE_COLORS.items():
            c = cur[cur.state == st]
            if len(c):
                ms = (9 if st in ("dividing", "dying") else 5) * (0.6 if big else 1)
                ax.plot(c.x, c.y, "o", mfc="none", mec=col, ms=ms, mew=1.5 if not big else 1.0)
        if not big:
            for r in cur.itertuples():
                ax.text(r.x + 6, r.y - 6, str(int(r.track_id)), color="white", fontsize=5, alpha=0.8)
        ax.text(5, 20, f"t={t}  ({t * cfg['acquisition']['frame_interval_min'] / 60:.2f} h)   n={len(cur)}",
                color="white", fontsize=10, backgroundcolor="black")
        ax.set_xlim(0, s); ax.set_ylim(s, 0); ax.axis("off")
        canvas.draw()
        frame = np.asarray(canvas.buffer_rgba())[..., :3]
        writer.append_data(frame)
    writer.close()
    plt.close(fig)
    log.info("Wrote %s", path)


def _track_color(tid: int):
    """Stable, well-separated colour per track id (golden-ratio hue), as (r, g, b) in 0..255."""
    import colorsys

    h = (tid * 0.618033988749895) % 1.0
    r, g, b = colorsys.hsv_to_rgb(h, 0.65, 1.0)
    return int(r * 255), int(g * 255), int(b * 255)


def cmd_export_web(args, cfg):
    """Write a dataset for the browser viewer in demo/web: JPEG frames, RGBA mask overlays coloured
    by track id, tracks.json (per track: lineage, parent, [frame, x, y, state, speed, confidence]),
    the transition matrices as JSON and meta.json. ~10-15 MB per 96-frame 768 px crop."""
    import json
    import shutil

    import imageio.v3 as iio

    crop = args.crop or (1200, 1200, 768)
    args.crop = crop
    img, lab, start, stop = _load_window(cfg, args)
    img, lab = np.asarray(img), np.asarray(lab)
    d = load_states(cfg, args.name, start, stop, crop, args.fraction, args.seed)
    tag = args.tag or f"{Path(cfg['paths']['raw_tif']).stem.replace('_stack', '')}_f{start}-{stop}"
    from .common import ROOT

    out = Path(args.out or (ROOT / "demo" / "web" / "data")) / tag
    for sub in ("frames", "masks"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    # frames + mask overlays (colour = track id so identities can be followed on the masks too)
    tid_of = {(int(f), int(l)): int(t) for f, l, t in d[["frame", "label", "track_id"]].itertuples(index=False)}
    for i in range(stop - start):
        iio.imwrite(out / "frames" / f"{i:03d}.jpg", img[i], quality=85)
        m = lab[i]
        labels = np.unique(m)
        labels = labels[labels > 0]
        lut = np.zeros((int(m.max()) + 1, 4), dtype=np.uint8)
        for l in labels:
            tid = tid_of.get((start + i, int(l)))
            lut[l] = (*_track_color(tid), 255) if tid is not None else (140, 140, 140, 255)
        iio.imwrite(out / "masks" / f"{i:03d}.png", lut[m], compress_level=6)

    # tracks
    codes = {s: k for k, s in enumerate(STATE_COLORS)}
    tid_track = d.set_index("node_id")["track_id"].to_dict()
    tracks = []
    for tid, g in d.sort_values("frame").groupby("track_id"):
        p = g.iloc[0]["parent_node"]
        parent = int(tid_track[p]) if p != -1 and p in tid_track and tid_track[p] != tid else None
        pts = [[int(r.frame_local), round(float(r.x), 1), round(float(r.y), 1), int(codes.get(r.state, 1)),
                round(float(r.speed_px_per_frame), 1) if pd.notna(r.speed_px_per_frame) else None,
                round(float(r.edge_conf_in), 2) if "edge_conf_in" in g and pd.notna(r.edge_conf_in) else None]
               for r in g.itertuples()]
        tracks.append(dict(id=int(tid), lineage=int(g.iloc[0]["lineage_id"]), parent=parent,
                           start=str(g.iloc[0].get("start_type", "")), end=str(g.iloc[-1].get("end_type", "")),
                           color=_track_color(int(tid)), pts=pts))
    ev = d[d.state.isin(["dividing", "dying"])]
    events = [[int(r.frame_local), round(float(r.x), 1), round(float(r.y), 1), str(r.state), int(r.track_id)] for r in ev.itertuples()]

    # transition matrices (whole movie) + summary numbers
    tr = Path(cfg["paths"]["analysis_dir"]) / f"transition_{args.name}"
    mats = {}
    for lag in (1, 4, 12):
        p = tr / f"probs_lag{lag}_all.csv"
        if p.exists():
            m = pd.read_csv(p, index_col=0)
            n = pd.read_csv(tr / f"counts_lag{lag}_all.csv", index_col=0).sum(axis=1) if (tr / f"counts_lag{lag}_all.csv").exists() else None
            mats[f"lag{lag}"] = dict(rows=list(m.index), cols=list(m.columns), probs=m.round(3).values.tolist(),
                                     n=[int(x) for x in n.reindex(m.index).fillna(0)] if n is not None else None)
    summary = json.loads((tr / "summary.json").read_text()) if (tr / "summary.json").exists() else {}
    meta = dict(tag=tag, name=args.name, source=Path(cfg["paths"]["raw_tif"]).name, start=start, stop=stop, crop=list(crop),
                size=int(crop[2]), n_frames=stop - start, frame_interval_min=cfg["acquisition"]["frame_interval_min"],
                pixel_size_um=cfg["acquisition"]["pixel_size_um"], n_tracks=len(tracks), n_cells=len(d),
                n_lineages=int(d.lineage_id.nunique()), n_divisions=int((d.state == "dividing").sum()),
                n_deaths=int((d.state == "dying").sum()), states=list(STATE_COLORS),
                state_colors={s: [int(255 * c) for c in col] for s, col in STATE_COLORS.items()},
                state_fraction=d.state.value_counts(normalize=True).round(3).to_dict(),
                speed_thresholds_px_per_frame=summary.get("speed_thresholds_px_per_frame"),
                dwell_minutes=summary.get("dwell_minutes_uncensored"), event_rates=summary.get("event_rates"),
                movie_state_fraction=summary.get("state_fraction"), transitions=mats)
    (out / "meta.json").write_text(json.dumps(meta, indent=1, default=float))
    (out / "tracks.json").write_text(json.dumps(dict(tracks=tracks, events=events), separators=(",", ":")))
    # register in the viewer's dataset index
    idx_p = out.parent / "index.json"
    idx = json.loads(idx_p.read_text()) if idx_p.exists() else []
    idx = [e for e in idx if e["tag"] != tag] + [dict(tag=tag, label=args.label or f"{meta['source']}  frames {start}-{stop - 1} "
                                                       f"({start * meta['frame_interval_min'] / 60:.0f}-{stop * meta['frame_interval_min'] / 60:.0f} h)")]
    idx_p.write_text(json.dumps(sorted(idx, key=lambda e: e["tag"]), indent=1))
    size = sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) / 1e6
    log.info("Wrote %s (%.1f MB): %d frames, %d tracks, %d events", out, size, stop - start, len(tracks), len(events))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["napari", "movie", "export-demo", "export-web"])
    ap.add_argument("--tag", default=None, help="export-web: dataset id (default <well>_f<start>-<stop>)")
    ap.add_argument("--label", default=None, help="export-web: label shown in the dataset menu")
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default=None, help="export-demo: output zip path")
    ap.add_argument("--name", default="hoct")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=10**9)
    ap.add_argument("--crop", type=int, nargs=3, default=None, metavar=("Y0", "X0", "SIZE"))
    ap.add_argument("--fps", type=int, default=6)
    ap.add_argument("--fraction", type=float, default=1.0,
                    help="randomly show only this fraction of lineages (0-1], e.g. 0.1 for 10%%; mothers and daughters are kept together")
    ap.add_argument("--seed", type=int, default=0, help="random seed for --fraction (change it to draw a different subset)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    {"napari": cmd_napari, "movie": cmd_movie, "export-demo": cmd_export_demo, "export-web": cmd_export_web}[args.cmd](args, cfg)


if __name__ == "__main__":
    main()
