#!/usr/bin/env python
"""Interactive demo of the celltrack-dynamics results - no raw data or GPU needed.

Setup (any laptop, Python >= 3.10):

    pip install "napari[pyqt6]" numpy pandas pyarrow

Run:

    python celltrack_demo.py                 # well A1, frames 0-95 (sparse, 0-4 h), 768x768 crop
    python celltrack_demo.py late            # well A1, frames 480-575 (dense, 20-24 h)
    python celltrack_demo.py --fraction 0.3  # random 30 % of lineages (less clutter)
    python celltrack_demo.py path/to/bundle.zip

The ~45 MB demo bundle is downloaded once from the GitHub release page and cached in
~/.cache/celltrack_demo. It contains the flattened phase-contrast crop, the Cellpose-SAM masks,
the per-cell state table (track id, lineage, state, speed, events, link confidence) and the
population transition matrices.

In napari:
    - play / scrub time with the slider at the bottom (2.5 min per frame)
    - scroll to zoom, drag to pan; toggle layers in the panel on the left
    - "tracks (colour = track id)": each trail is one track; a trail that changes colour is an
      identity switch, a trail that stops is a lost/left/merged/divided cell
    - "tracks (colour = state)": blue = stationary, yellow = slow, red = fast
    - "events": green circle = confirmed division, magenta = death
    - "masks": segmentation (hidden by default)
    - hover a cell with the tracks layer selected to see its track id, state and speed
The lag-1 transition matrix is shown in a dock on the right; all matrices are in the bundle's
transition/ folder (CSV + PNG).
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = "dushanw/celltrack-dynamics"
RELEASE_TAG = "demo-v1"
RELEASE = f"https://github.com/{REPO}/releases/download/{RELEASE_TAG}"
BUNDLES = {
    "early": f"{RELEASE}/celltrack_demo_A1_early.zip",
    "late": f"{RELEASE}/celltrack_demo_A1_late.zip",
}
CACHE = Path.home() / ".cache" / "celltrack_demo"

# keep in sync with celltrack/view.py
STATE_COLORS = {
    "stationary": (0.3, 0.6, 1.0),
    "slow": (1.0, 0.85, 0.2),
    "fast": (1.0, 0.3, 0.2),
    "dividing": (0.2, 1.0, 0.4),
    "dying": (0.9, 0.2, 0.9),
}


def download(url: str, dest: Path):
    """Plain HTTPS download. If the asset is not publicly readable (e.g. the repository was made
    private again), fall back to the GitHub CLI, which can read release assets you have access to."""
    import shutil
    import subprocess
    import urllib.error

    print(f"downloading {url}")

    def hook(n, bs, total):
        if total > 0:
            sys.stdout.write(f"\r  {n * bs / 1e6:6.1f} / {total / 1e6:.1f} MB")
            sys.stdout.flush()

    try:
        urllib.request.urlretrieve(url, dest, hook)
        print()
        return
    except urllib.error.HTTPError as e:
        if e.code not in (401, 403, 404):
            raise
        print(f"\n  direct download failed ({e.code}); trying the GitHub CLI")
    if shutil.which("gh") is None:
        sys.exit("The release asset is not publicly readable. Either install the GitHub CLI and run `gh auth login`\n"
                 "(you need access to " + REPO + "), or ask for the bundle zip and pass its path.")
    r = subprocess.run(["gh", "release", "download", RELEASE_TAG, "--repo", REPO, "--pattern", dest.name,
                        "--dir", str(dest.parent), "--clobber"])
    if r.returncode != 0 or not dest.exists():
        sys.exit("gh release download failed; run `gh auth login` and make sure you have access to " + REPO)


def fetch(name_or_path: str) -> Path:
    """Return a directory with the extracted bundle, downloading it if needed."""
    if name_or_path in BUNDLES:
        url = BUNDLES[name_or_path]
        zpath = CACHE / Path(url).name
        if not zpath.exists():
            CACHE.mkdir(parents=True, exist_ok=True)
            download(url, zpath)
    else:
        zpath = Path(name_or_path)
        if not zpath.exists():
            sys.exit(f"{zpath} not found (use 'early', 'late' or a path to a bundle zip)")
    ddir = zpath.with_suffix("")
    if not (ddir / "meta.json").exists():
        with zipfile.ZipFile(zpath) as z:
            z.extractall(ddir)
    return ddir


def load(ddir: Path, fraction: float, seed: int):
    meta = json.loads((ddir / "meta.json").read_text())
    img = np.load(ddir / "image.npy")
    lab = np.load(ddir / "masks.npy")
    d = pd.read_parquet(ddir / "states.parquet")
    if fraction < 1.0:
        ids = d["lineage_id"].unique()
        keep = np.random.default_rng(seed).choice(ids, size=max(1, int(round(fraction * len(ids)))), replace=False)
        d = d[d["lineage_id"].isin(keep)]
        print(f"showing {fraction:.0%} of lineages: {len(keep)} of {len(ids)} (seed {seed})")
    return meta, img, lab, d


def tracks_arrays(d: pd.DataFrame):
    data = d[["track_id", "frame_local", "y", "x"]].to_numpy(dtype=float)
    graph = {}
    tid_of = d.set_index("node_id")["track_id"].to_dict()
    for tid, g in d.groupby("track_id"):
        p = g.sort_values("frame_local").iloc[0]["parent_node"]
        if p != -1 and p in tid_of and tid_of[p] != tid:
            graph[int(tid)] = [int(tid_of[p])]
    codes = {s: i for i, s in enumerate(STATE_COLORS)}
    props = {
        "state": d["state"].astype(str).to_numpy(),
        "state_code": d["state"].map(codes).fillna(1).to_numpy(dtype=float),
        "speed_px_per_frame": d["speed_px_per_frame"].fillna(0).to_numpy(),
        "lineage_id": d["lineage_id"].to_numpy(),
    }
    return data, graph, props


def show(meta, img, lab, d, ddir: Path):
    import napari

    data, graph, props = tracks_arrays(d)
    t0, t1 = meta["start"], meta["stop"] - 1
    h0, h1 = t0 * meta["frame_interval_min"] / 60, t1 * meta["frame_interval_min"] / 60
    v = napari.Viewer(title=f"celltrack-dynamics demo - {meta['source']}  frames {t0}-{t1} ({h0:.1f}-{h1:.1f} h), "
                            f"{d.track_id.nunique()} tracks")
    v.add_image(img, name="phase contrast (flattened)", contrast_limits=(0, 255))
    v.add_labels(lab, name="masks", opacity=0.25, visible=False)
    v.add_tracks(data, graph=graph, properties=props, name="tracks (colour = track id)", tail_length=10, tail_width=3)
    v.add_tracks(data, graph=graph, properties=props, name="tracks (colour = state)", tail_length=10, tail_width=3,
                 color_by="state_code", colormap="turbo", visible=False)
    ev = d[d.state.isin(["dividing", "dying"])]
    if len(ev):
        v.add_points(ev[["frame_local", "y", "x"]].to_numpy(), size=30, face_color="transparent",
                     border_color=[STATE_COLORS[s] for s in ev.state], name="events (green=division, magenta=death)")

    # lag-1 transition matrix of the whole movie as a dock widget
    heat = ddir / "transition" / "heatmap_lag1_all.png"
    if heat.exists():
        try:
            from qtpy.QtGui import QPixmap
            from qtpy.QtWidgets import QLabel

            lbl = QLabel()
            lbl.setPixmap(QPixmap(str(heat)).scaledToWidth(520))
            lbl.setToolTip("P(state at t+1 | state at t), whole movie; CSVs in the bundle's transition/ folder")
            v.window.add_dock_widget(lbl, name="transition matrix (lag 1 = 2.5 min)", area="right")
        except Exception as e:  # the viewer is still useful without the dock
            print("could not add heatmap dock:", e)

    print(f"\n{meta['source']}: frames {t0}-{t1}, crop y={meta['crop'][0]} x={meta['crop'][1]} size={meta['crop'][2]} px")
    print(f"{len(d)} cell-frames, {d.track_id.nunique()} tracks, {d.lineage_id.nunique()} lineages, "
          f"{(d.state == 'dividing').sum()} divisions, {(d.state == 'dying').sum()} deaths")
    print("state fractions:", d.state.value_counts(normalize=True).round(3).to_dict())
    p = ddir / "transition" / "probs_lag1_all.csv"
    if p.exists():
        print("\nlag-1 transition probabilities (rows: state at t, columns: state at t+1):")
        print(pd.read_csv(p, index_col=0).round(2).to_string())
    napari.run()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bundle", nargs="?", default="early", help="'early' (default), 'late' or a path to a bundle zip")
    ap.add_argument("--fraction", type=float, default=1.0, help="show only this fraction of lineages (0-1]")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-gui", action="store_true", help="download + load only (smoke test)")
    args = ap.parse_args(argv)
    ddir = fetch(args.bundle)
    meta, img, lab, d = load(ddir, args.fraction, args.seed)
    if args.no_gui:
        print(meta, img.shape, lab.shape, len(d))
        return
    show(meta, img, lab, d, ddir)


if __name__ == "__main__":
    main()
