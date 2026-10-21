"""Phase 1: Cellpose-SAM instance segmentation of every frame.

Usage:
    # quick test on a handful of frames at two resolutions, writes QC overlays
    python -m celltrack.segment test --frames 0 80 160 240 320 400 480 575 --downsample 1 2

    # full run (resumable: already-written frames are skipped)
    python -m celltrack.segment run [--start 0 --stop 576] [--downsample 1]

    # post-process masks -> detections.parquet / debris.parquet (+ QC plots)
    python -m celltrack.segment features
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy import ndimage as ndi
from skimage import measure, segmentation, transform
from tqdm import tqdm

from .common import ensure_dir, load_config, pick_device, setup_logging

log = setup_logging("segment")


# --------------------------------------------------------------------------- model
def load_model(cfg):
    import torch
    from cellpose import models

    device = pick_device(cfg["segmentation"]["device"])
    # bf16 verified on MPS: IoU 0.993 vs fp32, ~18% faster. CUDA GPUs older than Ampere
    # (e.g. V100, compute capability 7.x) have no native bf16 -> fp32 there.
    bf16 = device == "mps" or (device == "cuda" and torch.cuda.get_device_capability()[0] >= 8)
    log.info("Loading Cellpose model %s on %s (bf16=%s)", cfg["segmentation"]["model"], device, bf16)
    model = models.CellposeModel(
        gpu=device != "cpu",
        pretrained_model=cfg["segmentation"]["model"],
        device=torch.device(device),
        use_bfloat16=bf16,
    )
    return model


def segment_frame(model, img: np.ndarray, cfg, downsample: int = 1) -> np.ndarray:
    seg = cfg["segmentation"]
    x = img
    if downsample > 1:
        x = transform.downscale_local_mean(img.astype(np.float32), (downsample, downsample))
    masks, _flows, _styles = model.eval(
        x,
        batch_size=16,
        flow_threshold=seg["flow_threshold"],
        cellprob_threshold=seg["cellprob_threshold"],
        min_size=max(15, seg["min_area_px"] // (4 * downsample * downsample)),
        normalize=True,
    )
    masks = masks.astype(np.uint16)
    if downsample > 1:
        masks = np.kron(masks, np.ones((downsample, downsample), dtype=np.uint16))
        masks = masks[: img.shape[0], : img.shape[1]]
    return masks


# --------------------------------------------------------------------------- post-processing
def split_concave_masks(masks: np.ndarray, min_area: int, area_factor: float = 1.3,
                        solidity_max: float = 0.88, min_peak_dist: int = 9) -> tuple[np.ndarray, int]:
    """Split under-segmented touching cells.

    Candidates are objects that are both large (area > area_factor * median area) and
    concave (solidity < solidity_max). They are split by a marker-controlled watershed on
    the distance transform (2-3 markers). A split is accepted only if every part is at
    least `min_area` px and at least 25% of the whole, so convex pre-mitotic cells
    (high solidity) are left intact. Returns (new_masks, n_split).
    """
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed

    props = measure.regionprops(masks)
    if not props:
        return masks, 0
    med_area = np.median([p.area for p in props])
    out = masks.copy()
    next_label = int(masks.max()) + 1
    n_split = 0
    for p in props:
        if p.area <= area_factor * med_area or p.solidity >= solidity_max:
            continue
        sl = p.slice
        sub = masks[sl] == p.label
        dist = ndi.distance_transform_edt(sub)
        dist = ndi.gaussian_filter(dist, 1.0)
        peaks = peak_local_max(dist, min_distance=min_peak_dist, labels=sub.astype(int), num_peaks=3, exclude_border=False)
        if len(peaks) < 2:
            continue
        markers = np.zeros(sub.shape, int)
        for k, (yy, xx) in enumerate(peaks, 1):
            markers[yy, xx] = k
        lab = watershed(-dist, markers, mask=sub)
        parts = [(lab == k).sum() for k in range(1, len(peaks) + 1)]
        if min(parts) < min_area or min(parts) < 0.25 * p.area:
            # try with only the two strongest peaks
            if len(peaks) > 2:
                markers[:] = 0
                for k, (yy, xx) in enumerate(peaks[:2], 1):
                    markers[yy, xx] = k
                lab = watershed(-dist, markers, mask=sub)
                parts = [(lab == k).sum() for k in (1, 2)]
                if min(parts) < min_area or min(parts) < 0.25 * p.area:
                    continue
            else:
                continue
        region = out[sl]
        for k in range(2, len(parts) + 1):
            region[(lab == k)] = next_label
            next_label += 1
        n_split += 1
    return out, n_split


def clean_masks(masks: np.ndarray, bubble: np.ndarray, cfg) -> tuple[np.ndarray, np.ndarray]:
    """Remove objects overlapping bubbles, fill holes, split concave merged cells, split off debris.

    Returns (cell_masks, debris_masks); both relabelled consecutively.
    """
    seg = cfg["segmentation"]
    masks = masks.copy()
    if seg.get("split_concave", True):
        masks, _ = split_concave_masks(masks, seg["min_area_px"], seg.get("split_area_factor", 1.3),
                                       seg.get("split_solidity_max", 0.88))
    # drop anything touching the bubble/artifact mask
    bad = np.unique(masks[bubble])
    bad = bad[bad > 0]
    if bad.size:
        masks[np.isin(masks, bad)] = 0
    # fill holes per object (cheap approximation: fill holes of the foreground, assign to nearest label)
    fg = masks > 0
    filled = ndi.binary_fill_holes(fg)
    holes = filled & ~fg
    if holes.any():
        # assign hole pixels the label of their enclosing object via nearest non-zero label
        idx = ndi.distance_transform_edt(~fg, return_distances=False, return_indices=True)
        masks[holes] = masks[idx[0][holes], idx[1][holes]]
    areas = np.bincount(masks.ravel())
    labels = np.arange(areas.size)
    small = labels[(areas > 0) & (areas < seg["min_area_px"])]
    small = small[small > 0]
    big = labels[areas > seg["max_area_px"]]
    big = big[big > 0]
    debris = np.where(np.isin(masks, small), masks, 0)
    cells = np.where(np.isin(masks, np.concatenate([small, big])), 0, masks)
    cells, _, _ = segmentation.relabel_sequential(cells)
    debris, _, _ = segmentation.relabel_sequential(debris)
    return cells.astype(np.uint16), debris.astype(np.uint16)


def frame_features(t: int, cells: np.ndarray, img: np.ndarray, cfg) -> pd.DataFrame:
    """Per-object regionprops plus neighbourhood descriptors."""
    H, W = cells.shape
    margin = cfg["segmentation"]["border_margin_px"]
    props = measure.regionprops_table(
        cells,
        intensity_image=img,
        properties=(
            "label",
            "centroid",
            "area",
            "perimeter",
            "eccentricity",
            "solidity",
            "major_axis_length",
            "minor_axis_length",
            "intensity_mean",
            "intensity_max",
            "intensity_std",
            "bbox",
        ),
    )
    df = pd.DataFrame(props)
    if df.empty:
        return df
    df = df.rename(columns={"centroid-0": "y", "centroid-1": "x"})
    df.insert(0, "frame", t)
    # touching neighbours: count distinct labels within a 2-px dilation ring
    dil = ndi.grey_dilation(cells, size=(5, 5))
    ero = ndi.grey_erosion(np.where(cells > 0, cells, np.iinfo(np.uint16).max), size=(5, 5))
    touch = {}
    boundary = (dil != cells) | (ero != cells)
    a = cells[boundary]
    b = dil[boundary]
    pairs = np.unique(np.stack([a, b], 1), axis=0)
    pairs = pairs[(pairs[:, 0] > 0) & (pairs[:, 1] > 0) & (pairs[:, 0] != pairs[:, 1])]
    for u, v in pairs:
        touch.setdefault(int(u), set()).add(int(v))
        touch.setdefault(int(v), set()).add(int(u))
    df["n_touching"] = df["label"].map(lambda l: len(touch.get(int(l), ()))).astype(int)
    # nearest-neighbour centroid distance
    from scipy.spatial import cKDTree

    xy = df[["x", "y"]].to_numpy()
    if len(xy) > 1:
        d, _ = cKDTree(xy).query(xy, k=2)
        df["nn_dist"] = d[:, 1]
    else:
        df["nn_dist"] = np.nan
    df["dist_to_border"] = np.minimum.reduce([df["x"], df["y"], W - 1 - df["x"], H - 1 - df["y"]])
    df["near_border"] = df["dist_to_border"] < margin
    df["circularity"] = 4 * np.pi * df["area"] / np.clip(df["perimeter"] ** 2, 1, None)
    return df


# --------------------------------------------------------------------------- commands
def cmd_test(args, cfg):
    P = cfg["paths"]
    flat = zarr.open_array(P["flat_zarr"], mode="r")
    bubble = np.load(P["bubble_mask"])
    qc = ensure_dir(Path(P["qc_dir"]) / "segment_test")
    model = load_model(cfg)
    results = []
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.color import label2rgb

    for ds in args.downsample:
        for t in args.frames:
            img = flat[t]
            t0 = time.time()
            masks = segment_frame(model, img, cfg, downsample=ds)
            dt = time.time() - t0
            cells, debris = clean_masks(masks, bubble, cfg)
            n = int(cells.max())
            areas = np.bincount(cells.ravel())[1:]
            results.append(
                dict(frame=t, downsample=ds, seconds=dt, n_cells=n, n_debris=int(debris.max()),
                     area_median=float(np.median(areas)) if n else np.nan)
            )
            log.info("ds=%d t=%d: %d cells, %d debris, %.1fs", ds, t, n, int(debris.max()), dt)
            # overlay of a central crop
            y0, x0, s = 1200, 1200, 768
            crop = img[y0:y0 + s, x0:x0 + s]
            lab = cells[y0:y0 + s, x0:x0 + s]
            fig, ax = plt.subplots(1, 2, figsize=(16, 8))
            ax[0].imshow(crop, cmap="gray", vmin=0, vmax=255)
            ax[0].set_title(f"t={t}")
            ax[1].imshow(label2rgb(lab, image=crop, bg_label=0, alpha=0.35))
            ax[1].contour(lab > 0, levels=[0.5], colors="yellow", linewidths=0.5)
            ax[1].set_title(f"ds={ds}  n={n} in frame")
            for a in ax:
                a.axis("off")
            plt.tight_layout()
            plt.savefig(qc / f"overlay_t{t:03d}_ds{ds}.png", dpi=70)
            plt.close(fig)
    df = pd.DataFrame(results)
    df.to_csv(qc / "test_summary.csv", index=False)
    print(df.to_string(index=False))
    print(df.groupby("downsample")[["seconds", "n_cells"]].mean())


def cmd_run(args, cfg):
    P = cfg["paths"]
    flat = zarr.open_array(P["flat_zarr"], mode="r")
    T, H, W = flat.shape
    ds = args.downsample or cfg["segmentation"]["downsample"]
    mpath = Path(P["masks_zarr"])
    if mpath.exists():
        masks = zarr.open_array(str(mpath), mode="r+")
    else:
        masks = zarr.create_array(
            store=str(mpath), shape=(T, H, W), chunks=(1, *cfg["ingest"]["chunk_yx"]),
            dtype=np.uint16, overwrite=False,
        )
        masks.attrs.update({"axes": "TYX", "model": cfg["segmentation"]["model"], "downsample": ds})
    done_path = mpath.parent / "masks_done.json"

    def read_done() -> set[int]:
        try:
            return set(json.load(open(done_path))) if done_path.exists() else set()
        except json.JSONDecodeError:  # concurrent writer mid-write
            time.sleep(0.2)
            return set(json.load(open(done_path)))

    def write_done(mine: set[int]):
        # merge with other concurrent workers before writing
        merged = read_done() | mine
        tmp = done_path.with_suffix(".tmp")
        json.dump(sorted(merged), open(tmp, "w"))
        tmp.replace(done_path)

    done = read_done()
    todo = [t for t in range(args.start, min(args.stop, T)) if t not in done]
    log.info("%d frames to segment (%d already done), downsample=%d", len(todo), len(done), ds)
    if not todo:
        return
    model = load_model(cfg)
    t_start = time.time()
    mine: set[int] = set()
    n_done = 0
    for i, t in enumerate(tqdm(todo, desc="segment")):
        if i % 10 == 0 and t in read_done():
            continue  # another worker got here first
        masks[t] = segment_frame(model, flat[t], cfg, downsample=ds)
        mine.add(t)
        n_done += 1
        if n_done % 3 == 0 or t == todo[-1]:
            write_done(mine)
    write_done(mine)
    log.info("Segmented %d frames in %.1f min", n_done, (time.time() - t_start) / 60)


def _features_one(t: int, cfg, clean: bool):
    """Worker: clean masks of frame t in place and compute features. Returns (det, debris, counts)."""
    P = cfg["paths"]
    flat = zarr.open_array(P["flat_zarr"], mode="r")
    masks = zarr.open_array(P["masks_zarr"], mode="r+")
    m = masks[t]
    if clean:
        bubble = np.load(P["bubble_mask"], mmap_mode="r")
        cells, debris = clean_masks(m, np.asarray(bubble), cfg)
        masks[t] = cells  # store the cleaned, relabelled version
    else:
        cells, debris = m, np.zeros_like(m)
    img = flat[t]
    df = frame_features(t, cells, img, cfg)
    dd = None
    if debris.max():
        dp = measure.regionprops_table(debris, intensity_image=img,
                                       properties=("label", "centroid", "area", "intensity_mean"))
        dd = pd.DataFrame(dp).rename(columns={"centroid-0": "y", "centroid-1": "x"})
        dd.insert(0, "frame", t)
    return df, dd, dict(frame=t, n_cells=int(cells.max()), n_debris=int(debris.max()))


def cmd_features(args, cfg):
    """Clean masks (once per frame, tracked in masks_cleaned.json so labels never change
    afterwards) and compute per-object features. Existing detections for frames not in
    this run are kept, so the table can be built incrementally."""
    from concurrent.futures import ProcessPoolExecutor

    P = cfg["paths"]
    masks = zarr.open_array(P["masks_zarr"], mode="r+")
    done_path = Path(P["masks_zarr"]).parent / "masks_done.json"
    cleaned_path = Path(P["masks_zarr"]).parent / "masks_cleaned.json"
    frames = sorted(json.load(open(done_path))) if done_path.exists() else list(range(masks.shape[0]))
    if args.frames:
        frames = args.frames
    cleaned = set(json.load(open(cleaned_path))) if cleaned_path.exists() else set()
    dets, debs, counts = [], [], []
    with ProcessPoolExecutor(args.workers) as ex:
        futs = {t: ex.submit(_features_one, t, cfg, t not in cleaned) for t in frames}
        for t, fut in tqdm(futs.items(), desc="features"):
            df, dd, c = fut.result()
            dets.append(df)
            if dd is not None:
                debs.append(dd)
            counts.append(c)
            cleaned.add(t)
    json.dump(sorted(cleaned), open(cleaned_path, "w"))
    det = pd.concat(dets, ignore_index=True)
    deb = pd.concat(debs, ignore_index=True) if debs else pd.DataFrame(columns=["frame", "label", "y", "x", "area", "intensity_mean"])
    cdf = pd.DataFrame(counts)
    # merge with previously computed frames
    det_p, deb_p = Path(P["detections"]), Path(P["debris"])
    if det_p.exists():
        old = pd.read_parquet(det_p)
        det = pd.concat([old[~old.frame.isin(frames)], det], ignore_index=True)
    if deb_p.exists():
        old = pd.read_parquet(deb_p)
        deb = pd.concat([old[~old.frame.isin(frames)], deb], ignore_index=True)
    cpath = Path(P["qc_dir"]) / "counts_per_frame.csv"
    if cpath.exists():
        old = pd.read_csv(cpath)
        cdf = pd.concat([old[~old.frame.isin(frames)], cdf], ignore_index=True)
    det = det.sort_values(["frame", "label"]).reset_index(drop=True)
    deb = deb.sort_values(["frame", "label"]).reset_index(drop=True)
    cdf = cdf.sort_values("frame").reset_index(drop=True)
    det.to_parquet(P["detections"], index=False)
    deb.to_parquet(P["debris"], index=False)
    qc = ensure_dir(P["qc_dir"])
    cdf.to_csv(qc / "counts_per_frame.csv", index=False)
    log.info("Wrote %d detections (%d frames in this run, %d total)", len(det), len(frames), det.frame.nunique())
    _qc_plots(det, cdf, qc, cfg)


def _qc_plots(det: pd.DataFrame, cdf: pd.DataFrame, qc: Path, cfg):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    hours = cdf["frame"] * cfg["acquisition"]["frame_interval_min"] / 60
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(hours, cdf["n_cells"], label="cells")
    ax[0, 0].plot(hours, cdf["n_debris"], label="debris")
    ax[0, 0].set_xlabel("time (h)"); ax[0, 0].set_ylabel("count"); ax[0, 0].legend()
    g = det.groupby("frame")
    gh = g.size().index.to_numpy() * cfg["acquisition"]["frame_interval_min"] / 60
    ax[0, 1].plot(gh, g["area"].median(), label="median area")
    ax[0, 1].fill_between(gh, g["area"].quantile(0.1), g["area"].quantile(0.9), alpha=0.3)
    ax[0, 1].set_xlabel("time (h)"); ax[0, 1].set_ylabel("area (px^2)")
    ax[1, 0].plot(gh, g["intensity_mean"].median())
    ax[1, 0].set_xlabel("time (h)"); ax[1, 0].set_ylabel("median mean intensity")
    ax[1, 1].plot(gh, g["n_touching"].apply(lambda s: (s >= 2).mean()))
    ax[1, 1].set_xlabel("time (h)"); ax[1, 1].set_ylabel("fraction with >=2 touching neighbours")
    plt.tight_layout()
    plt.savefig(qc / "segmentation_qc.png", dpi=80)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("test")
    s.add_argument("--frames", type=int, nargs="+", default=[0, 80, 160, 240, 320, 400, 480, 575])
    s.add_argument("--downsample", type=int, nargs="+", default=[1, 2])
    r = sub.add_parser("run")
    r.add_argument("--start", type=int, default=0)
    r.add_argument("--stop", type=int, default=10**9)
    r.add_argument("--downsample", type=int, default=None)
    r.add_argument("--device", default=None, help="override segmentation.device (e.g. cpu for an extra worker)")
    r.add_argument("--threads", type=int, default=None, help="torch CPU threads")
    f = sub.add_parser("features")
    f.add_argument("--frames", type=int, nargs="*", default=None)
    f.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    ensure_dir(cfg["paths"]["qc_dir"])
    if getattr(args, "device", None):
        cfg["segmentation"]["device"] = args.device
    if getattr(args, "threads", None):
        import torch

        torch.set_num_threads(args.threads)
    {"test": cmd_test, "run": cmd_run, "features": cmd_features}[args.cmd](args, cfg)


if __name__ == "__main__":
    main()
