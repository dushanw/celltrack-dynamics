"""Phase 0: TIFF -> chunked zarr, temporal-median background, flattened copy, bubble mask.

Usage:
    python -m celltrack.ingest [--config configs/default.yaml] [--skip-flat]
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import tifffile
import zarr
from scipy import ndimage as ndi
from skimage import measure
from tqdm import tqdm

from .common import ensure_dir, load_config, setup_logging

log = setup_logging("ingest")


def tif_to_zarr(tif_path: str, zarr_path: str, chunk_yx, n_workers: int = 6) -> zarr.Array:
    with tifffile.TiffFile(tif_path) as tif:
        series = tif.series[0]
        shape = series.shape
        dtype = series.dtype
        log.info("TIFF shape %s dtype %s axes %s", shape, dtype, series.axes)
        # Uncompressed, one strip per page: collect byte offsets so threads can read
        # with independent file handles (tifffile page access is not thread-safe).
        offsets, nbytes, compressed = [], [], False
        for page in tif.pages:
            if page.compression != 1 or len(page.dataoffsets) != 1:
                compressed = True
                break
            offsets.append(page.dataoffsets[0])
            nbytes.append(page.databytecounts[0])
    arr = zarr.create_array(
        store=zarr_path,
        shape=shape,
        chunks=(1, *chunk_yx),
        dtype=dtype,
        overwrite=True,
    )
    frame_shape = shape[1:]

    if compressed:
        log.info("Compressed/multistrip TIFF: falling back to single-threaded tifffile reads")
        with tifffile.TiffFile(tif_path) as tif:
            for t in tqdm(range(shape[0]), desc="tif->zarr"):
                arr[t] = tif.pages[t].asarray()
    else:

        def work(t):
            with open(tif_path, "rb") as fh:
                fh.seek(offsets[t])
                buf = fh.read(nbytes[t])
            arr[t] = np.frombuffer(buf, dtype=dtype).reshape(frame_shape)
            return t

        with ThreadPoolExecutor(n_workers) as ex:
            list(tqdm(ex.map(work, range(shape[0])), total=shape[0], desc="tif->zarr"))
    arr.attrs.update({"axes": "TYX", "source": tif_path})
    return arr


def temporal_median_background(raw: zarr.Array, n_frames: int) -> np.ndarray:
    T = raw.shape[0]
    idx = np.linspace(0, T - 1, n_frames).round().astype(int)
    log.info("Computing temporal median background from %d frames", len(idx))
    stack = np.stack([raw[t] for t in idx], axis=0)
    bg = np.median(stack, axis=0).astype(np.float32)
    # smooth lightly so that slow-moving cells that survive the median do not leave imprints
    bg = ndi.gaussian_filter(bg, 2.0)
    return bg


def bubble_mask_from_background(bg: np.ndarray, dilate_px: int) -> np.ndarray:
    """Static bright/dark artifacts (air bubbles, dust rings) show up in the background."""
    med = np.median(bg)
    mad = np.median(np.abs(bg - med)) + 1e-6
    hi = bg > med + 8 * mad
    lo = bg < med - 8 * mad
    cand = hi | lo
    cand = ndi.binary_closing(cand, iterations=3)
    lab = measure.label(cand)
    keep = np.zeros_like(cand)
    for r in measure.regionprops(lab):
        if r.area >= 1500:
            keep[lab == r.label] = True
    keep = ndi.binary_fill_holes(keep)
    if dilate_px > 0:
        keep = ndi.binary_dilation(keep, iterations=dilate_px)
    return keep


def write_flat(raw: zarr.Array, bg: np.ndarray, flat_path: str, chunk_yx, n_workers: int = 6) -> zarr.Array:
    """Divide by the background (multiplicative illumination model) and rescale to uint8."""
    gain = np.median(bg) / np.clip(bg, 1.0, None)
    flat = zarr.create_array(
        store=flat_path, shape=raw.shape, chunks=(1, *chunk_yx), dtype=np.uint8, overwrite=True
    )

    def work(t):
        f = raw[t].astype(np.float32) * gain
        flat[t] = np.clip(f, 0, 255).astype(np.uint8)
        return t

    with ThreadPoolExecutor(n_workers) as ex:
        list(tqdm(ex.map(work, range(raw.shape[0])), total=raw.shape[0], desc="flatten"))
    flat.attrs.update({"axes": "TYX", "note": "raw / background * median(background)"})
    return flat


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-flat", action="store_true")
    ap.add_argument("--skip-convert", action="store_true", help="raw.zarr already exists")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    P = cfg["paths"]
    ensure_dir(P["out_dir"])
    chunk_yx = tuple(cfg["ingest"]["chunk_yx"])

    if args.skip_convert:
        raw = zarr.open_array(P["raw_zarr"], mode="r")
    else:
        raw = tif_to_zarr(P["raw_tif"], P["raw_zarr"], chunk_yx)

    bg = temporal_median_background(raw, cfg["ingest"]["background_n_frames"])
    np.save(P["background"], bg)
    bubbles = bubble_mask_from_background(bg, cfg["segmentation"]["bubble_dilate_px"])
    np.save(P["bubble_mask"], bubbles)
    log.info("Bubble/artifact mask covers %.2f%% of the field", 100 * bubbles.mean())

    if not args.skip_flat:
        write_flat(raw, bg, P["flat_zarr"], chunk_yx)

    meta = {
        "n_frames": int(raw.shape[0]),
        "shape_yx": [int(raw.shape[1]), int(raw.shape[2])],
        "frame_interval_min": cfg["acquisition"]["frame_interval_min"],
        "pixel_size_um": cfg["acquisition"]["pixel_size_um"],
        "background_median": float(np.median(bg)),
        "bubble_fraction": float(bubbles.mean()),
    }
    with open(f"{P['out_dir']}/ingest_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    log.info("Done: %s", meta)


if __name__ == "__main__":
    main()
