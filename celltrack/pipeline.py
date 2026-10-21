"""Run every stage after segmentation in order, logging timings.

    python -m celltrack.pipeline [--skip features track ...] [--primary hoct] [--compare-frames 113]

Stages: features -> tracking (primary linker on the full movie in cached chunks; Trackastra
        ilp_nodiv / greedy and LapTrack on a comparison window) -> validation -> states
        -> transition matrices -> motility -> overlay movies.
Every stage is resumable/idempotent, so the pipeline can simply be re-run.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from .common import load_config, setup_logging

log = setup_logging("pipeline")


_CONFIG: list[str] = []  # ["--config", path] forwarded to every stage when a non-default config is used


def run(mod: str, *args: str, check=True) -> int:
    cmd = [sys.executable, "-m", f"celltrack.{mod}", *_CONFIG, *args]
    log.info(">>> %s", " ".join(cmd))
    t0 = time.time()
    r = subprocess.run(cmd)
    log.info("<<< %s finished in %.1f min (exit %d)", mod, (time.time() - t0) / 60, r.returncode)
    if check and r.returncode != 0:
        raise SystemExit(f"stage {mod} {args} failed")
    return r.returncode


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip", nargs="*", default=[], help="stages to skip: features track compare validate states transition motility view")
    ap.add_argument("--primary", default=None, help="linker for the full movie (default: tracking.primary in config)")
    ap.add_argument("--compare-frames", type=int, default=113, help="frames of the linker-comparison window (0 = off)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    if args.config:
        _CONFIG[:] = ["--config", str(Path(args.config).resolve())]
    tcfg = cfg["tracking"]
    primary = args.primary or tcfg.get("primary", "hoct")
    chunk = str(tcfg.get("chunk_frames", 24))
    T = cfg["acquisition"]["n_frames"]
    done_path = Path(cfg["paths"]["masks_zarr"]).parent / "masks_done.json"
    n_done = len(json.load(open(done_path))) if done_path.exists() else 0
    log.info("%d / %d frames segmented; primary linker = %s", n_done, T, primary)
    stop = str(n_done)
    prim_name = primary if n_done >= T else f"{primary}_f0-{n_done}"

    if "features" not in args.skip:
        run("segment", "features")
    if "track" not in args.skip:
        if primary == "hoct":
            run("track", "hoct", "--stop", stop, "--chunk", chunk, "--name", "hoct")
            # time-reversed run -> reverse-consistency flag for every edge (high-confidence matrices)
            run("track", "hoct", "--reverse", "--stop", stop, "--chunk", chunk, "--name", "hoct", check=False)
        else:
            run("track", "trackastra", "--mode", tcfg["trackastra"].get("mode", "ilp_nodiv"), "--stop", stop,
                "--chunk", str(tcfg["trackastra"].get("ilp_chunk", 12)), "--name", primary)
    if "compare" not in args.skip and args.compare_frames:
        cs = str(min(args.compare_frames, n_done))
        run("track", "trackastra", "--mode", "ilp_nodiv", "--stop", cs, "--chunk", str(tcfg["trackastra"].get("ilp_chunk", 12)),
            "--name", "ta_nodiv", check=False)
        run("track", "trackastra", "--mode", "greedy", "--stop", cs, "--name", "ta_greedy", check=False)
        run("track", "laptrack", "--stop", cs, "--name", "laptrack", check=False)
        if primary != "hoct":
            run("track", "hoct", "--stop", cs, "--chunk", chunk, "--name", "hoct", check=False)
    if "validate" not in args.skip:
        run("validate", "all", "--name", prim_name, check=False)
    if "states" not in args.skip:
        run("states", "--name", prim_name)
        for n in (f"ta_nodiv_f0-{min(args.compare_frames, n_done)}", f"laptrack_f0-{min(args.compare_frames, n_done)}"):
            if (Path(cfg["paths"]["tracks_dir"]) / f"{n}_tracks.parquet").exists():
                run("states", "--name", n, check=False)
    if "transition" not in args.skip:
        run("transition", "--name", prim_name)
    if "motility" not in args.skip:
        run("motility", "--name", prim_name)
        run("transition", "--name", prim_name)  # re-render the report with motility included
    if "view" not in args.skip:
        for a, b in [(0, 96), (n_done // 2 - 48, n_done // 2 + 48), (n_done - 96, n_done)]:
            run("view", "movie", "--name", prim_name, "--crop", "1200", "1200", "768", "--start", str(max(0, a)), "--stop", str(b), check=False)
        run("view", "movie", "--name", prim_name, "--crop", "0", "0", "3072", "--fps", "12", check=False)  # full-field overview
    log.info("Pipeline complete. Report: %s", Path(cfg["paths"]["analysis_dir"]) / f"transition_{prim_name}" / "report.html")


if __name__ == "__main__":
    main()
