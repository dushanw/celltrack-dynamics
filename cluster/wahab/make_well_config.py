#!/usr/bin/env python
"""Write a per-well config for the cluster: configs/wells/<well>.yaml

    python cluster/wahab/make_well_config.py well_B3 [--template configs/wahab.yaml] [--data data/gdrive]

Paths in a config are resolved relative to the config file, so the generated file points at
../../data/gdrive/<well>_stack.tif and ../../outputs/<well>/... . acquisition.n_frames is read
from the TIFF so movies of different length work unchanged.
"""
import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("well", help="e.g. well_B3 (file <well>_stack.tif)")
    ap.add_argument("--template", default=ROOT / "configs" / "wahab.yaml")
    ap.add_argument("--data", default="data/gdrive", help="directory with the <well>_stack.tif files, relative to repo root")
    args = ap.parse_args()

    tif = ROOT / args.data / f"{args.well}_stack.tif"
    if not tif.exists():
        raise SystemExit(f"{tif} not found")
    import tifffile

    with tifffile.TiffFile(tif) as t:
        n_frames = len(t.pages)
        shape = list(t.pages[0].shape)

    text = Path(args.template).read_text()
    out_rel = f"../../outputs/{args.well}"
    text = re.sub(r"^(  raw_tif: ).*$", rf"\1../../{args.data}/{args.well}_stack.tif", text, flags=re.M)
    text = re.sub(r"^(  [a-z_]+: )\.\./outputs", rf"\g<1>{out_rel}", text, flags=re.M)
    text = re.sub(r"^(  n_frames: )\d+", rf"\g<1>{n_frames}", text, flags=re.M)
    text = re.sub(r"^(  shape_yx: ).*$", rf"\g<1>[{shape[0]}, {shape[1]}]", text, flags=re.M)
    text = f"# Auto-generated for {args.well} by cluster/wahab/make_well_config.py\n" + text

    out = ROOT / "configs" / "wells" / f"{args.well}.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(out, n_frames, shape)


if __name__ == "__main__":
    main()
