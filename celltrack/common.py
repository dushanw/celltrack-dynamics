"""Shared helpers: config loading, unit conversion, zarr access, logging."""

from __future__ import annotations

import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = ROOT / "configs" / "default.yaml"


def setup_logging(name: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )
    return logging.getLogger(name)


def load_config(path: str | os.PathLike | None = None) -> dict[str, Any]:
    cfg_path = Path(path) if path else DEFAULT_CONFIG
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)
    # resolve paths relative to the config file directory
    base = cfg_path.parent
    cfg["paths"] = {k: str((base / v).resolve()) for k, v in cfg["paths"].items()}
    cfg["_root"] = str(base)
    return cfg


def ensure_dir(p: str | os.PathLike) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


@dataclass
class Units:
    """Convert between pixel/frame units and physical units.

    pixel_size_um is None until the user supplies it; conversions then return
    NaN and label strings fall back to px / frame.
    """

    frame_interval_min: float
    pixel_size_um: float | None = None

    @classmethod
    def from_config(cls, cfg: dict) -> "Units":
        acq = cfg["acquisition"]
        return cls(float(acq["frame_interval_min"]), acq.get("pixel_size_um"))

    @property
    def has_pixel_size(self) -> bool:
        return self.pixel_size_um is not None

    def frames_to_hours(self, frames) -> np.ndarray:
        return np.asarray(frames, dtype=float) * self.frame_interval_min / 60.0

    def frames_to_minutes(self, frames) -> np.ndarray:
        return np.asarray(frames, dtype=float) * self.frame_interval_min

    def px_to_um(self, px) -> np.ndarray:
        if not self.has_pixel_size:
            return np.full_like(np.asarray(px, dtype=float), np.nan)
        return np.asarray(px, dtype=float) * self.pixel_size_um

    def speed_px_per_frame_to_um_per_min(self, v) -> np.ndarray:
        if not self.has_pixel_size:
            return np.full_like(np.asarray(v, dtype=float), np.nan)
        return np.asarray(v, dtype=float) * self.pixel_size_um / self.frame_interval_min

    def length_label(self) -> str:
        return "um" if self.has_pixel_size else "px"

    def speed_label(self) -> str:
        return "um/min" if self.has_pixel_size else "px/frame"


def open_zarr(path, mode="r"):
    import zarr

    return zarr.open(path, mode=mode)


def pick_device(requested: str) -> str:
    import torch

    if requested == "mps" and torch.backends.mps.is_available():
        return "mps"
    if requested == "cuda" and torch.cuda.is_available():
        return "cuda"
    return "cpu"
