"""Low-frame-rate label-free cell tracking and state-transition analysis.

Pipeline stages (each a CLI module):
    ingest      TIFF -> zarr, temporal-median background, flattened copy, bubble mask
    segment     Cellpose-SAM instance segmentation -> masks.zarr + detections.parquet
    track       Trackastra / HOCT / LapTrack linking -> GEFF + tracks.parquet
    validate    ground-truth-free validation suite
    states      per-frame cell states and events
    transition  transition matrices, dwell times, event rates, report
    motility    MSD, speed, persistence
    view        napari viewer and overlay movies
"""

from .common import load_config, Units

__all__ = ["load_config", "Units"]
