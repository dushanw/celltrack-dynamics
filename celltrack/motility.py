"""Phase 5a: motility metrics.

Per-track: mean speed, net and total displacement, straightness, turning angles.
Ensemble: time-averaged and ensemble MSD with power-law exponent alpha and
short-lag diffusion coefficient (2D: MSD = 4 D tau), velocity autocorrelation and
persistence time. Split by density regime (time windows) and by cluster membership.

Usage:
    python -m celltrack.motility --name trackastra_ilp
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .common import Units, ensure_dir, load_config, setup_logging

log = setup_logging("motility")


def per_track_metrics(d: pd.DataFrame, units: Units, min_len: int = 5) -> pd.DataFrame:
    rows = []
    for tid, g in d.sort_values(["track_id", "frame"]).groupby("track_id"):
        if len(g) < min_len:
            continue
        xy = g[["x", "y"]].to_numpy()
        fr = g["frame"].to_numpy()
        steps = np.diff(xy, axis=0)
        dts = np.diff(fr)
        sl = np.hypot(*steps.T)
        total = sl.sum()
        net = np.hypot(*(xy[-1] - xy[0]))
        ang = np.arctan2(steps[:, 1], steps[:, 0])
        turn = np.angle(np.exp(1j * np.diff(ang)))
        rows.append(dict(
            track_id=tid, lineage_id=g.lineage_id.iloc[0], start_frame=fr[0], end_frame=fr[-1], n_frames=len(g),
            duration_frames=fr[-1] - fr[0],
            mean_speed_px_per_frame=float((sl / dts).mean()),
            net_displacement_px=float(net), total_path_px=float(total),
            straightness=float(net / total) if total > 0 else np.nan,
            mean_turning_angle_rad=float(np.abs(turn).mean()) if len(turn) else np.nan,
            directional_persistence=float(np.cos(turn).mean()) if len(turn) else np.nan,
            cluster_fraction=float(g["in_cluster"].mean()),
            mean_edge_conf=float(g["edge_conf_in"].mean()),
            start_hours=float(units.frames_to_hours(fr[0])),
        ))
    m = pd.DataFrame(rows)
    if len(m):
        m["mean_speed_um_per_min"] = units.speed_px_per_frame_to_um_per_min(m["mean_speed_px_per_frame"])
    return m


def ensemble_msd(d: pd.DataFrame, max_lag: int = 24, min_len: int = 5) -> pd.DataFrame:
    """Ensemble MSD over all (track, t, t+tau) pairs with contiguous frames."""
    sums = np.zeros(max_lag + 1)
    cnt = np.zeros(max_lag + 1)
    for _, g in d.sort_values(["track_id", "frame"]).groupby("track_id"):
        if len(g) < min_len:
            continue
        fr = g["frame"].to_numpy()
        xy = g[["x", "y"]].to_numpy()
        # only use tracks with no gaps for simplicity
        if (np.diff(fr) != 1).any():
            continue
        n = len(xy)
        for tau in range(1, min(max_lag, n - 1) + 1):
            dd = ((xy[tau:] - xy[:-tau]) ** 2).sum(1)
            sums[tau] += dd.sum()
            cnt[tau] += len(dd)
    tau = np.arange(1, max_lag + 1)
    msd = sums[1:] / np.clip(cnt[1:], 1, None)
    return pd.DataFrame(dict(lag_frames=tau, msd_px2=msd, n_pairs=cnt[1:]))


def msd_fit(m: pd.DataFrame, units: Units, fit_lags=(1, 10)) -> dict:
    sel = (m.lag_frames >= fit_lags[0]) & (m.lag_frames <= fit_lags[1]) & (m.n_pairs > 50)
    if sel.sum() < 3:
        return {}
    x, y = np.log(m.lag_frames[sel]), np.log(m.msd_px2[sel])
    alpha, logk = np.polyfit(x, y, 1)
    D_px2_per_frame = float(m.msd_px2.iloc[0] / 4)
    return dict(alpha=float(alpha), msd_lag1_px2=float(m.msd_px2.iloc[0]),
                D_px2_per_frame=D_px2_per_frame,
                D_um2_per_min=float(D_px2_per_frame * units.pixel_size_um ** 2 / units.frame_interval_min) if units.has_pixel_size else None,
                fit_lags=list(fit_lags))


def velocity_autocorrelation(d: pd.DataFrame, max_lag: int = 12, min_len: int = 6) -> pd.DataFrame:
    num = np.zeros(max_lag + 1)
    cnt = np.zeros(max_lag + 1)
    for _, g in d.sort_values(["track_id", "frame"]).groupby("track_id"):
        if len(g) < min_len:
            continue
        fr = g["frame"].to_numpy()
        if (np.diff(fr) != 1).any():
            continue
        v = np.diff(g[["x", "y"]].to_numpy(), axis=0)
        for tau in range(0, min(max_lag, len(v) - 1) + 1):
            dots = (v[tau:] * v[: len(v) - tau]).sum(1)
            num[tau] += dots.sum()
            cnt[tau] += len(dots)
    c = num / np.clip(cnt, 1, None)
    c = c / c[0] if c[0] > 0 else c
    return pd.DataFrame(dict(lag_frames=np.arange(max_lag + 1), vacf=c, n=cnt))


def persistence_time(vacf: pd.DataFrame, units: Units) -> dict:
    c = vacf.vacf.to_numpy()
    below = np.where(c < np.exp(-1))[0]
    if len(below) == 0:
        return dict(persistence_frames=None, note="VACF never drops below 1/e within max lag")
    i = below[0]
    if i == 0:
        return dict(persistence_frames=0.0)
    # linear interpolation between i-1 and i
    t = (i - 1) + (c[i - 1] - np.exp(-1)) / (c[i - 1] - c[i])
    return dict(persistence_frames=float(t), persistence_minutes=float(units.frames_to_minutes(t)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--name", default="hoct")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    units = Units.from_config(cfg)
    out = ensure_dir(Path(cfg["paths"]["analysis_dir"]) / f"motility_{args.name}")
    d = pd.read_parquet(Path(cfg["paths"]["analysis_dir"]) / f"states_{args.name}.parquet")
    conf_thr = cfg["tracking"]["high_confidence_prob"]
    win_h = cfg["states"]["time_window_hours"]

    per = per_track_metrics(d, units)
    per.to_parquet(out / "per_track.parquet", index=False)
    summary = dict(name=args.name, n_tracks_used=int(len(per)), units=units.speed_label())
    summary["per_track_median"] = per[["mean_speed_px_per_frame", "net_displacement_px", "total_path_px", "straightness",
                                       "directional_persistence", "n_frames"]].median().round(3).to_dict()

    # subsets: all / high-confidence tracks / isolated vs cluster / time windows
    d_hc = d[d.track_id.isin(per.loc[per.mean_edge_conf >= conf_thr, "track_id"])]
    iso_ids = per.loc[per.cluster_fraction < 0.5, "track_id"]
    clu_ids = per.loc[per.cluster_fraction >= 0.5, "track_id"]
    subsets = {"all": d, "high_conf_tracks": d_hc, "isolated_tracks": d[d.track_id.isin(iso_ids)],
               "cluster_tracks": d[d.track_id.isin(clu_ids)]}
    d["window"] = (d["hours"] // win_h).astype(int)
    for w, g in d.groupby("window"):
        subsets[f"window_{w*win_h:.0f}-{(w+1)*win_h:.0f}h"] = g
    msd_tables, vacf_tables = {}, {}
    summary["msd"] = {}
    summary["persistence"] = {}
    for name, g in subsets.items():
        if g.track_id.nunique() < 20:
            continue
        m = ensemble_msd(g)
        m.to_csv(out / f"msd_{name}.csv", index=False)
        msd_tables[name] = m
        summary["msd"][name] = msd_fit(m, units)
        v = velocity_autocorrelation(g)
        vacf_tables[name] = v
        summary["persistence"][name] = persistence_time(v, units)
        summary["msd"][name]["median_step_px"] = float(np.nanmedian(g["step_px"]))
        summary["msd"][name]["n_tracks"] = int(g.track_id.nunique())

    # speed vs time
    sp = d.groupby("window")["speed_px_per_frame"].agg(["median", "mean", "count"])
    sp.index = [f"{w*win_h:.0f}-{(w+1)*win_h:.0f}h" for w in sp.index]
    sp.to_csv(out / "speed_by_window.csv")
    summary["speed_by_window_px_per_frame"] = sp.round(3).to_dict()
    with open(out / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    _plots(out, msd_tables, vacf_tables, d, per, units, win_h)
    log.info("%s", json.dumps(summary, indent=1, default=float)[:3000])


def _plots(out, msd_tables, vacf_tables, d, per, units, win_h):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for name, m in msd_tables.items():
        ax.loglog(m.lag_frames * units.frame_interval_min, m.msd_px2, "o-", ms=3, label=name)
    ax.set_xlabel("lag (min)"); ax.set_ylabel("MSD (px^2)"); ax.legend(fontsize=7); ax.set_title("ensemble MSD")
    plt.tight_layout(); plt.savefig(out / "msd.png", dpi=90); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    for name, v in vacf_tables.items():
        if name in ("all", "high_conf_tracks", "isolated_tracks", "cluster_tracks"):
            ax.plot(v.lag_frames * units.frame_interval_min, v.vacf, "o-", ms=3, label=name)
    ax.axhline(np.exp(-1), ls="--", c="gray"); ax.set_xlabel("lag (min)"); ax.set_ylabel("velocity autocorrelation")
    ax.legend(fontsize=7); plt.tight_layout(); plt.savefig(out / "velocity_autocorrelation.png", dpi=90); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    hrs = d["hours"]
    bins = np.arange(0, hrs.max() + 1, 1.0)
    med = d.groupby(pd.cut(hrs, bins), observed=True)["speed_px_per_frame"].median()
    q1 = d.groupby(pd.cut(hrs, bins), observed=True)["speed_px_per_frame"].quantile(0.25)
    q3 = d.groupby(pd.cut(hrs, bins), observed=True)["speed_px_per_frame"].quantile(0.75)
    c = [i.mid for i in med.index]
    ax.plot(c, med.values, label="median step speed"); ax.fill_between(c, q1.values, q3.values, alpha=0.3, label="IQR")
    ax.set_xlabel("time (h)"); ax.set_ylabel("speed (px/frame)"); ax.legend()
    plt.tight_layout(); plt.savefig(out / "speed_by_time.png", dpi=90); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist(per["mean_turning_angle_rad"].dropna(), bins=36, range=(0, np.pi))
    ax.set_xlabel("mean |turning angle| per track (rad)"); ax.set_ylabel("tracks")
    plt.tight_layout(); plt.savefig(out / "turning_angles.png", dpi=90); plt.close(fig)


if __name__ == "__main__":
    main()
