#!/usr/bin/env python
"""Compare Gothic observed SW with NLDAS-2 and NLDAS-3 point values for 2022.

This script is intentionally point-only: it reads the nearest grid-cell SWdown
value for each requested hour instead of loading full NLDAS grids. It writes a
cached comparison table, monthly bias metrics, a coverage report, and a figure
with model-vs-observed KDE panels plus monthly bias.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from build_nldas_point_eval import (
    DEFAULT_NLDAS2_DIR,
    DEFAULT_NLDAS3_DIR,
    DEFAULT_OBS_DIR,
    GOTHIC_LAT,
    GOTHIC_LON,
    build_times,
    ensure_nldas3_files,
    find_nldas2_file,
    find_nldas3_file,
    metsim_clear_sky_one_time,
    parse_hours,
    read_observed_sw,
    read_point_sw,
)


OUT_DIR = Path("analysis/output_17_nldas_point_eval")
NLDAS2_RE = re.compile(r"A(\d{8})\.(\d{2})00")
NLDAS3_RE = re.compile(r"(\d{8})_(\d{2})00")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2022)
    parser.add_argument("--hours", default="14-23,0", help="UTC hours to evaluate; use 0-23 for all hours")
    parser.add_argument("--lat", type=float, default=GOTHIC_LAT)
    parser.add_argument("--lon", type=float, default=GOTHIC_LON)
    parser.add_argument("--nldas2-dir", type=Path, default=DEFAULT_NLDAS2_DIR)
    parser.add_argument("--nldas3-dir", type=Path, default=DEFAULT_NLDAS3_DIR)
    parser.add_argument("--obs-dir", type=Path, default=DEFAULT_OBS_DIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--download-nldas3", action="store_true", help="Download missing NLDAS3 subset files")
    parser.add_argument("--force-nldas3", action="store_true", help="Force re-download of NLDAS3 subset files")
    parser.add_argument("--require-both-models", action="store_true", help="Keep only rows with NLDAS2 and NLDAS3")
    return parser.parse_args()


def available_nldas2_times(nldas2_dir: Path) -> set[pd.Timestamp]:
    out = set()
    for path in nldas2_dir.glob("*"):
        match = NLDAS2_RE.search(path.name)
        if match:
            out.add(pd.Timestamp(f"{match.group(1)} {match.group(2)}:00"))
    return out


def available_nldas3_times(nldas3_dir: Path) -> set[pd.Timestamp]:
    out = set()
    for path in nldas3_dir.glob("*.nc"):
        match = NLDAS3_RE.search(path.name)
        if match:
            out.add(pd.Timestamp(f"{match.group(1)} {match.group(2)}:00"))
    return out


def source_coverage(times: pd.DatetimeIndex, nldas2_dir: Path, nldas3_dir: Path) -> pd.DataFrame:
    rows = []
    nldas2_times = available_nldas2_times(nldas2_dir)
    nldas3_times = available_nldas3_times(nldas3_dir)
    df = pd.DataFrame({"time_utc": times})
    df["month"] = df["time_utc"].dt.to_period("M").astype(str)
    df["nldas2_available"] = df["time_utc"].isin(nldas2_times)
    df["nldas3_available"] = df["time_utc"].isin(nldas3_times)
    for month, group in df.groupby("month"):
        rows.append(
            {
                "month": month,
                "requested_hours": len(group),
                "nldas2_available_hours": int(group["nldas2_available"].sum()),
                "nldas3_available_hours": int(group["nldas3_available"].sum()),
                "both_available_hours": int((group["nldas2_available"] & group["nldas3_available"]).sum()),
            }
        )
    return pd.DataFrame(rows)


def read_model_rows(times: pd.DatetimeIndex, args: argparse.Namespace) -> pd.DataFrame:
    rows = []

    nldas3_files = {}
    if args.download_nldas3 or args.force_nldas3:
        # ensure_nldas3_files expects the downloader-related attributes used by
        # build_nldas_point_eval.py.
        args.nldas3_python = Path("/glade/work/cdalden/conda-envs/goes_downloading/bin/python")
        args.nldas3_downloader = Path(__file__).resolve().parent / "download_nldas3_sw_subset.py"
        args.point_box_deg = 0.05
        nldas3_files = ensure_nldas3_files(times, args)

    for dt in times:
        dt = pd.Timestamp(dt)
        try:
            rows.append(read_point_sw(find_nldas2_file(dt, args.nldas2_dir), dt, args.lat, args.lon, "nldas2"))
        except FileNotFoundError:
            pass

        try:
            nldas3_path = nldas3_files.get(dt) if nldas3_files else find_nldas3_file(dt, args.nldas3_dir)
            rows.append(read_point_sw(nldas3_path, dt, args.lat, args.lon, "nldas3"))
        except FileNotFoundError:
            pass

    if not rows:
        raise RuntimeError("No NLDAS rows were read. Check local NLDAS paths and cache coverage.")

    long = pd.DataFrame(rows)
    return (
        long.pivot_table(index="time_utc", columns="source", values="sw_wm2", aggfunc="first")
        .reset_index()
        .rename_axis(None, axis=1)
    )


def build_comparison(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    start = f"{args.year}-01-01"
    end = f"{args.year}-12-31"
    times = build_times(start, end, parse_hours(args.hours))
    coverage = source_coverage(times, args.nldas2_dir, args.nldas3_dir)

    compare = read_model_rows(times, args)
    obs = read_observed_sw(pd.DatetimeIndex(compare["time_utc"]), args.obs_dir)
    compare = compare.merge(obs, on="time_utc", how="left")
    compare["clear_sky_sw"] = [metsim_clear_sky_one_time(t, args.lat, args.lon) for t in compare["time_utc"]]
    compare["month"] = compare["time_utc"].dt.to_period("M").astype(str)

    if args.require_both_models:
        compare = compare.dropna(subset=["nldas2", "nldas3"])

    metric_rows = []
    for source in ("nldas2", "nldas3"):
        if source not in compare:
            continue
        for month, group in compare.groupby("month"):
            valid = group[[source, "gothic_obs_sw"]].dropna()
            if valid.empty:
                continue
            err = valid[source] - valid["gothic_obs_sw"]
            metric_rows.append(
                {
                    "source": source,
                    "month": month,
                    "n": len(valid),
                    "bias_wm2": err.mean(),
                    "mae_wm2": err.abs().mean(),
                    "rmse_wm2": np.sqrt((err**2).mean()),
                    "corr": valid[source].corr(valid["gothic_obs_sw"]),
                }
            )
    metrics = pd.DataFrame(metric_rows)
    return compare, metrics, coverage


def plot_density_and_bias(compare: pd.DataFrame, metrics: pd.DataFrame, out_png: Path, hours: str) -> None:
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
    source_info = [("nldas2", "NLDAS-2", "tab:blue"), ("nldas3", "NLDAS-3", "tab:green")]

    for ax, (source, label, color) in zip(axes[0], source_info):
        if source not in compare:
            ax.set_axis_off()
            continue
        valid = compare[["gothic_obs_sw", source]].dropna()
        if len(valid) >= 5:
            sns.kdeplot(
                data=valid,
                x="gothic_obs_sw",
                y=source,
                fill=True,
                cmap="Blues" if source == "nldas2" else "Greens",
                levels=8,
                thresh=0.02,
                ax=ax,
            )
            ax.scatter(valid["gothic_obs_sw"], valid[source], s=8, alpha=0.25, color=color, edgecolor="none")
        max_sw = np.nanmax(valid.to_numpy()) if not valid.empty else 1000.0
        ax.plot([0, max_sw], [0, max_sw], color="0.2", lw=1.5, ls="--")
        ax.set_xlim(0, max_sw * 1.05)
        ax.set_ylim(0, max_sw * 1.05)
        ax.set_xlabel("Gothic observed SW (W m$^{-2}$)")
        ax.set_ylabel(f"{label} SW (W m$^{-2}$)")
        ax.set_title(f"{label} vs Gothic obs, 2022 hours {hours}")

    for ax, (source, label, color) in zip(axes[1], source_info):
        metric = metrics[metrics["source"] == source].copy()
        if metric.empty:
            ax.set_axis_off()
            continue
        metric["month_start"] = pd.PeriodIndex(metric["month"], freq="M").to_timestamp()
        ax.axhline(0, color="0.25", lw=1)
        ax.bar(metric["month_start"], metric["bias_wm2"], width=25, color=color, alpha=0.75)
        for _, row in metric.iterrows():
            ax.text(row["month_start"], row["bias_wm2"], f"n={int(row['n'])}", ha="center", va="bottom", fontsize=8)
        ax.set_title(f"{label} monthly bias")
        ax.set_xlabel("Month")
        ax.set_ylabel("Model - observed SW (W m$^{-2}$)")
        ax.tick_params(axis="x", rotation=45)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    compare, metrics, coverage = build_comparison(args)

    tag = f"{args.year}_h{args.hours.replace(',', '-')}"
    compare_csv = args.out_dir / f"gothic_nldas_point_sw_compare_{tag}.csv"
    metrics_csv = args.out_dir / f"gothic_nldas_point_sw_monthly_bias_{tag}.csv"
    coverage_csv = args.out_dir / f"gothic_nldas_point_sw_coverage_{tag}.csv"
    figure_png = args.out_dir / f"gothic_nldas_point_sw_density_bias_{tag}.png"

    compare.to_csv(compare_csv, index=False)
    metrics.to_csv(metrics_csv, index=False)
    coverage.to_csv(coverage_csv, index=False)
    plot_density_and_bias(compare, metrics, figure_png, args.hours)

    print(f"Wrote {compare_csv}")
    print(f"Wrote {metrics_csv}")
    print(f"Wrote {coverage_csv}")
    print(f"Wrote {figure_png}")
    print("Rows:", len(compare))
    print("Non-null counts:")
    print(compare[["nldas2", "nldas3", "gothic_obs_sw"]].notna().sum())
    print("Monthly metrics:")
    print(metrics)
    print("Coverage:")
    print(coverage)


if __name__ == "__main__":
    main()
