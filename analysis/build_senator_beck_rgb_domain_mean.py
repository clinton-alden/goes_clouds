#!/usr/bin/env python3
"""Build Senator Beck domain-mean RGB time series from daily GOES RGB NetCDF files.

Run this on a compute node (for example via PBS), then load the output CSV in
Notebook 11c instead of processing RGB files interactively on a login node.
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def goes16_xy_to_latlon(x_1d: np.ndarray, y_1d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert GOES fixed-grid x/y coordinates to lat/lon."""
    lon_origin = -75.0
    perspective_point_height = 35786023.0
    semi_major_axis = 6378137.0
    semi_minor_axis = 6356752.31414

    h = perspective_point_height + semi_major_axis
    r_eq = semi_major_axis
    r_pol = semi_minor_axis

    x2d, y2d = np.meshgrid(x_1d, y_1d)
    lambda_0 = np.deg2rad(lon_origin)

    a = np.sin(x2d) ** 2 + np.cos(x2d) ** 2 * (
        np.cos(y2d) ** 2 + (r_eq**2 / r_pol**2) * np.sin(y2d) ** 2
    )
    b = -2.0 * h * np.cos(x2d) * np.cos(y2d)
    c = h**2 - r_eq**2

    rs = (-b - np.sqrt(b**2 - 4.0 * a * c)) / (2.0 * a)
    sx = rs * np.cos(x2d) * np.cos(y2d)
    sy = -rs * np.sin(x2d)
    sz = rs * np.cos(x2d) * np.sin(y2d)

    lat = np.rad2deg(np.arctan((r_eq**2 / r_pol**2) * (sz / np.sqrt((h - sx) ** 2 + sy**2))))
    lon = np.rad2deg(lambda_0 - np.arctan(sy / (h - sx)))
    return lat, lon


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build domain-mean RGB CSV for Senator Beck.")
    parser.add_argument(
        "--rgb-dir",
        default="/glade/derecho/scratch/cdalden/senator_beck/goes16/rgb_composite",
        help="Directory containing daily RGB NetCDF files.",
    )
    parser.add_argument(
        "--pattern",
        default="goes16_C02_C05_C13_rgb_senator_beck_*.nc",
        help="Glob pattern for RGB files within --rgb-dir.",
    )
    parser.add_argument(
        "--output-csv",
        default="/glade/u/home/cdalden/goes_work/analysis/output_11c_senator_beck/senator_beck_rgb_domain_mean_all.csv",
        help="Output CSV path.",
    )
    parser.add_argument("--lat-min", type=float, default=37.86)
    parser.add_argument("--lat-max", type=float, default=37.96)
    parser.add_argument("--lon-min", type=float, default=-107.80)
    parser.add_argument("--lon-max", type=float, default=-107.64)
    parser.add_argument(
        "--start-time",
        default=None,
        help="Optional inclusive UTC start time filter, e.g. 2022-01-01.",
    )
    parser.add_argument(
        "--end-time",
        default=None,
        help="Optional inclusive UTC end time filter, e.g. 2024-12-31.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="Print progress every N files.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    rgb_files = sorted(glob.glob(os.path.join(args.rgb_dir, args.pattern)))
    if not rgb_files:
        raise FileNotFoundError(f"No RGB files found in {args.rgb_dir} with pattern {args.pattern}")

    print(f"RGB files found: {len(rgb_files)}")
    print(f"First file: {rgb_files[0]}")
    print(f"Last file : {rgb_files[-1]}")

    with xr.open_dataset(rgb_files[0]) as ds0:
        lat2d, lon2d = goes16_xy_to_latlon(ds0["x"].values, ds0["y"].values)

    domain_mask = (
        (lat2d >= args.lat_min)
        & (lat2d <= args.lat_max)
        & (lon2d >= args.lon_min)
        & (lon2d <= args.lon_max)
    )

    pixel_count = int(domain_mask.sum())
    print(f"Pixels in requested domain: {pixel_count}")
    if pixel_count == 0:
        raise ValueError("No pixels found in requested lat/lon bounds.")

    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp_out.exists():
        tmp_out.unlink()

    wrote_header = False
    total_rows = 0

    for i, file_path in enumerate(rgb_files, start=1):
        with xr.open_dataset(file_path) as ds:
            red_mean = ds["red"].where(domain_mask).mean(dim=("y", "x"), skipna=True)
            green_mean = ds["green"].where(domain_mask).mean(dim=("y", "x"), skipna=True)
            blue_mean = ds["blue"].where(domain_mask).mean(dim=("y", "x"), skipna=True)

            chunk = pd.DataFrame(
                {
                    "time": pd.to_datetime(ds["t"].values),
                    "red": red_mean.values,
                    "green": green_mean.values,
                    "blue": blue_mean.values,
                }
            )

        chunk = chunk.dropna(subset=["red", "green", "blue"])
        if args.start_time is not None:
            chunk = chunk[chunk["time"] >= pd.Timestamp(args.start_time)]
        if args.end_time is not None:
            chunk = chunk[chunk["time"] <= pd.Timestamp(args.end_time)]

        if not chunk.empty:
            chunk = chunk.sort_values("time")
            chunk.to_csv(tmp_out, mode="a", header=not wrote_header, index=False)
            wrote_header = True
            total_rows += len(chunk)

        if i % max(1, args.progress_every) == 0 or i == len(rgb_files):
            print(f"Processed {i}/{len(rgb_files)} RGB files")

    if not wrote_header:
        raise ValueError("No valid RGB rows found after filtering.")

    tmp_out.replace(out_path)

    print(f"Saved: {out_path}")
    print(f"Rows: {total_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
