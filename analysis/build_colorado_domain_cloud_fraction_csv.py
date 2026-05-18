#!/usr/bin/env python3
"""Build Colorado domain cloud-fraction time series from ACM and RGB-mask files."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
import xarray as xr

DATE_RE = re.compile(r"(\d{8})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a CSV of ACM and RGB-mask cloud fraction over a fixed Colorado domain."
    )
    parser.add_argument(
        "--acm-dir",
        default="/glade/u/home/cdalden/scratch/colorado_acm/goes16/daily_nc_latlon",
        help="Directory containing daily Colorado ACM lat/lon NetCDF files.",
    )
    parser.add_argument(
        "--rgb-dir",
        default="/glade/u/home/cdalden/scratch/colorado_rgbmasks/cloud_mask_tempbin_10c",
        help="Directory containing daily Colorado RGB-mask NetCDF files.",
    )
    parser.add_argument(
        "--output-csv",
        default="/glade/u/home/cdalden/goes_work/analysis/output_08c/colorado_domain_cloud_fraction_14z_00z.csv",
        help="Output CSV path.",
    )
    parser.add_argument("--lat-min", type=float, default=38.90)
    parser.add_argument("--lat-max", type=float, default=39.03)
    parser.add_argument("--lon-min", type=float, default=-107.05)
    parser.add_argument("--lon-max", type=float, default=-106.92)
    parser.add_argument(
        "--start-hour",
        type=int,
        default=14,
        help="Keep timestamps with hour >= start-hour in UTC.",
    )
    parser.add_argument(
        "--end-hour",
        type=int,
        default=24,
        help="Keep timestamps with hour < end-hour in UTC.",
    )
    parser.add_argument(
        "--match-tolerance-min",
        type=float,
        default=3.0,
        help="Nearest-time merge tolerance in minutes.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=10,
        help="Print progress every N processed dates.",
    )
    return parser.parse_args()


def extract_date_key(path: Path) -> str | None:
    match = DATE_RE.search(path.name)
    return match.group(1) if match else None


def subset_latlon(ds: xr.Dataset, lat_min: float, lat_max: float, lon_min: float, lon_max: float) -> xr.Dataset:
    lat_vals = ds["latitude"].values
    lon_vals = ds["longitude"].values

    lat_slice = slice(lat_max, lat_min) if lat_vals[0] > lat_vals[-1] else slice(lat_min, lat_max)
    lon_slice = slice(lon_max, lon_min) if lon_vals[0] > lon_vals[-1] else slice(lon_min, lon_max)
    return ds.sel(latitude=lat_slice, longitude=lon_slice)


def cloud_fraction_df(da: xr.DataArray, cloudy_value: int, time_col: str, value_name: str) -> pd.DataFrame:
    cloud_mask = xr.where(da == cloudy_value, 1.0, 0.0)
    cf = cloud_mask.mean(dim=("latitude", "longitude"), skipna=True)
    df = cf.to_dataframe(name=value_name).reset_index()
    df = df[[time_col, value_name]].rename(columns={time_col: "t"})
    df["t"] = pd.to_datetime(df["t"])
    return df.sort_values("t").reset_index(drop=True)


def filter_utc_window(df: pd.DataFrame, start_hour: int, end_hour: int) -> pd.DataFrame:
    hour = df["t"].dt.hour
    keep = (hour >= start_hour) & (hour < end_hour)
    return df.loc[keep].copy()


def list_date_files(directory: Path, pattern: str) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for path in sorted(directory.glob(pattern)):
        date_key = extract_date_key(path)
        if date_key is not None:
            out[date_key] = path
    return out


def process_day(
    date_key: str,
    acm_path: Path,
    rgb_path: Path,
    *,
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
    start_hour: int,
    end_hour: int,
    tolerance: pd.Timedelta,
) -> pd.DataFrame:
    with xr.open_dataset(acm_path) as acm_ds, xr.open_dataset(rgb_path) as rgb_ds:
        acm_sub = subset_latlon(acm_ds, lat_min, lat_max, lon_min, lon_max)
        rgb_sub = subset_latlon(rgb_ds, lat_min, lat_max, lon_min, lon_max)

        acm_df = cloud_fraction_df(acm_sub["BCM"], cloudy_value=1, time_col="time", value_name="acm_cloud_frac")
        rgb_df = cloud_fraction_df(
            rgb_sub["cloud_binary"],
            cloudy_value=1,
            time_col="t",
            value_name="rgb_cloud_frac",
        )

        acm_df = filter_utc_window(acm_df, start_hour, end_hour)
        rgb_df = filter_utc_window(rgb_df, start_hour, end_hour)

        if acm_df.empty or rgb_df.empty:
            return pd.DataFrame()

        merged = pd.merge_asof(
            rgb_df,
            acm_df,
            on="t",
            direction="nearest",
            tolerance=tolerance,
        )

        merged["date"] = pd.to_datetime(date_key, format="%Y%m%d")
        merged["date_key"] = date_key
        merged["acm_n_pixels"] = acm_sub.sizes["latitude"] * acm_sub.sizes["longitude"]
        merged["rgb_n_pixels"] = rgb_sub.sizes["latitude"] * rgb_sub.sizes["longitude"]
        merged["acm_path"] = str(acm_path)
        merged["rgb_path"] = str(rgb_path)

        cols = [
            "date",
            "date_key",
            "t",
            "acm_cloud_frac",
            "rgb_cloud_frac",
            "acm_n_pixels",
            "rgb_n_pixels",
            "acm_path",
            "rgb_path",
        ]
        return merged[cols]


def main() -> None:
    args = parse_args()

    acm_dir = Path(args.acm_dir)
    rgb_dir = Path(args.rgb_dir)
    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    acm_files = list_date_files(acm_dir, "goes16_acm_colorado_*.nc")
    rgb_files = list_date_files(rgb_dir, "goes16_C02_C05_C13_rgb_colorado_*_cloud_binary_tempbin10c.nc")

    common_dates = sorted(set(acm_files) & set(rgb_files))
    if not common_dates:
        raise RuntimeError("No overlapping ACM and RGB-mask dates were found.")

    tolerance = pd.Timedelta(minutes=args.match_tolerance_min)
    frames: list[pd.DataFrame] = []

    for i, date_key in enumerate(common_dates, start=1):
        try:
            day_df = process_day(
                date_key,
                acm_files[date_key],
                rgb_files[date_key],
                lat_min=args.lat_min,
                lat_max=args.lat_max,
                lon_min=args.lon_min,
                lon_max=args.lon_max,
                start_hour=args.start_hour,
                end_hour=args.end_hour,
                tolerance=tolerance,
            )
        except Exception as exc:
            print(f"Skipping {date_key}: {exc}")
            continue

        if not day_df.empty:
            frames.append(day_df)

        if i % max(args.progress_every, 1) == 0 or i == len(common_dates):
            rows = sum(len(frame) for frame in frames)
            print(f"Processed {i}/{len(common_dates)} dates; rows collected={rows}")

    if not frames:
        raise RuntimeError("No cloud-fraction rows were produced.")

    out_df = pd.concat(frames, ignore_index=True).sort_values("t").reset_index(drop=True)
    out_df.to_csv(output_csv, index=False)

    print(f"Wrote {len(out_df)} rows to {output_csv}")


if __name__ == "__main__":
    main()
