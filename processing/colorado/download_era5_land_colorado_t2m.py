#!/usr/bin/env python3
"""Download monthly ERA5-Land 2m temperature for the Colorado RGB domain."""

from __future__ import annotations

import argparse
import calendar
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import xarray as xr

DEFAULT_ERA5_PADDING_DEG = 0.2


@dataclass(frozen=True)
class YearMonth:
    year: int
    month: int


def iter_months(start_year: int, start_month: int, end_year: int, end_month: int):
    start = date(start_year, start_month, 1)
    end = date(end_year, end_month, 1)
    if end < start:
        raise ValueError("End year-month must be on or after start year-month")

    current_year = start.year
    current_month = start.month
    while (current_year, current_month) <= (end.year, end.month):
        yield YearMonth(current_year, current_month)
        current_month += 1
        if current_month == 13:
            current_month = 1
            current_year += 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download monthly ERA5-Land 2m temperature for the Colorado RGB domain"
    )
    parser.add_argument(
        "--rgb-file",
        default=(
            "/glade/u/home/cdalden/scratch/colorado/goes16/rgb_composite/"
            "goes16_C02_C05_C13_rgb_colorado_20220325.nc"
        ),
        help="Reference RGB file used to infer the Colorado domain bounds",
    )
    parser.add_argument("--start-year", type=int, default=2021)
    parser.add_argument("--start-month", type=int, default=10)
    parser.add_argument("--end-year", type=int, default=2023)
    parser.add_argument("--end-month", type=int, default=6)
    parser.add_argument(
        "--padding-deg",
        type=float,
        default=DEFAULT_ERA5_PADDING_DEG,
        help="Padding added around the RGB domain before downloading ERA5-Land",
    )
    parser.add_argument(
        "--output-dir",
        default="/glade/derecho/scratch/cdalden/tmp/colorado/era5_land/t2m_hourly",
        help="Directory for monthly ERA5-Land NetCDF files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload files even if they already exist",
    )
    return parser


def derive_bounds(rgb_path: Path, padding_deg: float) -> list[float]:
    with xr.open_dataset(rgb_path) as ds:
        lon = np.asarray(ds["longitude"].values, dtype=float)
        lat = np.asarray(ds["latitude"].values, dtype=float)

    return [
        float(lat.max() + padding_deg),
        float(lon.min() - padding_deg),
        float(lat.min() - padding_deg),
        float(lon.max() + padding_deg),
    ]


def main() -> int:
    args = build_parser().parse_args()

    rgb_path = Path(args.rgb_file)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not rgb_path.exists():
        raise FileNotFoundError(f"Missing RGB file: {rgb_path}")

    try:
        import cdsapi
    except ImportError as exc:
        raise SystemExit("cdsapi is required. Install it, then rerun.") from exc

    area = derive_bounds(rgb_path, padding_deg=args.padding_deg)
    client = cdsapi.Client()
    hours = [f"{hour:02d}:00" for hour in range(24)]

    print("Downloading ERA5-Land 2m temperature for Colorado RGB domain")
    print(f"Reference RGB file: {rgb_path}")
    print(f"Bounds with padding (N,W,S,E): {area}")
    print(
        f"Period: {args.start_year}-{args.start_month:02d} to "
        f"{args.end_year}-{args.end_month:02d}"
    )
    print(f"Output directory: {output_dir}")

    for ym in iter_months(args.start_year, args.start_month, args.end_year, args.end_month):
        out_path = output_dir / f"era5land_t2m_colorado_{ym.year}{ym.month:02d}.nc"
        if out_path.exists() and not args.overwrite:
            print(f"[skip] {out_path} already exists")
            continue

        ndays = calendar.monthrange(ym.year, ym.month)[1]
        request = {
            "variable": ["2m_temperature"],
            "year": [f"{ym.year}"],
            "month": [f"{ym.month:02d}"],
            "day": [f"{day:02d}" for day in range(1, ndays + 1)],
            "time": hours,
            "data_format": "netcdf",
            "download_format": "unarchived",
            "area": area,
        }
        print(f"[download] {ym.year}-{ym.month:02d} -> {out_path}")
        client.retrieve("reanalysis-era5-land", request, str(out_path))

    print("ERA5-Land download complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
