#!/usr/bin/env python3
"""Download ERA5 hourly 2m temperature for the Gothic domain by month."""

from __future__ import annotations

import argparse
import calendar
from dataclasses import dataclass
from datetime import date
from pathlib import Path


GOTHIC_BOUNDS = {
    "lon_min": -107.08,
    "lat_min": 38.89,
    "lon_max": -106.90,
    "lat_max": 39.03,
}


@dataclass(frozen=True)
class YearMonth:
    year: int
    month: int


def iter_months(start_year: int, start_month: int, end_year: int, end_month: int):
    start = date(start_year, start_month, 1)
    end = date(end_year, end_month, 1)
    if end < start:
        raise ValueError("End year-month must be on/after start year-month")

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
        description="Download ERA5 hourly 2m temperature for Gothic domain"
    )
    parser.add_argument("--start-year", type=int, default=2021)
    parser.add_argument("--start-month", type=int, default=10)
    parser.add_argument("--end-year", type=int, default=2023)
    parser.add_argument("--end-month", type=int, default=6)
    parser.add_argument(
        "--output-dir",
        default="/glade/derecho/scratch/cdalden/gothic/era5/t2m_hourly",
        help="Directory for monthly ERA5 NetCDF files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload files even if they already exist",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    try:
        import cdsapi
    except ImportError as exc:
        raise SystemExit(
            "cdsapi is required. Install it in your environment, then rerun."
        ) from exc

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    area = [
        GOTHIC_BOUNDS["lat_max"],
        GOTHIC_BOUNDS["lon_min"],
        GOTHIC_BOUNDS["lat_min"],
        GOTHIC_BOUNDS["lon_max"],
    ]

    print("Downloading ERA5 2m temperature for Gothic")
    print(f"Bounds (N,W,S,E): {area}")
    print(
        f"Period: {args.start_year}-{args.start_month:02d} to "
        f"{args.end_year}-{args.end_month:02d}"
    )
    print(f"Output: {output_dir}")

    client = cdsapi.Client()
    hours = [f"{h:02d}:00" for h in range(24)]

    for ym in iter_months(
        args.start_year, args.start_month, args.end_year, args.end_month
    ):
        ndays = calendar.monthrange(ym.year, ym.month)[1]
        out_path = output_dir / f"era5_t2m_gothic_{ym.year}{ym.month:02d}.nc"

        if out_path.exists() and not args.overwrite:
            print(f"[skip] {out_path} already exists")
            continue

        request = {
            "product_type": "reanalysis",
            "variable": ["2m_temperature"],
            "year": f"{ym.year}",
            "month": f"{ym.month:02d}",
            "day": [f"{d:02d}" for d in range(1, ndays + 1)],
            "time": hours,
            "data_format": "netcdf",
            "download_format": "unarchived",
            "area": area,
        }

        print(f"[download] {ym.year}-{ym.month:02d} -> {out_path}")
        client.retrieve("reanalysis-era5-single-levels", request, str(out_path))

    print("ERA5 download complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
