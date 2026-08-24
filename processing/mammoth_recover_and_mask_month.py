#!/usr/bin/env python3
"""Build all Mammoth RGB days for a month, then exact-tree masks."""

from __future__ import annotations

import argparse
import calendar
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("/glade/derecho/scratch/cdalden/mammoth")
PYTHON = "/glade/work/cdalden/conda-envs/goes_downloading/bin/python"
BOUNDS = (-119.53390, 37.13740, -118.53390, 38.13740)
DOMAIN = "mammoth_1deg"
ERA_DIR = BASE / "era5_land/t2m_hourly_mammoth_1deg"
RULES = ROOT / "analysis/output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv"


def segments(year: int, month: int) -> list[tuple[str, int, int]]:
    end = calendar.monthrange(year, month)[1]
    if (year, month) < (2023, 1):
        return [("goes17", 1, end)]
    if (year, month) == (2023, 1):
        return [("goes17", 1, 2), ("goes18", 3, end)]
    return [("goes18", 1, end)]


def validate_rgb(args, satellite: str, start_day: int, end_day: int) -> bool:
    rgb_dir = BASE / satellite / "rgb_composite"
    first = f"{args.year:04d}-{args.month:02d}-{start_day:02d}"
    last = f"{args.year:04d}-{args.month:02d}-{end_day:02d}"
    command = [
        args.python, str(ROOT / "processing/goes_rgb_domain/validate_rgb_dates.py"),
        str(rgb_dir), "--start", first, "--end", last, "--goes", satellite,
        "--domain", DOMAIN, "--bounds", *map(str, BOUNDS),
    ]
    return subprocess.run(command).returncode == 0


def recover_segment(args, satellite: str, start_day: int, end_day: int) -> None:
    if validate_rgb(args, satellite, start_day, end_day):
        print(f"RGB segment complete: {satellite} days {start_day}-{end_day}", flush=True)
        return
    for day in range(start_day, end_day + 1):
        subprocess.run([
            args.python, str(ROOT / "processing/recover_goes_rgb_day.py"),
            "--base", str(BASE), "--goes", satellite, "--domain", DOMAIN,
            "--year", str(args.year), "--month", str(args.month), "--day", str(day),
            "--bounds", *map(str, BOUNDS),
            "--dem", str(BASE / "static/SRTMGL3_mammoth_1deg.tif"),
            "--python", args.python, "--attempts", "3",
        ], check=True)
    if not validate_rgb(args, satellite, start_day, end_day):
        raise RuntimeError(f"RGB validation failed after recovery: {satellite} {args.year}-{args.month:02d}")


def mask_segment(args, satellite: str, start_day: int, end_day: int) -> None:
    subprocess.run([
        args.python, str(ROOT / "processing/run_tree_mask_month.py"),
        "--year", str(args.year), "--month", str(args.month),
        "--day-start", str(start_day), "--day-end", str(end_day),
        "--satellite", satellite, "--domain", DOMAIN, "--base", str(BASE),
        "--era-dir", str(ERA_DIR), "--rules", str(RULES),
        "--workers", str(args.workers), "--python", args.python,
    ], check=True)


def validate_all_masks(year: int, month: int, month_segments) -> None:
    expected_hours = {0, 1, *range(14, 24)}
    for satellite, start_day, end_day in month_segments:
        mask_dir = BASE / satellite / "cloud_mask_tempbin_tree"
        for day in range(start_day, end_day + 1):
            token = f"{year}{month:02d}{day:02d}"
            path = mask_dir / f"{satellite}_C02_C05_C13_rgb_{DOMAIN}_{token}_cloud_binary_tempbin_tree.nc"
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f"Cannot delete ERA5: missing mask {path}")
            with xr.open_dataset(path) as ds:
                times = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
                lat = np.asarray(ds.latitude, dtype=float)
                lon = np.asarray(ds.longitude, dtype=float)
                if len(times) != 144 or set(times.hour) != expected_hours:
                    raise RuntimeError(f"Cannot delete ERA5: invalid timestamps in {path}")
                if not np.all(np.diff(lat) > 0) or not np.all(np.diff(lon) > 0):
                    raise RuntimeError(f"Cannot delete ERA5: invalid orientation in {path}")
                if ds.attrs.get("tree_logic") != "AND within each leaf; OR across cloudy leaves":
                    raise RuntimeError(f"Cannot delete ERA5: wrong rule provenance in {path}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("year", type=int)
    p.add_argument("month", type=int)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--python", default=PYTHON)
    args = p.parse_args()
    if not ((2019, 3) <= (args.year, args.month) <= (2025, 12)):
        raise ValueError("Mammoth range is March 2019 through December 2025")

    month_segments = segments(args.year, args.month)
    for satellite, start_day, end_day in month_segments:
        recover_segment(args, satellite, start_day, end_day)
    for satellite, start_day, end_day in month_segments:
        mask_segment(args, satellite, start_day, end_day)

    validate_all_masks(args.year, args.month, month_segments)
    era = ERA_DIR / f"era5land_t2m_{DOMAIN}_{args.year}{args.month:02d}.nc"
    if era.exists():
        era.unlink()
        print(f"DELETED disposable ERA5 month after complete mask validation: {era}", flush=True)
    print(f"COMPLETE Mammoth {args.year}-{args.month:02d}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
