#!/usr/bin/env python3
"""Download/validate one ERA5-Land month and build East River cloud masks."""

from __future__ import annotations

import argparse
import calendar
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from processing.colorado.apply_tempbin_thresholds_colorado import build_cloud_mask

WEST, SOUTH, EAST, NORTH = -107.25, 38.70, -106.75, 39.25
ERA5_AREA = [39.45, -107.45, 38.50, -106.55]
HOURS = {0, 1, *range(14, 24)}


def era5_path(root: Path, year: int, month: int) -> Path:
    return root / f"era5land_t2m_east_river_{year}{month:02d}.nc"


def validate_era5(path: Path, year: int, month: int) -> None:
    expected = calendar.monthrange(year, month)[1] * 24
    with xr.open_dataset(path) as ds:
        if "t2m" not in ds:
            raise ValueError(f"{path}: missing t2m")
        time_name = "valid_time" if "valid_time" in ds.coords else "time"
        times = pd.DatetimeIndex(pd.to_datetime(ds[time_name].values))
        if len(times) != expected or times.nunique() != expected:
            raise ValueError(f"{path}: expected {expected} unique hourly times, found {len(times)}")
        if times.min() != pd.Timestamp(year, month, 1):
            raise ValueError(f"{path}: unexpected first time {times.min()}")
        lat = np.asarray(ds.latitude.values, dtype=float)
        lon = np.asarray(ds.longitude.values, dtype=float)
        if lat.min() > SOUTH + 0.051 or lat.max() < NORTH - 0.051:
            raise ValueError(f"{path}: latitude does not cover East River")
        if lon.min() > WEST + 0.051 or lon.max() < EAST - 0.051:
            raise ValueError(f"{path}: longitude does not cover East River")


def ensure_era5(root: Path, year: int, month: int) -> Path:
    path = era5_path(root, year, month)
    if path.is_file() and path.stat().st_size:
        validate_era5(path, year, month)
        return path
    import cdsapi

    root.mkdir(parents=True, exist_ok=True)
    ndays = calendar.monthrange(year, month)[1]
    request = {
        "variable": ["2m_temperature"],
        "year": [str(year)],
        "month": [f"{month:02d}"],
        "day": [f"{day:02d}" for day in range(1, ndays + 1)],
        "time": [f"{hour:02d}:00" for hour in range(24)],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": ERA5_AREA,
    }
    temporary = path.with_suffix(".nc.part")
    if temporary.exists():
        temporary.unlink()
    cdsapi.Client().retrieve("reanalysis-era5-land", request, str(temporary))
    validate_era5(temporary, year, month)
    temporary.replace(path)
    return path


def validate_mask(path: Path) -> None:
    with xr.open_dataset(path) as ds:
        if "cloud_binary" not in ds:
            raise ValueError(f"{path}: missing cloud_binary")
        times = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
        if len(times) != 144 or set(times.hour) != HOURS:
            raise ValueError(f"{path}: expected 144 times in 00,01,14-23Z")
        lat = np.asarray(ds.latitude.values, dtype=float)
        lon = np.asarray(ds.longitude.values, dtype=float)
        if not np.all(np.diff(lat) > 0) or not np.all(np.diff(lon) > 0):
            raise ValueError(f"{path}: coordinate orientation is invalid")
        if lat.min() > SOUTH + 0.02 or lat.max() < NORTH - 0.02:
            raise ValueError(f"{path}: latitude does not cover East River")
        if lon.min() > WEST + 0.02 or lon.max() < EAST - 0.02:
            raise ValueError(f"{path}: longitude does not cover East River")
        values = np.unique(ds.cloud_binary.values)
        if not set(values.tolist()).issubset({0, 1}):
            raise ValueError(f"{path}: non-binary values {values}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("year", type=int)
    parser.add_argument("month", type=int)
    parser.add_argument("--base", type=Path,
                        default=Path("/glade/derecho/scratch/cdalden/east_river_goes"))
    parser.add_argument("--include-diagnostics", action="store_true")
    args = parser.parse_args()
    rgb_dir = args.base / "goes16/rgb_composite"
    mask_dir = args.base / "goes16/cloud_mask_tempbin_10c"
    era5_dir = args.base / "era5_land/t2m_hourly"
    threshold_csv = REPO / "analysis/output_12_rgb_threshold_transfer/gothic_temp_bin_rgb_thresholds_10c.csv"
    ndays = calendar.monthrange(args.year, args.month)[1]

    # Skip partial months before making a CDS request. A later rerun of the
    # array will process them once all daily RGB inputs are present.
    rgb_paths = []
    missing = []
    for day in range(1, ndays + 1):
        date = f"{args.year}{args.month:02d}{day:02d}"
        path = rgb_dir / f"goes16_C02_C05_C13_rgb_east_river_{date}.nc"
        if not path.is_file() or not path.stat().st_size:
            missing.append(path.name)
        rgb_paths.append(path)
    if missing:
        print(
            f"SKIPPED {args.year}-{args.month:02d}: "
            f"{len(missing)}/{ndays} daily RGB inputs missing or empty",
            flush=True,
        )
        return 0

    era5 = ensure_era5(era5_dir, args.year, args.month)
    mask_dir.mkdir(parents=True, exist_ok=True)
    for rgb in rgb_paths:
        output = mask_dir / f"{rgb.stem}_cloud_binary_tempbin10c.nc"
        if output.is_file() and output.stat().st_size:
            validate_mask(output)
            print(f"VALID {output}", flush=True)
            continue
        temporary = output.with_suffix(".nc.part")
        build_cloud_mask(rgb, era5, threshold_csv, temporary, target_hours=HOURS,
                         include_diagnostics=args.include_diagnostics)
        validate_mask(temporary)
        temporary.replace(output)
        print(f"WROTE {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
