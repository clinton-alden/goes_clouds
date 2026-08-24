#!/usr/bin/env python3
"""Build one RGB temperature-bin mask, downloading a validated ERA5-Land month."""

import argparse
import calendar
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from processing.colorado.apply_tempbin_thresholds_colorado import build_cloud_mask, derive_bounds


def validate_era(path: Path, rgb: Path) -> None:
    with xr.open_dataset(path) as e, xr.open_dataset(rgb) as r:
        if "t2m" not in e:
            raise ValueError(f"{path}: missing t2m")
        lat, lon = np.asarray(e.latitude), np.asarray(e.longitude)
        if lat.min() > float(r.latitude.min()) + 0.051 or lat.max() < float(r.latitude.max()) - 0.051:
            raise ValueError(f"{path}: latitude coverage is insufficient")
        if lon.min() > float(r.longitude.min()) + 0.051 or lon.max() < float(r.longitude.max()) - 0.051:
            raise ValueError(f"{path}: longitude coverage is insufficient")


def ensure_era(rgb: Path, era_dir: Path) -> Path:
    token = rgb.stem.rsplit("_", 1)[-1]
    year, month = int(token[:4]), int(token[4:6])
    try:
        domain = rgb.stem.split("_rgb_", 1)[1].rsplit("_", 1)[0]
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Could not derive domain from RGB filename: {rgb.name}") from exc
    path = era_dir / f"era5land_t2m_{domain}_{year}{month:02d}.nc"
    if path.is_file() and path.stat().st_size:
        validate_era(path, rgb)
        return path
    import cdsapi
    era_dir.mkdir(parents=True, exist_ok=True)
    with xr.open_dataset(rgb) as ds:
        area = derive_bounds(ds, padding_deg=0.2)
    request = {
        "variable": ["2m_temperature"], "year": [str(year)], "month": [f"{month:02d}"],
        "day": [f"{d:02d}" for d in range(1, calendar.monthrange(year, month)[1] + 1)],
        "time": [f"{h:02d}:00" for h in range(24)], "data_format": "netcdf",
        "download_format": "unarchived", "area": area,
    }
    temporary = path.with_suffix(".nc.part")
    cdsapi.Client().retrieve("reanalysis-era5-land", request, str(temporary))
    validate_era(temporary, rgb)
    temporary.replace(path)
    return path


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("rgb", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--era-dir", type=Path, required=True)
    p.add_argument("--thresholds", type=Path, default=ROOT / "analysis/output_12_rgb_threshold_transfer/gothic_temp_bin_rgb_thresholds_10c.csv")
    a = p.parse_args()
    era = ensure_era(a.rgb, a.era_dir)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = a.output.with_suffix(".nc.part")
    build_cloud_mask(a.rgb, era, a.thresholds, temporary,
                     target_hours={0, 1, *range(14, 24)}, include_diagnostics=False)
    temporary.replace(a.output)
    print(f"Wrote mask: {a.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
