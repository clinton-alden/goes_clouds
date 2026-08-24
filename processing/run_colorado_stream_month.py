#!/usr/bin/env python3
"""Run and validate one optimized full-resolution Colorado GOES month."""

from __future__ import annotations

import argparse
import calendar
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "processing"))
from stream_goes_rgb_mask_day import outputs_are_valid  # noqa: E402

BASE = Path("/glade/derecho/scratch/cdalden/colorado_goes")
BOUNDS = (-109.0, 37.0, -104.0, 41.0)
DOMAIN = "colorado_5x4"
RULES = ROOT / "analysis/output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("year", type=int)
    p.add_argument("month", type=int)
    p.add_argument("--base", type=Path, default=BASE)
    p.add_argument("--satellite", default="goes16")
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--workers", type=int, default=8)
    a = p.parse_args()
    if not ((2017, 1) <= (a.year, a.month) <= (2024, 12)):
        raise ValueError("Supported Colorado range is January 2017–December 2024")

    end = calendar.monthrange(a.year, a.month)[1]
    day_script = ROOT / "processing/stream_goes_rgb_mask_day.py"
    for day in range(1, end + 1):
        date = f"{a.year:04d}-{a.month:02d}-{day:02d}"
        subprocess.run([
            a.python, str(day_script), "--date", date, "--base", str(a.base),
            "--satellite", a.satellite, "--domain", DOMAIN,
            "--bounds", *map(str, BOUNDS),
            "--dem", str(a.base / "static/SRTMGL3_colorado_5x4.tif"),
            "--era-dir", str(a.base / "era5_land/t2m_hourly"),
            "--rules", str(RULES), "--python", a.python, "--workers", str(a.workers),
        ], check=True)

    rgb_dir = a.base / a.satellite / "rgb_composite_packed"
    mask_dir = a.base / a.satellite / "cloud_mask_tempbin_tree_packed"
    invalid = []
    for day in range(1, end + 1):
        token = f"{a.year}{a.month:02d}{day:02d}"
        rgb = rgb_dir / f"{a.satellite}_C02_C05_C13_rgb_{DOMAIN}_{token}.nc"
        mask = mask_dir / f"{a.satellite}_C02_C05_C13_rgb_{DOMAIN}_{token}_cloud_binary_tempbin_tree.nc"
        if not outputs_are_valid(rgb, mask, BOUNDS):
            invalid.append(token)
    if invalid:
        raise RuntimeError(f"Month validation failed for {len(invalid)} days: {','.join(invalid)}")

    era = a.base / "era5_land/t2m_hourly" / f"era5land_t2m_{DOMAIN}_{a.year}{a.month:02d}.nc"
    era.unlink(missing_ok=True)
    print(f"COMPLETE Colorado {a.year}-{a.month:02d}: {end} packed RGB/mask days; ERA5 removed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
