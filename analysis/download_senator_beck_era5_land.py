#!/usr/bin/env python3
"""Download missing monthly ERA5-Land t2m files covering Senator Beck."""

from __future__ import annotations

import argparse
import calendar
from pathlib import Path

import cdsapi
import pandas as pd


OUT_DIR = Path("/glade/derecho/scratch/cdalden/tmp/colorado/era5_land/t2m_hourly")
SITE_LAT = 37.90
SITE_LON = -107.72


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", default="2023-07")
    parser.add_argument("--end", default="2024-07")
    args = parser.parse_args()
    periods = pd.period_range(args.start, args.end, freq="M")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    client = cdsapi.Client()
    for period in periods:
        target = OUT_DIR / (
            f"era5land_t2m_colorado_{period.year}{period.month:02d}.nc"
        )
        if target.exists() and target.stat().st_size > 0:
            print("Exists:", target, flush=True)
            continue
        n_days = calendar.monthrange(period.year, period.month)[1]
        request = {
            "variable": ["2m_temperature"],
            "year": [f"{period.year}"],
            "month": [f"{period.month:02d}"],
            "day": [f"{day:02d}" for day in range(1, n_days + 1)],
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "data_format": "netcdf",
            "download_format": "unarchived",
            # Small bounding box surrounding Senator Beck; ERA5-Land is 0.1°.
            "area": [SITE_LAT + 0.11, SITE_LON - 0.11,
                     SITE_LAT - 0.11, SITE_LON + 0.11],
        }
        temporary = target.with_suffix(".nc.part")
        print("Downloading:", target, flush=True)
        client.retrieve("reanalysis-era5-land", request, str(temporary))
        temporary.replace(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
