#!/usr/bin/env python3
"""Download a monthly ERA5-Land t2m grid covering the Mammoth GOES domain."""

from __future__ import annotations

import argparse
import calendar
import os
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / f"era5land_t2m_mammoth_{args.year}{args.month:02d}.nc"
    if out.exists() and out.stat().st_size > 0:
        print(f"[skip] {out}")
        return 0
    import cdsapi

    request = {
        "variable": ["2m_temperature"],
        "year": [str(args.year)],
        "month": [f"{args.month:02d}"],
        "day": [f"{day:02d}" for day in range(1, calendar.monthrange(args.year, args.month)[1] + 1)],
        "time": [f"{hour:02d}:00" for hour in range(24)],
        "data_format": "netcdf",
        "download_format": "unarchived",
        # Padding guarantees a nearest ERA5-Land cell beyond every GOES edge.
        "area": [38.0, -119.6, 37.3, -118.5],
    }
    temporary = out.with_suffix(f".nc.part.{os.getpid()}")
    print(f"[download] {args.year}-{args.month:02d} -> {out}", flush=True)
    cdsapi.Client().retrieve("reanalysis-era5-land", request, str(temporary))
    os.replace(temporary, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
