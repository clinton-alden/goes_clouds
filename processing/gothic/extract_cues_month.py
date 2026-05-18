#!/usr/bin/env python
"""Extract a month-range slice from a CUES netCDF file."""

import argparse
import os
import re
from calendar import monthrange
from datetime import datetime

import xarray as xr


def _coalesce_station_slug(ds):
    station_name = ds.attrs.get("station_name", "cues")
    station_name = str(station_name).strip().lower()
    slug = re.sub(r"[^a-z0-9]+", "_", station_name)
    slug = re.sub(r"^_|_$", "", slug)
    return slug or "cues"


def extract_month(source, year, month, output_dir, prefix, overwrite):
    start = datetime(year, month, 1)
    _, ndays = monthrange(year, month)
    end = datetime(year, month, ndays, 23, 59, 59)
    os.makedirs(output_dir, exist_ok=True)

    with xr.open_dataset(source) as ds:
        if "datetime" not in ds:
            raise KeyError("Expected CUES variable 'datetime' in source dataset.")

        # Keep only requested month.
        sliced = ds.sel(datetime=slice(start, end))

        if len(sliced["datetime"]) == 0 and not overwrite:
            print(
                f"No CUES rows for {year:04d}-{month:02d} in {source}; skipping"
            )
            return None

        station_slug = _coalesce_station_slug(sliced)
        month_s = f"{year:04d}{month:02d}"
        out_name = f"{prefix}_{station_slug}_{month_s}.nc"
        out_path = os.path.join(output_dir, out_name)

        if not overwrite and os.path.exists(out_path):
            print(f"Output exists, skipping: {out_path}")
            return out_path

        sliced.attrs["cues_subset"] = f"month={year:04d}-{month:02d}"
        sliced.attrs["source_path"] = source
        sliced.to_netcdf(out_path)
        return out_path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Subset a CUES file to a single YYYY-MM month and save one NetCDF."
        )
    )
    parser.add_argument("--source", required=True, help="Input CUES NetCDF file.")
    parser.add_argument("--year", required=True, type=int)
    parser.add_argument("--month", required=True, type=int)
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument(
        "--prefix",
        default="cues_month",
        help="Output filename prefix (default: cues_month).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output file if it already exists.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    out_path = extract_month(
        source=args.source,
        year=args.year,
        month=args.month,
        output_dir=args.output_dir,
        prefix=args.prefix,
        overwrite=args.overwrite,
    )
    if out_path:
        print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
