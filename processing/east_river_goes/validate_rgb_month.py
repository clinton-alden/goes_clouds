#!/usr/bin/env python3
"""Validate completeness and geographic orientation of daily RGB files."""

import argparse
import calendar
from pathlib import Path

import numpy as np
import xarray as xr


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("rgb_dir", type=Path)
    parser.add_argument("year", type=int)
    parser.add_argument("month", type=int)
    parser.add_argument("--goes", default="goes16")
    parser.add_argument("--domain", default="east_river")
    parser.add_argument("--bounds", type=float, nargs=4, required=True,
                        metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    args = parser.parse_args()
    west, south, east, north = args.bounds
    failures = []
    days = calendar.monthrange(args.year, args.month)[1]

    for day in range(1, days + 1):
        date = f"{args.year}{args.month:02d}{day:02d}"
        path = args.rgb_dir / f"{args.goes}_C02_C05_C13_rgb_{args.domain}_{date}.nc"
        if not path.is_file() or path.stat().st_size == 0:
            failures.append(f"{date}: missing or empty")
            continue
        try:
            with xr.open_dataset(path) as ds:
                lat = np.asarray(ds.latitude.values)
                lon = np.asarray(ds.longitude.values)
                if lat.size < 2 or lon.size < 2 or ds.sizes.get("t", 0) == 0:
                    failures.append(f"{date}: empty coordinate/time dimension")
                    continue
                if not np.all(np.diff(lat) > 0):
                    failures.append(f"{date}: latitude is not strictly south-to-north")
                if not np.all(np.diff(lon) > 0):
                    failures.append(f"{date}: longitude is not strictly west-to-east")
                lat_tol = max(0.02, 2 * abs(float(np.median(np.diff(lat)))))
                lon_tol = max(0.02, 2 * abs(float(np.median(np.diff(lon)))))
                if abs(float(lat[0]) - south) > lat_tol or abs(float(lat[-1]) - north) > lat_tol:
                    failures.append(f"{date}: latitude coverage {lat[0]:.4f}..{lat[-1]:.4f}")
                if abs(float(lon[0]) - west) > lon_tol or abs(float(lon[-1]) - east) > lon_tol:
                    failures.append(f"{date}: longitude coverage {lon[0]:.4f}..{lon[-1]:.4f}")
                for name in ("red", "green", "blue"):
                    if name not in ds or ds[name].sizes.get("latitude", 0) == 0 or ds[name].sizes.get("longitude", 0) == 0:
                        failures.append(f"{date}: invalid {name} channel")
        except Exception as exc:
            failures.append(f"{date}: unreadable ({exc})")

    if failures:
        print("RGB validation failed:")
        print("\n".join(failures))
        return 1
    print(f"Validated {days} RGB days: complete, west-to-east and south-to-north")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
