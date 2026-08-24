#!/usr/bin/env python3
"""Validate daily RGB files for an inclusive date slice and domain bounds."""

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("rgb_dir", type=Path)
    p.add_argument("--start", type=dt.date.fromisoformat, required=True)
    p.add_argument("--end", type=dt.date.fromisoformat, required=True)
    p.add_argument("--goes", required=True)
    p.add_argument("--domain", required=True)
    p.add_argument("--bounds", type=float, nargs=4, required=True,
                   metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    a = p.parse_args()
    west, south, east, north = a.bounds
    failures = []
    day = a.start
    count = 0
    while day <= a.end:
        count += 1
        token = day.strftime("%Y%m%d")
        path = a.rgb_dir / f"{a.goes}_C02_C05_C13_rgb_{a.domain}_{token}.nc"
        if not path.is_file() or not path.stat().st_size:
            failures.append(f"{token}: missing or empty")
            day += dt.timedelta(days=1)
            continue
        try:
            with xr.open_dataset(path) as ds:
                lat = np.asarray(ds.latitude.values, dtype=float)
                lon = np.asarray(ds.longitude.values, dtype=float)
                times = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
                requested_hours = {0, 1, *range(14, 24)}
                requested_count = int(np.isin(times.hour, sorted(requested_hours)).sum())
                if requested_count != 144:
                    failures.append(
                        f"{token}: expected 144 requested-hour times, found {requested_count}"
                    )
                if not times.is_unique:
                    failures.append(f"{token}: duplicate timestamps")
                if lat.size < 2 or not np.all(np.diff(lat) > 0):
                    failures.append(f"{token}: latitude orientation invalid")
                if lon.size < 2 or not np.all(np.diff(lon) > 0):
                    failures.append(f"{token}: longitude orientation invalid")
                if lat.size > 1 and (lat.min() > south + 0.02 or lat.max() < north - 0.02):
                    failures.append(f"{token}: latitude coverage invalid")
                if lon.size > 1 and (lon.min() > west + 0.02 or lon.max() < east - 0.02):
                    failures.append(f"{token}: longitude coverage invalid")
                for name in ("red", "green", "blue"):
                    if name not in ds:
                        failures.append(f"{token}: missing {name}")
        except Exception as exc:
            failures.append(f"{token}: unreadable ({exc})")
        day += dt.timedelta(days=1)
    if failures:
        print("RGB validation failed:\n" + "\n".join(failures))
        return 1
    print(f"Validated {count} RGB days from {a.start} through {a.end}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
