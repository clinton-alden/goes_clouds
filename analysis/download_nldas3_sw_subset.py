#!/usr/bin/env python
"""Download small NLDAS-3 SWdown subsets from the public NASA WaterInsight S3 bucket.

The NLDAS-3 beta daily forcing files are roughly 13 GB each. This script uses
anonymous S3 range reads and writes one small hourly NetCDF per requested time,
containing only SWdown over either one nearest grid cell or a requested lat/lon box.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


SCRATCH_DEPS = Path("/glade/derecho/scratch/cdalden/tmp/python_deps_h5py")
if SCRATCH_DEPS.exists():
    sys.path.insert(0, str(SCRATCH_DEPS))

import pandas as pd
import s3fs
import xarray as xr


S3_TEMPLATE = "s3://nasa-waterinsight/NLDAS3/forcing/hourly/{ym}/NLDAS_FOR0010_H.A{ymd}.030.beta.nc"


def output_path(out_dir: Path, dt: pd.Timestamp) -> Path:
    return out_dir / f"nldas3_sw_{dt:%Y%m%d_%H}00.nc"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", required=True, help="UTC start time, e.g. 2022-05-08T14:00")
    parser.add_argument("--end", required=True, help="UTC end time, inclusive, e.g. 2022-05-09T00:00")
    parser.add_argument("--hours", help="Optional UTC hours to download for each date, e.g. '14-23,0'")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--point-lat", type=float, help="If set with --point-lon, write only the nearest grid cell")
    parser.add_argument("--point-lon", type=float, help="If set with --point-lat, write only the nearest grid cell")
    parser.add_argument("--lat-min", type=float)
    parser.add_argument("--lat-max", type=float)
    parser.add_argument("--lon-min", type=float)
    parser.add_argument("--lon-max", type=float)
    parser.add_argument("--force", action="store_true", help="Overwrite existing hourly subset files")
    return parser.parse_args()


def parse_hours(spec: str) -> list[int]:
    hours: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = [int(x) for x in part.split("-", 1)]
            if start <= end:
                hours.extend(range(start, end + 1))
            else:
                hours.extend(list(range(start, 24)) + list(range(0, end + 1)))
        else:
            hours.append(int(part))
    unique = []
    for hour in hours:
        if hour < 0 or hour > 23:
            raise ValueError(f"Invalid hour {hour}; expected 0-23")
        if hour not in unique:
            unique.append(hour)
    return unique


def build_times(start: str, end: str, hours: list[int]) -> pd.DatetimeIndex:
    times = []
    for day in pd.date_range(pd.Timestamp(start).normalize(), pd.Timestamp(end).normalize(), freq="D"):
        for hour in hours:
            offset_day = day + pd.Timedelta(days=1) if hour == 0 and any(h > 0 for h in hours) else day
            times.append(offset_day + pd.Timedelta(hours=hour))
    return pd.DatetimeIndex(times).drop_duplicates().sort_values()


def subset_day(fs: s3fs.S3FileSystem, day: pd.Timestamp, times: pd.DatetimeIndex, args: argparse.Namespace) -> None:
    uri = S3_TEMPLATE.format(ym=day.strftime("%Y%m"), ymd=day.strftime("%Y%m%d"))
    print(f"Opening {uri}", flush=True)

    with fs.open(uri, "rb", block_size=16_777_216, cache_type="none") as remote_file:
        ds = xr.open_dataset(remote_file, engine="h5netcdf", phony_dims="sort")
        if args.point_lat is not None or args.point_lon is not None:
            if args.point_lat is None or args.point_lon is None:
                raise ValueError("--point-lat and --point-lon must be provided together")
            point = ds[["SWdown"]].sel(lat=args.point_lat, lon=args.point_lon, method="nearest")
            # Keep lat/lon as length-1 dimensions so downstream code can use the
            # same nearest-point reader for both point and box subset files.
            sw = point.expand_dims(lat=[float(point.lat)], lon=[float(point.lon)])
        else:
            missing = [
                name
                for name in ("lat_min", "lat_max", "lon_min", "lon_max")
                if getattr(args, name) is None
            ]
            if missing:
                raise ValueError(
                    "Provide either --point-lat/--point-lon or all box bounds "
                    "(--lat-min --lat-max --lon-min --lon-max)"
                )
            sw = ds[["SWdown"]].sel(
                lat=slice(args.lat_min, args.lat_max),
                lon=slice(args.lon_min, args.lon_max),
            )

        for dt in times:
            out = output_path(args.out_dir, dt)
            if out.exists() and not args.force:
                print(f"Exists {out}", flush=True)
                continue

            hour = sw.sel(time=dt.to_datetime64())
            hour.load()
            hour.attrs.update(
                {
                    "source": uri,
                    "subset_note": "SWdown point subset from NLDAS-3 beta public S3 forcing file"
                    if args.point_lat is not None
                    else "SWdown box subset from NLDAS-3 beta public S3 forcing file",
                }
            )
            tmp_out = out.with_name(f".{out.name}.{os.getpid()}.tmp")
            if tmp_out.exists():
                tmp_out.unlink()
            try:
                hour.to_netcdf(tmp_out)
                tmp_out.replace(out)
            except Exception:
                if tmp_out.exists():
                    tmp_out.unlink()
                raise
            print(f"Wrote {out}", flush=True)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.hours:
        times = build_times(args.start, args.end, parse_hours(args.hours))
    else:
        times = pd.date_range(pd.Timestamp(args.start), pd.Timestamp(args.end), freq="1h")
    if not args.force:
        times = pd.DatetimeIndex([dt for dt in times if not output_path(args.out_dir, pd.Timestamp(dt)).exists()])
        if len(times) == 0:
            print("All requested NLDAS-3 SW subset files already exist.", flush=True)
            return
    fs = s3fs.S3FileSystem(anon=True, default_fill_cache=False, default_cache_type="none")

    for day in pd.DatetimeIndex(times.normalize().unique()):
        day_times = times[times.normalize() == day]
        subset_day(fs, day, day_times, args)


if __name__ == "__main__":
    main()
