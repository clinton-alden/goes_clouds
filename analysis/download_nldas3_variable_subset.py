#!/usr/bin/env python
"""Download small NLDAS-3 variable subsets from the public WaterInsight S3 bucket."""

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


def output_path(out_dir: Path, dt: pd.Timestamp, variables: list[str]) -> Path:
    tag = "_".join(variable.lower() for variable in variables)
    return out_dir / f"nldas3_{tag}_{dt:%Y%m%d_%H}00.nc"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", required=True, help="UTC start time, e.g. 2022-05-08T14:00")
    parser.add_argument("--end", required=True, help="UTC end time, inclusive, e.g. 2022-05-09T00:00")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--variables", default="LWdown", help="Comma-separated NLDAS3 variables")
    parser.add_argument("--point-lat", type=float, help="If set with --point-lon, write only the nearest grid cell")
    parser.add_argument("--point-lon", type=float, help="If set with --point-lat, write only the nearest grid cell")
    parser.add_argument("--lat-min", type=float)
    parser.add_argument("--lat-max", type=float)
    parser.add_argument("--lon-min", type=float)
    parser.add_argument("--lon-max", type=float)
    parser.add_argument("--force", action="store_true", help="Overwrite existing hourly subset files")
    return parser.parse_args()


def parse_variables(spec: str) -> list[str]:
    variables = [part.strip() for part in spec.split(",") if part.strip()]
    if not variables:
        raise ValueError("At least one variable is required")
    return variables


def subset_day(fs: s3fs.S3FileSystem, day: pd.Timestamp, times: pd.DatetimeIndex, args: argparse.Namespace) -> None:
    variables = parse_variables(args.variables)
    uri = S3_TEMPLATE.format(ym=day.strftime("%Y%m"), ymd=day.strftime("%Y%m%d"))
    print(f"Opening {uri}", flush=True)

    with fs.open(uri, "rb", block_size=16_777_216, cache_type="none") as remote_file:
        ds = xr.open_dataset(remote_file, engine="h5netcdf", phony_dims="sort")
        missing_vars = [variable for variable in variables if variable not in ds.data_vars]
        if missing_vars:
            raise KeyError(f"{uri} is missing requested variables {missing_vars}; has {list(ds.data_vars)}")

        if args.point_lat is not None or args.point_lon is not None:
            if args.point_lat is None or args.point_lon is None:
                raise ValueError("--point-lat and --point-lon must be provided together")
            point = ds[variables].sel(lat=args.point_lat, lon=args.point_lon, method="nearest")
            subset = point.expand_dims(lat=[float(point.lat)], lon=[float(point.lon)])
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
            subset = ds[variables].sel(
                lat=slice(args.lat_min, args.lat_max),
                lon=slice(args.lon_min, args.lon_max),
            )

        for dt in times:
            out = output_path(args.out_dir, dt, variables)
            if out.exists() and not args.force:
                print(f"Exists {out}", flush=True)
                continue

            hour = subset.sel(time=dt.to_datetime64())
            hour.load()
            hour.attrs.update(
                {
                    "source": uri,
                    "subset_note": f"{','.join(variables)} point subset from NLDAS-3 beta public S3 forcing file"
                    if args.point_lat is not None
                    else f"{','.join(variables)} box subset from NLDAS-3 beta public S3 forcing file",
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

    times = pd.date_range(pd.Timestamp(args.start), pd.Timestamp(args.end), freq="1h")
    fs = s3fs.S3FileSystem(anon=True, default_fill_cache=False, default_cache_type="none")

    for day in pd.DatetimeIndex(times.normalize().unique()):
        day_times = times[times.normalize() == day]
        subset_day(fs, day, day_times, args)


if __name__ == "__main__":
    main()
