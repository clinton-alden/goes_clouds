#!/usr/bin/env python3
"""Combine clipped ACM files into one daily NetCDF with a time dimension."""

from __future__ import annotations

import argparse
import calendar
from datetime import datetime
from pathlib import Path
import re

import numpy as np
import xarray as xr


SCAN_RE = re.compile(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write daily NetCDF files from clipped GOES ACM inputs."
    )
    parser.add_argument("--base-dir", required=True)
    parser.add_argument("--year", required=True, type=int)
    parser.add_argument("--month", required=True, type=int)
    parser.add_argument("--goes", default="goes16")
    parser.add_argument("--domain", default="colorado")
    parser.add_argument("--product", default="ABI-L2-ACMC")
    parser.add_argument("--input-subdir", default="clipped")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--mask-var", default="AUTO")
    return parser.parse_args()


def resolve_mask_var(ds: xr.Dataset, requested: str) -> str:
    if requested != "AUTO":
        if requested not in ds:
            raise KeyError(f"{requested} missing from dataset")
        return requested
    for candidate in ("BCM", "ACM"):
        if candidate in ds:
            return candidate
    raise KeyError("Neither BCM nor ACM found in dataset")


def infer_time(path: Path) -> np.datetime64:
    match = SCAN_RE.search(path.name)
    if match is None:
        raise ValueError(f"Could not parse scan time from {path.name}")
    year, jday, hour, minute, second = match.groups()
    dt = datetime.strptime(
        f"{year} {jday} {hour} {minute} {second}",
        "%Y %j %H %M %S",
    )
    return np.datetime64(dt, "s")


def find_clipped_files(base_dir: Path, goes: str, year: int, month: int, day: int, product: str, input_subdir: str) -> list[Path]:
    day_dir = base_dir / goes / str(year) / str(month) / str(day) / product / input_subdir
    if not day_dir.exists():
        return []
    return sorted(day_dir.rglob("*.nc"))


def load_one(path: Path, mask_var: str) -> xr.Dataset:
    ds = xr.open_dataset(path)
    mask_var = resolve_mask_var(ds, mask_var)

    if "time" in ds.coords:
        time_value = ds["time"].values
    elif "t" in ds:
        time_value = ds["t"].values
    else:
        time_value = infer_time(path)

    time_value = np.asarray(time_value).reshape(-1)[0]
    keep_names = [mask_var]
    if "DQF" in ds:
        keep_names.append("DQF")

    out = ds[keep_names].assign_coords(
        latitude=ds["latitude"],
        longitude=ds["longitude"],
    ).expand_dims(time=[time_value]).load()
    ds.close()
    return out


def output_path(base_dir: Path, goes: str, domain: str, year: int, month: int, day: int, output_dir: str | None) -> Path:
    if output_dir is None:
        out_dir = base_dir / goes / "daily_nc_latlon"
    else:
        out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"{goes}_acm_{domain}_{year:04d}{month:02d}{day:02d}.nc"


def write_day(base_dir: Path, goes: str, domain: str, year: int, month: int, day: int, product: str, input_subdir: str, output_dir: str | None, mask_var: str) -> bool:
    files = find_clipped_files(base_dir, goes, year, month, day, product, input_subdir)
    if not files:
        print(f"{year:04d}-{month:02d}-{day:02d}: no clipped files found")
        return False

    datasets = []
    for path in files:
        try:
            datasets.append(load_one(path, mask_var))
        except Exception as exc:
            print(f"Skipping {path}: {exc}")

    if not datasets:
        print(f"{year:04d}-{month:02d}-{day:02d}: no valid clipped files")
        return False

    combined = xr.concat(datasets, dim="time").sortby("time")
    combined.attrs["source_product"] = "GOES ABI-L2-ACMC clipped in native fixed-grid coordinates"
    combined.attrs["mask_variable"] = next(iter(combined.data_vars))

    enc = {}
    for name in combined.data_vars:
        enc[name] = {"zlib": True, "complevel": 4}
    out_path = output_path(base_dir, goes, domain, year, month, day, output_dir)
    combined.to_netcdf(out_path, encoding=enc)
    combined.close()
    for ds in datasets:
        ds.close()
    print(f"Wrote {out_path}")
    return True


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir)
    num_days = calendar.monthrange(args.year, args.month)[1]
    wrote = 0
    for day in range(1, num_days + 1):
        wrote += int(
            write_day(
                base_dir,
                args.goes,
                args.domain,
                args.year,
                args.month,
                day,
                args.product,
                args.input_subdir,
                args.output_dir,
                args.mask_var,
            )
        )
    print(f"Finished {args.year:04d}-{args.month:02d}: wrote {wrote}/{num_days} daily files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
