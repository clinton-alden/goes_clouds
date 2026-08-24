#!/usr/bin/env python3
"""Save a byte-range spatial subset of one day of GOES ACM from public S3."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import s3fs
import xarray as xr
from pyproj import CRS, Transformer


HOURS = (0, 1, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)


def fixed_grid_slice(ds, bounds):
    lon_min, lat_min, lon_max, lat_max = bounds
    attrs = ds.goes_imager_projection.attrs
    height = float(attrs["perspective_point_height"])
    crs = CRS.from_proj4(
        f"+proj=geos +h={height} +lon_0={attrs['longitude_of_projection_origin']} "
        f"+sweep={attrs['sweep_angle_axis']} +a={attrs['semi_major_axis']} "
        f"+b={attrs['semi_minor_axis']} +units=m +no_defs"
    )
    transformer = Transformer.from_crs(4326, crs, always_xy=True)
    lons = np.array([lon_min, lon_min, lon_max, lon_max])
    lats = np.array([lat_min, lat_max, lat_min, lat_max])
    projected_x, projected_y = transformer.transform(lons, lats)
    x = np.asarray(ds.x, dtype=float)
    y = np.asarray(ds.y, dtype=float)
    dx = abs(float(np.median(np.diff(x))))
    dy = abs(float(np.median(np.diff(y))))
    # Native-pixel padding retains source pixels displaced into the final box
    # by terrain/parallax during orthorectification.
    ix = np.where((x >= projected_x.min() / height - 4 * dx) & (x <= projected_x.max() / height + 4 * dx))[0]
    iy = np.where((y >= projected_y.min() / height - 4 * dy) & (y <= projected_y.max() / height + 4 * dy))[0]
    if not len(ix) or not len(iy):
        raise ValueError("Mammoth bounds selected no native ACM pixels")
    return slice(int(iy.min()), int(iy.max()) + 1), slice(int(ix.min()), int(ix.max()) + 1)


def save_remote_subset(fs, key, output, bounds):
    """Read HDF5 chunks via S3 ranges; never stage the full CONUS object."""
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and output.stat().st_size > 0:
        return
    temporary = output.with_suffix(".nc.part")
    with fs.open(key, "rb", block_size=4 * 1024**2, cache_type="blockcache") as remote:
        # Some valid GOES-17 objects carry fill values in the unused
        # time_bounds variable. CF decoding tries to interpret those fills as
        # enormous timestamps and aborts before BCM can be read, so retain the
        # source numeric time coordinates and let the local subset decode them.
        with xr.open_dataset(remote, engine="h5netcdf", decode_times=False) as ds:
            iy, ix = fixed_grid_slice(ds, bounds)
            if "BCM" not in ds:
                raise KeyError(f"BCM is absent from {key}")
            subset = ds[["BCM", "goes_imager_projection"]].isel(y=iy, x=ix).load()
            # Decoded scan angles are already floating point; discard the
            # source packed-coordinate encoding before writing the subset.
            subset.x.encoding = {}
            subset.y.encoding = {}
            if "t" in subset.coords:
                subset.t.encoding = {}
            subset.attrs.update({
                "remote_source": f"s3://{key}",
                "spatial_subset_bounds": ",".join(map(str, bounds)),
                "download_method": "S3 byte-range reads; full CONUS object not staged",
            })
            subset.to_netcdf(temporary, engine="h5netcdf")
    temporary.replace(output)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", required=True, help="YYYY-MM-DD")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--satellite", default="goes18", choices=("goes16", "goes17", "goes18"))
    parser.add_argument("--bounds", nargs=4, type=float, default=(-119.32, 37.41, -118.75, 37.86))
    parser.add_argument("--workers", type=int, default=8, help="Concurrent S3 subset reads")
    parser.add_argument("--max-files", type=int, help="Testing only")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    date = pd.Timestamp(args.date)
    fs = s3fs.S3FileSystem(anon=True)
    keys = []
    for hour in HOURS:
        for scan_mode in (3, 6):
            prefix = (
                f"noaa-{args.satellite}/ABI-L2-ACMC/{date.year}/{date.dayofyear:03d}/{hour:02d}/"
                f"OR_ABI-L2-ACMC-M{scan_mode}_G{args.satellite[-2:]}_*.nc"
            )
            keys.extend(fs.glob(prefix))
    keys = sorted(set(keys))
    if args.max_files:
        keys = keys[:args.max_files]
    if not keys:
        print(f"No {args.satellite.upper()} ACM objects for {date.date()}")
        return 0
    bounds = tuple(args.bounds)
    worker_count = min(args.workers, len(keys))
    with ThreadPoolExecutor(max_workers=worker_count) as pool:
        futures = {
            pool.submit(save_remote_subset, fs, key, args.output_dir / Path(key).name, bounds): key
            for key in keys
        }
        for number, future in enumerate(as_completed(futures), start=1):
            key = futures[future]
            future.result()
            print(f"Saved spatial subset {number}/{len(keys)}: {Path(key).name}", flush=True)
    print(f"Saved {len(keys)} spatial ACM subsets for {date.date()} to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
