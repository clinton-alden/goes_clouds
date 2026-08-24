#!/usr/bin/env python3
"""Extract Mammoth/CUES-south hourly cloud fraction from native GOES ACM files."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS, Transformer


CUES_LAT = 37.6431
CUES_LON = -119.0291
LAT_MIN = CUES_LAT - 5.0 / 111.32
LAT_MAX = CUES_LAT
LON_HALF_WIDTH = 5.0 / (111.32 * np.cos(np.deg2rad(CUES_LAT)))
LON_MIN = CUES_LON - LON_HALF_WIDTH
LON_MAX = CUES_LON + LON_HALF_WIDTH


def projection(ds: xr.Dataset) -> CRS:
    attrs = ds.goes_imager_projection.attrs
    return CRS.from_proj4(
        f"+proj=geos +h={attrs['perspective_point_height']} "
        f"+lon_0={attrs['longitude_of_projection_origin']} "
        f"+sweep={attrs['sweep_angle_axis']} +a={attrs['semi_major_axis']} "
        f"+b={attrs['semi_minor_axis']} +units=m +no_defs"
    )


def native_domain(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    crs = projection(ds)
    transformer = Transformer.from_crs(4326, crs, always_xy=True)
    corner_lon = np.array([LON_MIN, LON_MIN, LON_MAX, LON_MAX])
    corner_lat = np.array([LAT_MIN, LAT_MAX, LAT_MIN, LAT_MAX])
    corner_x, corner_y = transformer.transform(corner_lon, corner_lat)
    height = float(ds.goes_imager_projection.attrs["perspective_point_height"])
    x = np.asarray(ds.x, dtype=float)
    y = np.asarray(ds.y, dtype=float)
    dx = abs(float(np.median(np.diff(x))))
    dy = abs(float(np.median(np.diff(y))))
    # Pad transformed-corner bounds because lines of constant latitude/longitude
    # are curved and skewed in GOES fixed-grid coordinates at Mammoth.
    ix = np.where(
        (x >= np.min(corner_x) / height - 2 * dx)
        & (x <= np.max(corner_x) / height + 2 * dx)
    )[0]
    iy = np.where(
        (y >= np.min(corner_y) / height - 2 * dy)
        & (y <= np.max(corner_y) / height + 2 * dy)
    )[0]
    if not len(ix) or not len(iy):
        raise ValueError("Requested CUES-south domain selected no ACM pixels")
    # Include every half-bordered pixel: transform all four half-grid corners
    # and retain any native pixel whose footprint intersects the geographic box.
    xx, yy = np.meshgrid(x[ix], y[iy])
    inverse = Transformer.from_crs(crs, 4326, always_xy=True)
    footprint_lons = []
    footprint_lats = []
    for x_sign, y_sign in ((-0.5, -0.5), (-0.5, 0.5), (0.5, -0.5), (0.5, 0.5)):
        lon, lat = inverse.transform(
            (xx + x_sign * dx) * height,
            (yy + y_sign * dy) * height,
        )
        footprint_lons.append(lon)
        footprint_lats.append(lat)
    pixel_lon_min = np.min(footprint_lons, axis=0)
    pixel_lon_max = np.max(footprint_lons, axis=0)
    pixel_lat_min = np.min(footprint_lats, axis=0)
    pixel_lat_max = np.max(footprint_lats, axis=0)
    geographic_mask = (
        (pixel_lat_max >= LAT_MIN) & (pixel_lat_min <= LAT_MAX)
        & (pixel_lon_max >= LON_MIN) & (pixel_lon_min <= LON_MAX)
    )
    return iy, ix, geographic_mask


def process_files(files: list[Path]) -> pd.DataFrame:
    rows = []
    for number, path in enumerate(files, start=1):
        try:
            with xr.open_dataset(path, mask_and_scale=True) as ds:
                iy, ix, geographic_mask = native_domain(ds)
                bcm = np.asarray(ds.BCM.isel(y=iy, x=ix), dtype=float)
                valid = np.isfinite(bcm) & geographic_mask
                rows.append({
                    "time": pd.Timestamp(ds.t.values),
                    "acm_cloud_pixels": int(((bcm == 1) & valid).sum()),
                    "acm_valid_pixels": int(valid.sum()),
                    "acm_grid_pixels": int(geographic_mask.sum()),
                    "source_file": str(path),
                })
        except (OSError, ValueError, KeyError) as exc:
            print(f"Skipping unreadable ACM file {path}: {exc}", flush=True)
        if number % 100 == 0 or number == len(files):
            print(f"Processed {number}/{len(files)} ACM scans", flush=True)
    if not rows:
        raise RuntimeError("No readable ACM scans were found")
    scans = pd.DataFrame(rows).sort_values("time").drop_duplicates("time")
    scans["hour"] = scans.time.dt.floor("h")
    hourly = scans.groupby("hour", as_index=False).agg(
        acm_cloud_pixels=("acm_cloud_pixels", "sum"),
        acm_valid_pixels=("acm_valid_pixels", "sum"),
        acm_n_scans=("time", "size"),
        acm_grid_pixels=("acm_grid_pixels", "median"),
    )
    hourly["acm_cloud_fraction"] = hourly.acm_cloud_pixels / hourly.acm_valid_pixels
    return hourly


def process_ortho_zarr(stores: list[Path]) -> pd.DataFrame:
    """Extract hourly BCM fractions from daily orthorectified Mammoth stores."""
    rows = []
    for number, path in enumerate(stores, start=1):
        try:
            with xr.open_zarr(path) as ds:
                variable = "BCM" if "BCM" in ds else "ACM"
                if variable != "BCM":
                    raise KeyError(f"{path} has ACM but not the required binary BCM field")
                latitude = np.asarray(ds.latitude, dtype=float)
                longitude = np.asarray(ds.longitude, dtype=float)
                geographic_mask = (
                    (latitude[:, None] >= LAT_MIN) & (latitude[:, None] <= LAT_MAX)
                    & (longitude[None, :] >= LON_MIN) & (longitude[None, :] <= LON_MAX)
                )
                if not geographic_mask.any():
                    raise ValueError("Requested CUES-south domain selected no orthorectified ACM pixels")
                times = np.atleast_1d(pd.to_datetime(ds.t.values))
                for scan_index, scan_time in enumerate(times):
                    bcm = np.asarray(ds[variable].isel(t=scan_index), dtype=float)
                    valid = np.isfinite(bcm) & geographic_mask
                    rows.append({
                        "time": pd.Timestamp(scan_time),
                        "acm_cloud_pixels": int(((bcm == 1) & valid).sum()),
                        "acm_valid_pixels": int(valid.sum()),
                        "acm_grid_pixels": int(geographic_mask.sum()),
                        "source_file": str(path),
                    })
        except (OSError, ValueError, KeyError) as exc:
            print(f"Skipping unreadable orthorectified ACM store {path}: {exc}", flush=True)
        if number % 10 == 0 or number == len(stores):
            print(f"Processed {number}/{len(stores)} orthorectified ACM days", flush=True)
    if not rows:
        raise RuntimeError("No readable orthorectified ACM scans were found")
    scans = pd.DataFrame(rows).sort_values("time").drop_duplicates("time")
    scans["hour"] = scans.time.dt.floor("h")
    hourly = scans.groupby("hour", as_index=False).agg(
        acm_cloud_pixels=("acm_cloud_pixels", "sum"),
        acm_valid_pixels=("acm_valid_pixels", "sum"),
        acm_n_scans=("time", "size"),
        acm_grid_pixels=("acm_grid_pixels", "median"),
    )
    hourly["acm_cloud_fraction"] = hourly.acm_cloud_pixels / hourly.acm_valid_pixels
    return hourly


def process_ortho_files(files: list[Path]) -> pd.DataFrame:
    """Extract hourly BCM fractions directly from orthorectified scan files."""
    rows = []
    for number, path in enumerate(files, start=1):
        try:
            with xr.open_dataset(path) as ds:
                if "BCM" not in ds:
                    raise KeyError(f"BCM is absent from {path}")
                latitude = np.asarray(ds.latitude, dtype=float)
                longitude = np.asarray(ds.longitude, dtype=float)
                geographic_mask = (
                    (latitude[:, None] >= LAT_MIN) & (latitude[:, None] <= LAT_MAX)
                    & (longitude[None, :] >= LON_MIN) & (longitude[None, :] <= LON_MAX)
                )
                bcm = np.asarray(ds.BCM, dtype=float)
                valid = np.isfinite(bcm) & geographic_mask
                rows.append({
                    "time": pd.Timestamp(ds.t.values),
                    "acm_cloud_pixels": int(((bcm == 1) & valid).sum()),
                    "acm_valid_pixels": int(valid.sum()),
                    "acm_grid_pixels": int(geographic_mask.sum()),
                    "source_file": str(path),
                })
        except (OSError, ValueError, KeyError) as exc:
            print(f"Skipping unreadable orthorectified ACM file {path}: {exc}", flush=True)
        if number % 100 == 0 or number == len(files):
            print(f"Processed {number}/{len(files)} orthorectified ACM scans", flush=True)
    if not rows:
        raise RuntimeError("No readable orthorectified ACM scans were found")
    scans = pd.DataFrame(rows).sort_values("time").drop_duplicates("time")
    scans["hour"] = scans.time.dt.floor("h")
    hourly = scans.groupby("hour", as_index=False).agg(
        acm_cloud_pixels=("acm_cloud_pixels", "sum"),
        acm_valid_pixels=("acm_valid_pixels", "sum"),
        acm_n_scans=("time", "size"),
        acm_grid_pixels=("acm_grid_pixels", "median"),
    )
    hourly["acm_cloud_fraction"] = hourly.acm_cloud_pixels / hourly.acm_valid_pixels
    return hourly


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--combine-dir", type=Path)
    parser.add_argument("--combine-output", type=Path)
    parser.add_argument("--ortho-zarr-dir", type=Path)
    parser.add_argument("--ortho-dir", type=Path)
    parser.add_argument("--year-month", help="YYYYMM filter for --ortho-zarr-dir")
    args = parser.parse_args()
    if args.ortho_dir:
        files = sorted(args.ortho_dir.rglob("*_ortho.nc"))
        if not files:
            raise FileNotFoundError(f"No orthorectified ACM files under {args.ortho_dir}")
        hourly = process_ortho_files(files)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        hourly.to_csv(args.output, index=False)
        print(f"Saved {len(hourly)} hourly rows from {len(files)} orthorectified scans: {args.output}")
        return 0
    if args.ortho_zarr_dir:
        if not args.year_month or len(args.year_month) != 6:
            parser.error("--ortho-zarr-dir requires --year-month YYYYMM")
        stores = sorted(args.ortho_zarr_dir.glob(f"*_{args.year_month}??.zarr"))
        if not stores:
            raise FileNotFoundError(f"No {args.year_month} orthorectified ACM stores in {args.ortho_zarr_dir}")
        hourly = process_ortho_zarr(stores)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        hourly.to_csv(args.output, index=False)
        print(f"Saved {len(hourly)} hourly rows from {len(stores)} orthorectified days: {args.output}")
        return 0
    if args.combine_dir:
        files = sorted(args.combine_dir.glob("mammoth_acm_hourly_????????.csv"))
        if not files:
            raise FileNotFoundError(f"No daily ACM tables in {args.combine_dir}")
        combined = pd.concat((pd.read_csv(p, parse_dates=["hour"]) for p in files), ignore_index=True)
        combined = combined.sort_values("hour").drop_duplicates("hour")
        args.combine_output.parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(args.combine_output, index=False)
        print(f"Combined {len(files)} days and {len(combined):,} hours: {args.combine_output}")
        return 0
    if args.input_dir is None or args.output is None:
        parser.error("provide native input, orthorectified Zarr input, or combine input")
    files = sorted(args.input_dir.rglob("OR_ABI-L2-ACMC-*.nc"))
    if not files:
        raise FileNotFoundError(f"No ACM files under {args.input_dir}")
    hourly = process_files(files)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    hourly.to_csv(args.output, index=False)
    print(f"Saved {len(hourly)} hourly rows: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
