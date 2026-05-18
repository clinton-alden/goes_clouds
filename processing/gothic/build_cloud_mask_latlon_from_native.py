#!/usr/bin/env python3
"""Remap daily Gothic cloud masks from native GOES x/y grid to lat/lon grid."""

from __future__ import annotations

import argparse
import calendar
from pathlib import Path

import xarray as xr


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create daily cloud mask NetCDF files on lat/lon grid while preserving "
            "5-minute timesteps (t dimension)."
        )
    )
    parser.add_argument(
        "--native-mask-dir",
        default="/glade/derecho/scratch/cdalden/gothic/goes16/cloud_mask_tempbin_10c",
        help="Directory containing native-grid daily cloud mask files",
    )
    parser.add_argument(
        "--ortho-map-path",
        default="/glade/u/home/cdalden/.ortho_cache/ortho_map_ce81dc15_00bab0b9.nc",
        help="Path to ortho map file containing dem_px_angle_x/y and lat/lon grid",
    )
    parser.add_argument(
        "--output-dir",
        default="/glade/derecho/scratch/cdalden/gothic/goes16/cloud_mask_tempbin_10c_latlon",
        help="Output directory for daily lat/lon cloud mask NetCDF files",
    )
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing output files"
    )
    return parser


def remap_day(native_path: Path, out_path: Path, dem_x: xr.DataArray, dem_y: xr.DataArray) -> None:
    ds = xr.open_dataset(native_path)

    cloud_native = ds["cloud_binary"]
    # Drop non-dimension coords (for example scalar t from ortho map) that can
    # conflict with the daily time coordinate during interpolation.
    dem_x_clean = dem_x.reset_coords(drop=True)
    dem_y_clean = dem_y.reset_coords(drop=True)
    cloud_ll = cloud_native.interp(x=dem_x_clean, y=dem_y_clean, method="nearest")
    cloud_ll = cloud_ll.fillna(0).astype("uint8")

    out_ds = xr.Dataset(
        data_vars={
            "cloud_binary": cloud_ll,
            "air_temp_c": ds["air_temp_c"],
            "temp_bin_index": ds["temp_bin_index"],
        },
        attrs={
            "title": "Gothic cloud mask (lat/lon grid) from ERA5-temperature-selected RGB thresholds",
            "source_native_mask": str(native_path),
            "remap_method": "nearest",
        },
    )

    out_ds["cloud_binary"].attrs.update(
        {
            "long_name": "cloud mask (0=clear, 1=cloudy)",
            "grid": "latlon",
        }
    )

    encoding = {
        "cloud_binary": {"zlib": True, "complevel": 4, "dtype": "uint8"},
        "air_temp_c": {"zlib": True, "complevel": 4, "dtype": "float32"},
        "temp_bin_index": {"zlib": True, "complevel": 4, "dtype": "int16"},
    }
    out_ds.to_netcdf(out_path, encoding=encoding)

    ds.close()
    out_ds.close()


def main() -> int:
    args = build_parser().parse_args()

    native_mask_dir = Path(args.native_mask_dir)
    ortho_map_path = Path(args.ortho_map_path)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not ortho_map_path.exists():
        raise FileNotFoundError(f"Missing ortho map file: {ortho_map_path}")

    map_ds = xr.open_dataset(ortho_map_path)
    dem_x = map_ds["dem_px_angle_x"]
    dem_y = map_ds["dem_px_angle_y"]

    ndays = calendar.monthrange(args.year, args.month)[1]
    print(f"Building lat/lon cloud masks for {args.year}-{args.month:02d} ({ndays} days)")
    print(f"Native input dir: {native_mask_dir}")
    print(f"Output dir: {output_dir}")

    for day in range(1, ndays + 1):
        date = f"{args.year}{args.month:02d}{day:02d}"
        native_path = native_mask_dir / f"goes16_cloud_binary_tempbin10c_gothic_{date}.nc"
        out_path = output_dir / f"goes16_cloud_binary_tempbin10c_gothic_latlon_{date}.nc"

        if not native_path.exists():
            print(f"[skip] missing native mask {native_path}")
            continue
        if out_path.exists() and not args.overwrite:
            print(f"[skip] already exists {out_path}")
            continue

        print(f"[remap] {date}")
        remap_day(native_path, out_path, dem_x, dem_y)

    map_ds.close()
    print("Lat/lon cloud mask generation complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
