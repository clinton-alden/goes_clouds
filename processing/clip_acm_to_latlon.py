#!/usr/bin/env python3
"""Resample GOES ACM files onto an exact lat/lon grid over a target domain."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from pyproj import CRS, Transformer
from scipy.interpolate import griddata
import xarray as xr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Clip a GOES ABI-L2-ACMC file to an exact lat/lon grid."
    )
    parser.add_argument("input", help="Input GOES ACM NetCDF path")
    parser.add_argument("output", help="Output clipped NetCDF path")
    parser.add_argument("--lon-min", type=float, required=True)
    parser.add_argument("--lon-max", type=float, required=True)
    parser.add_argument("--lat-min", type=float, required=True)
    parser.add_argument("--lat-max", type=float, required=True)
    parser.add_argument("--mask-var", default="AUTO")
    parser.add_argument("--target-width", type=int, default=241)
    parser.add_argument("--target-height", type=int, default=161)
    return parser.parse_args()


def resolve_mask_var(ds: xr.Dataset, requested: str) -> str:
    if requested != "AUTO":
        if requested not in ds:
            raise KeyError(f"Variable '{requested}' not found in dataset.")
        return requested
    for candidate in ("BCM", "ACM"):
        if candidate in ds:
            return candidate
    raise KeyError("Neither BCM nor ACM found in dataset.")


def calculate_lat_lon(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    projection_info = ds["goes_imager_projection"]
    crs_geos = CRS.from_cf(
        {
            "grid_mapping_name": projection_info.attrs["grid_mapping_name"],
            "perspective_point_height": projection_info.attrs["perspective_point_height"],
            "semi_major_axis": projection_info.attrs["semi_major_axis"],
            "semi_minor_axis": projection_info.attrs["semi_minor_axis"],
            "latitude_of_projection_origin": projection_info.attrs["latitude_of_projection_origin"],
            "longitude_of_projection_origin": projection_info.attrs["longitude_of_projection_origin"],
            "sweep_angle_axis": projection_info.attrs["sweep_angle_axis"],
        }
    )
    transformer = Transformer.from_crs(crs_geos, CRS.from_epsg(4326), always_xy=True)

    x_1d = ds["x"].values * projection_info.attrs["perspective_point_height"]
    y_1d = ds["y"].values * projection_info.attrs["perspective_point_height"]
    x_2d, y_2d = np.meshgrid(x_1d, y_1d)
    lon_2d, lat_2d = transformer.transform(x_2d, y_2d)
    return lat_2d, lon_2d


def build_target_grid(lon_min: float, lon_max: float, lat_min: float, lat_max: float, target_width: int, target_height: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    target_lon = np.linspace(lon_min, lon_max, target_width)
    target_lat = np.linspace(lat_min, lat_max, target_height)
    target_lon_2d, target_lat_2d = np.meshgrid(target_lon, target_lat)
    return target_lon, target_lat, target_lon_2d, target_lat_2d


def resample_to_latlon_grid(ds: xr.Dataset, mask_var: str, lon_min: float, lon_max: float, lat_min: float, lat_max: float, target_width: int, target_height: int) -> xr.Dataset:
    mask_var = resolve_mask_var(ds, mask_var)
    lat_2d, lon_2d = calculate_lat_lon(ds)
    data = ds[mask_var].values.astype(float)

    valid = (
        np.isfinite(lat_2d)
        & np.isfinite(lon_2d)
        & np.isfinite(data)
        & (lon_2d >= lon_min)
        & (lon_2d <= lon_max)
        & (lat_2d >= lat_min)
        & (lat_2d <= lat_max)
    )
    if not np.any(valid):
        raise ValueError("No valid mask pixels found inside the requested bounds.")

    target_lon, target_lat, target_lon_2d, target_lat_2d = build_target_grid(
        lon_min, lon_max, lat_min, lat_max, target_width, target_height
    )
    resampled = griddata(
        np.column_stack((lon_2d[valid], lat_2d[valid])),
        data[valid],
        (target_lon_2d, target_lat_2d),
        method="nearest",
    )

    out = xr.Dataset(
        data_vars={
            mask_var: (("latitude", "longitude"), resampled.astype(np.float32)),
        },
        coords={
            "latitude": target_lat.astype(np.float32),
            "longitude": target_lon.astype(np.float32),
        },
        attrs=dict(ds.attrs),
    )

    if "DQF" in ds:
        dqf = ds["DQF"].values.astype(float)
        dqf_valid = np.isfinite(dqf) & valid
        out["DQF"] = (
            ("latitude", "longitude"),
            griddata(
                np.column_stack((lon_2d[dqf_valid], lat_2d[dqf_valid])),
                dqf[dqf_valid],
                (target_lon_2d, target_lat_2d),
                method="nearest",
            ).astype(np.float32),
        )

    if "t" in ds.coords:
        out = out.assign_coords(time=ds["t"])

    out.attrs["domain_bounds"] = f"{lon_min},{lon_max},{lat_min},{lat_max}"
    out.attrs["clipping_method"] = "nearest-neighbor resample to exact lat/lon grid"
    out.attrs["mask_variable"] = mask_var
    return out


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(input_path) as ds:
        clipped = resample_to_latlon_grid(
            ds,
            mask_var=args.mask_var,
            lon_min=args.lon_min,
            lon_max=args.lon_max,
            lat_min=args.lat_min,
            lat_max=args.lat_max,
            target_width=args.target_width,
            target_height=args.target_height,
        )
        encoding = {name: {"zlib": True, "complevel": 4} for name in clipped.data_vars}
        clipped.to_netcdf(output_path, encoding=encoding)

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
