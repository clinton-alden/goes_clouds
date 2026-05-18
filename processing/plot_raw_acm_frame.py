#!/usr/bin/env python3
"""Plot one raw GOES ACM file in geographic coordinates."""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
from pyproj import CRS, Transformer
from scipy.interpolate import griddata
import xarray as xr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a raw, non-orthorectified GOES ACM frame in lon/lat space."
    )
    parser.add_argument("--input", required=True, help="Raw GOES ACM NetCDF file")
    parser.add_argument("--output", required=True, help="Output image path")
    parser.add_argument("--lon-min", type=float, required=True)
    parser.add_argument("--lon-max", type=float, required=True)
    parser.add_argument("--lat-min", type=float, required=True)
    parser.add_argument("--lat-max", type=float, required=True)
    parser.add_argument(
        "--target-width",
        type=int,
        default=1250,
        help="Output longitude grid size for resampled plot.",
    )
    parser.add_argument(
        "--target-height",
        type=int,
        default=750,
        help="Output latitude grid size for resampled plot.",
    )
    return parser.parse_args()


def calculate_degrees(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
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

    x_coordinate_1d = ds["x"].values * projection_info.attrs["perspective_point_height"]
    y_coordinate_1d = ds["y"].values * projection_info.attrs["perspective_point_height"]
    x_coordinate_2d, y_coordinate_2d = np.meshgrid(x_coordinate_1d, y_coordinate_1d)

    abi_lon, abi_lat = transformer.transform(x_coordinate_2d, y_coordinate_2d)
    return abi_lat, abi_lon


def main() -> int:
    args = parse_args()
    ds = xr.open_dataset(args.input)
    try:
        acm = ds["ACM"].values.astype(float)
        abi_lat, abi_lon = calculate_degrees(ds)

        bounds_mask = (
            (abi_lon >= args.lon_min)
            & (abi_lon <= args.lon_max)
            & (abi_lat >= args.lat_min)
            & (abi_lat <= args.lat_max)
        )
        valid_mask = np.isfinite(abi_lon) & np.isfinite(abi_lat) & np.isfinite(acm) & bounds_mask
        if not np.any(valid_mask):
            raise ValueError("No valid ACM pixels found inside the requested bounds.")
        source_lon = abi_lon[valid_mask]
        source_lat = abi_lat[valid_mask]
        source_acm = acm[valid_mask]

        target_lon = np.linspace(args.lon_min, args.lon_max, args.target_width)
        target_lat = np.linspace(args.lat_min, args.lat_max, args.target_height)
        target_lon_2d, target_lat_2d = np.meshgrid(target_lon, target_lat)
        resampled_acm = griddata(
            np.column_stack((source_lon, source_lat)),
            source_acm,
            (target_lon_2d, target_lat_2d),
            method="nearest",
        )

        cmap = ListedColormap(["#1b9e77", "#66a61e", "#fdb863", "#d73027"])
        norm = BoundaryNorm(boundaries=np.arange(-0.5, 4.5, 1), ncolors=cmap.N)
        cmap.set_bad(color=(1, 1, 1, 0))

        fig, ax = plt.subplots(figsize=(9, 7), dpi=180)
        mesh = ax.imshow(
            resampled_acm,
            extent=[args.lon_min, args.lon_max, args.lat_min, args.lat_max],
            origin="lower",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            aspect="auto",
        )
        cbar = fig.colorbar(mesh, ax=ax, shrink=0.82, ticks=[0, 1, 2, 3])
        cbar.ax.set_yticklabels(
            ["clear", "probably clear", "probably cloudy", "cloudy"]
        )

        time_str = np.datetime_as_string(ds["t"].values, unit="m")
        ax.set_xlim(args.lon_min, args.lon_max)
        ax.set_ylim(args.lat_min, args.lat_max)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"Raw GOES-16 ACM (non-ortho)\n{time_str} UTC")
        fig.tight_layout()
        fig.savefig(args.output, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {args.output}")
    finally:
        ds.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
