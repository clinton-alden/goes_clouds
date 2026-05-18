#!/usr/bin/env python3
"""Render a GIF loop from a daily ACM NetCDF file."""

from __future__ import annotations

import argparse
from pathlib import Path

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import xarray as xr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a GIF from a daily ACM file.")
    parser.add_argument("--input", required=True, help="Daily NetCDF path")
    parser.add_argument("--output", required=True, help="Output GIF path")
    parser.add_argument("--domain", default="colorado")
    parser.add_argument("--mask-var", default="AUTO")
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--frame-step", type=int, default=1)
    return parser.parse_args()


def acm_style():
    cmap = ListedColormap(
        ["#1b9e77", "#66a61e", "#fdb863", "#d73027"]
    )
    norm = BoundaryNorm(boundaries=np.arange(-0.5, 4.5, 1), ncolors=cmap.N)
    cmap.set_bad(color=(1, 1, 1, 0))
    return cmap, norm


def bcm_style():
    cmap = ListedColormap(["#2166ac", "#ffffff"])
    norm = BoundaryNorm(boundaries=np.arange(-0.5, 2.5, 1), ncolors=cmap.N)
    cmap.set_bad(color=(1, 1, 1, 0))
    return cmap, norm


def resolve_mask_var(ds: xr.Dataset, requested: str) -> str:
    if requested != "AUTO":
        if requested not in ds:
            raise KeyError(f"{requested} not found in dataset")
        return requested
    for candidate in ("BCM", "ACM"):
        if candidate in ds:
            return candidate
    raise KeyError("Neither BCM nor ACM found in dataset")


def main() -> int:
    args = parse_args()
    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    frame_paths = []

    with xr.open_dataset(args.input) as ds:
        mask_var = resolve_mask_var(ds, args.mask_var)
        if mask_var == "BCM":
            cmap, norm = bcm_style()
            tick_values = [0, 1]
            tick_labels = ["clear", "cloudy"]
        else:
            cmap, norm = acm_style()
            tick_values = [0, 1, 2, 3]
            tick_labels = ["clear", "probably clear", "probably cloudy", "cloudy"]

        lon = ds["longitude"].values
        lat = ds["latitude"].values
        time_values = ds["time"].values
        data = ds[mask_var].values

        for idx in range(0, data.shape[0], args.frame_step):
            ts = np.datetime_as_string(np.datetime64(time_values[idx]), unit="m")
            frame_path = tmp_dir / f"frame_{idx:03d}.png"

            fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
            mesh = ax.pcolormesh(lon, lat, data[idx], cmap=cmap, norm=norm, shading="nearest")
            cbar = fig.colorbar(mesh, ax=ax, shrink=0.82, ticks=tick_values)
            cbar.ax.set_yticklabels(tick_labels)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.set_title(f"GOES-16 {mask_var} {args.domain} {ts} UTC")
            fig.tight_layout()
            fig.savefig(frame_path, bbox_inches="tight")
            plt.close(fig)
            frame_paths.append(frame_path)

    if not frame_paths:
        raise SystemExit("No frames were generated.")

    with imageio.get_writer(args.output, mode="I", duration=0.4) as writer:
        for frame_path in frame_paths:
            writer.append_data(imageio.imread(frame_path))

    for frame_path in frame_paths:
        frame_path.unlink(missing_ok=True)

    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
