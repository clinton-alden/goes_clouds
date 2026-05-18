#!/usr/bin/env python3
"""Render a GIF loop from orthorectified GOES ACM NetCDF files."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import xarray as xr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a GIF from orthorectified GOES ACM NetCDF files."
    )
    parser.add_argument("--input-dir", required=True, help="Day directory to scan")
    parser.add_argument("--output", required=True, help="Output GIF path")
    parser.add_argument("--date", required=True, help="YYYYMMDD date token")
    parser.add_argument("--start-hour", type=int, default=14, help="Inclusive UTC hour")
    parser.add_argument("--end-hour", type=int, default=23, help="Inclusive UTC hour")
    parser.add_argument(
        "--frame-step",
        type=int,
        default=1,
        help="Keep every Nth frame after time filtering.",
    )
    parser.add_argument("--domain", default="colorado", help="Domain label for title")
    parser.add_argument(
        "--tmp-dir",
        required=True,
        help="Scratch directory for temporary PNG frames",
    )
    return parser.parse_args()


def iter_ortho_files(day_dir: str) -> list[str]:
    files = []
    for root, _, names in os.walk(day_dir):
        for name in names:
            if name.endswith("_ortho.nc"):
                files.append(os.path.join(root, name))
    return sorted(files)


def pick_mask_variable(ds: xr.Dataset) -> str | None:
    for candidate in ("BCM", "ACM"):
        if candidate in ds:
            return candidate
    return None


def timestamp_parts(ds: xr.Dataset, nc_path: str) -> tuple[int, int]:
    if "t" in ds:
        try:
            timestamp = ds["t"].values
            ts = np.asarray(timestamp).reshape(-1)[0]
            ts64 = np.datetime64(ts)
            return int(str(ts64)[11:13]), int(str(ts64)[14:16])
        except Exception:
            pass

    match = re.search(r"_s\d{7}(\d{2})(\d{2})\d{2}", Path(nc_path).name)
    if match:
        return int(match.group(1)), int(match.group(2))

    raise KeyError(f"Could not determine timestamp for {nc_path}")


def main() -> int:
    args = parse_args()
    os.makedirs(args.tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    frame_paths: list[str] = []

    kept_idx = 0
    for idx, nc_path in enumerate(iter_ortho_files(args.input_dir)):
        with xr.open_dataset(nc_path) as ds:
            var_name = pick_mask_variable(ds)
            if var_name is None:
                continue

            hour, minute = timestamp_parts(ds, nc_path)
            if hour < args.start_hour or hour > args.end_hour:
                continue

            if kept_idx % args.frame_step != 0:
                kept_idx += 1
                continue

            acm = ds[var_name].values
            lon = ds["longitude"].values
            lat = ds["latitude"].values
            kept_idx += 1

        if var_name == "BCM":
            cmap = ListedColormap(
                [
                    "#1b9e77",  # clear
                    "#d73027",  # cloudy
                ]
            )
            norm = BoundaryNorm(boundaries=np.arange(-0.5, 2.5, 1), ncolors=cmap.N)
            tick_values = [0, 1]
            tick_labels = ["clear", "cloudy"]
            title_prefix = "BCM"
        else:
            cmap = ListedColormap(
                [
                    "#1b9e77",  # clear
                    "#66a61e",  # probably clear
                    "#fdb863",  # probably cloudy
                    "#d73027",  # cloudy
                ]
            )
            norm = BoundaryNorm(boundaries=np.arange(-0.5, 4.5, 1), ncolors=cmap.N)
            tick_values = [0, 1, 2, 3]
            tick_labels = ["clear", "probably clear", "probably cloudy", "cloudy"]
            title_prefix = "ACM"
        cmap.set_bad(color=(1, 1, 1, 0))

        hhmm = f"{hour:02d}{minute:02d}"
        frame_path = os.path.join(args.tmp_dir, f"{args.date}_{idx:03d}_{hhmm}.png")

        fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
        im = ax.imshow(
            acm,
            extent=[float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())],
            origin="lower",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
        )
        cbar = fig.colorbar(im, ax=ax, shrink=0.82, ticks=tick_values)
        cbar.ax.set_yticklabels(tick_labels)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"GOES-16 {title_prefix} {args.domain} {args.date} {hour:02d}:{minute:02d} UTC")
        fig.tight_layout()
        fig.savefig(frame_path, bbox_inches="tight")
        plt.close(fig)
        frame_paths.append(frame_path)

    if not frame_paths:
        raise SystemExit("No orthorectified ACM frames matched the requested time window.")

    with imageio.get_writer(args.output, mode="I", duration=0.5) as writer:
        for frame_path in frame_paths:
            writer.append_data(imageio.imread(frame_path))

    for frame_path in frame_paths:
        try:
            os.remove(frame_path)
        except FileNotFoundError:
            pass

    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
