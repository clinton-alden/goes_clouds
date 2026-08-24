#!/usr/bin/env python3
"""Create a timestamp-aligned, two-panel GOES RGB/cloud-mask GIF."""

import argparse
from pathlib import Path

import imageio.v2 as imageio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import xarray as xr


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("rgb", type=Path)
    p.add_argument("mask", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--fps", type=float, default=8)
    p.add_argument("--region", help="Display name for the mapped region")
    a = p.parse_args()
    a.output.parent.mkdir(parents=True, exist_ok=True)

    satellite = a.rgb.name.split("_", 1)[0].upper().replace("GOES", "GOES-")
    if a.region:
        region = a.region
    elif "mammoth" in a.rgb.name.lower():
        region = "Mammoth"
    elif "east_river" in a.rgb.name.lower():
        region = "East River"
    else:
        region = "GOES domain"

    with xr.open_dataset(a.rgb) as rgb_ds, xr.open_dataset(a.mask) as mask_ds:
        rgb_ds, mask_ds = xr.align(rgb_ds, mask_ds, join="inner")
        if rgb_ds.sizes.get("t", 0) == 0:
            raise ValueError("RGB and mask have no matching timestamps")
        rgb = np.stack([rgb_ds[c].values for c in ("red", "green", "blue")], axis=-1)
        rgb = np.clip(rgb, 0, 1)
        mask = mask_ds.cloud_binary.values
        times = pd.to_datetime(rgb_ds.t.values)
        lon = rgb_ds.longitude.values
        lat = rgb_ds.latitude.values
        if not np.all(np.diff(lon) > 0) or not np.all(np.diff(lat) > 0):
            raise ValueError("Expected west-to-east longitude and south-to-north latitude")
        extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]

        with imageio.get_writer(a.output, mode="I", duration=1 / a.fps, loop=0) as writer:
            for i, timestamp in enumerate(times):
                fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
                axes[0].imshow(rgb[i], origin="lower", extent=extent, aspect="auto")
                axes[0].set_title(f"{satellite} RGB")
                axes[1].imshow(mask[i], origin="lower", extent=extent, aspect="auto",
                               cmap=ListedColormap(["#2455a4", "#ffffff"]), vmin=0, vmax=1,
                               interpolation="nearest")
                axes[1].set_title("10°C-bin cloud mask\nblue = clear, white = cloudy")
                for ax in axes:
                    ax.set_xlabel("Longitude")
                    ax.set_ylabel("Latitude")
                    ax.grid(color="black", alpha=0.18, linewidth=0.5)
                fig.suptitle(f"{region} — {timestamp:%Y-%m-%d %H:%M UTC}", fontsize=14)
                fig.canvas.draw()
                frame = np.asarray(fig.canvas.buffer_rgba())[..., :3]
                writer.append_data(frame)
                plt.close(fig)
    print(f"Wrote {len(times)} frames: {a.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
