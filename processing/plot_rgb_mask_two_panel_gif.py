#!/usr/bin/env python3
"""Render a coordinate-aware two-panel RGB/cloud-mask GIF."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import ListedColormap


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("rgb", type=Path)
    p.add_argument("mask", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--start-hour", type=int, default=16)
    p.add_argument("--end-hour", type=int, default=24)
    p.add_argument("--duration", type=float, default=0.18)
    a = p.parse_args()
    frames = a.output.parent / f".{a.output.stem}_frames"
    shutil.rmtree(frames, ignore_errors=True)
    frames.mkdir(parents=True)
    paths = []
    try:
        with xr.open_dataset(a.rgb) as r, xr.open_dataset(a.mask) as m:
            times = pd.DatetimeIndex(pd.to_datetime(r.t.values))
            hour = times.hour + times.minute / 60
            keep = np.where((hour >= a.start_hour) & (hour < a.end_hour))[0]
            if not len(keep):
                raise ValueError("No RGB frames fall in the requested time range")
            extent = [float(r.longitude.min()), float(r.longitude.max()), float(r.latitude.min()), float(r.latitude.max())]
            cmap = ListedColormap(["#2563eb", "#ffffff"])
            for frame_number, i in enumerate(keep):
                rgb = np.dstack([np.asarray(r[name].isel(t=i)) for name in ("red", "green", "blue")])
                mask = np.asarray(m.cloud_binary.sel(t=r.t.isel(t=i), method="nearest"))
                fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)
                axes[0].imshow(rgb, origin="lower", extent=extent, aspect="auto", interpolation="nearest")
                axes[1].imshow(mask, origin="lower", extent=extent, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
                for ax, title in zip(axes, ("Day Cloud Phase RGB", "Decision-tree cloud mask")):
                    ax.set_title(title, fontweight="bold")
                    ax.set_xlabel("Longitude")
                    ax.set_ylabel("Latitude")
                fig.suptitle(times[i].strftime("%Y-%m-%d %H:%M UTC"), fontweight="bold")
                path = frames / f"frame_{frame_number:04d}.png"
                fig.savefig(path, dpi=100)
                plt.close(fig)
                paths.append(path)
        a.output.parent.mkdir(parents=True, exist_ok=True)
        with imageio.get_writer(a.output, mode="I", duration=a.duration, loop=0) as writer:
            for path in paths:
                writer.append_data(imageio.imread(path))
    finally:
        shutil.rmtree(frames, ignore_errors=True)
    print(f"Wrote {len(paths)} frames: {a.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
