#!/usr/bin/env python3
"""Build a 5-panel GOES/GFS GIF for March 25, 2022."""

from __future__ import annotations

import argparse
from pathlib import Path

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd
import xarray as xr


LON_MIN = -109.0
LON_MAX = -104.0
LAT_MIN = 37.0
LAT_MAX = 41.0

ROOT = Path("/glade/u/home/cdalden/goes_work")
THRESHOLD_PATH = (
    ROOT
    / "analysis/output_12_rgb_threshold_transfer/"
    "gothic_temp_bin_rgb_thresholds_10c.csv"
)

OUT_DIR = ROOT / "analysis/output_15_goes_surface_temp_investigation"


MASK_CLEAR_CLOUD_CMAP = ListedColormap(["#2166ac", "#ffffff"])
BCM_NORM = BoundaryNorm([-0.5, 0.5, 1.5], MASK_CLEAR_CLOUD_CMAP.N)
MASK_NORM = BoundaryNorm([-0.5, 0.5, 1.5], MASK_CLEAR_CLOUD_CMAP.N)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a 16Z-00Z GOES/GFS five-panel GIF for one YYYYMMDD date."
    )
    parser.add_argument("date", nargs="?", default="20220325", help="Date as YYYYMMDD")
    return parser.parse_args()


def date_paths(date: str) -> dict[str, Path]:
    day = pd.Timestamp(date)
    return {
        "rgb": Path(
            "/glade/u/home/cdalden/scratch/colorado/goes16/rgb_composite/"
            f"goes16_C02_C05_C13_rgb_colorado_{date}.nc"
        ),
        "acm": Path(
            "/glade/u/home/cdalden/scratch/colorado_acm/goes16/"
            f"{day.year}/{day.month}/{day.day}/ABI-L2-ACMC/clipped"
        ),
        "gfs": (
            ROOT
            / f"data_download/gfs_surface_{date}/"
            f"gfs_{date}_00z_f000_f023_surface_temperature_snow_mask_109W_104W_37N_41N.nc"
        ),
        "frames": OUT_DIR / f"five_panel_frames_{date}_16z_00z_rgbmask_10c",
        "gif": OUT_DIR / f"goes_gfs_five_panel_{date}_16z_00z_rgbmask_10c.gif",
    }


def ensure_inputs(paths: dict[str, Path]) -> None:
    missing = [
        path
        for path in [paths["rgb"], paths["acm"], paths["gfs"], THRESHOLD_PATH]
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError("Missing required input(s):\n" + "\n".join(map(str, missing)))


def acm_time_table(acm_dir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(acm_dir.glob("*_clip.nc")):
        with xr.open_dataset(path) as ds:
            time_name = "time" if "time" in ds.coords else "t"
            rows.append({"acm_time": pd.Timestamp(ds[time_name].item()), "path": path})
    out = pd.DataFrame(rows).sort_values("acm_time").reset_index(drop=True)
    if out.empty:
        raise FileNotFoundError(f"No clipped ACM files found in {acm_dir}")
    return out


def matched_goes_times(rgb_times: np.ndarray, acm_df: pd.DataFrame) -> pd.DataFrame:
    rgb_df = pd.DataFrame(
        {
            "frame_idx": np.arange(len(rgb_times), dtype=int),
            "rgb_time": pd.to_datetime(rgb_times),
        }
    )
    merged = pd.merge_asof(
        rgb_df.sort_values("rgb_time"),
        acm_df,
        left_on="rgb_time",
        right_on="acm_time",
        direction="nearest",
        tolerance=pd.Timedelta("90s"),
    )
    if merged["path"].isna().any():
        missing = merged.loc[merged["path"].isna(), "rgb_time"].head(5).to_list()
        raise RuntimeError(f"No ACM match within 90 s for RGB time(s), e.g. {missing}")
    return merged.sort_values("frame_idx").reset_index(drop=True)


def load_rgb_thresholds() -> pd.DataFrame:
    thresholds = pd.read_csv(THRESHOLD_PATH)
    thresholds = thresholds.loc[thresholds["status"].eq("trained")].copy()
    required = [
        "temp_left_c",
        "temp_right_c",
        "rule",
        "red_threshold",
        "red_direction",
        "green_threshold",
        "green_direction",
        "blue_threshold",
        "blue_direction",
    ]
    missing = [column for column in required if column not in thresholds.columns]
    if missing:
        raise ValueError(f"Threshold table is missing required column(s): {missing}")
    if thresholds.empty:
        raise ValueError(f"No trained threshold rows found in {THRESHOLD_PATH}")
    return thresholds.sort_values("temp_left_c").reset_index(drop=True)


def threshold_vote(channel: np.ndarray, threshold: float, direction: str) -> np.ndarray:
    if direction == ">":
        return channel > threshold
    if direction == "<":
        return channel < threshold
    raise ValueError(f"Unsupported RGB threshold direction: {direction!r}")


def threshold_cloud_mask(
    *,
    red: np.ndarray,
    green: np.ndarray,
    blue: np.ndarray,
    surface_temperature_c: np.ndarray,
    thresholds: pd.DataFrame,
) -> np.ndarray:
    mask = np.zeros(red.shape, dtype=np.int8)
    for row_idx, row in thresholds.iterrows():
        temp_left = float(row["temp_left_c"])
        temp_right = float(row["temp_right_c"])
        if row_idx == 0:
            in_bin = surface_temperature_c < temp_right
        elif row_idx == len(thresholds) - 1:
            in_bin = surface_temperature_c >= temp_left
        else:
            in_bin = (surface_temperature_c >= temp_left) & (surface_temperature_c < temp_right)

        votes = (
            threshold_vote(red, float(row["red_threshold"]), row["red_direction"]).astype(np.int8)
            + threshold_vote(green, float(row["green_threshold"]), row["green_direction"]).astype(np.int8)
            + threshold_vote(blue, float(row["blue_threshold"]), row["blue_direction"]).astype(np.int8)
        )

        rule = row["rule"]
        if rule == "majority":
            cloudy = votes >= 2
        elif rule == "union":
            cloudy = votes >= 1
        elif rule == "intersection":
            cloudy = votes == 3
        else:
            raise ValueError(f"Unsupported RGB threshold rule: {rule!r}")
        mask[in_bin & cloudy] = 1
    return mask


def gfs_temperature_on_rgb_grid(gfs_hour: xr.Dataset, rgb_ds: xr.Dataset) -> np.ndarray:
    target_lat = xr.DataArray(rgb_ds["latitude"].values, dims="latitude")
    target_lon = xr.DataArray(rgb_ds["longitude"].values, dims="longitude")
    return (
        gfs_hour["surface_temperature_c"]
        .sortby("latitude")
        .sortby("longitude")
        .interp(latitude=target_lat, longitude=target_lon, method="nearest")
        .values
    )


def as_rgb_image(red: np.ndarray, green: np.ndarray, blue: np.ndarray) -> np.ndarray:
    rgb = np.stack([red, green, blue], axis=-1)
    return np.clip(np.nan_to_num(rgb, nan=0.0), 0.0, 1.0)


def imshow_latlon(ax, data, lon, lat, *, cmap=None, norm=None, vmin=None, vmax=None, title=""):
    lat_values = np.asarray(lat)
    origin = "lower" if lat_values[0] < lat_values[-1] else "upper"
    return ax.imshow(
        data,
        extent=[float(np.min(lon)), float(np.max(lon)), float(np.min(lat)), float(np.max(lat))],
        origin=origin,
        cmap=cmap,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        aspect="auto",
    )


def label_geo_axis(ax) -> None:
    ax.set_xlim(LON_MIN, LON_MAX)
    ax.set_ylim(LAT_MIN, LAT_MAX)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")


def render_frame(
    *,
    rgb_ds: xr.Dataset,
    gfs_ds: xr.Dataset,
    match_row: pd.Series,
    thresholds: pd.DataFrame,
    frame_dir: Path,
    temp_vmin: float,
    temp_vmax: float,
) -> Path:
    frame_idx = int(match_row["frame_idx"])
    rgb_time = pd.Timestamp(match_row["rgb_time"])
    acm_time = pd.Timestamp(match_row["acm_time"])
    gfs_time = rgb_time.floor("h")

    rgb_slice = rgb_ds.isel(t=frame_idx)[["red", "green", "blue"]].load()
    red = rgb_slice["red"].values
    green = rgb_slice["green"].values
    blue = rgb_slice["blue"].values
    rgb_image = as_rgb_image(red, green, blue)

    with xr.open_dataset(match_row["path"]) as acm_ds:
        bcm = acm_ds["BCM"].values
        acm_lon = acm_ds["longitude"].values
        acm_lat = acm_ds["latitude"].values

    gfs_hour = gfs_ds.sel(valid_time=np.datetime64(gfs_time))
    gfs_temp_rgb = gfs_temperature_on_rgb_grid(gfs_hour, rgb_ds)
    rgb_mask = threshold_cloud_mask(
        red=red,
        green=green,
        blue=blue,
        surface_temperature_c=gfs_temp_rgb,
        thresholds=thresholds,
    )

    fig = plt.figure(figsize=(11, 12), dpi=125, constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05])
    ax_rgb = fig.add_subplot(gs[0, 0])
    ax_bcm = fig.add_subplot(gs[0, 1])
    ax_rgb_mask = fig.add_subplot(gs[1, 0])
    ax_snow = fig.add_subplot(gs[1, 1])
    ax_temp = fig.add_subplot(gs[2, 0])
    ax_blank = fig.add_subplot(gs[2, 1])

    imshow_latlon(
        ax_rgb,
        rgb_image,
        rgb_ds["longitude"].values,
        rgb_ds["latitude"].values,
        title="RGB composite",
    )
    ax_rgb.set_title("RGB composite")
    label_geo_axis(ax_rgb)

    bcm_im = imshow_latlon(
        ax_bcm,
        bcm,
        acm_lon,
        acm_lat,
        cmap=MASK_CLEAR_CLOUD_CMAP,
        norm=BCM_NORM,
    )
    ax_bcm.set_title(f"GOES BCM cloud mask\nACM nearest {acm_time:%H:%M:%S}Z")
    label_geo_axis(ax_bcm)
    cbar = fig.colorbar(bcm_im, ax=ax_bcm, ticks=[0, 1], shrink=0.78)
    cbar.ax.set_yticklabels(["clear", "cloudy"])

    rgb_mask_im = imshow_latlon(
        ax_rgb_mask,
        rgb_mask,
        rgb_ds["longitude"].values,
        rgb_ds["latitude"].values,
        cmap=MASK_CLEAR_CLOUD_CMAP,
        norm=MASK_NORM,
    )
    ax_rgb_mask.set_title("RGB cloud mask\n10C temp-bin thresholds")
    label_geo_axis(ax_rgb_mask)
    cbar = fig.colorbar(rgb_mask_im, ax=ax_rgb_mask, ticks=[0, 1], shrink=0.78)
    cbar.ax.set_yticklabels(["clear", "cloudy"])

    snow_im = imshow_latlon(
        ax_snow,
        gfs_hour["snow_mask"].astype("int8").values,
        gfs_hour["longitude"].values,
        gfs_hour["latitude"].values,
        cmap=MASK_CLEAR_CLOUD_CMAP,
        norm=MASK_NORM,
    )
    ax_snow.set_title(f"GFS snow mask\nusing {gfs_time:%H}:00Z")
    label_geo_axis(ax_snow)
    cbar = fig.colorbar(snow_im, ax=ax_snow, ticks=[0, 1], shrink=0.78)
    cbar.ax.set_yticklabels(["no snow", "snow"])

    temp_im = imshow_latlon(
        ax_temp,
        gfs_hour["surface_temperature_c"].values,
        gfs_hour["longitude"].values,
        gfs_hour["latitude"].values,
        cmap="coolwarm",
        vmin=temp_vmin,
        vmax=temp_vmax,
    )
    ax_temp.set_title(f"GFS surface temperature, degC | using {gfs_time:%H}:00Z")
    label_geo_axis(ax_temp)
    fig.colorbar(temp_im, ax=ax_temp, label="degC", shrink=0.78)
    ax_blank.axis("off")

    fig.suptitle(f"GOES/GFS Colorado domain | {rgb_time:%Y-%m-%d %H:%M}Z")

    frame_path = frame_dir / f"five_panel_{frame_idx:03d}_{rgb_time:%H%M}.png"
    fig.savefig(frame_path)
    plt.close(fig)
    return frame_path


def main() -> int:
    args = parse_args()
    date = args.date
    day = pd.Timestamp(date)
    start_time = day + pd.Timedelta(hours=16)
    end_time = day + pd.Timedelta(days=1)
    paths = date_paths(date)

    ensure_inputs(paths)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths["frames"].mkdir(parents=True, exist_ok=True)

    acm_df = acm_time_table(paths["acm"])
    thresholds = load_rgb_thresholds()
    frame_paths: list[Path] = []

    with xr.open_dataset(paths["rgb"]) as rgb_ds, xr.open_dataset(paths["gfs"]) as gfs_ds:
        matches = matched_goes_times(rgb_ds["t"].values, acm_df)
        matches = matches.loc[
            (matches["rgb_time"] >= start_time) & (matches["rgb_time"] < end_time)
        ].reset_index(drop=True)
        if matches.empty:
            raise RuntimeError(f"No matched GOES frames found from {start_time} to {end_time}")
        temp_vmin = float(gfs_ds["surface_temperature_c"].min())
        temp_vmax = float(gfs_ds["surface_temperature_c"].max())

        print(f"RGB frames in file: {rgb_ds.sizes['t']}")
        print(f"Matched frames to render: {len(matches)} ({start_time:%HZ} to {end_time:%HZ})")
        print(f"ACM frames: {len(acm_df)}")
        print(f"GFS hourly frames: {gfs_ds.sizes['valid_time']}")
        print(f"RGB mask thresholds: {THRESHOLD_PATH}")
        print(f"Writing frames to {paths['frames']}")

        for row_idx, row in matches.iterrows():
            frame_path = render_frame(
                rgb_ds=rgb_ds,
                gfs_ds=gfs_ds,
                match_row=row,
                thresholds=thresholds,
                frame_dir=paths["frames"],
                temp_vmin=temp_vmin,
                temp_vmax=temp_vmax,
            )
            frame_paths.append(frame_path)
            if (row_idx + 1) % 24 == 0:
                print(f"Rendered {row_idx + 1}/{len(matches)} frames")

    with imageio.get_writer(paths["gif"], mode="I", duration=0.18, loop=0) as writer:
        for frame_path in frame_paths:
            writer.append_data(imageio.imread(frame_path))

    print(f"Wrote {paths['gif']}")
    print(f"Frames: {len(frame_paths)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
