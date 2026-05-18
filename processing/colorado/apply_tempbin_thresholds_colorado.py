#!/usr/bin/env python3
"""Apply temperature-bin RGB thresholds to a Colorado RGB file using ERA5-Land."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import xarray as xr

START_HOUR_UTC = 14
END_HOUR_UTC = 24
DEFAULT_ERA5_PADDING_DEG = 0.2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download ERA5-Land 2m air temperature for the RGB domain if needed, "
            "apply temperature-selected RGB thresholds, and render a binary mask GIF."
        )
    )
    parser.add_argument(
        "--rgb-file",
        default=(
            "/glade/u/home/cdalden/scratch/colorado/goes16/rgb_composite/"
            "goes16_C02_C05_C13_rgb_colorado_20220325.nc"
        ),
        help="Input Colorado RGB NetCDF file",
    )
    parser.add_argument(
        "--threshold-csv",
        default=(
            "/glade/u/home/cdalden/goes_work/analysis/output_12_rgb_threshold_transfer/"
            "gothic_temp_bin_rgb_thresholds_10c.csv"
        ),
        help="CSV containing the temperature-bin RGB thresholds",
    )
    parser.add_argument(
        "--era5-dir",
        default="/glade/derecho/scratch/cdalden/tmp/colorado/era5_land/t2m_hourly",
        help="Directory for monthly ERA5-Land files",
    )
    parser.add_argument(
        "--mask-dir",
        default="/glade/derecho/scratch/cdalden/tmp/colorado/goes16/cloud_mask_tempbin_10c",
        help="Directory for output binary mask NetCDF files",
    )
    parser.add_argument(
        "--gif-dir",
        default="/glade/derecho/scratch/cdalden/tmp/colorado/goes16/gif_loops_tempbin_10c",
        help="Directory for output GIF files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing ERA5, mask, and GIF outputs",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Do not download ERA5-Land even if the monthly file is missing",
    )
    parser.add_argument(
        "--frame-duration",
        type=float,
        default=0.15,
        help="GIF frame duration in seconds",
    )
    parser.add_argument(
        "--era5-padding-deg",
        type=float,
        default=DEFAULT_ERA5_PADDING_DEG,
        help="Padding added around RGB bounds when downloading ERA5-Land",
    )
    return parser


def load_thresholds(path: Path) -> pd.DataFrame:
    thresholds = pd.read_csv(path)
    thresholds = thresholds.loc[thresholds["status"] == "trained"].copy()
    thresholds = thresholds.sort_values("temp_left_c").reset_index(drop=True)
    if thresholds.empty:
        raise ValueError(f"No trained thresholds found in {path}")
    return thresholds


def infer_output_paths(rgb_path: Path, mask_dir: Path, gif_dir: Path) -> tuple[Path, Path]:
    stem = rgb_path.stem
    mask_path = mask_dir / f"{stem}_cloud_binary_tempbin10c.nc"
    gif_path = gif_dir / f"{stem}_cloud_binary_tempbin10c.gif"
    return mask_path, gif_path


def parse_rgb_date(rgb_path: Path) -> tuple[int, int, int]:
    token = rgb_path.stem.rsplit("_", 1)[-1]
    if len(token) != 8 or not token.isdigit():
        raise ValueError(f"Could not parse YYYYMMDD from RGB filename: {rgb_path.name}")
    return int(token[:4]), int(token[4:6]), int(token[6:8])


def derive_bounds(ds: xr.Dataset, padding_deg: float = 0.0) -> list[float]:
    lon = np.asarray(ds["longitude"].values, dtype=float)
    lat = np.asarray(ds["latitude"].values, dtype=float)
    return [
        float(lat.max() + padding_deg),
        float(lon.min() - padding_deg),
        float(lat.min() - padding_deg),
        float(lon.max() + padding_deg),
    ]


def ensure_era5_land_month(
    rgb_path: Path,
    era5_dir: Path,
    overwrite_download: bool,
    skip_download: bool,
    padding_deg: float,
) -> Path:
    year, month, _ = parse_rgb_date(rgb_path)
    out_path = era5_dir / f"era5land_t2m_colorado_{year}{month:02d}.nc"
    if out_path.exists() and (skip_download or not overwrite_download):
        return out_path
    if skip_download and not out_path.exists():
        raise FileNotFoundError(
            f"Missing ERA5-Land file {out_path} and --skip-download was requested"
        )

    try:
        import cdsapi
    except ImportError as exc:
        raise RuntimeError("cdsapi is required to download ERA5-Land data") from exc

    era5_dir.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(rgb_path) as ds:
        area = derive_bounds(ds, padding_deg=padding_deg)

    days = [f"{day:02d}" for day in range(1, 32)]
    hours = [f"{hour:02d}:00" for hour in range(24)]
    request = {
        "variable": ["2m_temperature"],
        "year": [f"{year}"],
        "month": [f"{month:02d}"],
        "day": days,
        "time": hours,
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": area,
    }

    client = cdsapi.Client()
    client.retrieve("reanalysis-era5-land", request, str(out_path))
    return out_path


def resolve_time_coord_name(da: xr.DataArray) -> str:
    if "time" in da.coords:
        return "time"
    if "valid_time" in da.coords:
        return "valid_time"
    for dim_name in da.dims:
        if "time" in dim_name:
            return dim_name
    raise KeyError("Could not find a time coordinate in the ERA5-Land data")


def load_era5_temp_field(era5_path: Path) -> xr.DataArray:
    ds = xr.open_dataset(era5_path)
    if "t2m" not in ds:
        raise KeyError(f"Expected variable 't2m' in {era5_path}")
    return ds["t2m"] - 273.15


def choose_bin(temp_c: np.ndarray, left_edges: np.ndarray, right_edges: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(right_edges, temp_c, side="right")
    idx = np.clip(idx, 0, len(left_edges) - 1)
    idx[temp_c < left_edges[0]] = 0
    idx[temp_c >= right_edges[-1]] = len(left_edges) - 1
    return idx


def select_target_hours(ds: xr.Dataset) -> xr.Dataset:
    times = pd.DatetimeIndex(pd.to_datetime(ds["t"].values))
    hour_float = times.hour + times.minute / 60.0 + times.second / 3600.0
    keep = (hour_float >= START_HOUR_UTC) & (hour_float < END_HOUR_UTC)
    if not np.any(keep):
        raise ValueError("No GOES timesteps found in the requested 14Z-00Z window")
    return ds.isel(t=np.where(keep)[0])


def band_condition(values: np.ndarray, threshold: float, direction: str) -> np.ndarray:
    if direction == ">":
        return values > threshold
    if direction == "<=":
        return values <= threshold
    raise ValueError(f"Unsupported threshold direction: {direction}")


def apply_rule(c1: np.ndarray, c2: np.ndarray, c3: np.ndarray, rule: str) -> np.ndarray:
    if rule == "union":
        return c1 | c2 | c3
    if rule == "intersection":
        return c1 & c2 & c3
    if rule == "majority":
        return (c1.astype(np.uint8) + c2.astype(np.uint8) + c3.astype(np.uint8)) >= 2
    raise ValueError(f"Unsupported combine rule: {rule}")


def interpolate_temp_to_goes_grid(
    t2m_field: xr.DataArray,
    goes_times: pd.DatetimeIndex,
    goes_lat: np.ndarray,
    goes_lon: np.ndarray,
) -> xr.DataArray:
    era5_time_name = resolve_time_coord_name(t2m_field)
    target = {
        era5_time_name: xr.DataArray(goes_times, dims=("t",)),
        "latitude": xr.DataArray(goes_lat, dims=("latitude",)),
        "longitude": xr.DataArray(goes_lon, dims=("longitude",)),
    }

    # Use linear interpolation in time and space, then patch edge NaNs with nearest-neighbor.
    linear = t2m_field.interp(target, method="linear")
    nearest = t2m_field.interp(target, method="nearest")
    temp_on_goes = linear.where(np.isfinite(linear), nearest)

    if temp_on_goes.isnull().any():
        temp_on_goes = temp_on_goes.ffill("t").bfill("t")
    return temp_on_goes.rename({era5_time_name: "t"}) if era5_time_name in temp_on_goes.dims else temp_on_goes


def build_cloud_mask(
    rgb_path: Path,
    era5_path: Path,
    threshold_csv: Path,
    mask_path: Path,
) -> tuple[pd.DataFrame, pd.Series]:
    thresholds = load_thresholds(threshold_csv)
    t2m_field = load_era5_temp_field(era5_path)

    with xr.open_dataset(rgb_path) as ds:
        ds = select_target_hours(ds)
        goes_times = pd.DatetimeIndex(pd.to_datetime(ds["t"].values))
        goes_lat = np.asarray(ds["latitude"].values, dtype=np.float64)
        goes_lon = np.asarray(ds["longitude"].values, dtype=np.float64)
        temp_at_goes = interpolate_temp_to_goes_grid(
            t2m_field=t2m_field,
            goes_times=goes_times,
            goes_lat=goes_lat,
            goes_lon=goes_lon,
        )

        left_edges = thresholds["temp_left_c"].to_numpy()
        right_edges = thresholds["temp_right_c"].to_numpy()
        temp_values = np.asarray(temp_at_goes.values, dtype=np.float32)
        bin_idx = choose_bin(temp_values.astype(float), left_edges, right_edges)

        red = np.asarray(ds["red"].values, dtype=np.float32)
        green = np.asarray(ds["green"].values, dtype=np.float32)
        blue = np.asarray(ds["blue"].values, dtype=np.float32)
        cloud_mask = np.zeros(red.shape, dtype=np.uint8)

        for idx in np.unique(bin_idx):
            row = thresholds.iloc[int(idx)]
            sel = bin_idx == idx
            c1 = band_condition(red[sel], float(row["red_threshold"]), str(row["red_direction"]))
            c2 = band_condition(
                green[sel], float(row["green_threshold"]), str(row["green_direction"])
            )
            c3 = band_condition(
                blue[sel], float(row["blue_threshold"]), str(row["blue_direction"])
            )
            cloud_mask[sel] = apply_rule(c1, c2, c3, str(row["rule"])).astype(np.uint8)

        out_ds = xr.Dataset(
            data_vars={
                "cloud_binary": (("t", "latitude", "longitude"), cloud_mask),
                "air_temp_c": (("t", "latitude", "longitude"), temp_values.astype(np.float32)),
                "air_temp_domain_mean_c": (
                    ("t",),
                    temp_values.mean(axis=(1, 2)).astype(np.float32),
                ),
                "temp_bin_index": (("t", "latitude", "longitude"), bin_idx.astype(np.int16)),
            },
            coords={
                "t": ds["t"],
                "latitude": ds["latitude"],
                "longitude": ds["longitude"],
            },
            attrs={
                "title": "Colorado cloud mask from Gothic RGB temperature-bin thresholds",
                "rgb_source": str(rgb_path),
                "era5_land_source": str(era5_path),
                "threshold_csv": str(threshold_csv),
            },
        )
        out_ds["cloud_binary"].attrs["long_name"] = "cloud mask (0=clear, 1=cloudy)"
        out_ds["air_temp_c"].attrs["long_name"] = "ERA5-Land 2m air temperature interpolated to GOES grid"
        out_ds["air_temp_c"].attrs["units"] = "degC"
        out_ds["air_temp_domain_mean_c"].attrs["long_name"] = (
            "domain-mean ERA5-Land 2m air temperature after interpolation to GOES grid"
        )
        out_ds["air_temp_domain_mean_c"].attrs["units"] = "degC"
        out_ds["temp_bin_index"].attrs["long_name"] = (
            "row index in threshold CSV used per pixel and timestep"
        )

        encoding = {
            "cloud_binary": {"zlib": True, "complevel": 4, "dtype": "uint8"},
            "air_temp_c": {"zlib": True, "complevel": 4, "dtype": "float32"},
            "air_temp_domain_mean_c": {"zlib": True, "complevel": 4, "dtype": "float32"},
            "temp_bin_index": {"zlib": True, "complevel": 4, "dtype": "int16"},
        }
        out_ds.to_netcdf(mask_path, encoding=encoding)
        out_ds.close()

    applied_bins = np.unique(bin_idx)
    return thresholds.iloc[applied_bins].copy(), pd.Series(
        temp_values.mean(axis=(1, 2)),
        index=goes_times,
    )


def render_gif(mask_path: Path, gif_path: Path, frame_duration: float) -> None:
    frames_dir = gif_path.parent / f"{gif_path.stem}_frames"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    frame_paths: list[Path] = []
    cloud_cmap = ListedColormap(["#2563eb", "#ffffff"])
    with xr.open_dataset(mask_path) as ds:
        mask = ds["cloud_binary"].values
        times = pd.to_datetime(ds["t"].values)
        lon = np.asarray(ds["longitude"].values, dtype=float)
        lat = np.asarray(ds["latitude"].values, dtype=float)
        extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]

        for i, timestamp in enumerate(times):
            fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
            ax.imshow(
                mask[i],
                origin="lower",
                cmap=cloud_cmap,
                vmin=0,
                vmax=1,
                extent=extent,
                interpolation="nearest",
                aspect="auto",
            )
            ax.set_title(f"RGB Cloud Mask\n{timestamp.strftime('%Y-%m-%d %H:%M UTC')}")
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

            frame_path = frames_dir / f"frame_{i:04d}.png"
            fig.savefig(frame_path, dpi=120)
            plt.close(fig)
            frame_paths.append(frame_path)

    with imageio.get_writer(gif_path, mode="I", duration=frame_duration) as writer:
        for frame_path in frame_paths:
            writer.append_data(imageio.imread(frame_path))
    shutil.rmtree(frames_dir, ignore_errors=True)


def main() -> int:
    args = build_parser().parse_args()

    rgb_path = Path(args.rgb_file)
    threshold_csv = Path(args.threshold_csv)
    era5_dir = Path(args.era5_dir)
    mask_dir = Path(args.mask_dir)
    gif_dir = Path(args.gif_dir)

    if not rgb_path.exists():
        raise FileNotFoundError(f"Missing RGB file: {rgb_path}")
    if not threshold_csv.exists():
        raise FileNotFoundError(f"Missing threshold CSV: {threshold_csv}")

    mask_dir.mkdir(parents=True, exist_ok=True)
    gif_dir.mkdir(parents=True, exist_ok=True)
    mask_path, gif_path = infer_output_paths(rgb_path, mask_dir, gif_dir)

    if args.overwrite or not mask_path.exists():
        era5_path = ensure_era5_land_month(
            rgb_path=rgb_path,
            era5_dir=era5_dir,
            overwrite_download=False,
            skip_download=args.skip_download,
            padding_deg=args.era5_padding_deg,
        )
        thresholds, temp_at_goes = build_cloud_mask(
            rgb_path=rgb_path,
            era5_path=era5_path,
            threshold_csv=threshold_csv,
            mask_path=mask_path,
        )
        print(f"Wrote cloud mask: {mask_path}")
        print(f"ERA5-Land monthly file: {era5_path}")
        print(
            "Applied temperature bins: "
            + ", ".join(thresholds["temp_bin"].astype(str).tolist())
        )
    else:
        print(f"Cloud mask exists, reusing: {mask_path}")

    if args.overwrite or not gif_path.exists():
        render_gif(mask_path=mask_path, gif_path=gif_path, frame_duration=args.frame_duration)
        print(f"Wrote GIF: {gif_path}")
    else:
        print(f"GIF exists, reusing: {gif_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
