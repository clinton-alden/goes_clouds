#!/usr/bin/env python3
"""Build Gothic cloud masks from RGB files using ERA5-temperature-selected thresholds."""

from __future__ import annotations

import argparse
import calendar
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create cloud binary masks at native Gothic RGB grid using 10C threshold bins "
            "selected by ERA5 2m temperature at each timestep."
        )
    )
    parser.add_argument(
        "--rgb-dir",
        default="/glade/derecho/scratch/cdalden/gothic/goes16/rgb_composite",
        help="Directory containing daily Gothic RGB NetCDF files",
    )
    parser.add_argument(
        "--era5-dir",
        default="/glade/derecho/scratch/cdalden/gothic/era5/t2m_hourly",
        help="Directory with monthly ERA5 t2m files from download_era5_gothic_t2m.py",
    )
    parser.add_argument(
        "--threshold-csv",
        default="/glade/u/home/cdalden/goes_work/analysis/output_12_rgb_threshold_transfer/gothic_temp_bin_rgb_thresholds_10c.csv",
        help="CSV of trained 10C thresholds",
    )
    parser.add_argument(
        "--output-dir",
        default="/glade/derecho/scratch/cdalden/gothic/goes16/cloud_mask_tempbin_10c",
        help="Output directory for daily cloud mask NetCDF files",
    )
    parser.add_argument("--year", type=int, default=2022)
    parser.add_argument("--month", type=int, default=7)
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing output files"
    )
    return parser


def load_thresholds(path: Path) -> pd.DataFrame:
    thresholds = pd.read_csv(path)
    thresholds = thresholds.loc[thresholds["status"] == "trained"].copy()
    thresholds = thresholds.sort_values("temp_left_c").reset_index(drop=True)
    if thresholds.empty:
        raise ValueError(f"No trained thresholds found in {path}")
    return thresholds


def choose_bin(temp_c: np.ndarray, left_edges: np.ndarray, right_edges: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(right_edges, temp_c, side="right")
    idx = np.clip(idx, 0, len(left_edges) - 1)
    idx[temp_c < left_edges[0]] = 0
    idx[temp_c >= right_edges[-1]] = len(left_edges) - 1
    return idx


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


def load_era5_month(era5_dir: Path, year: int, month: int) -> xr.DataArray:
    era5_path = era5_dir / f"era5_t2m_gothic_{year}{month:02d}.nc"
    if not era5_path.exists():
        raise FileNotFoundError(f"Missing ERA5 file: {era5_path}")

    ds = xr.open_dataset(era5_path)
    if "t2m" not in ds:
        raise KeyError(f"Expected variable 't2m' in {era5_path}")

    # Domain mean temp (K -> C) since thresholding model is keyed to scalar air temp.
    return ds["t2m"].mean(dim=["latitude", "longitude"]) - 273.15


def resolve_time_coord_name(da: xr.DataArray) -> str:
    if "time" in da.coords:
        return "time"
    if "valid_time" in da.coords:
        return "valid_time"
    for dim_name in da.dims:
        if "time" in dim_name:
            return dim_name
    raise KeyError("Could not find time coordinate in ERA5 temperature series")


def build_day_mask(
    rgb_path: Path,
    t2m_series: xr.DataArray,
    thresholds: pd.DataFrame,
    out_path: Path,
) -> None:
    ds = xr.open_dataset(rgb_path)

    goes_times = pd.DatetimeIndex(pd.to_datetime(ds["t"].values))
    era5_time_name = resolve_time_coord_name(t2m_series)
    era5_times = pd.DatetimeIndex(pd.to_datetime(t2m_series[era5_time_name].values))
    temp_series = pd.Series(t2m_series.values.astype(float), index=era5_times)

    # Time interpolation avoids abrupt hour-step jumps at 5-minute GOES cadence.
    temp_at_goes = (
        temp_series.reindex(temp_series.index.union(goes_times))
        .sort_index()
        .interpolate(method="time")
        .reindex(goes_times)
    )
    if temp_at_goes.isna().any():
        temp_at_goes = temp_at_goes.ffill().bfill()

    left_edges = thresholds["temp_left_c"].to_numpy()
    right_edges = thresholds["temp_right_c"].to_numpy()
    bin_idx = choose_bin(temp_at_goes.values.astype(float), left_edges, right_edges)

    red = ds["red"].values
    green = ds["green"].values
    blue = ds["blue"].values
    cloud_mask = np.zeros(red.shape, dtype=np.uint8)

    for idx in np.unique(bin_idx):
        row = thresholds.iloc[int(idx)]
        sel = bin_idx == idx

        c1 = band_condition(red[sel], float(row["red_threshold"]), str(row["red_direction"]))
        c2 = band_condition(
            green[sel], float(row["green_threshold"]), str(row["green_direction"])
        )
        c3 = band_condition(blue[sel], float(row["blue_threshold"]), str(row["blue_direction"]))
        cloud_mask[sel] = apply_rule(c1, c2, c3, str(row["rule"])).astype(np.uint8)

    out_ds = xr.Dataset(
        data_vars={
            "cloud_binary": (("t", "y", "x"), cloud_mask),
            "air_temp_c": (("t",), temp_at_goes.values.astype(np.float32)),
            "temp_bin_index": (("t",), bin_idx.astype(np.int16)),
        },
        coords={"t": ds["t"], "y": ds["y"], "x": ds["x"]},
        attrs={
            "title": "Gothic cloud mask from RGB thresholds selected by ERA5 2m temperature",
            "rgb_source": str(rgb_path),
        },
    )
    out_ds["cloud_binary"].attrs["long_name"] = "cloud mask (0=clear, 1=cloudy)"
    out_ds["air_temp_c"].attrs["long_name"] = "ERA5 domain-mean 2m air temperature"
    out_ds["air_temp_c"].attrs["units"] = "degC"
    out_ds["temp_bin_index"].attrs["long_name"] = "row index in threshold CSV used per timestep"

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

    rgb_dir = Path(args.rgb_dir)
    era5_dir = Path(args.era5_dir)
    threshold_csv = Path(args.threshold_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    thresholds = load_thresholds(threshold_csv)
    t2m_series = load_era5_month(era5_dir, args.year, args.month)

    ndays = calendar.monthrange(args.year, args.month)[1]
    print(
        f"Building cloud masks for {args.year}-{args.month:02d} ({ndays} days) from {rgb_dir}"
    )
    print(f"ERA5 source month file: {era5_dir / f'era5_t2m_gothic_{args.year}{args.month:02d}.nc'}")
    print(f"Output dir: {output_dir}")

    for day in range(1, ndays + 1):
        date = f"{args.year}{args.month:02d}{day:02d}"
        rgb_path = rgb_dir / f"goes16_C02_C05_C13_rgb_gothic_{date}.nc"
        out_path = output_dir / f"goes16_cloud_binary_tempbin10c_gothic_{date}.nc"

        if not rgb_path.exists():
            print(f"[skip] missing RGB file {rgb_path}")
            continue
        if out_path.exists() and not args.overwrite:
            print(f"[skip] already exists {out_path}")
            continue

        print(f"[build] {date}")
        build_day_mask(rgb_path, t2m_series, thresholds, out_path)

    print("Cloud mask generation complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
