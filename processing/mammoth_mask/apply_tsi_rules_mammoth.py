#!/usr/bin/env python3
"""Apply Notebook 8d TSI RGB/ERA5-Land rules to native Mammoth RGB pixels."""

from __future__ import annotations

import argparse
import calendar
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS, Transformer


GOES18_PROJ = CRS.from_proj4(
    "+proj=geos +h=35786023 +lon_0=-137 +sweep=x "
    "+a=6378137 +b=6356752.31414 +units=m +no_defs"
)
COND_RE = re.compile(r"(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)")
FILL_VALUE = np.uint8(255)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, required=True)
    p.add_argument(
        "--rgb-dir",
        type=Path,
        default=Path("/glade/derecho/scratch/cdalden/mammoth/goes18/rgb_composite"),
    )
    p.add_argument(
        "--era5-dir",
        type=Path,
        default=Path("/glade/derecho/scratch/cdalden/mammoth/era5_land/t2m_hourly"),
    )
    p.add_argument(
        "--rules-csv",
        type=Path,
        default=Path("/glade/u/home/cdalden/goes_work/analysis/tsi_rgb_decision_tree_era5_t2m_10c_rules.csv"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/glade/derecho/scratch/cdalden/mammoth/rgb_mask"),
    )
    p.add_argument("--overwrite", action="store_true")
    return p


def rgb_path(rgb_dir: Path, day: pd.Timestamp) -> Path:
    return rgb_dir / f"goes18_C02_C05_C13_rgb_mammoth_{day:%Y%m%d}.nc"


def era5_path(era5_dir: Path, day: pd.Timestamp) -> Path:
    return era5_dir / f"era5land_t2m_mammoth_{day:%Y%m}.nc"


def geolocate(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xx, yy = np.meshgrid(x * 35786023.0, y * 35786023.0)
    transformer = Transformer.from_crs(GOES18_PROJ, 4326, always_xy=True)
    lon, lat = transformer.transform(xx, yy)
    return np.asarray(lat, np.float32), np.asarray(lon, np.float32)


def load_rules(path: Path) -> tuple[dict[str, list[str]], dict[str, float]]:
    rules = pd.read_csv(path)
    rules = rules.loc[rules["predicted_class"].eq(1)]
    cloudy = {
        str(label): group["rule"].astype(str).tolist()
        for label, group in rules.groupby("temp_bin", sort=False)
    }
    if not cloudy:
        raise ValueError(f"No cloudy leaves in {path}")
    midpoints = {}
    for label in cloudy:
        left, right = re.findall(r"-?\d+(?:\.\d+)?", label)
        midpoints[label] = (float(left) + float(right)) / 2.0
    return cloudy, midpoints


def route_bins(temp: np.ndarray, midpoints: dict[str, float]) -> np.ndarray:
    labels = np.asarray(list(midpoints), dtype=object)
    mids = np.asarray([midpoints[label] for label in labels])
    raw_left = np.floor(temp / 10.0) * 10.0
    raw = np.full(temp.shape, "", dtype=object)
    for left in np.unique(raw_left[np.isfinite(raw_left)]):
        raw[raw_left == left] = f"[{int(left)}, {int(left + 10)})"
    unknown = ~np.isin(raw, labels) & np.isfinite(temp)
    if unknown.any():
        nearest = np.abs(temp[unknown, None] - mids[None, :]).argmin(axis=1)
        raw[unknown] = labels[nearest]
    raw[~np.isfinite(temp)] = ""
    return raw


def apply_condition(values: dict[str, np.ndarray], rule: str) -> np.ndarray:
    result = np.ones(values["red"].shape, dtype=bool)
    conditions = COND_RE.findall(rule)
    if not conditions:
        raise ValueError(f"Could not parse rule: {rule}")
    for feature, operator, threshold_text in conditions:
        threshold = float(threshold_text)
        if operator == "<=":
            result &= values[feature] <= threshold
        else:
            result &= values[feature] > threshold
    return result


def classify(red, green, blue, bins, cloudy_rules) -> np.ndarray:
    out = np.zeros(red.shape, dtype=np.uint8)
    for label in np.unique(bins):
        if not label:
            continue
        selected = bins == label
        values = {"red": red[selected], "green": green[selected], "blue": blue[selected]}
        cloudy = np.zeros(values["red"].shape, dtype=bool)
        for rule in cloudy_rules[str(label)]:
            cloudy |= apply_condition(values, rule)
        out[selected] = cloudy.astype(np.uint8)
    invalid = ~np.isfinite(red) | ~np.isfinite(green) | ~np.isfinite(blue) | (bins == "")
    out[invalid] = FILL_VALUE
    return out


def time_name(da: xr.DataArray) -> str:
    for name in ("valid_time", "time"):
        if name in da.dims:
            return name
    raise KeyError("ERA5 t2m has neither valid_time nor time dimension")


def load_rgb_window(rgb_dir: Path, day: pd.Timestamp) -> xr.Dataset | None:
    start = day + pd.Timedelta(hours=15)
    end = day + pd.Timedelta(days=1, hours=1)
    primary = rgb_path(rgb_dir, day)
    if not primary.exists() or primary.stat().st_size == 0:
        return None
    pieces = []
    for source_day in (day, day + pd.Timedelta(days=1)):
        path = rgb_path(rgb_dir, source_day)
        if not path.exists() or path.stat().st_size == 0:
            continue
        with xr.open_dataset(path) as ds:
            times = pd.DatetimeIndex(pd.to_datetime(ds["t"].values))
            keep = np.flatnonzero((times >= start) & (times < end))
            if keep.size:
                pieces.append(ds[["red", "green", "blue"]].isel(t=keep).load())
    if not pieces:
        return None
    return xr.concat(pieces, dim="t").sortby("t")


def load_temperature(
    era5_dir: Path, day: pd.Timestamp, times: xr.DataArray, lat: np.ndarray, lon: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    paths = []
    for timestamp in pd.DatetimeIndex(pd.to_datetime(times.values)):
        path = era5_path(era5_dir, timestamp)
        if path not in paths:
            paths.append(path)
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing ERA5-Land file(s): " + ", ".join(missing))
    arrays = []
    for path in paths:
        with xr.open_dataset(path) as ds:
            arrays.append(ds["t2m"].load())
    t2m = xr.concat(arrays, dim=time_name(arrays[0])).sortby(time_name(arrays[0]))
    tdim = time_name(t2m)
    if float(t2m.longitude.max()) > 180:
        target_lon = np.mod(lon, 360.0)
    else:
        target_lon = lon
    target_lat_da = xr.DataArray(lat, dims=("y", "x"))
    target_lon_da = xr.DataArray(target_lon, dims=("y", "x"))
    # Spatial selection is nearest ERA5-Land pixel for every GOES pixel.
    pixel_series = t2m.sel(
        latitude=target_lat_da, longitude=target_lon_da, method="nearest"
    )
    # Temperature evolves continuously between the hourly ERA5-Land analyses.
    selected = pixel_series.interp({tdim: xr.DataArray(times.values, dims="t")}) - 273.15
    selected = selected.transpose("t", "y", "x")
    era_lat = np.broadcast_to(np.asarray(pixel_series.latitude), lat.shape).astype(np.float32)
    era_lon = np.broadcast_to(np.asarray(pixel_series.longitude), lon.shape).astype(np.float32)
    return np.asarray(selected, np.float32), era_lat, era_lon


def build_day(args, day, cloudy_rules, midpoints) -> bool:
    ds = load_rgb_window(args.rgb_dir, day)
    if ds is None:
        print(f"[skip] {day:%Y-%m-%d}: no nonempty RGB data in 15Z-01Z window", flush=True)
        return False
    lat, lon = geolocate(np.asarray(ds.x), np.asarray(ds.y))
    temp, era_lat, era_lon = load_temperature(args.era5_dir, day, ds.t, lat, lon)
    bins = route_bins(temp, midpoints)
    red = np.asarray(ds.red, np.float32)
    green = np.asarray(ds.green, np.float32)
    blue = np.asarray(ds.blue, np.float32)
    mask = classify(red, green, blue, bins, cloudy_rules)

    labels = list(midpoints)
    bin_index = np.full(bins.shape, -1, dtype=np.int8)
    for index, label in enumerate(labels):
        bin_index[bins == label] = index
    out = xr.Dataset(
        data_vars={
            "cloud_binary": (("t", "y", "x"), mask),
            "air_temp_c": (("t", "y", "x"), temp),
            "temp_bin_index": (("t", "y", "x"), bin_index),
            "latitude": (("y", "x"), lat),
            "longitude": (("y", "x"), lon),
            "era5_land_nearest_latitude": (("y", "x"), era_lat),
            "era5_land_nearest_longitude": (("y", "x"), era_lon),
        },
        coords={"t": ds.t, "y": ds.y, "x": ds.x},
        attrs={
            "title": "Mammoth pixel-wise RGB binary cloud mask using Notebook 8d rules",
            "rules_csv": str(args.rules_csv),
            "era5_land_mapping": "nearest spatial pixel per GOES pixel; linear interpolation in time",
            "time_window_utc": "15:00 on nominal date through 00:59:59 on following date",
            "temp_bin_labels": ";".join(f"{i}:{label}" for i, label in enumerate(labels)),
        },
    )
    out.cloud_binary.attrs.update(long_name="binary cloud mask", flag_values=[0, 1, 255], flag_meanings="clear cloud missing")
    out.air_temp_c.attrs.update(long_name="ERA5-Land 2 m air temperature", units="degC")
    out.temp_bin_index.attrs.update(long_name="trained temperature-bin index; see global temp_bin_labels", missing_value=np.int8(-1))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = args.output_dir / f"goes18_rgb_cloud_binary_mammoth_{day:%Y%m%d}.nc"
    encoding = {
        "cloud_binary": {"zlib": True, "complevel": 4, "dtype": "uint8", "_FillValue": FILL_VALUE},
        "air_temp_c": {"zlib": True, "complevel": 4, "dtype": "float32"},
        "temp_bin_index": {"zlib": True, "complevel": 4, "dtype": "int8", "_FillValue": np.int8(-1)},
    }
    for name in ("latitude", "longitude", "era5_land_nearest_latitude", "era5_land_nearest_longitude"):
        encoding[name] = {"zlib": True, "complevel": 4, "dtype": "float32"}
    out.to_netcdf(path, encoding=encoding)
    ds.close()
    out.close()
    print(f"[wrote] {path} ({mask.shape[0]} times; cloud={np.mean(mask == 1):.3f})", flush=True)
    return True


def main() -> int:
    args = parser().parse_args()
    cloudy_rules, midpoints = load_rules(args.rules_csv)
    count = 0
    for day_number in range(1, calendar.monthrange(args.year, args.month)[1] + 1):
        day = pd.Timestamp(args.year, args.month, day_number)
        path = args.output_dir / f"goes18_rgb_cloud_binary_mammoth_{day:%Y%m%d}.nc"
        if path.exists() and not args.overwrite:
            print(f"[skip] exists: {path}", flush=True)
            continue
        count += build_day(args, day, cloudy_rules, midpoints)
    print(f"Completed {args.year}-{args.month:02d}: wrote {count} daily files", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
