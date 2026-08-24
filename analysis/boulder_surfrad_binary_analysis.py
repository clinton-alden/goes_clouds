#!/usr/bin/env python3
"""Build Table Mountain SURFRAD and Colorado GOES RGB binary time series.

The workflow mirrors the preprocessing used by Notebook 13:

* spatially average RGB over a 10-km east-west by 5-km north-south box
  extending south from the site;
* select the notebook-5b TSI-trained RGB thresholds with ERA5-Land air temperature;
* average one-minute SURFRAD observations onto the five-minute GOES clock;
* classify shortwave using the same unambiguous k_t cutoffs as Notebook 13.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import urllib.error
import urllib.request

import numpy as np
import pandas as pd
import xarray as xr

try:
    from mammoth_cues_binary_analysis import apply_rgb_rules
except ModuleNotFoundError:  # Imported as analysis.boulder_surfrad_binary_analysis.
    from analysis.mammoth_cues_binary_analysis import apply_rgb_rules


SITE_LAT = 40.12498
SITE_LON = -105.23680
SITE_ALTITUDE_KM = 1.689
RGB_DIR = Path("/glade/derecho/scratch/cdalden/colorado/goes16/rgb_composite")
SURFRAD_DIR = Path(
    "/glade/derecho/scratch/cdalden/surface_obs/colorado/surfrad_boulder"
)
THRESHOLDS = Path(__file__).resolve().parent / (
    "tsi_rgb_decision_tree_era5_t2m_10c_rules.csv"
)
ERA5_DIR = Path("/glade/derecho/scratch/cdalden/tmp/colorado/era5_land/t2m_hourly")
OUT_DIR = Path(__file__).resolve().parent / "output_13b_boulder_domain_cloud_binary"
SURFRAD_URL = "https://gml.noaa.gov/aftp/data/radiation/surfrad/Boulder_CO"

# Columns after the two-line daily-file header. Every measured value is followed
# by its integer QC flag.
SURFRAD_COLUMNS = [
    "year", "jday", "month", "day", "hour", "minute", "decimal_time", "zenith",
    "dw_solar", "qc_dw_solar", "uw_solar", "qc_uw_solar",
    "direct_n", "qc_direct_n", "diffuse", "qc_diffuse",
    "dw_ir", "qc_dw_ir", "dw_case_temp", "qc_dw_case_temp",
    "dw_dome_temp", "qc_dw_dome_temp", "uw_ir", "qc_uw_ir",
    "uw_case_temp", "qc_uw_case_temp", "uw_dome_temp", "qc_uw_dome_temp",
    "uvb", "qc_uvb", "par", "qc_par", "net_solar", "qc_net_solar",
    "net_ir", "qc_net_ir", "total_net", "qc_total_net",
    "air_temp_c", "qc_air_temp", "rh", "qc_rh", "wind_speed", "qc_wind_speed",
    "wind_direction", "qc_wind_direction", "pressure", "qc_pressure",
]


def rgb_date(path: Path) -> pd.Timestamp:
    match = re.search(r"(20\d{6})", path.name)
    if match is None:
        raise ValueError(f"Cannot find YYYYMMDD in {path.name}")
    return pd.to_datetime(match.group(1), format="%Y%m%d")


def inventory_rgb(rgb_dir: Path, max_files: int | None = None) -> list[Path]:
    paths = sorted(path for path in rgb_dir.glob("*.nc") if path.stat().st_size > 0)
    if max_files is not None:
        paths = paths[:max_files]
    if not paths:
        raise FileNotFoundError(f"No nonempty RGB NetCDF files found in {rgb_dir}")
    return paths


def surfrad_path(day: pd.Timestamp, cache_dir: Path) -> Path:
    return cache_dir / str(day.year) / f"tbl{day:%y}{day.dayofyear:03d}.dat"


def download_one_surfrad(
    day: pd.Timestamp, cache_dir: Path, overwrite: bool = False
) -> Path | None:
    target = surfrad_path(day, cache_dir)
    if target.exists() and target.stat().st_size > 0 and not overwrite:
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    url = f"{SURFRAD_URL}/{day.year}/{target.name}"
    temporary = target.with_suffix(".dat.part")
    try:
        urllib.request.urlretrieve(url, temporary)
        temporary.replace(target)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            print(f"SURFRAD file unavailable (404): {url}", flush=True)
            return None
        raise
    finally:
        temporary.unlink(missing_ok=True)
    return target


def download_surfrad(
    days: list[pd.Timestamp], cache_dir: Path, workers: int = 8, overwrite: bool = False
) -> list[Path]:
    with ThreadPoolExecutor(max_workers=workers) as pool:
        downloaded = list(
            pool.map(lambda day: download_one_surfrad(day, cache_dir, overwrite), days)
        )
    paths = [path for path in downloaded if path is not None]
    if not paths:
        raise FileNotFoundError("No matching SURFRAD files were available")
    print(f"SURFRAD inventory: {len(paths)}/{len(days)} requested daily files", flush=True)
    return paths


def read_surfrad(paths: list[Path]) -> pd.DataFrame:
    frames = []
    wanted = [
        "year", "month", "day", "hour", "minute", "zenith",
        "dw_solar", "qc_dw_solar", "direct_n", "qc_direct_n",
        "diffuse", "qc_diffuse", "air_temp_c", "qc_air_temp",
    ]
    for path in paths:
        frame = pd.read_csv(
            path, sep=r"\s+", skiprows=2, names=SURFRAD_COLUMNS,
            usecols=wanted, na_values=-9999.9,
        )
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    data["time"] = pd.to_datetime(
        data[["year", "month", "day", "hour", "minute"]]
    )
    cos_sza = np.clip(np.cos(np.deg2rad(data.zenith.to_numpy(float))), 0, 1)
    component_good = (
        data.qc_direct_n.eq(0) & data.qc_diffuse.eq(0)
        & data.direct_n.notna() & data.diffuse.notna()
    )
    global_good = data.qc_dw_solar.eq(0) & data.dw_solar.notna()
    reconstructed = data.direct_n * cos_sza + data.diffuse
    data["sw_obs"] = reconstructed.where(component_good)
    data.loc[data.sw_obs.isna() & global_good, "sw_obs"] = data.dw_solar
    data["air_temp_c"] = data.air_temp_c.where(
        data.qc_air_temp.eq(0) & data.air_temp_c.between(-60, 60)
    )
    # Input paths and records are already chronological. Avoid datetime argsort
    # here because the current Casper pandas/numpy build fails when sorting a
    # large DatetimeArray containing no missing values.
    data = data.dropna(subset=["time"]).drop_duplicates("time")
    indexed = data.set_index("time")
    data["sw_obs_5min"] = (
        indexed.sw_obs.rolling("5min", center=True, min_periods=3).mean().to_numpy()
    )
    data["air_temp_5min_c"] = (
        indexed.air_temp_c.rolling("5min", center=True, min_periods=3).mean().to_numpy()
    )
    return data[
        ["time", "zenith", "sw_obs", "sw_obs_5min", "air_temp_c", "air_temp_5min_c"]
    ].reset_index(drop=True)


def site_box_indices(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Return indices for a 10-km E/W box extending 5 km south of the site."""
    lat_min, lat_max = SITE_LAT - 5.0 / 111.32, SITE_LAT
    lon_half = 5.0 / (111.32 * np.cos(np.deg2rad(SITE_LAT)))
    lon_min, lon_max = SITE_LON - lon_half, SITE_LON + lon_half
    lat = np.asarray(ds.latitude, float)
    lon = np.asarray(ds.longitude, float)
    lat_idx = np.flatnonzero((lat >= lat_min) & (lat <= lat_max))
    lon_idx = np.flatnonzero((lon >= lon_min) & (lon <= lon_max))
    if not len(lat_idx) or not len(lon_idx):
        raise ValueError("The Table Mountain 10-km E/W by 5-km south box selected no GOES pixels")
    return lat_idx, lon_idx


def rebuild_rgb_means(paths: list[Path], output_csv: Path) -> pd.DataFrame:
    frames = []
    pixel_counts = set()
    for number, path in enumerate(paths, start=1):
        with xr.open_dataset(path) as ds:
            lat_idx, lon_idx = site_box_indices(ds)
            pixel_counts.add(len(lat_idx) * len(lon_idx))
            subset = ds[["red", "green", "blue"]].isel(
                latitude=lat_idx, longitude=lon_idx
            )
            # Some historical files carry a misleading date in their internal
            # time coordinate. Retain time of day but trust the daily filename.
            internal = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
            time = rgb_date(path) + (internal - internal.normalize())
            values = {
                band: np.asarray(
                    subset[band].mean(("latitude", "longitude"), skipna=True), float
                )
                for band in ("red", "green", "blue")
            }
        frames.append(pd.DataFrame({"time": time, **values, "rgb_file": str(path)}))
        if number % 50 == 0 or number == len(paths):
            print(f"Processed RGB files: {number}/{len(paths)}", flush=True)
    rgb = pd.concat(frames, ignore_index=True)
    rgb = (
        rgb.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["red", "green", "blue"])
        .sort_values("time").drop_duplicates("time").reset_index(drop=True)
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    rgb.to_csv(output_csv, index=False)
    pd.DataFrame(
        [{
            "source_file_count": len(paths),
            "start_date": min(map(rgb_date, paths)),
            "end_date": max(map(rgb_date, paths)),
            "selected_pixel_count_min": min(pixel_counts),
            "selected_pixel_count_max": max(pixel_counts),
            "site_lat": SITE_LAT, "site_lon": SITE_LON,
        }]
    ).to_csv(output_csv.with_name("boulder_rgb_rebuild_metadata.csv"), index=False)
    return rgb


def load_thresholds(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    table = table.loc[table.status.eq("trained")].sort_values("temp_left_c")
    if table.empty:
        raise ValueError(f"No trained thresholds in {path}")
    return table.reset_index(drop=True)


def threshold_rows(temperature: np.ndarray, rules: pd.DataFrame) -> np.ndarray:
    right = rules.temp_right_c.to_numpy(float)
    indices = np.searchsorted(right, temperature, side="right")
    return np.clip(indices, 0, len(rules) - 1)


def apply_rgb_thresholds(rgb: pd.DataFrame, rules: pd.DataFrame) -> pd.DataFrame:
    out = rgb.copy()
    indices = threshold_rows(out.air_temp_5min_c.to_numpy(float), rules)
    out["temp_bin"] = rules.iloc[indices].temp_bin.to_numpy()
    votes = np.zeros((len(out), 3), dtype=bool)
    for column, band in enumerate(("red", "green", "blue")):
        thresholds = rules.iloc[indices][f"{band}_threshold"].to_numpy(float)
        directions = rules.iloc[indices][f"{band}_direction"].to_numpy(str)
        values = out[band].to_numpy(float)
        votes[:, column] = np.where(directions == ">", values > thresholds, values <= thresholds)
    # The transferred 10-C table currently uses majority in every trained bin.
    out["rgb_cloud_binary"] = (votes.sum(axis=1) >= 2).astype(np.int8)
    return out


def load_era5_land_site_temperature(
    target_times: pd.Series,
    era5_dir: Path = ERA5_DIR,
    lat: float = SITE_LAT,
    lon: float = SITE_LON,
) -> np.ndarray:
    """Interpolate nearest-pixel hourly ERA5-Land t2m to GOES timestamps."""
    periods = pd.PeriodIndex(pd.to_datetime(target_times).dt.to_period("M").unique())
    hourly_times, hourly_values = [], []
    for period in periods:
        path = era5_dir / f"era5land_t2m_colorado_{period.year}{period.month:02d}.nc"
        if not path.exists():
            raise FileNotFoundError(f"Missing ERA5-Land month: {path}")
        with xr.open_dataset(path) as ds:
            time_name = "valid_time" if "valid_time" in ds.t2m.dims else "time"
            target_lon = lon % 360 if float(ds.longitude.max()) > 180 else lon
            field = ds.t2m.sel(latitude=lat, longitude=target_lon, method="nearest")
            hourly_times.append(pd.DatetimeIndex(pd.to_datetime(field[time_name].values)))
            hourly_values.append(np.asarray(field.values, float) - 273.15)
    source_time = pd.DatetimeIndex(np.concatenate([x.values for x in hourly_times]))
    source_temp = np.concatenate(hourly_values)
    target = pd.DatetimeIndex(pd.to_datetime(target_times))
    return np.interp(
        target.asi8.astype(float), source_time.asi8.astype(float), source_temp
    )


def match_and_classify(
    rgb: pd.DataFrame, observations: pd.DataFrame,
    cloudy_kt: float = 0.55, clear_kt: float = 0.85,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observations_for_merge = observations[
        ["time", "sw_obs_5min", "air_temp_5min_c"]
    ].rename(columns={"air_temp_5min_c": "surfrad_air_temp_5min_c"})
    matched = pd.merge_asof(
        rgb,
        observations_for_merge,
        on="time", direction="nearest", tolerance=pd.Timedelta("45s"),
    )
    matched = matched.dropna(subset=["era5_temp_c"]).reset_index(drop=True)
    rgb_binary = apply_rgb_rules(matched, THRESHOLDS)
    cos_sza = solar_cosine_zenith(rgb_binary.time)
    rgb_binary["cos_sza"] = cos_sza
    rgb_binary["sw_clear"] = 1361.0 * cos_sza * (0.78 + 0.04 * SITE_ALTITUDE_KM)
    rgb_binary["k_t"] = rgb_binary.sw_obs_5min / rgb_binary.sw_clear
    valid = (cos_sza >= 0.25) & rgb_binary.sw_obs_5min.notna()
    rgb_binary["sw_cloud_binary"] = pd.Series(pd.NA, index=rgb_binary.index, dtype="Int8")
    rgb_binary.loc[valid & (rgb_binary.k_t < cloudy_kt), "sw_cloud_binary"] = 1
    rgb_binary.loc[valid & (rgb_binary.k_t > clear_kt), "sw_cloud_binary"] = 0
    sw_binary = rgb_binary.loc[
        rgb_binary.sw_cloud_binary.notna(),
        ["time", "sw_obs_5min", "sw_clear", "cos_sza", "k_t", "sw_cloud_binary"],
    ].copy()
    sw_binary["sw_cloud_binary"] = sw_binary.sw_cloud_binary.astype(int)
    return rgb_binary, sw_binary


def solar_cosine_zenith(times_utc: pd.Series) -> np.ndarray:
    times = pd.DatetimeIndex(pd.to_datetime(times_utc))
    doy = times.dayofyear.to_numpy()
    hour = times.hour + times.minute / 60 + times.second / 3600
    gamma = 2 * np.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eqtime = 229.18 * (
        0.000075 + 0.001868*np.cos(gamma) - 0.032077*np.sin(gamma)
        - 0.014615*np.cos(2*gamma) - 0.040849*np.sin(2*gamma)
    )
    decl = (
        0.006918 - 0.399912*np.cos(gamma) + 0.070257*np.sin(gamma)
        - 0.006758*np.cos(2*gamma) + 0.000907*np.sin(2*gamma)
        - 0.002697*np.cos(3*gamma) + 0.00148*np.sin(3*gamma)
    )
    solar_minutes = hour * 60 + eqtime + 4 * SITE_LON
    hour_angle = np.deg2rad(solar_minutes / 4 - 180)
    lat_rad = np.deg2rad(SITE_LAT)
    return np.clip(
        np.sin(lat_rad)*np.sin(decl)
        + np.cos(lat_rad)*np.cos(decl)*np.cos(hour_angle), 0, 1
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rgb-dir", type=Path, default=RGB_DIR)
    parser.add_argument("--surfrad-dir", type=Path, default=SURFRAD_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--max-files", type=int)
    parser.add_argument("--rebuild-rgb", action="store_true")
    parser.add_argument("--overwrite-surfrad", action="store_true")
    parser.add_argument("--cloudy-kt", type=float, default=0.55)
    parser.add_argument("--clear-kt", type=float, default=0.85)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = inventory_rgb(args.rgb_dir, args.max_files)
    days = sorted(set(map(rgb_date, paths)))
    print(f"GOES inventory: {len(paths)} daily files, {days[0].date()} to {days[-1].date()}")
    daily_obs = download_surfrad(
        days, args.surfrad_dir, workers=args.workers, overwrite=args.overwrite_surfrad
    )
    observations = read_surfrad(daily_obs)
    observations.to_csv(args.output_dir / "boulder_surfrad_observations.csv", index=False)

    rgb_means_path = args.output_dir / "boulder_rgb_means_all_available.csv"
    if rgb_means_path.exists() and not args.rebuild_rgb and args.max_files is None:
        rgb = pd.read_csv(rgb_means_path, parse_dates=["time"])
        print(f"Reusing RGB means: {rgb_means_path}")
    else:
        rgb = rebuild_rgb_means(paths, rgb_means_path)
    rgb["era5_temp_c"] = load_era5_land_site_temperature(rgb.time)

    rgb_binary, sw_binary = match_and_classify(
        rgb, observations, cloudy_kt=args.cloudy_kt, clear_kt=args.clear_kt
    )
    rgb_path = args.output_dir / "boulder_rgb_binary_all_available.csv"
    sw_path = args.output_dir / "boulder_sw_binary_all_available.csv"
    rgb_binary.to_csv(rgb_path, index=False)
    sw_binary.to_csv(sw_path, index=False)
    print(f"RGB binary: {rgb_path} ({len(rgb_binary):,} rows)")
    print(f"SW binary:  {sw_path} ({len(sw_binary):,} unambiguous daytime rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
