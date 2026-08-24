#!/usr/bin/env python3
"""Build hourly CUES k_t and Mammoth GOES RGB-mask cloud fractions.

GOES cloud fraction is calculated from the completed exact-tree binary masks
over a 10-km E/W by 5-km southward box whose north edge is centered on CUES.
GOES-17 is used through 2023-01-02 and GOES-18 afterward. The hourly value is
the mean of all valid binary pixel predictions from all scans in the UTC hour.

CUES cloud fraction is the fraction of the twelve five-minute samples in each
UTC hour for which k_t < 0.65. Hours without all twelve samples are retained in
the intermediate columns but are not assigned an SW cloud fraction.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from mammoth_cues_binary_analysis import (
    goes_latlon,
    load_clean_cues_sw,
    load_era5_cues,
    solar_cosine_zenith,
)


HERE = Path(__file__).resolve().parent
MAMMOTH_BASE = Path("/glade/derecho/scratch/cdalden/mammoth")
MASK_DIRS = {
    "goes17": MAMMOTH_BASE / "goes17/cloud_mask_tempbin_tree",
    "goes18": MAMMOTH_BASE / "goes18/cloud_mask_tempbin_tree",
}
ERA5_DIR = Path("/glade/derecho/scratch/cdalden/mammoth/era5_land/t2m_hourly")
CUES_PATH = HERE / "CUES_1min_data_atmos_radiation_soiltemp_precip_2015to2025.nc"
RULES_PATH = HERE / "output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv"
OUT_DIR = HERE / "output_19_sw_temporal_cf_eval"
CUES_LAT = 37.6431
CUES_LON = -119.0291
COND_RE = re.compile(r"(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)")


def rule_mask(values: dict[str, np.ndarray], rule: str) -> np.ndarray:
    selected = np.ones_like(values["red"], dtype=bool)
    conditions = COND_RE.findall(rule)
    if not conditions:
        raise ValueError(f"Rule contains no recognized conditions: {rule}")
    for feature, operator, threshold_text in conditions:
        threshold = float(threshold_text)
        if operator == "<=":
            selected &= values[feature] <= threshold
        else:
            selected &= values[feature] > threshold
    return selected


def load_rules() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    rules = pd.read_csv(RULES_PATH)
    rules = rules.loc[rules.status.eq("trained")].copy()
    cloudy = {
        str(label): group.loc[group.prediction.eq(1), "rule"].astype(str).tolist()
        for label, group in rules.groupby("temp_bin", sort=False)
    }
    return rules, cloudy


def assign_temp_bins(temperature: np.ndarray, rules: pd.DataFrame) -> np.ndarray:
    labels = np.full(len(temperature), None, dtype=object)
    for label, group in rules.groupby("temp_bin", sort=False):
        left = float(group.temp_left_c.iloc[0])
        right = float(group.temp_right_c.iloc[0])
        labels[(temperature >= left) & (temperature < right)] = str(label)
    return labels


def domain_mask(ds: xr.Dataset) -> tuple[np.ndarray, dict[str, float | int | str]]:
    lat, lon, projection_lon = goes_latlon(np.asarray(ds.x), np.asarray(ds.y))
    lat_min = CUES_LAT - 5.0 / 111.32
    lon_half_width = 5.0 / (111.32 * np.cos(np.deg2rad(CUES_LAT)))
    mask = (
        (lat >= lat_min) & (lat <= CUES_LAT)
        & (lon >= CUES_LON - lon_half_width)
        & (lon <= CUES_LON + lon_half_width)
    )
    if not mask.any():
        raise ValueError("The requested CUES-south domain selected no GOES pixels")
    metadata = {
        "domain_lat_min": lat_min,
        "domain_lat_max": CUES_LAT,
        "domain_lon_min": CUES_LON - lon_half_width,
        "domain_lon_max": CUES_LON + lon_half_width,
        "domain_grid_pixels": int(mask.sum()),
        "projection_origin_lon": projection_lon,
    }
    return mask, metadata


def corrected_times(path: Path, ds: xr.Dataset) -> pd.DatetimeIndex:
    date = pd.to_datetime(path.stem.rsplit("_", 1)[-1], format="%Y%m%d")
    internal = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
    return date + (internal - internal.normalize())


def expected_satellite(path: Path) -> bool:
    match = re.search(r"_(\d{8})_cloud_binary", path.name)
    if not match:
        return False
    day = pd.to_datetime(match.group(1), format="%Y%m%d")
    return (path.name.startswith("goes17_") and day <= pd.Timestamp("2023-01-02")) or (
        path.name.startswith("goes18_") and day >= pd.Timestamp("2023-01-03")
    )


def build_goes(year: int, max_files: int | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    files = sorted(
        path
        for directory in MASK_DIRS.values()
        for path in directory.glob(f"goes*_C02_C05_C13_rgb_mammoth_1deg_{year}*_cloud_binary_tempbin_tree.nc")
        if path.stat().st_size > 0 and expected_satellite(path)
    )
    if max_files is not None:
        files = files[:max_files]
    if not files:
        raise FileNotFoundError(f"No nonempty {year} exact-tree mask files in {MASK_DIRS}")

    frames: list[pd.DataFrame] = []
    metadata_rows: list[dict] = []
    for number, path in enumerate(files, start=1):
        with xr.open_dataset(path) as ds:
            lat_min = CUES_LAT - 5.0 / 111.32
            lon_half_width = 5.0 / (111.32 * np.cos(np.deg2rad(CUES_LAT)))
            subset = ds.cloud_binary.sel(
                latitude=slice(lat_min, CUES_LAT),
                longitude=slice(CUES_LON - lon_half_width, CUES_LON + lon_half_width),
            )
            if subset.sizes.get("latitude", 0) == 0 or subset.sizes.get("longitude", 0) == 0:
                raise ValueError(f"CUES-south domain selected no pixels in {path}")
            values = np.asarray(subset, dtype=np.float32)
            times = pd.DatetimeIndex(pd.to_datetime(subset.t.values))
            valid_pixels = np.isfinite(values)
            cloud_count = np.nansum(values == 1, axis=(1, 2)).astype(np.int64)
            valid_count = valid_pixels.sum(axis=(1, 2)).astype(np.int64)
            scan_fraction = np.divide(
                cloud_count, valid_count,
                out=np.full(len(times), np.nan), where=valid_count > 0,
            )
            frames.append(pd.DataFrame({
                "time": times,
                "goes_cloud_pixels": cloud_count,
                "goes_valid_pixels": valid_count,
                "goes_scan_cloud_fraction": scan_fraction,
                "source_file": str(path),
            }))
            metadata_rows.append({
                "source_file": str(path),
                "satellite": path.name.split("_", 1)[0],
                "domain_lat_min": lat_min, "domain_lat_max": CUES_LAT,
                "domain_lon_min": CUES_LON - lon_half_width,
                "domain_lon_max": CUES_LON + lon_half_width,
                "domain_grid_pixels": int(subset.sizes["latitude"] * subset.sizes["longitude"]),
                "mask_tree_logic": ds.attrs.get("tree_logic", ""),
            })
        if number % 25 == 0 or number == len(files):
            print(f"Processed {number}/{len(files)} mask files for {year}", flush=True)

    scans = pd.concat(frames, ignore_index=True).sort_values("time").drop_duplicates("time")
    scans["hour"] = scans.time.dt.floor("h")
    hourly = scans.groupby("hour", as_index=False).agg(
        goes_cloud_pixels=("goes_cloud_pixels", "sum"),
        goes_valid_pixels=("goes_valid_pixels", "sum"),
        goes_n_scans=("time", "size"),
        goes_mean_scan_cloud_fraction=("goes_scan_cloud_fraction", "mean"),
    )
    hourly["goes_cloud_fraction"] = (
        hourly.goes_cloud_pixels / hourly.goes_valid_pixels
    )
    metadata = pd.DataFrame(metadata_rows).drop_duplicates(
        ["satellite", "domain_grid_pixels", "mask_tree_logic"]
    )
    return hourly, metadata


def build_sw(year: int, kt_threshold: float) -> pd.DataFrame:
    start = pd.Timestamp(year=year, month=1, day=1)
    end = pd.Timestamp(year=year + 1, month=1, day=1) - pd.Timedelta(minutes=1)
    sw, _ = load_clean_cues_sw(CUES_PATH, start=start, end=end)
    five_minute = (
        sw.set_index("time")[["sw_obs_5min"]]
        .resample("5min").first()
        .loc[start:end]
        .reset_index()
    )
    five_minute["cos_sza"] = solar_cosine_zenith(
        five_minute.time, lat=CUES_LAT, lon=CUES_LON
    )
    five_minute["sw_clear"] = 1361.0 * five_minute.cos_sza * (0.78 + 0.04 * 2.9)
    five_minute["k_t"] = five_minute.sw_obs_5min / five_minute.sw_clear
    five_minute["valid"] = (
        five_minute.sw_obs_5min.notna() & (five_minute.cos_sza >= 0.25)
        & np.isfinite(five_minute.k_t)
    )
    five_minute["sw_cloud"] = np.where(
        five_minute.valid, five_minute.k_t < kt_threshold, np.nan
    )
    five_minute["hour"] = five_minute.time.dt.floor("h")
    hourly = five_minute.groupby("hour", as_index=False).agg(
        sw_cloud_count=("sw_cloud", "sum"),
        sw_n_valid=("valid", "sum"),
        sw_mean_k_t=("k_t", "mean"),
    )
    hourly["sw_cloud_fraction"] = hourly.sw_cloud_count / 12.0
    hourly.loc[hourly.sw_n_valid.ne(12), "sw_cloud_fraction"] = np.nan
    hourly["kt_threshold"] = kt_threshold
    return hourly


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--kt-threshold", type=float, default=0.65)
    parser.add_argument("--max-files", type=int, help="Testing only: process the first N daily files")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.year}" if args.max_files is None else f"_{args.year}_test"
    output = OUT_DIR / f"mammoth_cues_hourly_cloud_fraction{suffix}.csv"
    if output.exists() and not args.force:
        print(f"Output exists; use --force to replace it: {output}")
        return 0

    goes, metadata = build_goes(args.year, max_files=args.max_files)
    sw = build_sw(args.year, args.kt_threshold)
    comparison = goes.merge(sw, on="hour", how="left", validate="one_to_one")
    comparison["residual_goes_minus_sw"] = (
        comparison.goes_cloud_fraction - comparison.sw_cloud_fraction
    )
    comparison.to_csv(output, index=False)
    metadata.assign(year=args.year).to_csv(
        OUT_DIR / f"mammoth_cues_hourly_cloud_fraction_metadata{suffix}.csv", index=False
    )
    print(f"Saved {len(comparison):,} GOES hours to {output}")
    print(f"Complete SW hours: {comparison.sw_cloud_fraction.notna().sum():,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
