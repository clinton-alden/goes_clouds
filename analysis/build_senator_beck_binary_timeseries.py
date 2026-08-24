#!/usr/bin/env python3
"""Build Senator Beck SWin and notebook-5b TSI/ERA5-Land RGB binaries."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

try:
    from mammoth_cues_binary_analysis import apply_rgb_rules
except ModuleNotFoundError:
    from analysis.mammoth_cues_binary_analysis import apply_rgb_rules


ROOT = Path(__file__).resolve().parent
RGB_CSV = ROOT / "output_11c_senator_beck/senator_beck_rgb_domain_mean_all.csv"
SBSP_CSV = ROOT / "SBSP_1hr_2010-2024.csv"
RULES = ROOT / "tsi_rgb_decision_tree_era5_t2m_10c_rules.csv"
ERA5_DIR = Path("/glade/derecho/scratch/cdalden/tmp/colorado/era5_land/t2m_hourly")
OUT = ROOT / "output_13c_senator_beck_domain_cloud_binary"

SITE_LAT = 37.90
SITE_LON = -107.72
SITE_ALTITUDE_KM = 3.714


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


def load_station() -> pd.DataFrame:
    station = pd.read_csv(SBSP_CSV)
    hour_code = station.Hour.astype(int)
    local_hour = hour_code // 100
    local_minute = hour_code % 100
    local_standard_time = (
        pd.to_datetime(
            station.Year.astype(int).astype(str)
            + "-" + station.DOY.astype(int).astype(str),
            format="%Y-%j",
        )
        + pd.to_timedelta(local_hour, unit="h")
        + pd.to_timedelta(local_minute, unit="m")
    )
    # The SBSP logger uses fixed Mountain Standard Time (UTC-7).
    station["time"] = local_standard_time + pd.Timedelta(hours=7)
    # "Up" denotes the upper Senator Beck station, not upwelling radiation.
    station["sw_obs"] = pd.to_numeric(station.PyUp_Unfilt_W, errors="coerce")
    station = station.loc[station.sw_obs.between(0, 1400), ["time", "sw_obs"]]
    order = np.argsort(station.time.to_numpy(dtype="datetime64[ns]"))
    return station.iloc[order].drop_duplicates("time").reset_index(drop=True)


def load_rgb() -> pd.DataFrame:
    rgb = pd.read_csv(RGB_CSV, parse_dates=["time"])
    rgb = rgb.dropna(subset=["time", "red", "green", "blue"]).copy()
    order = np.argsort(rgb.time.to_numpy(dtype="datetime64[ns]"))
    return rgb.iloc[order].drop_duplicates("time").reset_index(drop=True)


def available_era5_periods() -> dict[pd.Period, Path]:
    result = {}
    for path in ERA5_DIR.glob("era5land_t2m_colorado_*.nc"):
        token = path.stem.rsplit("_", 1)[-1]
        result[pd.Period(token, freq="M")] = path
    if not result:
        raise FileNotFoundError(f"No Colorado ERA5-Land files in {ERA5_DIR}")
    return result


def load_era5_temperature(target_times: pd.Series) -> np.ndarray:
    files = available_era5_periods()
    needed = pd.PeriodIndex(pd.to_datetime(target_times).dt.to_period("M").unique())
    frames = []
    for period in needed:
        path = files.get(period)
        if path is None:
            raise FileNotFoundError(f"Missing ERA5-Land month {period}")
        with xr.open_dataset(path) as ds:
            time_name = "valid_time" if "valid_time" in ds.t2m.dims else "time"
            lon = SITE_LON % 360 if float(ds.longitude.max()) > 180 else SITE_LON
            field = ds.t2m.sel(latitude=SITE_LAT, longitude=lon, method="nearest")
            frames.append(pd.DataFrame({
                "time": pd.to_datetime(field[time_name].values),
                "temperature": np.asarray(field.values, float) - 273.15,
            }))
    hourly = pd.concat(frames, ignore_index=True).drop_duplicates("time")
    target = pd.DatetimeIndex(pd.to_datetime(target_times))
    source = pd.DatetimeIndex(hourly.time)
    return np.interp(
        target.asi8.astype(float),
        source.asi8.astype(float),
        hourly.temperature.to_numpy(float),
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rgb = load_rgb()
    station = load_station()
    era5_files = available_era5_periods()
    first_month, last_month = min(era5_files), max(era5_files)
    start = first_month.start_time
    end = last_month.end_time
    rgb = rgb.loc[rgb.time.between(start, end)].reset_index(drop=True)
    station = station.loc[station.time.between(start, end)].reset_index(drop=True)

    # One independent satellite sample per hourly station observation. GOES
    # scans are centered at hh:02:30, so a 3-minute nearest match is sufficient.
    matched = pd.merge_asof(
        station, rgb, on="time", direction="nearest",
        tolerance=pd.Timedelta(minutes=3),
    ).dropna(subset=["red", "green", "blue"]).reset_index(drop=True)
    matched["era5_temp_c"] = load_era5_temperature(matched.time)
    matched = apply_rgb_rules(matched, RULES)

    matched["cos_sza"] = solar_cosine_zenith(matched.time)
    matched["sw_clear"] = (
        1361.0 * matched.cos_sza * (0.78 + 0.04 * SITE_ALTITUDE_KM)
    )
    matched["k_t"] = matched.sw_obs / matched.sw_clear
    valid = (matched.cos_sza >= 0.25) & matched.sw_obs.notna()
    matched["sw_cloud_binary"] = pd.Series(pd.NA, index=matched.index, dtype="Int8")
    matched.loc[valid & (matched.k_t < 0.55), "sw_cloud_binary"] = 1
    matched.loc[valid & (matched.k_t > 0.85), "sw_cloud_binary"] = 0

    rgb_out = OUT / "senator_beck_rgb_binary_all_available.csv"
    sw_out = OUT / "senator_beck_sw_binary_all_available.csv"
    matched.to_csv(rgb_out, index=False)
    sw = matched.loc[
        matched.sw_cloud_binary.notna(),
        ["time", "sw_obs", "sw_clear", "cos_sza", "k_t", "sw_cloud_binary"],
    ].copy()
    sw["sw_cloud_binary"] = sw.sw_cloud_binary.astype(int)
    sw.to_csv(sw_out, index=False)
    pd.DataFrame([{
        "rgb_source": str(RGB_CSV),
        "station_source": str(SBSP_CSV),
        "sw_variable": "PyUp_Unfilt_W (upper-station downwelling pyranometer)",
        "station_time_basis": "fixed MST (UTC-7)",
        "era5_first_month": str(first_month),
        "era5_last_month": str(last_month),
        "matched_hourly_rows": len(matched),
        "unambiguous_daytime_rows": len(sw),
    }]).to_csv(OUT / "senator_beck_build_metadata.csv", index=False)
    print(f"Hourly RGB/station/ERA5 matches: {len(matched):,}")
    print(f"Unambiguous daytime SWin rows: {len(sw):,}")
    print("Saved to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
