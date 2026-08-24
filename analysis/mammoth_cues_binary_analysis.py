"""Reusable processing functions for Notebook 13 Mammoth/CUES binary comparison."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS, Transformer


COND_RE = re.compile(r"(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)")


def load_rgb_means(path: Path) -> pd.DataFrame:
    rgb = pd.read_csv(path, parse_dates=["time"])
    required = {"time", "red", "green", "blue"}
    missing = required.difference(rgb.columns)
    if missing:
        raise KeyError(f"Missing RGB columns: {sorted(missing)}")
    rgb = rgb.dropna(subset=list(required)).sort_values("time").drop_duplicates("time")
    daytime = (rgb.time.dt.hour >= 15) | (rgb.time.dt.hour < 1)
    return rgb.loc[(rgb.time.dt.year == 2024) & daytime].reset_index(drop=True)


def goes_latlon(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    # The archive switches satellite locations and the RGB files do not retain
    # their projection variable. Choose the known origin that maps the grid
    # center nearest Mammoth; GOES-18 was near -89.5 during 2022 checkout and
    # at -137 after becoming GOES-West.
    center_x = float(np.nanmedian(x)) * 35786023.0
    center_y = float(np.nanmedian(y)) * 35786023.0
    candidates = (-137.0, -89.5, -75.0)
    distances = []
    for candidate in candidates:
        candidate_crs = CRS.from_proj4(
            f"+proj=geos +h=35786023 +lon_0={candidate} +sweep=x "
            "+a=6378137 +b=6356752.31414 +units=m +no_defs"
        )
        center_lon, center_lat = Transformer.from_crs(
            candidate_crs, 4326, always_xy=True
        ).transform(center_x, center_y)
        distances.append((center_lat - 37.6431) ** 2 + (center_lon + 119.0291) ** 2)
    projection_lon = candidates[int(np.argmin(distances))]
    projection = CRS.from_proj4(
        f"+proj=geos +h=35786023 +lon_0={projection_lon} +sweep=x "
        "+a=6378137 +b=6356752.31414 +units=m +no_defs"
    )
    xx, yy = np.meshgrid(x * 35786023.0, y * 35786023.0)
    lon, lat = Transformer.from_crs(projection, 4326, always_xy=True).transform(xx, yy)
    return np.asarray(lat), np.asarray(lon), projection_lon


def rebuild_rgb_means(
    rgb_dir: Path,
    output_csv: Path,
    year: int = 2024,
    cues_lat: float = 37.6431,
    cues_lon: float = -119.0291,
) -> tuple[pd.DataFrame, dict]:
    """Recompute band means directly from raw RGB files; no prior notebook outputs."""
    files = sorted(rgb_dir.glob(f"goes18_C02_C05_C13_rgb_mammoth_{year}*.nc"))
    files = [path for path in files if path.stat().st_size > 0]
    if not files:
        raise FileNotFoundError(f"No nonempty {year} RGB files in {rgb_dir}")

    lat_min = cues_lat - 5.0 / 111.32
    lat_max = cues_lat
    lon_half_width = 5.0 / (111.32 * np.cos(np.deg2rad(cues_lat)))
    lon_min, lon_max = cues_lon - lon_half_width, cues_lon + lon_half_width

    frames = []
    mask_cache = {}
    selected_pixel_counts = set()
    projection_origins = set()
    for index, path in enumerate(files, start=1):
        with xr.open_dataset(path) as ds:
            x_values, y_values = np.asarray(ds.x), np.asarray(ds.y)
            grid_key = (
                len(x_values), len(y_values), float(x_values[0]), float(x_values[-1]),
                float(y_values[0]), float(y_values[-1]),
            )
            if grid_key not in mask_cache:
                lat, lon, projection_lon = goes_latlon(np.asarray(ds.x), np.asarray(ds.y))
                pixel_mask = ((lat >= lat_min) & (lat <= lat_max)
                              & (lon >= lon_min) & (lon <= lon_max))
                if not pixel_mask.any():
                    raise ValueError("CUES-north box selected no GOES pixels")
                mask_cache[grid_key] = pixel_mask
                selected_pixel_counts.add(int(pixel_mask.sum()))
                projection_origins.add(projection_lon)
            pixel_mask = mask_cache[grid_key]
            date = pd.to_datetime(path.stem.rsplit("_", 1)[-1], format="%Y%m%d")
            internal = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
            corrected_time = date + (internal - internal.normalize())
            values = {}
            for band in ("red", "green", "blue"):
                array = np.asarray(ds[band], dtype=np.float32)
                selected = array[:, pixel_mask]
                valid_count = np.isfinite(selected).sum(axis=1)
                values[band] = np.divide(
                    np.nansum(selected, axis=1), valid_count,
                    out=np.full(valid_count.shape, np.nan, dtype=np.float32),
                    where=valid_count > 0,
                )
            frames.append(pd.DataFrame({"time": corrected_time, **values}))
        if index % 50 == 0 or index == len(files):
            print(f"Processed raw RGB files: {index}/{len(files)}", flush=True)

    rgb = pd.concat(frames, ignore_index=True)
    rgb = rgb.replace([np.inf, -np.inf], np.nan).dropna(subset=["red", "green", "blue"])
    rgb = rgb.sort_values("time").drop_duplicates("time").reset_index(drop=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    rgb.to_csv(output_csv, index=False)
    metadata = {
        "source_file_count": len(files),
        "selected_pixel_count_min": min(selected_pixel_counts),
        "selected_pixel_count_max": max(selected_pixel_counts),
        "lat_min": lat_min, "lat_max": lat_max, "lon_min": lon_min, "lon_max": lon_max,
        "projection_origin_lon": ";".join(str(value) for value in sorted(projection_origins)),
    }
    return rgb, metadata


def load_era5_cues(era5_dir: Path, target_times: pd.Series, lat=37.6431, lon=-119.0291):
    frames = []
    months = list(pd.PeriodIndex(target_times.dt.to_period("M").unique()))
    next_month = max(months) + 1
    next_path = era5_dir / f"era5land_t2m_mammoth_{next_month.year}{next_month.month:02d}.nc"
    if next_path.exists():
        months.append(next_month)
    for month in months:
        path = era5_dir / f"era5land_t2m_mammoth_{month.year}{month.month:02d}.nc"
        if not path.exists():
            raise FileNotFoundError(f"ERA5-Land month is not ready: {path}")
        with xr.open_dataset(path) as ds:
            tname = "valid_time" if "valid_time" in ds.t2m.dims else "time"
            target_lon = lon % 360 if float(ds.longitude.max()) > 180 else lon
            da = ds.t2m.sel(latitude=lat, longitude=target_lon, method="nearest") - 273.15
            frames.append(pd.DataFrame({"time": pd.to_datetime(da[tname].values), "era5_temp_c": da.values}))
    hourly = pd.concat(frames, ignore_index=True).drop_duplicates("time").set_index("time")
    source_time = pd.DatetimeIndex(hourly.index)
    target_time = pd.DatetimeIndex(pd.to_datetime(target_times))
    interpolated = np.interp(
        target_time.asi8.astype(float),
        source_time.asi8.astype(float),
        hourly.era5_temp_c.to_numpy(float),
    )
    return interpolated, hourly


def load_cloudy_rules(path: Path):
    rules = pd.read_csv(path)
    rules = rules.loc[rules.predicted_class.eq(1)]
    cloudy = {str(k): g.rule.astype(str).tolist() for k, g in rules.groupby("temp_bin", sort=False)}
    midpoints = {}
    for label in cloudy:
        left, right = re.findall(r"-?\d+(?:\.\d+)?", label)
        midpoints[label] = (float(left) + float(right)) / 2
    return cloudy, midpoints


def temperature_bin(temp_c: float, midpoints: dict[str, float]) -> str:
    nominal_left = np.floor(temp_c / 10) * 10
    nominal = f"[{int(nominal_left)}, {int(nominal_left + 10)})"
    if nominal in midpoints:
        return nominal
    return min(midpoints, key=lambda label: abs(temp_c - midpoints[label]))


def rule_matches(row: pd.Series, rule: str) -> bool:
    for feature, operator, threshold_text in COND_RE.findall(rule):
        value, threshold = float(row[feature]), float(threshold_text)
        if operator == "<=" and not value <= threshold:
            return False
        if operator == ">" and not value > threshold:
            return False
    return True


def apply_rgb_rules(rgb: pd.DataFrame, rules_path: Path) -> pd.DataFrame:
    cloudy, midpoints = load_cloudy_rules(rules_path)
    out = rgb.copy()
    out["temp_bin"] = out.era5_temp_c.map(lambda value: temperature_bin(value, midpoints))
    out["rgb_cloud_binary"] = [
        int(any(rule_matches(row, rule) for rule in cloudy[row.temp_bin]))
        for _, row in out.iterrows()
    ]
    return out


def apply_gothic_tree_rules(
    rgb: pd.DataFrame,
    rules_path: Path,
    temp_col: str,
) -> pd.DataFrame:
    """Apply the latest Gothic temperature-bin decision-tree leaf rules."""
    rules = pd.read_csv(rules_path)
    required = {
        "temp_bin", "temp_left_c", "temp_right_c", "status",
        "rule", "prediction",
    }
    missing = required.difference(rules.columns)
    if missing:
        raise KeyError(f"Missing Gothic rule columns: {sorted(missing)}")
    rules = rules.loc[rules.status.eq("trained")].copy()
    if rules.empty:
        raise ValueError(f"No trained Gothic rules in {rules_path}")

    out = rgb.copy()
    out["temp_bin"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out["rgb_cloud_binary"] = pd.Series(pd.NA, index=out.index, dtype="Int8")
    temperature = pd.to_numeric(out[temp_col], errors="coerce")
    for label, group in rules.groupby("temp_bin", sort=False):
        left = float(group.temp_left_c.iloc[0])
        right = float(group.temp_right_c.iloc[0])
        selected = temperature.ge(left) & temperature.lt(right)
        if not selected.any():
            continue
        out.loc[selected, "temp_bin"] = str(label)
        rows = out.loc[selected]
        predictions = []
        for _, row in rows.iterrows():
            matching = group.loc[group.rule.map(lambda rule: rule_matches(row, str(rule)))]
            if len(matching) != 1:
                raise ValueError(
                    f"Expected one Gothic tree leaf for {label}, found {len(matching)}"
                )
            predictions.append(int(matching.prediction.iloc[0]))
        out.loc[selected, "rgb_cloud_binary"] = predictions
    return out.dropna(subset=["rgb_cloud_binary"]).assign(
        rgb_cloud_binary=lambda frame: frame.rgb_cloud_binary.astype(np.int8)
    )


def solar_cosine_zenith(times_utc, lat=37.6431, lon=-119.0291):
    times = pd.DatetimeIndex(pd.to_datetime(times_utc))
    doy = times.dayofyear.to_numpy()
    hour = times.hour + times.minute / 60 + times.second / 3600
    gamma = 2 * np.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eqtime = 229.18 * (0.000075 + 0.001868*np.cos(gamma) - 0.032077*np.sin(gamma)
                       - 0.014615*np.cos(2*gamma) - 0.040849*np.sin(2*gamma))
    decl = (0.006918 - 0.399912*np.cos(gamma) + 0.070257*np.sin(gamma)
            - 0.006758*np.cos(2*gamma) + 0.000907*np.sin(2*gamma)
            - 0.002697*np.cos(3*gamma) + 0.00148*np.sin(3*gamma))
    solar_minutes = hour * 60 + eqtime + 4 * lon
    hour_angle = np.deg2rad(solar_minutes / 4 - 180)
    lat_rad = np.deg2rad(lat)
    return np.clip(np.sin(lat_rad)*np.sin(decl) + np.cos(lat_rad)*np.cos(decl)*np.cos(hour_angle), 0, 1)


def load_clean_cues_sw(
    path: Path, start: pd.Timestamp | None = None, end: pd.Timestamp | None = None
) -> tuple[pd.DataFrame, dict]:
    with xr.open_dataset(path) as ds:
        attrs = dict(ds.attrs)
        sw = pd.DataFrame({
            "time_pst": pd.to_datetime(ds.datetime.values),
            "sw_obs": np.asarray(ds.downwelling_shortwave_platform, float),
            "instrument_flag": np.asarray(ds.downwelling_shortwave_platform_instrument_flag, int),
            "air_temp_c": np.asarray(ds.air_temperature, float),
        })
    # The file declares fixed America/PST, UTC-8 (not daylight-saving local time).
    sw["time"] = sw.time_pst + pd.Timedelta(hours=8)
    sw = sw.loc[sw.sw_obs.between(0, 1400) & (sw.instrument_flag > 0)]
    if start is not None:
        sw = sw.loc[sw.time >= pd.Timestamp(start)]
    if end is not None:
        sw = sw.loc[sw.time <= pd.Timestamp(end)]
    sw = sw.dropna(subset=["sw_obs"]).drop_duplicates("time")
    # Centered five-minute average matches the spatially averaged 5-minute GOES sample.
    sw["sw_obs_5min"] = sw.set_index("time").sw_obs.rolling("5min", center=True, min_periods=3).mean().to_numpy()
    sw["air_temp_5min_c"] = (
        sw.set_index("time").air_temp_c.rolling("5min", center=True, min_periods=3).mean().to_numpy()
    )
    return sw.reset_index(drop=True), attrs


def attach_shortwave(
    rgb: pd.DataFrame, sw: pd.DataFrame, altitude_km=2.9,
    cloudy_kt=0.55, clear_kt=0.85,
) -> pd.DataFrame:
    rgb = rgb.copy()
    sw = sw.copy()
    rgb["time"] = pd.to_datetime(rgb.time).astype("datetime64[ns]")
    sw["time"] = pd.to_datetime(sw.time).astype("datetime64[ns]")
    matched = pd.merge_asof(
        rgb,
        sw[["time", "sw_obs_5min", "air_temp_5min_c"]],
        on="time", direction="nearest", tolerance=pd.Timedelta("45s")
    )
    matched["cos_sza"] = solar_cosine_zenith(matched.time)
    matched["sw_clear"] = 1361.0 * matched.cos_sza * (0.78 + 0.04 * altitude_km)
    matched["k_t"] = matched.sw_obs_5min / matched.sw_clear
    valid = (matched.cos_sza >= 0.25) & matched.sw_obs_5min.notna() & (matched.sw_clear > 0)
    matched["sw_cloud_binary"] = pd.Series(pd.NA, index=matched.index, dtype="Int8")
    matched.loc[valid & (matched.k_t < cloudy_kt), "sw_cloud_binary"] = 1
    matched.loc[valid & (matched.k_t > clear_kt), "sw_cloud_binary"] = 0
    return matched
