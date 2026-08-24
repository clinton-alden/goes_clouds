#!/usr/bin/env python
"""Build a point SW comparison table for NLDAS-2, NLDAS-3, and Gothic obs.

The expensive part of this workflow is opening many small files, not the
calculation itself. This script reads only the nearest grid-cell SWdown value
for each requested hour, computes a Metsim clear-sky reference at that grid
cell, merges Gothic ARM shortwave observations, and writes one cached CSV.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import types
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path(__file__).resolve().parent
DEFAULT_NLDAS2_DIR = Path("/glade/u/home/cdalden/scratch/nldas")
DEFAULT_NLDAS3_DIR = Path("/glade/derecho/scratch/cdalden/tmp/nldas3_sw")
DEFAULT_OBS_DIR = Path("/glade/u/home/cdalden/scratch/surface_obs/colorado/gucsebsM1.b1")
DEFAULT_OUT_DIR = ROOT / "output_17_nldas_point_eval"
DEFAULT_NLDAS3_PYTHON = Path("/glade/work/cdalden/conda-envs/goes_downloading/bin/python")
DEFAULT_NLDAS3_DOWNLOADER = ROOT / "download_nldas3_sw_subset.py"

GOTHIC_LAT = 38.96
GOTHIC_LON = -106.99
GOTHIC_ELEV_M = 2900.0
METSIM_LAPSE_RATE = 0.0065
METSIM_TBASE = 0.90
LOCAL_TZ = "America/Denver"
OBS_VAR = "down_short_hemisp"

SW_VAR_CANDIDATES = (
    "SWdown",
    "SWDOWN",
    "DSWRF",
    "DSWRF_surface",
    "dswrf",
    "surface_downwelling_shortwave_flux",
    "shortwave_down",
)


try:
    from metsim import constants as metsim_constants
    from metsim.disaggregate import shortwave as metsim_shortwave
    from metsim.physics import solar_geom
except ModuleNotFoundError:
    metsim_pkg = Path("/glade/work/cdalden/.conda/pkgs/metsim-2.4.4-pyhd8ed1ab_0/site-packages")
    if metsim_pkg.exists():
        sys.path.append(str(metsim_pkg))
    if "numba" not in sys.modules:
        numba_stub = types.ModuleType("numba")

        def _jit(*args, **kwargs):
            if args and callable(args[0]) and len(args) == 1 and not kwargs:
                return args[0]
            return lambda func: func

        numba_stub.jit = _jit
        sys.modules["numba"] = numba_stub
    from metsim import constants as metsim_constants
    from metsim.disaggregate import shortwave as metsim_shortwave
    from metsim.physics import solar_geom


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", default="2022-05-01", help="First local/UTC date, YYYY-MM-DD")
    parser.add_argument("--end", default="2022-05-31", help="Last date, YYYY-MM-DD")
    parser.add_argument("--hours", default="14-23,0", help="UTC hours, e.g. '14-23,0' or '0-23'")
    parser.add_argument("--lat", type=float, default=GOTHIC_LAT)
    parser.add_argument("--lon", type=float, default=GOTHIC_LON)
    parser.add_argument("--nldas2-dir", type=Path, default=DEFAULT_NLDAS2_DIR)
    parser.add_argument("--nldas3-dir", type=Path, default=DEFAULT_NLDAS3_DIR)
    parser.add_argument("--obs-dir", type=Path, default=DEFAULT_OBS_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--download-nldas3", action="store_true", help="Download missing NLDAS3 point subsets")
    parser.add_argument("--skip-nldas2", action="store_true", help="Skip NLDAS-2 reads and write nldas2 as NaN")
    parser.add_argument("--force-nldas3", action="store_true", help="Force re-download of NLDAS3 subset files")
    parser.add_argument("--nldas3-python", type=Path, default=DEFAULT_NLDAS3_PYTHON)
    parser.add_argument("--nldas3-downloader", type=Path, default=DEFAULT_NLDAS3_DOWNLOADER)
    parser.add_argument("--point-box-deg", type=float, default=0.05, help="Half-width for NLDAS3 subset download")
    return parser.parse_args()


def parse_hours(spec: str) -> list[int]:
    hours: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = [int(x) for x in part.split("-", 1)]
            if start <= end:
                hours.extend(range(start, end + 1))
            else:
                hours.extend(list(range(start, 24)) + list(range(0, end + 1)))
        else:
            hours.append(int(part))
    unique = []
    for hour in hours:
        if hour < 0 or hour > 23:
            raise ValueError(f"Invalid hour {hour}; expected 0-23")
        if hour not in unique:
            unique.append(hour)
    return unique


def build_times(start: str, end: str, hours: list[int]) -> pd.DatetimeIndex:
    times = []
    for day in pd.date_range(pd.Timestamp(start), pd.Timestamp(end), freq="D"):
        for hour in hours:
            offset_day = day + pd.Timedelta(days=1) if hour == 0 and any(h > 0 for h in hours) else day
            times.append(offset_day + pd.Timedelta(hours=hour))
    return pd.DatetimeIndex(times).drop_duplicates().sort_values()


def find_nldas2_file(dt: pd.Timestamp, nldas2_dir: Path) -> Path:
    pattern = f"*A{dt:%Y%m%d}.{dt:%H}00*"
    matches = sorted(nldas2_dir.glob(pattern))
    if len(matches) > 1:
        clean_matches = [path for path in matches if not path.name.startswith("HTTP_services.cgi")]
        if len(clean_matches) == 1:
            return clean_matches[0]
        nc4_matches = [path for path in clean_matches if path.suffix in (".nc", ".nc4") or ".nc" in path.name]
        if len(nc4_matches) == 1:
            return nc4_matches[0]
        raise FileExistsError(f"Multiple NLDAS2 matches for {pattern}; found {matches}")
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one NLDAS2 match for {pattern}; found {len(matches)}")
    return matches[0]


def nldas3_local_name(dt: pd.Timestamp, nldas3_dir: Path) -> Path:
    return nldas3_dir / f"nldas3_sw_{dt:%Y%m%d_%H}00.nc"


def find_nldas3_file(dt: pd.Timestamp, nldas3_dir: Path) -> Path:
    preferred = nldas3_local_name(dt, nldas3_dir)
    if preferred.exists():
        return preferred

    patterns = [
        f"*A{dt:%Y%m%d}.{dt:%H}00*",
        f"*{dt:%Y%m%d}*{dt:%H}00*",
        f"*{dt:%Y%m%d}_{dt:%H}00*",
    ]
    matches = []
    for pattern in patterns:
        matches.extend(nldas3_dir.glob(pattern))
    matches = sorted({m for m in matches if m.is_file() and not m.name.startswith(".")})
    if len(matches) == 1:
        return matches[0]
    if preferred in matches:
        return preferred
    if len(matches) > 1:
        raise FileExistsError(f"Multiple NLDAS3 matches for {dt}: {matches}")
    raise FileNotFoundError(f"No cached NLDAS3 file for {dt} in {nldas3_dir}")


def run_nldas3_downloader(times: pd.DatetimeIndex, args: argparse.Namespace) -> None:
    if len(times) == 0:
        return
    if not args.nldas3_python.exists():
        raise FileNotFoundError(f"NLDAS3 Python not found: {args.nldas3_python}")
    if not args.nldas3_downloader.exists():
        raise FileNotFoundError(f"NLDAS3 downloader not found: {args.nldas3_downloader}")

    cmd = [
        str(args.nldas3_python),
        str(args.nldas3_downloader),
        "--start",
        times.min().strftime("%Y-%m-%dT%H:%M"),
        "--end",
        times.max().strftime("%Y-%m-%dT%H:%M"),
        "--out-dir",
        str(args.nldas3_dir),
        "--point-lat",
        str(args.lat),
        "--point-lon",
        str(args.lon),
    ]
    if args.force_nldas3:
        cmd.append("--force")
    result = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"NLDAS3 subset download failed.\nSTDERR:\n{result.stderr}")


def ensure_nldas3_files(times: pd.DatetimeIndex, args: argparse.Namespace) -> dict[pd.Timestamp, Path]:
    args.nldas3_dir.mkdir(parents=True, exist_ok=True)

    missing = []
    found: dict[pd.Timestamp, Path] = {}
    if args.force_nldas3:
        missing = [pd.Timestamp(dt) for dt in times]
    else:
        for dt in times:
            try:
                found[pd.Timestamp(dt)] = find_nldas3_file(pd.Timestamp(dt), args.nldas3_dir)
            except FileNotFoundError:
                missing.append(pd.Timestamp(dt))

    if missing and (args.download_nldas3 or args.force_nldas3):
        # Download only requested contiguous hour blocks. Calling the downloader
        # once for a whole month would fill every intervening nighttime hour.
        missing_idx = pd.DatetimeIndex(sorted(missing))
        block_start = missing_idx[0]
        prev = missing_idx[0]
        for dt in missing_idx[1:]:
            if dt - prev != pd.Timedelta(hours=1):
                run_nldas3_downloader(pd.date_range(block_start, prev, freq="1h"), args)
                block_start = dt
            prev = dt
        run_nldas3_downloader(pd.date_range(block_start, prev, freq="1h"), args)
        for dt in missing:
            found[dt] = find_nldas3_file(dt, args.nldas3_dir)

    return found


def find_coord_name(ds: xr.Dataset, candidates: tuple[str, ...]) -> str:
    names = set(ds.coords) | set(ds.dims) | set(ds.variables)
    for name in candidates:
        if name in names:
            return name
    lower_lookup = {name.lower(): name for name in names}
    for name in candidates:
        if name.lower() in lower_lookup:
            return lower_lookup[name.lower()]
    raise KeyError(f"Could not find coordinate among {candidates}; available names: {sorted(names)[:30]}")


def find_sw_var(ds: xr.Dataset) -> str:
    for name in SW_VAR_CANDIDATES:
        if name in ds.data_vars:
            return name
    for name, da in ds.data_vars.items():
        attrs = " ".join(str(da.attrs.get(k, "")) for k in ("standard_name", "long_name", "name")).lower()
        if "shortwave" in attrs and ("down" in attrs or "surface_downwelling" in attrs):
            return name
    raise KeyError(f"No downwelling shortwave variable found in {list(ds.data_vars)}")


def read_point_sw(path: Path, time_utc: pd.Timestamp, lat: float, lon: float, source: str) -> dict[str, object]:
    with xr.open_dataset(path) as ds:
        sw_name = find_sw_var(ds)
        lat_name = find_coord_name(ds, ("lat", "latitude", "y"))
        lon_name = find_coord_name(ds, ("lon", "longitude", "x"))
        da = ds[sw_name]

        for time_name in ("time", "valid_time", "time_utc"):
            if time_name in da.dims:
                try:
                    da = da.sel({time_name: time_utc}, method="nearest")
                except (KeyError, TypeError, ValueError):
                    da = da.isel({time_name: 0})
                break

        point_lon = lon
        if float(ds[lon_name].max()) > 180 and lon < 0:
            point_lon = lon % 360
        point = da.sel({lat_name: lat, lon_name: point_lon}, method="nearest")
        value = float(point)
        if np.isclose(value, -9999.0):
            value = np.nan

        return {
            "time_utc": pd.Timestamp(time_utc),
            "source": source,
            "sw_wm2": value,
            "grid_lat": float(point[lat_name]),
            "grid_lon": float(point[lon_name]),
            "sw_var": sw_name,
            "path": str(path),
        }


def daylight_saving_offset_hours(days: pd.DatetimeIndex, local_tz: str) -> np.ndarray:
    noons = pd.DatetimeIndex(days) + pd.Timedelta(hours=12)
    return np.array([noon.tz_localize(local_tz).dst().total_seconds() / 3600 for noon in noons])


@lru_cache(maxsize=32)
def metsim_solar_arrays(lat: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    metsim_constants.TBASE = METSIM_TBASE
    try:
        solar_geom.recompile()
    except AttributeError:
        pass
    tiny_rad_fract, daylength, flat_potrad, tt_max0 = solar_geom(GOTHIC_ELEV_M, float(lat), METSIM_LAPSE_RATE)
    return tiny_rad_fract, daylength, flat_potrad * tt_max0


def metsim_clear_sky_one_time(time_utc: pd.Timestamp, lat: float, lon: float) -> float:
    time_local = pd.Timestamp(time_utc).tz_localize("UTC").tz_convert(LOCAL_TZ).tz_localize(None)
    day = time_local.normalize()
    yday0 = day.dayofyear - 1
    tiny_rad_fract, daylength, clear_daily_by_day = metsim_solar_arrays(round(float(lat), 6))
    clear_daily = np.array([clear_daily_by_day[yday0]])
    params = {"time_step": 60, "method": "mtclim", "utc_offset": False, "lon": float(lon)}
    clear_subdaily = metsim_shortwave(
        clear_daily,
        np.array([daylength[yday0]]),
        np.array([yday0 + 1]),
        tiny_rad_fract,
        params,
    )
    clear_times_solar = pd.date_range(day, periods=24, freq="1h")
    dst_offset = daylight_saving_offset_hours(pd.DatetimeIndex([day]), LOCAL_TZ)[0]
    clear_times_local = clear_times_solar + pd.to_timedelta(dst_offset, unit="h")
    nearest_idx = int(np.argmin(np.abs(clear_times_local - time_local)))
    return float(clear_subdaily[nearest_idx])


def find_obs_file(day: pd.Timestamp, obs_dir: Path) -> Path:
    path = obs_dir / f"gucsebsM1.b1.{day:%Y%m%d}.000000.custom.cdf"
    if not path.exists():
        raise FileNotFoundError(f"Missing Gothic obs file: {path}")
    return path


def read_observed_sw(times: pd.DatetimeIndex, obs_dir: Path) -> pd.DataFrame:
    days = pd.DatetimeIndex(pd.to_datetime(times.normalize().unique()))
    rows = []
    missing_days = []
    for day in days:
        try:
            path = find_obs_file(pd.Timestamp(day), obs_dir)
        except FileNotFoundError:
            missing_days.append(pd.Timestamp(day))
            continue
        with xr.open_dataset(path) as ds:
            values = np.asarray(ds[OBS_VAR].to_numpy(), dtype=float)
        time_utc = pd.date_range(pd.Timestamp(day), periods=len(values), freq="30min")
        rows.append(pd.DataFrame({"time_utc": time_utc, "gothic_obs_sw": values}))

    if missing_days:
        by_month = Counter(day.strftime("%Y-%m") for day in missing_days)
        summary = ", ".join(f"{month}: {count}" for month, count in sorted(by_month.items()))
        print(f"Missing Gothic obs days: {summary}", flush=True)
    if not rows:
        return pd.DataFrame({"time_utc": times, "gothic_obs_sw": np.nan})

    obs = pd.concat(rows, ignore_index=True)
    obs.loc[np.isclose(obs["gothic_obs_sw"], -9999.0), "gothic_obs_sw"] = np.nan
    obs_hourly = (
        obs.dropna(subset=["gothic_obs_sw"])
        .set_index("time_utc")["gothic_obs_sw"]
        .reindex(obs.set_index("time_utc").index.union(times))
        .sort_index()
        .interpolate(method="time")
        .reindex(times)
        .rename("gothic_obs_sw")
        .reset_index()
        .rename(columns={"index": "time_utc"})
    )
    if missing_days:
        missing_norm = {day.normalize() for day in missing_days}
        obs_hourly.loc[obs_hourly["time_utc"].dt.normalize().isin(missing_norm), "gothic_obs_sw"] = np.nan
    return obs_hourly


def build_point_eval(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, Path, Path]:
    times = build_times(args.start, args.end, parse_hours(args.hours))
    nldas3_files = ensure_nldas3_files(times, args)

    rows = []
    missing_nldas2 = []
    for dt in times:
        dt = pd.Timestamp(dt)
        if not args.skip_nldas2:
            try:
                rows.append(read_point_sw(find_nldas2_file(dt, args.nldas2_dir), dt, args.lat, args.lon, "nldas2"))
            except FileNotFoundError as exc:
                missing_nldas2.append((dt, str(exc)))

        if dt in nldas3_files:
            rows.append(read_point_sw(nldas3_files[dt], dt, args.lat, args.lon, "nldas3"))

    if missing_nldas2:
        by_month = Counter(dt.strftime("%Y-%m") for dt, _ in missing_nldas2)
        summary = ", ".join(f"{month}: {count}" for month, count in sorted(by_month.items()))
        print(f"Missing NLDAS2 requested hours: {summary}", flush=True)
    missing_nldas3 = [pd.Timestamp(dt) for dt in times if pd.Timestamp(dt) not in nldas3_files]
    if missing_nldas3:
        by_month = Counter(dt.strftime("%Y-%m") for dt in missing_nldas3)
        summary = ", ".join(f"{month}: {count}" for month, count in sorted(by_month.items()))
        print(f"Missing NLDAS3 requested hours: {summary}", flush=True)

    long = pd.DataFrame(rows)
    if long.empty:
        raise RuntimeError("No NLDAS point rows were read.")

    table = long.pivot_table(index="time_utc", columns="source", values="sw_wm2", aggfunc="first").reset_index()
    table = table.rename_axis(None, axis=1)
    for source in ("nldas2", "nldas3"):
        if source not in table:
            table[source] = np.nan
    obs = read_observed_sw(pd.DatetimeIndex(table["time_utc"]), args.obs_dir)
    table = table.merge(obs, on="time_utc", how="left")

    clear_lat = args.lat
    clear_lon = args.lon
    if (long["source"] == "nldas2").any():
        first_v2 = long[long["source"] == "nldas2"].iloc[0]
        clear_lat = float(first_v2["grid_lat"])
        clear_lon = float(first_v2["grid_lon"])
    table["clear_sky_sw"] = [metsim_clear_sky_one_time(t, clear_lat, clear_lon) for t in table["time_utc"]]
    for source in ("nldas2", "nldas3", "gothic_obs_sw"):
        if source in table:
            table[f"{source}_kt"] = table[source] / table["clear_sky_sw"]

    metric_rows = []
    for source in ("nldas2", "nldas3"):
        if source not in table:
            continue
        valid = table[[source, "gothic_obs_sw"]].dropna()
        if valid.empty:
            continue
        err = valid[source] - valid["gothic_obs_sw"]
        metric_rows.append(
            {
                "source": source,
                "n": len(valid),
                "bias_wm2": err.mean(),
                "mae_wm2": err.abs().mean(),
                "rmse_wm2": np.sqrt((err**2).mean()),
                "corr": valid[source].corr(valid["gothic_obs_sw"]),
            }
        )
    metrics = pd.DataFrame(metric_rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{pd.Timestamp(args.start):%Y%m%d}_{pd.Timestamp(args.end):%Y%m%d}_h{args.hours.replace(',', '-')}"
    csv_path = args.out_dir / f"gothic_nldas_point_sw_{tag}.csv"
    metrics_path = args.out_dir / f"gothic_nldas_point_sw_metrics_{tag}.csv"
    table.to_csv(csv_path, index=False)
    metrics.to_csv(metrics_path, index=False)
    return table, metrics, csv_path, metrics_path


def main() -> None:
    args = parse_args()
    table, metrics, csv_path, metrics_path = build_point_eval(args)
    print(f"Wrote {csv_path}")
    print(f"Wrote {metrics_path}")
    print(f"Rows: {len(table)}")
    print(f"Columns: {list(table.columns)}")
    print(metrics)


if __name__ == "__main__":
    main()
