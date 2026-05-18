#!/usr/bin/env python3
"""Download selected GFS surface fields and save a cropped hourly NetCDF."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import requests
import xarray as xr


ROOT = Path("/glade/u/home/cdalden/goes_work")
CYCLE = "00"
FORECAST_HOURS = list(range(24))
PRODUCT_TEMPLATE = "pgrb2.0p25.f{forecast_hour:03d}"
LON_MIN = -109.0
LON_MAX = -104.0
LAT_MIN = 37.0
LAT_MAX = 41.0
GFS_LON_MIN = LON_MIN % 360
GFS_LON_MAX = LON_MAX % 360
FIELDS = {
    ("TMP", "surface"),
    ("WEASD", "surface"),
    ("SNOD", "surface"),
    ("CSNOW", "surface"),
    ("LAND", "surface"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build cropped GFS f000-f023 surface temperature/snow NetCDFs."
    )
    parser.add_argument("dates", nargs="+", help="Date(s) as YYYYMMDD")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing GRIB/NetCDF files")
    return parser.parse_args()


def out_dir(date: str) -> Path:
    return ROOT / f"data_download/gfs_surface_{date}"


def netcdf_path(date: str) -> Path:
    return (
        out_dir(date)
        / f"gfs_{date}_00z_f000_f023_surface_temperature_snow_mask_109W_104W_37N_41N.nc"
    )


def gfs_urls(date: str, forecast_hour: int) -> tuple[str, str]:
    product = PRODUCT_TEMPLATE.format(forecast_hour=forecast_hour)
    base = f"https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{date}/{CYCLE}/atmos"
    stem = f"gfs.t{CYCLE}z.{product}"
    return f"{base}/{stem}", f"{base}/{stem}.idx"


def output_paths(date: str, forecast_hour: int) -> tuple[Path, Path]:
    grib_path = out_dir(date) / (
        f"gfs_{date}_{CYCLE}z_f{forecast_hour:03d}_surface_snow_tmp_selected.grib2"
    )
    return grib_path, grib_path.with_suffix(".idx.txt")


def read_idx(idx_url: str) -> list[dict[str, object]]:
    response = requests.get(idx_url, timeout=60)
    response.raise_for_status()
    lines = response.text.splitlines()

    records = []
    for line_idx, line in enumerate(lines):
        parts = line.split(":")
        if len(parts) < 6:
            continue
        end = int(lines[line_idx + 1].split(":")[1]) - 1 if line_idx + 1 < len(lines) else None
        records.append(
            {
                "record": int(parts[0]),
                "start": int(parts[1]),
                "end": end,
                "var": parts[3],
                "level": parts[4],
                "forecast": parts[5],
                "line": line,
            }
        )
    return records


def selected_records(records: list[dict[str, object]]) -> list[dict[str, object]]:
    selected = []
    found = set()
    for record in records:
        key = (record["var"], record["level"])
        if key not in FIELDS or key in found:
            continue
        selected.append(record)
        found.add(key)

    missing = FIELDS - found
    if missing:
        raise ValueError(f"Missing requested GFS fields: {sorted(missing)}")
    return selected


def download_selected_grib(date: str, forecast_hour: int, overwrite: bool) -> Path:
    grib_path, idx_path = output_paths(date, forecast_hour)
    if grib_path.exists() and grib_path.stat().st_size > 0 and not overwrite:
        print(f"Skip existing {grib_path.name}")
        return grib_path

    grib_url, idx_url = gfs_urls(date, forecast_hour)
    records = selected_records(read_idx(idx_url))
    idx_path.write_text("\n".join(str(record["line"]) for record in records) + "\n")

    with grib_path.open("wb") as handle:
        for record in records:
            end = record["end"]
            byte_range = f"bytes={record['start']}-" if end is None else f"bytes={record['start']}-{end}"
            response = requests.get(grib_url, headers={"Range": byte_range}, timeout=120)
            response.raise_for_status()
            handle.write(response.content)

    print(f"Downloaded {grib_path.name} ({grib_path.stat().st_size / 1024**2:.2f} MB)")
    return grib_path


def open_and_subset_grib(path: Path, forecast_hour: int) -> xr.Dataset:
    ds = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"typeOfLevel": "surface"},
        },
    )
    analysis_time = pd.Timestamp(ds["time"].item())
    valid_time = pd.Timestamp(ds["valid_time"].item())

    ds = ds.sel(
        latitude=slice(LAT_MAX, LAT_MIN),
        longitude=slice(GFS_LON_MIN, GFS_LON_MAX),
    ).copy()
    ds = ds.drop_vars(["time", "valid_time"], errors="ignore")
    ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
    ds = ds.sortby("longitude")
    ds = ds.expand_dims(valid_time=[valid_time])
    ds = ds.assign_coords(
        cycle=("valid_time", [CYCLE]),
        forecast_hour=("valid_time", [forecast_hour]),
        analysis_time=("valid_time", [analysis_time]),
    )

    ds["surface_temperature_c"] = ds["t"] - 273.15
    ds["surface_temperature_c"].attrs.update(
        long_name="GFS surface temperature",
        units="degC",
    )
    ds["snow_mask"] = (ds["sdwe"] > 0) | (ds["sde"] > 0) | (ds["csnow"] > 0)
    ds["snow_mask"].attrs.update(
        long_name="Derived GFS snow mask",
        description="True where WEASD > 0, SNOD > 0, or CSNOW > 0",
        units="1",
    )
    return ds


def build_date(date: str, overwrite: bool) -> Path:
    out_dir(date).mkdir(parents=True, exist_ok=True)
    target = netcdf_path(date)
    if target.exists() and target.stat().st_size > 0 and not overwrite:
        print(f"Skip existing {target}")
        return target

    grib_paths = [
        download_selected_grib(date, forecast_hour, overwrite=overwrite)
        for forecast_hour in FORECAST_HOURS
    ]
    domain_ds = xr.concat(
        [
            open_and_subset_grib(path, forecast_hour)
            for path, forecast_hour in zip(grib_paths, FORECAST_HOURS)
        ],
        dim="valid_time",
        coords="minimal",
        compat="override",
    )
    domain_ds.attrs.update(
        title="GFS 00Z surface temperature and snow mask 24-hour forecast over Colorado domain",
        source="NOAA GFS gfs.t00z.pgrb2.0p25.f000-f023 forecast fields",
        date=date,
        cycle=f"{CYCLE}Z",
        forecast_hours=f"{FORECAST_HOURS[0]}-{FORECAST_HOURS[-1]}",
        domain=f"{abs(LON_MIN):.1f}W-{abs(LON_MAX):.1f}W, {LAT_MIN:.1f}N-{LAT_MAX:.1f}N",
        snow_mask_logic="WEASD > 0 or SNOD > 0 or CSNOW > 0",
    )
    domain_ds.to_netcdf(target)
    print(f"Saved {target}")
    print(f"Dimensions: {dict(domain_ds.sizes)}")
    return target


def main() -> int:
    args = parse_args()
    for date in args.dates:
        build_date(date, overwrite=args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
