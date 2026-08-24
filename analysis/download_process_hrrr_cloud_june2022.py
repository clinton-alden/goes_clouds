#!/usr/bin/env python3
"""Download byte-ranged HRRR TCDC records and subset them over Colorado."""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path
import tempfile
import time

import numpy as np
import pandas as pd
import requests
import xarray as xr


BASE_URL = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
VALID_HOURS = (0, 1, *range(14, 24))
COMPARISONS = ("analysis_f00", "day1_00z", "day2_00z")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", required=True, help="Valid date YYYYMMDD")
    parser.add_argument("--lon-min", type=float, default=-109.0)
    parser.add_argument("--lon-max", type=float, default=-104.0)
    parser.add_argument("--lat-min", type=float, default=37.0)
    parser.add_argument("--lat-max", type=float, default=41.0)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/glade/derecho/scratch/cdalden/hrrr_goes_june2022/hrrr_daily"),
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--retries", type=int, default=4)
    return parser.parse_args()


def source_for(valid: datetime, comparison: str) -> tuple[datetime, int]:
    if comparison == "analysis_f00":
        return valid.replace(minute=0, second=0, microsecond=0), 0
    if comparison == "day1_00z":
        init = valid.replace(hour=0, minute=0, second=0, microsecond=0)
        return init, valid.hour
    if comparison == "day2_00z":
        init = valid.replace(hour=0, minute=0, second=0, microsecond=0) - timedelta(days=1)
        return init, valid.hour + 24
    raise ValueError(comparison)


def object_url(init: datetime, lead: int) -> str:
    key = (
        f"hrrr.{init:%Y%m%d}/conus/"
        f"hrrr.t{init:%H}z.wrfsfcf{lead:02d}.grib2"
    )
    return f"{BASE_URL}/{key}"


def request_with_retries(session: requests.Session, url: str, *, retries: int, headers=None) -> requests.Response:
    for attempt in range(1, retries + 1):
        try:
            response = session.get(url, headers=headers, timeout=(20, 180))
            response.raise_for_status()
            return response
        except requests.RequestException:
            if attempt == retries:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def tcdc_byte_range(session: requests.Session, grib_url: str, retries: int) -> tuple[int, int | None]:
    response = request_with_retries(session, grib_url + ".idx", retries=retries)
    lines = response.text.splitlines()
    for index, line in enumerate(lines):
        parts = line.split(":")
        if len(parts) >= 5 and parts[3] == "TCDC" and parts[4] == "entire atmosphere":
            start = int(parts[1])
            end = int(lines[index + 1].split(":", 2)[1]) - 1 if index + 1 < len(lines) else None
            return start, end
    raise RuntimeError(f"TCDC entire-atmosphere record not found: {grib_url}.idx")


def load_tcdc_subset(
    session: requests.Session,
    init: datetime,
    lead: int,
    bounds: tuple[float, float, float, float],
    retries: int,
) -> xr.DataArray:
    url = object_url(init, lead)
    start, end = tcdc_byte_range(session, url, retries)
    range_value = f"bytes={start}-{'' if end is None else end}"
    response = request_with_retries(session, url, retries=retries, headers={"Range": range_value})
    if response.status_code != 206:
        raise RuntimeError(f"Server ignored byte-range request for {url}: HTTP {response.status_code}")

    with tempfile.NamedTemporaryFile(suffix=".grib2") as temp:
        temp.write(response.content)
        temp.flush()
        with xr.open_dataset(temp.name, engine="cfgrib", backend_kwargs={"indexpath": ""}) as ds:
            variable = "tcc" if "tcc" in ds else next(iter(ds.data_vars))
            da = ds[variable].load()
            lat = np.asarray(ds["latitude"].values)
            lon = np.asarray(ds["longitude"].values)

    lon = np.where(lon > 180.0, lon - 360.0, lon)
    lon_min, lon_max, lat_min, lat_max = bounds
    inside = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    yy, xx = np.where(inside)
    if not len(yy):
        raise RuntimeError("Requested bounds do not intersect the HRRR grid")
    y_slice = slice(max(0, yy.min() - 1), min(da.sizes["y"], yy.max() + 2))
    x_slice = slice(max(0, xx.min() - 1), min(da.sizes["x"], xx.max() + 2))
    values = np.asarray(da.isel(y=y_slice, x=x_slice).values, dtype="float32")
    da = xr.DataArray(
        values,
        dims=("y", "x"),
        coords={
            "latitude": (("y", "x"), lat[y_slice, x_slice].astype("float32")),
            "longitude": (("y", "x"), lon[y_slice, x_slice].astype("float32")),
        },
        name="hrrr_tcc",
    )
    da.attrs.update(long_name="HRRR total cloud cover", units="percent")
    return da


def main() -> int:
    args = parse_args()
    valid_date = datetime.strptime(args.date, "%Y%m%d")
    output = args.output_dir / f"hrrr_tcc_colorado_{args.date}.nc"
    if output.exists() and not args.overwrite:
        print(f"[skip] {output}")
        return 0
    args.output_dir.mkdir(parents=True, exist_ok=True)
    bounds = (args.lon_min, args.lon_max, args.lat_min, args.lat_max)

    session = requests.Session()
    by_comparison = []
    init_matrix = []
    lead_matrix = []
    valid_times = [valid_date + timedelta(hours=hour) for hour in VALID_HOURS]
    for comparison in COMPARISONS:
        hourly = []
        inits = []
        leads = []
        for valid in valid_times:
            init, lead = source_for(valid, comparison)
            print(f"[{args.date}] {comparison} valid={valid:%Y-%m-%dT%H} init={init:%Y-%m-%dT%H} f{lead:02d}", flush=True)
            hourly.append(load_tcdc_subset(session, init, lead, bounds, args.retries))
            inits.append(np.datetime64(init, "ns"))
            leads.append(lead)
        by_comparison.append(xr.concat(hourly, dim=pd.Index(valid_times, name="valid_time")))
        init_matrix.append(inits)
        lead_matrix.append(leads)

    tcc = xr.concat(by_comparison, dim=pd.Index(COMPARISONS, name="comparison"))
    dataset = tcc.to_dataset(name="hrrr_tcc")
    dataset["init_time"] = (("comparison", "valid_time"), np.asarray(init_matrix, dtype="datetime64[ns]"))
    dataset["lead_hour"] = (("comparison", "valid_time"), np.asarray(lead_matrix, dtype="int16"))
    dataset.attrs.update(
        title="HRRR total-cloud-cover comparisons over the Colorado GOES domain",
        source="NOAA HRRR AWS open-data archive; byte-ranged TCDC records",
        bounds_wsen=f"{args.lon_min},{args.lat_min},{args.lon_max},{args.lat_max}",
    )
    encoding = {
        "hrrr_tcc": {"zlib": True, "complevel": 4, "dtype": "float32"},
        "lead_hour": {"zlib": True, "complevel": 4, "dtype": "int16"},
    }
    temp_output = output.with_suffix(".tmp.nc")
    dataset.to_netcdf(temp_output, encoding=encoding)
    temp_output.replace(output)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
