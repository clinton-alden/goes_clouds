#!/usr/bin/env python3
"""Build MODIS daytime cloud-fraction time series over the Colorado domain.

This reads the downloaded MOD35_L2 and MYD35_L2 HDF4 granules under the local
`colorado_modis` archive, filters to daytime passes using the CMR manifest,
and writes one CSV row per daytime granule with the domain-mean cloud fraction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from pyhdf.SD import SD, SDC
except ImportError as exc:  # pragma: no cover - clearer failure in user environments
    raise SystemExit(
        "pyhdf is required to read MODIS HDF4 granules. Install it in the workspace environment first."
    ) from exc


LAT_MIN = 38.91617092558588
LAT_MAX = 39.00624907441412
LON_MIN = -107.05032765812774
LON_MAX = -106.93495234187225
MIN_PIXEL_COUNT = 50


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--modis-root",
        default="/glade/u/home/cdalden/scratch/colorado_modis",
        help="Root directory containing the downloaded MOD35_L2 and MYD35_L2 archives.",
    )
    parser.add_argument(
        "--manifest",
        default="/glade/u/home/cdalden/scratch/colorado_modis/cmr_manifest.jsonl",
        help="CMR manifest written by the downloader.",
    )
    parser.add_argument(
        "--output-csv",
        default="/glade/u/home/cdalden/goes_work/analysis/output_08c/colorado_modis_daytime_cloud_fraction.csv",
        help="Output CSV path.",
    )
    parser.add_argument("--lat-min", type=float, default=LAT_MIN)
    parser.add_argument("--lat-max", type=float, default=LAT_MAX)
    parser.add_argument("--lon-min", type=float, default=LON_MIN)
    parser.add_argument("--lon-max", type=float, default=LON_MAX)
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="Print progress every N processed daytime granules.",
    )
    return parser.parse_args()


def load_manifest(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as stream:
        for line in stream:
            if not line.strip():
                continue
            entry = json.loads(line)
            if entry.get("day_night_flag") != "DAY":
                continue
            granule_id = entry.get("producer_granule_id", "")
            if not granule_id.startswith(("MOD35_L2", "MYD35_L2")):
                continue
            rows.append(entry)
    if not rows:
        raise RuntimeError(f"No daytime MODIS granules found in {path}")
    df = pd.DataFrame(rows)
    df["time_start"] = pd.to_datetime(df["time_start"], utc=True)
    df["time_end"] = pd.to_datetime(df["time_end"], utc=True)
    df = df.sort_values("time_start").reset_index(drop=True)
    return df


def archive_path(root: Path, granule_id: str) -> Path:
    product, acquisition = granule_id.split(".")[:2]
    stamp = pd.to_datetime(acquisition[1:], format="%Y%j")
    return root / product / f"{stamp.year:04d}" / f"{stamp.dayofyear:03d}" / granule_id


def read_sds(sd: SD, name: str) -> np.ndarray:
    return sd.select(name).get()


def cloud_fraction_for_file(
    path: Path,
    *,
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
) -> dict[str, object]:
    sd = SD(str(path), SDC.READ)
    try:
        lat = read_sds(sd, "Latitude").astype(float)
        lon = read_sds(sd, "Longitude").astype(float)
        cloud_mask = read_sds(sd, "Cloud_Mask").astype(np.uint8)
    finally:
        sd.end()

    if cloud_mask.ndim != 3 or cloud_mask.shape[0] < 1:
        raise ValueError(f"Unexpected Cloud_Mask shape in {path}: {cloud_mask.shape}")

    along_5km, across_5km = lat.shape
    expected_along = along_5km * 5
    expected_across = across_5km * 5
    if cloud_mask.shape[1] < expected_along or cloud_mask.shape[2] < expected_across:
        raise ValueError(
            f"Cloud_Mask shape {cloud_mask.shape} does not cover the 5 km geolocation grid {lat.shape} in {path}"
        )

    # Byte 1 bits 1-2 encode the unobstructed FOV quality flag.
    # 00 = cloudy, 01 = uncertain, 10 = probably clear, 11 = confident clear.
    quality_flag = (cloud_mask[0, :expected_along, :expected_across] >> 1) & 0b11
    cloudy_1km = quality_flag <= 1

    # Aggregate the 1 km cloud mask onto the 5 km geolocation grid.
    cloudy_5km = cloudy_1km.reshape(along_5km, 5, across_5km, 5).mean(axis=(1, 3))
    valid_geo = np.isfinite(lat) & np.isfinite(lon)
    domain = valid_geo & (lat >= lat_min) & (lat <= lat_max) & (lon >= lon_min) & (lon <= lon_max)

    if not domain.any():
        return {
            "modis_cloud_fraction": np.nan,
            "modis_n_pixels": 0,
            "modis_n_cells": 0,
        }

    return {
        "modis_cloud_fraction": float(cloudy_5km[domain].mean()),
        "modis_n_pixels": int(domain.sum()) * 25,
        "modis_n_cells": int(domain.sum()),
    }


def main() -> int:
    args = parse_args()
    root = Path(args.modis_root)
    manifest = Path(args.manifest)
    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    entries = load_manifest(manifest)
    rows: list[dict[str, object]] = []

    for index, entry in enumerate(entries.to_dict("records"), start=1):
        granule_id = entry["producer_granule_id"]
        path = archive_path(root, granule_id)
        if not path.exists():
            print(f"Skipping missing file: {path}")
            continue

        try:
            stats = cloud_fraction_for_file(
                path,
                lat_min=args.lat_min,
                lat_max=args.lat_max,
                lon_min=args.lon_min,
                lon_max=args.lon_max,
            )
        except Exception as exc:
            print(f"Skipping {granule_id}: {exc}")
            continue

        if stats["modis_n_pixels"] < MIN_PIXEL_COUNT:
            continue

        rows.append(
            {
                "date": pd.to_datetime(entry["time_start"], utc=True).tz_convert(None).normalize(),
                "t": pd.to_datetime(entry["time_start"], utc=True).tz_convert(None),
                "time_start": pd.to_datetime(entry["time_start"], utc=True).tz_convert(None),
                "time_end": pd.to_datetime(entry["time_end"], utc=True).tz_convert(None),
                "day_night_flag": entry["day_night_flag"],
                "product": granule_id.split(".")[0],
                "granule_id": granule_id,
                "modis_path": str(path),
                **stats,
            }
        )

        if index % max(args.progress_every, 1) == 0 or index == len(entries):
            print(f"Processed {index}/{len(entries)} daytime granules; rows collected={len(rows)}")

    if not rows:
        raise RuntimeError("No MODIS cloud-fraction rows were produced.")

    out_df = pd.DataFrame(rows).sort_values("t").reset_index(drop=True)
    out_df.to_csv(output_csv, index=False)
    print(f"Wrote {len(out_df)} rows to {output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())