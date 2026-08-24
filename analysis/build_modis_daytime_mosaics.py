#!/usr/bin/env python3
"""Mosaic adjacent daytime MOD35/MYD35 granules over the Gothic TSI domain.

One compressed NPZ is written per satellite overpass. The companion CSV has
one row per mosaic, including cloud fraction, contributing granules, and a
geometric full-domain-coverage flag. Granules from the same product whose
start times are no more than 10 minutes apart are treated as one overpass.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC


LAT_MIN = 38.91617092558588
LAT_MAX = 39.00624907441412
LON_MIN = -107.05032765812774
LON_MAX = -106.93495234187225


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--modis-root", default="/glade/u/home/cdalden/scratch/colorado_modis"
    )
    parser.add_argument(
        "--manifest",
        default="/glade/u/home/cdalden/scratch/colorado_modis/cmr_manifest.jsonl",
    )
    parser.add_argument(
        "--output-dir",
        default="/glade/u/home/cdalden/scratch/colorado_modis/daytime_mosaics",
    )
    parser.add_argument(
        "--output-csv",
        default="/glade/u/home/cdalden/goes_work/analysis/output_08c/"
        "colorado_modis_daytime_mosaics.csv",
    )
    parser.add_argument("--lat-min", type=float, default=LAT_MIN)
    parser.add_argument("--lat-max", type=float, default=LAT_MAX)
    parser.add_argument("--lon-min", type=float, default=LON_MIN)
    parser.add_argument("--lon-max", type=float, default=LON_MAX)
    parser.add_argument("--max-gap-minutes", type=float, default=10.0)
    parser.add_argument("--progress-every", type=int, default=50)
    return parser.parse_args()


def load_daytime_manifest(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as stream:
        for line in stream:
            entry = json.loads(line)
            granule_id = entry.get("producer_granule_id", "")
            if entry.get("day_night_flag") == "DAY" and granule_id.startswith(
                ("MOD35_L2", "MYD35_L2")
            ):
                rows.append(entry)
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise RuntimeError(f"No daytime MODIS entries found in {path}")
    frame["time_start"] = pd.to_datetime(frame["time_start"], utc=True)
    frame["time_end"] = pd.to_datetime(frame["time_end"], utc=True)
    frame["product"] = frame["producer_granule_id"].str.split(".").str[0]
    return frame.sort_values(["product", "time_start"]).reset_index(drop=True)


def assign_overpasses(frame: pd.DataFrame, max_gap_minutes: float) -> pd.DataFrame:
    pieces = []
    max_gap = pd.Timedelta(minutes=max_gap_minutes)
    next_id = 0
    for _, group in frame.groupby("product", sort=True):
        group = group.sort_values("time_start").copy()
        new_pass = group["time_start"].diff().isna() | (
            group["time_start"].diff() > max_gap
        )
        local_id = new_pass.cumsum() - 1 + next_id
        group["overpass_id"] = local_id.astype(int)
        next_id = int(local_id.max()) + 1
        pieces.append(group)
    return pd.concat(pieces).sort_values("time_start").reset_index(drop=True)


def archive_path(root: Path, granule_id: str) -> Path:
    product, acquisition = granule_id.split(".")[:2]
    stamp = pd.to_datetime(acquisition[1:], format="%Y%j")
    return root / product / f"{stamp.year:04d}" / f"{stamp.dayofyear:03d}" / granule_id


def footprint(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    boundary_lat = np.concatenate(
        [lat[0, :], lat[1:, -1], lat[-1, -2::-1], lat[-2:0:-1, 0]]
    )
    boundary_lon = np.concatenate(
        [lon[0, :], lon[1:, -1], lon[-1, -2::-1], lon[-2:0:-1, 0]]
    )
    valid = (
        np.isfinite(boundary_lat)
        & np.isfinite(boundary_lon)
        & (np.abs(boundary_lat) <= 90)
        & (np.abs(boundary_lon) <= 180)
    )
    return np.column_stack([boundary_lon[valid], boundary_lat[valid]])


def point_in_polygon(x: float, y: float, polygon: np.ndarray) -> bool:
    if len(polygon) < 3:
        return False
    inside = False
    x1, y1 = polygon[-1]
    for x2, y2 in polygon:
        if (y1 > y) != (y2 > y):
            crossing_x = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < crossing_x:
                inside = not inside
        x1, y1 = x2, y2
    return inside


def read_domain_cells(
    path: Path, lat_min: float, lat_max: float, lon_min: float, lon_max: float
) -> dict[str, object]:
    sd = SD(str(path), SDC.READ)
    try:
        lat = sd.select("Latitude").get().astype(float)
        lon = sd.select("Longitude").get().astype(float)
        cloud_mask = sd.select("Cloud_Mask").get().astype(np.uint8)
    finally:
        sd.end()

    valid_geo = (
        np.isfinite(lat)
        & np.isfinite(lon)
        & (np.abs(lat) <= 90)
        & (np.abs(lon) <= 180)
    )
    # Find a small coarse-grid window around Gothic, then bilinearly interpolate
    # geolocation to the native 1 km mask pixels. MOD35 geolocation samples are
    # centered at 1 km pixel indices 2, 7, 12, ... (HDF sampling 3, 8, 13, ...).
    nearby = (
        valid_geo
        & (lat >= lat_min - 0.5)
        & (lat <= lat_max + 0.5)
        & (lon >= lon_min - 0.5)
        & (lon <= lon_max + 0.5)
    )
    if not nearby.any():
        return {
            "latitude": np.array([], dtype=float),
            "longitude": np.array([], dtype=float),
            "cloudy_count": np.array([], dtype=np.uint8),
            "determined_count": np.array([], dtype=np.uint8),
            "footprint": footprint(lat, lon),
        }

    coarse_rows, coarse_cols = np.where(nearby)
    along_5km, across_5km = lat.shape
    row_min = max(int(coarse_rows.min()) - 2, 0)
    row_max = min(int(coarse_rows.max()) + 2, along_5km - 1)
    col_min = max(int(coarse_cols.min()) - 2, 0)
    col_max = min(int(coarse_cols.max()) + 2, across_5km - 1)

    pixel_rows = np.arange(
        max(0, 2 + 5 * row_min), min(cloud_mask.shape[1], 3 + 5 * row_max)
    )
    pixel_cols = np.arange(
        max(0, 2 + 5 * col_min), min(cloud_mask.shape[2], 3 + 5 * col_max)
    )
    row_position = np.clip((pixel_rows - 2) / 5.0, 0, along_5km - 1)
    col_position = np.clip((pixel_cols - 2) / 5.0, 0, across_5km - 1)
    row0 = np.minimum(np.floor(row_position).astype(int), along_5km - 2)
    col0 = np.minimum(np.floor(col_position).astype(int), across_5km - 2)
    row_weight = (row_position - row0)[:, None]
    col_weight = (col_position - col0)[None, :]

    def bilinear(values: np.ndarray) -> np.ndarray:
        v00 = values[np.ix_(row0, col0)]
        v01 = values[np.ix_(row0, col0 + 1)]
        v10 = values[np.ix_(row0 + 1, col0)]
        v11 = values[np.ix_(row0 + 1, col0 + 1)]
        return (
            v00 * (1 - row_weight) * (1 - col_weight)
            + v01 * (1 - row_weight) * col_weight
            + v10 * row_weight * (1 - col_weight)
            + v11 * row_weight * col_weight
        )

    latitude_1km = bilinear(lat)
    longitude_1km = bilinear(lon)
    in_domain = (
        np.isfinite(latitude_1km)
        & np.isfinite(longitude_1km)
        & (latitude_1km >= lat_min)
        & (latitude_1km <= lat_max)
        & (longitude_1km >= lon_min)
        & (longitude_1km <= lon_max)
    )
    byte1 = cloud_mask[0][np.ix_(pixel_rows, pixel_cols)]
    determined = (byte1 & 0b1).astype(bool)
    quality = (byte1 >> 1) & 0b11
    cloudy = (quality <= 1) & determined

    return {
        "latitude": latitude_1km[in_domain],
        "longitude": longitude_1km[in_domain],
        "cloudy_count": cloudy[in_domain].astype(np.uint8),
        "determined_count": determined[in_domain].astype(np.uint8),
        "footprint": footprint(lat, lon),
    }


def mosaic_group(
    group: pd.DataFrame,
    root: Path,
    output_dir: Path,
    bounds: tuple[float, float, float, float],
) -> dict[str, object]:
    lat_min, lat_max, lon_min, lon_max = bounds
    product = group["product"].iloc[0]
    start = group["time_start"].min()
    end = group["time_end"].max()
    granule_ids = group["producer_granule_id"].tolist()

    cell_parts = []
    polygons = []
    source_paths = []
    for granule_id in granule_ids:
        path = archive_path(root, granule_id)
        if not path.exists():
            raise FileNotFoundError(path)
        cells = read_domain_cells(path, lat_min, lat_max, lon_min, lon_max)
        polygons.append(cells.pop("footprint"))
        cell_parts.append(cells)
        source_paths.append(str(path))

    if cell_parts:
        latitude = np.concatenate([part["latitude"] for part in cell_parts])
        longitude = np.concatenate([part["longitude"] for part in cell_parts])
        cloudy_count = np.concatenate([part["cloudy_count"] for part in cell_parts])
        determined_count = np.concatenate(
            [part["determined_count"] for part in cell_parts]
        )
    else:
        latitude = longitude = np.array([], dtype=float)
        cloudy_count = determined_count = np.array([], dtype=np.int16)

    # Remove duplicate 1 km pixel centers if neighboring granules overlap.
    if latitude.size:
        keys = np.round(np.column_stack([latitude, longitude]), 4)
        _, unique_index = np.unique(keys, axis=0, return_index=True)
        unique_index.sort()
        latitude = latitude[unique_index]
        longitude = longitude[unique_index]
        cloudy_count = cloudy_count[unique_index]
        determined_count = determined_count[unique_index]

    corners = [
        (lon_min, lat_min),
        (lon_min, lat_max),
        (lon_max, lat_min),
        (lon_max, lat_max),
    ]
    full_geo = all(
        any(point_in_polygon(x, y, polygon) for polygon in polygons)
        for x, y in corners
    )
    possible_pixels = int(latitude.size)
    valid_pixels = int(determined_count.sum())
    cloudy_pixels = int(cloudy_count.sum())
    cloud_fraction = cloudy_pixels / valid_pixels if valid_pixels else np.nan
    valid_fraction = valid_pixels / possible_pixels if possible_pixels else 0.0

    platform = "Terra" if product == "MOD35_L2" else "Aqua"
    filename = (
        f"{product}_{start.strftime('%Y%m%dT%H%M%S')}_"
        f"{end.strftime('%H%M%S')}_gothic_mosaic.npz"
    )
    output_path = output_dir / product / f"{start.year:04d}" / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        latitude=latitude.astype(np.float32),
        longitude=longitude.astype(np.float32),
        cloudy_count=cloudy_count.astype(np.uint8),
        determined_count=determined_count.astype(np.uint8),
        source_granules=np.asarray(granule_ids),
        time_start=np.asarray(start.isoformat()),
        time_end=np.asarray(end.isoformat()),
        bounds=np.asarray(bounds, dtype=np.float64),
    )

    midpoint = start + (end - start) / 2
    return {
        "date": start.tz_convert(None).normalize(),
        "t": midpoint.tz_convert(None),
        "time_start": start.tz_convert(None),
        "time_end": end.tz_convert(None),
        "product": product,
        "platform": platform,
        "n_granules": len(granule_ids),
        "granule_ids": ";".join(granule_ids),
        "source_paths": ";".join(source_paths),
        "mosaic_path": str(output_path),
        "modis_cloud_fraction": cloud_fraction,
        "modis_cloudy_pixels": cloudy_pixels,
        "modis_valid_pixels": valid_pixels,
        "modis_possible_pixels": possible_pixels,
        "modis_valid_fraction": valid_fraction,
        "modis_n_cells": int(latitude.size),
        "full_geolocation_coverage": bool(full_geo),
    }


def main() -> int:
    args = parse_args()
    root = Path(args.modis_root)
    output_dir = Path(args.output_dir)
    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    entries = assign_overpasses(
        load_daytime_manifest(Path(args.manifest)), args.max_gap_minutes
    )
    groups = list(entries.groupby("overpass_id", sort=False))
    bounds = (args.lat_min, args.lat_max, args.lon_min, args.lon_max)

    rows = []
    for index, (_, group) in enumerate(groups, 1):
        rows.append(mosaic_group(group, root, output_dir, bounds))
        if index % max(args.progress_every, 1) == 0 or index == len(groups):
            print(f"Built {index}/{len(groups)} daytime overpass mosaics", flush=True)

    frame = pd.DataFrame(rows).sort_values("t").reset_index(drop=True)
    frame.to_csv(output_csv, index=False)
    print(f"Wrote {len(frame)} mosaic rows to {output_csv}")
    print(
        "Coverage: "
        f"{int(frame.full_geolocation_coverage.sum())}/{len(frame)} full-domain; "
        f"{int((frame.n_granules > 1).sum())} used adjacent granules"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
