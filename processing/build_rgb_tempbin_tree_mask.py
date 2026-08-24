#!/usr/bin/env python3
"""Build a pixelwise cloud mask from temperature-binned decision-tree leaves."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from build_rgb_tempbin_mask import ensure_era
from colorado.apply_tempbin_thresholds_colorado import (
    interpolate_temp_to_goes_grid,
    load_era5_temp_field,
    select_target_hours,
)


COND_RE = re.compile(r"(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)")


def apply_leaf(values: dict[str, np.ndarray], rule: str) -> np.ndarray:
    conditions = COND_RE.findall(rule)
    if not conditions:
        raise ValueError(f"No conditions parsed from tree rule: {rule}")
    result = np.ones(values["red"].shape, dtype=bool)
    for feature, operator, threshold_text in conditions:
        threshold = float(threshold_text)
        if operator == ">":
            result &= values[feature] > threshold
        else:
            result &= values[feature] <= threshold
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("rgb", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--era-dir", type=Path, required=True)
    parser.add_argument("--rules", type=Path, required=True)
    parser.add_argument(
        "--binary-only", action="store_true",
        help="Store only the packed cloud mask, omitting reproducible diagnostics",
    )
    args = parser.parse_args()

    rules = pd.read_csv(args.rules)
    rules = rules.loc[(rules["status"] == "trained") & (rules["prediction"] == 1)].copy()
    bins = rules[["temp_bin", "temp_left_c", "temp_right_c"]].drop_duplicates()
    bins = bins.sort_values("temp_left_c").reset_index(drop=True)
    if bins.empty:
        raise ValueError(f"No trained cloudy leaves in {args.rules}")

    era_path = ensure_era(args.rgb, args.era_dir)
    t2m = load_era5_temp_field(era_path)
    with xr.open_dataset(args.rgb) as source:
        ds = select_target_hours(source, target_hours={0, 1, *range(14, 24)}).load()

    times = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
    lat = np.asarray(ds.latitude, dtype=np.float64)
    lon = np.asarray(ds.longitude, dtype=np.float64)
    temp = np.asarray(
        interpolate_temp_to_goes_grid(t2m, times, lat, lon).values,
        dtype=np.float32,
    )
    values = {name: np.asarray(ds[name], dtype=np.float32) for name in ("red", "green", "blue")}
    cloud = np.zeros(values["red"].shape, dtype=np.uint8)
    bin_index = np.full(cloud.shape, -1, dtype=np.int8)

    for index, row in bins.iterrows():
        selected = (temp >= float(row.temp_left_c)) & (temp < float(row.temp_right_c))
        bin_index[selected] = index
        if not selected.any():
            continue
        selected_values = {name: array[selected] for name, array in values.items()}
        cloudy = np.zeros(selected_values["red"].shape, dtype=bool)
        leaves = rules.loc[rules["temp_bin"] == row.temp_bin, "rule"]
        for leaf in leaves:
            cloudy |= apply_leaf(selected_values, str(leaf))
        cloud[selected] = cloudy.astype(np.uint8)

    invalid = ~np.isfinite(temp)
    for array in values.values():
        invalid |= ~np.isfinite(array)
    cloud[invalid] = np.uint8(255)

    data_vars = {"cloud_binary": (("t", "latitude", "longitude"), cloud)}
    if not args.binary_only:
        data_vars.update({
            "air_temp_c": (("t", "latitude", "longitude"), temp),
            "temp_bin_index": (("t", "latitude", "longitude"), bin_index),
        })
    out = xr.Dataset(
        data_vars=data_vars,
        coords={"t": ds.t, "latitude": ds.latitude, "longitude": ds.longitude},
        attrs={
            "title": "Cloud mask from exact Gothic temperature-bin decision-tree paths",
            "rgb_source": str(args.rgb),
            "era5_land_source": str(era_path),
            "rules_csv": str(args.rules),
            "tree_logic": "AND within each leaf; OR across cloudy leaves",
            "temp_bin_labels": ";".join(f"{i}:{label}" for i, label in enumerate(bins.temp_bin)),
        },
    )
    out.cloud_binary.attrs.update(
        long_name="binary cloud mask", flag_values=[0, 1, 255],
        flag_meanings="clear cloudy missing",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".nc.part")
    out.to_netcdf(
        temporary,
        encoding={name: encoding for name, encoding in {
            "cloud_binary": {"zlib": True, "complevel": 4, "dtype": "uint8", "_FillValue": np.uint8(255)},
            "air_temp_c": {"zlib": True, "complevel": 4, "dtype": "float32"},
            "temp_bin_index": {"zlib": True, "complevel": 4, "dtype": "int8", "_FillValue": np.int8(-1)},
        }.items() if name in out},
    )
    temporary.replace(args.output)
    valid = cloud != 255
    print(f"Wrote exact-tree mask: {args.output} (cloud fraction={cloud[valid].mean():.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
