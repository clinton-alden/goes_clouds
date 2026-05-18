#!/usr/bin/env python3
"""Build Mammoth domain cloud-fraction time series from GOES RGB files.

This script trains temperature-binned RGB cloud models from the Notebook 11b
training table, then applies the appropriate model to every pixel for each GOES
time slice and writes a CSV of domain cloud fraction.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from sklearn.tree import DecisionTreeClassifier

TEMP_BIN_WIDTH_C = 5.0
MIN_ROWS_PER_TEMP_BIN = 40
MIN_CLASS_COUNT_PER_TEMP_BIN = 10
DT_PARAMS = dict(max_depth=3, min_samples_split=20, class_weight="balanced", random_state=42)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Mammoth domain cloud-fraction CSV from GOES RGB files")
    parser.add_argument(
        "--train-table",
        default="/glade/u/home/cdalden/goes_work/analysis/output_11b_mammoth/mammoth_rgb_sw_training_table_2024.csv",
        help="Path to Notebook 11b training table CSV",
    )
    parser.add_argument(
        "--rgb-glob",
        default="/glade/u/home/cdalden/scratch/mammoth/goes18/rgb_composite/goes18_C02_C05_C13_rgb_mammoth_2024*.nc",
        help="Glob pattern for Mammoth RGB netCDF files",
    )
    parser.add_argument(
        "--output-csv",
        default="/glade/u/home/cdalden/goes_work/analysis/output_13_mammoth_domain_cloud_fraction/mammoth_domain_cloud_fraction_timeseries.csv",
        help="Output CSV path",
    )
    parser.add_argument(
        "--progress-every-files",
        type=int,
        default=30,
        help="Print progress every N files",
    )
    return parser.parse_args()


def build_temp_bin_models(train_df: pd.DataFrame):
    bin_df = train_df.copy()
    tmin = np.floor(min(bin_df["air_temp_c"].min(), 0.0) / TEMP_BIN_WIDTH_C) * TEMP_BIN_WIDTH_C
    tmax = np.ceil(max(bin_df["air_temp_c"].max(), 0.0) / TEMP_BIN_WIDTH_C) * TEMP_BIN_WIDTH_C
    edges = np.arange(tmin, tmax + TEMP_BIN_WIDTH_C, TEMP_BIN_WIDTH_C)
    bin_df["temp_bin"] = pd.cut(bin_df["air_temp_c"], bins=edges, right=False, include_lowest=True)

    models = {}
    for temp_bin, grp in bin_df.groupby("temp_bin", observed=False):
        if pd.isna(temp_bin):
            continue
        n_rows = len(grp)
        counts = grp["cloud_binary"].value_counts()
        n_clear = int(counts.get(0, 0))
        n_cloudy = int(counts.get(1, 0))
        if n_rows < MIN_ROWS_PER_TEMP_BIN or min(n_clear, n_cloudy) < MIN_CLASS_COUNT_PER_TEMP_BIN:
            continue

        clf = DecisionTreeClassifier(**DT_PARAMS)
        # Fit on numpy arrays so predict() on numpy features does not emit feature-name warnings.
        clf.fit(grp[["red", "green", "blue"]].to_numpy(), grp["cloud_binary"].to_numpy())
        models[str(temp_bin)] = {"bin": temp_bin, "model": clf}

    return models


def pick_model_key(air_temp_c: float, temp_bin_models: dict) -> str | None:
    for key, payload in temp_bin_models.items():
        temp_bin = payload["bin"]
        if (air_temp_c >= temp_bin.left) and (air_temp_c < temp_bin.right):
            return key
    return None


def predict_mask_for_frame(red2d: np.ndarray, green2d: np.ndarray, blue2d: np.ndarray, model):
    valid = np.isfinite(red2d) & np.isfinite(green2d) & np.isfinite(blue2d)
    mask = np.full(red2d.shape, np.nan, dtype=float)
    if valid.any():
        feats = np.column_stack([red2d[valid], green2d[valid], blue2d[valid]])
        pred = model.predict(feats)
        mask[valid] = pred.astype(float)
    return mask, valid


def nearest_air_temp(temp_lookup: pd.DataFrame, tstamp: pd.Timestamp) -> float:
    j = np.searchsorted(temp_lookup["time"].values, np.datetime64(tstamp), side="left")
    if j == 0:
        return float(temp_lookup.iloc[0]["air_temp_c"])
    if j >= len(temp_lookup):
        return float(temp_lookup.iloc[-1]["air_temp_c"])

    t0 = temp_lookup.iloc[j - 1]["time"]
    t1 = temp_lookup.iloc[j]["time"]
    a0 = temp_lookup.iloc[j - 1]["air_temp_c"]
    a1 = temp_lookup.iloc[j]["air_temp_c"]
    return float(a0 if abs(tstamp - t0) <= abs(t1 - tstamp) else a1)


def main() -> None:
    args = parse_args()

    train_table = Path(args.train_table)
    output_csv = Path(args.output_csv)
    rgb_files = sorted(glob.glob(args.rgb_glob))

    if not train_table.exists():
        raise FileNotFoundError(f"Training table not found: {train_table}")
    if not rgb_files:
        raise FileNotFoundError(f"No RGB files found for pattern: {args.rgb_glob}")

    # Skip zero-byte placeholders/corrupt outputs to avoid immediate job failure.
    rgb_files_nonempty = [fp for fp in rgb_files if Path(fp).stat().st_size > 0]
    n_empty = len(rgb_files) - len(rgb_files_nonempty)
    if not rgb_files_nonempty:
        raise RuntimeError("All matched RGB files are empty (0 bytes).")
    if n_empty:
        print(f"Skipping {n_empty} empty RGB files.")
    rgb_files = rgb_files_nonempty

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    train_df = pd.read_csv(train_table, parse_dates=["time"])
    train_df["time"] = pd.to_datetime(train_df["time"], errors="coerce").astype("datetime64[ns]")
    train_df = train_df.dropna(subset=["time", "red", "green", "blue", "cloud_binary", "air_temp_c"]).copy()
    train_df["cloud_binary"] = train_df["cloud_binary"].astype(int)

    temp_bin_models = build_temp_bin_models(train_df)
    if not temp_bin_models:
        raise RuntimeError("No trained temperature-bin models were created.")

    temp_lookup = (
        train_df[["time", "air_temp_c"]]
        .dropna()
        .sort_values("time")
        .drop_duplicates("time")
        .reset_index(drop=True)
    )

    header = ["time", "air_temp_c", "model_bin", "n_valid", "n_cloudy", "cloud_fraction_pct"]
    n_rows_written = 0

    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()

        for i, fp in enumerate(rgb_files, start=1):
            try:
                with xr.open_dataset(fp, engine="netcdf4") as ds:
                    if "t" not in ds or "red" not in ds or "green" not in ds or "blue" not in ds:
                        print(f"Skipping file with missing vars: {fp}")
                        continue

                    tvals = pd.to_datetime(np.asarray(ds["t"].values), errors="coerce")
                    for k, tstamp in enumerate(tvals):
                        if pd.isna(tstamp):
                            continue

                        air_temp_c = nearest_air_temp(temp_lookup, tstamp)
                        model_key = pick_model_key(air_temp_c, temp_bin_models)
                        if model_key is None:
                            writer.writerow(
                                {
                                    "time": tstamp,
                                    "air_temp_c": air_temp_c,
                                    "model_bin": "",
                                    "n_valid": 0,
                                    "n_cloudy": 0,
                                    "cloud_fraction_pct": "",
                                }
                            )
                            n_rows_written += 1
                            continue

                        model = temp_bin_models[model_key]["model"]
                        red2d = np.asarray(ds["red"].isel(t=k).values, dtype=float)
                        green2d = np.asarray(ds["green"].isel(t=k).values, dtype=float)
                        blue2d = np.asarray(ds["blue"].isel(t=k).values, dtype=float)

                        cloud_mask, valid = predict_mask_for_frame(red2d, green2d, blue2d, model)
                        n_valid = int(valid.sum())
                        n_cloudy = int(np.nansum(cloud_mask == 1)) if n_valid > 0 else 0
                        cloud_fraction_pct = (100.0 * n_cloudy / n_valid) if n_valid > 0 else np.nan

                        writer.writerow(
                            {
                                "time": tstamp,
                                "air_temp_c": air_temp_c,
                                "model_bin": model_key,
                                "n_valid": n_valid,
                                "n_cloudy": n_cloudy,
                                "cloud_fraction_pct": cloud_fraction_pct,
                            }
                        )
                        n_rows_written += 1
            except Exception as exc:
                print(f"Skipping unreadable file: {fp} ({exc})")
                continue

            if i % max(args.progress_every_files, 1) == 0 or i == len(rgb_files):
                f.flush()
                os.fsync(f.fileno())
                print(f"Processed {i}/{len(rgb_files)} files; rows written={n_rows_written}")

    print(f"Wrote cloud-fraction CSV: {output_csv}")


if __name__ == "__main__":
    main()
