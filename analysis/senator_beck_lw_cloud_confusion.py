#!/usr/bin/env python3
"""Build a Senator Beck longwave cloud flag from Pyrgeom_W and compare to GOES RGB.

Method notes
------------
- Follows the Word doc's longwave clear-sky-index (CSI) framing:
    CSI = LW_obs / (LW_clear + a + b * LW_obs)
  with cloud = CSI > 1.
- The Word doc in this workspace states the tuned constants were
  ``a = 14 W m^-2`` and ``b = 0``. Those fixed values are applied here.
- The doc excerpt does not include the explicit 17-method clear-sky longwave
  ensemble formulas, so this script uses a Brutsaert-style clear-sky emissivity
  estimate as the clear-sky longwave baseline before applying the fixed doc
  constants.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score
from sklearn.tree import DecisionTreeClassifier


ROOT = Path("/glade/u/home/cdalden/goes_work")
SBSP_CSV = ROOT / "analysis/SBSP_1hr_2010-2024.csv"
RGB_CSV = ROOT / "analysis/output_11c_senator_beck/senator_beck_rgb_domain_mean_all.csv"
GOTHIC_THRESH_CSV = ROOT / "analysis/output_12_rgb_threshold_transfer/gothic_temp_bin_rgb_thresholds.csv"
OUT_DIR = ROOT / "analysis/output_11c_senator_beck_lw_pyrgeom"

SITE_LAT = 37.90
SITE_LON = -107.72
SITE_ALT_M = 3714.0
MERGE_TOLERANCE = pd.Timedelta("30min")

RNG_SEED = 42
TEMP_BIN_WIDTH_C = 5.0
MIN_ROWS_PER_TEMP_BIN = 40
MIN_CLASS_COUNT_PER_TEMP_BIN = 10

DOC_A = 14.0
DOC_B = 0.0


def load_station() -> pd.DataFrame:
    df = pd.read_csv(SBSP_CSV)
    hour_code = df["Hour"].astype(int)
    hours = (hour_code // 100).astype(int)
    minutes = (hour_code % 100).astype(int)

    df["time"] = (
        pd.to_datetime(
            df["Year"].astype(int).astype(str) + "-" + df["DOY"].astype(int).astype(str),
            format="%Y-%j",
        )
        + pd.to_timedelta(hours, unit="h")
        + pd.to_timedelta(minutes, unit="m")
    )

    df["air_temp_c"] = (df["LoAir_Min_C"] + df["LoAir_Max_C"]) / 2.0
    df["rh_frac"] = np.clip(df["Lo_RH"] / 100.0, 0.0, 1.0)
    return df.sort_values("time").reset_index(drop=True)


def load_rgb() -> pd.DataFrame:
    rgb = pd.read_csv(RGB_CSV)
    rgb["time"] = pd.to_datetime(rgb["time"], errors="coerce")
    rgb = rgb.dropna(subset=["time", "red", "green", "blue"]).copy()
    rgb = rgb.iloc[np.argsort(rgb["time"].to_numpy(dtype="datetime64[ns]"))].reset_index(drop=True)

    # Collapse 5-minute GOES samples to one hourly RGB value so the evaluation
    # matches the hourly station Pyrgeom record one-to-one.
    rgb_hourly = (
        rgb.set_index("time")[["red", "green", "blue"]]
        .resample("1h")
        .mean()
        .dropna()
        .reset_index()
    )
    return rgb_hourly


def sort_by_time(df: pd.DataFrame, time_col: str = "time") -> pd.DataFrame:
    out = df.copy()
    out[time_col] = pd.to_datetime(out[time_col], errors="coerce")
    out = out.dropna(subset=[time_col]).copy()
    order = np.argsort(out[time_col].to_numpy(dtype="datetime64[ns]"))
    return out.iloc[order].reset_index(drop=True)


def calculate_solar_position(time: pd.Series, lat: float, lon: float) -> np.ndarray:
    time = pd.to_datetime(time)
    doy = time.dt.dayofyear.to_numpy()
    hour = time.dt.hour.to_numpy() + time.dt.minute.to_numpy() / 60.0 + time.dt.second.to_numpy() / 3600.0

    gamma = 2.0 * np.pi / 365.0 * (doy - 1 + (hour - 12.0) / 24.0)
    eqtime = 229.18 * (
        0.000075
        + 0.001868 * np.cos(gamma)
        - 0.032077 * np.sin(gamma)
        - 0.014615 * np.cos(2.0 * gamma)
        - 0.040849 * np.sin(2.0 * gamma)
    )
    decl = (
        0.006918
        - 0.399912 * np.cos(gamma)
        + 0.070257 * np.sin(gamma)
        - 0.006758 * np.cos(2.0 * gamma)
        + 0.000907 * np.sin(2.0 * gamma)
        - 0.002697 * np.cos(3.0 * gamma)
        + 0.00148 * np.sin(3.0 * gamma)
    )

    tst = hour * 60.0 + (eqtime + 4.0 * lon)
    ha = np.deg2rad((tst / 4.0) - 180.0)
    lat_rad = np.deg2rad(lat)
    cos_sza = np.sin(lat_rad) * np.sin(decl) + np.cos(lat_rad) * np.cos(decl) * np.cos(ha)
    return np.clip(cos_sza, 0.0, None)


def add_clear_sky_longwave(station: pd.DataFrame) -> pd.DataFrame:
    ta_c = station["air_temp_c"].to_numpy(dtype=float)
    ta_k = ta_c + 273.15
    rh = station["rh_frac"].to_numpy(dtype=float)

    es_kpa = 0.6112 * np.exp((17.67 * ta_c) / (ta_c + 243.5))
    ea_hpa = rh * es_kpa * 10.0

    sigma = 5.670374419e-8
    emissivity = 1.24 * np.power(np.clip(ea_hpa / np.clip(ta_k, 1.0, None), 1e-9, None), 1.0 / 7.0)

    station = station.copy()
    station["lw_clear_sky"] = emissivity * sigma * np.power(ta_k, 4.0)
    station["obs_lw"] = pd.to_numeric(station["Pyrgeom_W"], errors="coerce")
    station["cos_sza"] = calculate_solar_position(station["time"], SITE_LAT, SITE_LON)
    return station


def apply_lw_flag(station: pd.DataFrame, a: float, b: float) -> pd.DataFrame:
    station = station.copy()
    denom = station["lw_clear_sky"] + a + b * station["obs_lw"]
    station["lw_csi"] = np.where(
        np.isfinite(station["obs_lw"]) & np.isfinite(denom) & (denom > 1e-6),
        station["obs_lw"] / denom,
        np.nan,
    )
    station["cloud_binary"] = np.where(np.isfinite(station["lw_csi"]), (station["lw_csi"] > 1.0).astype(int), np.nan)
    return station


def assign_temp_bins(df: pd.DataFrame, temp_col: str = "air_temp_c", width_c: float = TEMP_BIN_WIDTH_C) -> pd.DataFrame:
    out = df.copy()
    tmin = math.floor(min(out[temp_col].min(), 0.0) / width_c) * width_c
    tmax = math.ceil(max(out[temp_col].max(), 0.0) / width_c) * width_c
    edges = np.arange(tmin, tmax + width_c, width_c)
    out["temp_bin"] = pd.cut(out[temp_col], bins=edges, right=False, include_lowest=True)
    return out


def confusion_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, object]:
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    return {
        "cm": cm,
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "clear_f1": float(f1_score(y_true, y_pred, pos_label=0, zero_division=0)),
        "cloudy_f1": float(f1_score(y_true, y_pred, pos_label=1, zero_division=0)),
    }


def extract_tree_thresholds(clf: DecisionTreeClassifier) -> dict[str, float]:
    feature_names = ["red", "green", "blue"]
    thresholds: dict[str, list[float]] = {name: [] for name in feature_names}
    tree = clf.tree_
    for node_id in range(tree.node_count):
        fid = tree.feature[node_id]
        if fid >= 0:
            thresholds[feature_names[fid]].append(float(tree.threshold[node_id]))
    return {
        f"{name}_threshold": (float(np.median(vals)) if vals else np.nan)
        for name, vals in thresholds.items()
    }


def train_sb_per_bin(train_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    sb_pred = pd.Series(index=train_df.index, dtype="float64")
    records: list[dict[str, object]] = []

    for temp_bin, grp in train_df.groupby("temp_bin", observed=False):
        if pd.isna(temp_bin):
            continue

        counts = grp["cloud_binary"].value_counts()
        n_rows = len(grp)
        n_clear = int(counts.get(0, 0))
        n_cloudy = int(counts.get(1, 0))

        if n_rows < MIN_ROWS_PER_TEMP_BIN or min(n_clear, n_cloudy) < MIN_CLASS_COUNT_PER_TEMP_BIN:
            records.append(
                {
                    "temp_bin": str(temp_bin),
                    "temp_left_c": float(temp_bin.left),
                    "temp_right_c": float(temp_bin.right),
                    "n_rows": n_rows,
                    "n_clear": n_clear,
                    "n_cloudy": n_cloudy,
                    "status": "skipped_low_sample",
                    "red_threshold": np.nan,
                    "green_threshold": np.nan,
                    "blue_threshold": np.nan,
                }
            )
            continue

        clf = DecisionTreeClassifier(
            max_depth=3,
            min_samples_split=20,
            class_weight="balanced",
            random_state=RNG_SEED,
        )
        clf.fit(grp[["red", "green", "blue"]], grp["cloud_binary"].astype(int))
        sb_pred.loc[grp.index] = clf.predict(grp[["red", "green", "blue"]])

        records.append(
            {
                "temp_bin": str(temp_bin),
                "temp_left_c": float(temp_bin.left),
                "temp_right_c": float(temp_bin.right),
                "n_rows": n_rows,
                "n_clear": n_clear,
                "n_cloudy": n_cloudy,
                "status": "trained",
                **extract_tree_thresholds(clf),
            }
        )

    eval_df = train_df.loc[sb_pred.notna()].copy()
    eval_df["sb_dt_pred"] = sb_pred.loc[eval_df.index].astype(int)
    metrics = confusion_metrics(eval_df["cloud_binary"].to_numpy(), eval_df["sb_dt_pred"].to_numpy())
    thresholds = pd.DataFrame(records).sort_values("temp_left_c").reset_index(drop=True)
    return eval_df, thresholds, metrics


def cond_for_feature(x: float, threshold: float, direction: str) -> bool:
    if pd.isna(threshold) or pd.isna(x):
        return False
    direction = str(direction).strip()
    if direction == ">":
        return bool(x > threshold)
    if direction == "<":
        return bool(x < threshold)
    return False


def apply_gothic_thresholds(train_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    gothic = pd.read_csv(GOTHIC_THRESH_CSV)
    gothic = gothic[gothic["status"] == "trained"].copy()
    thr_map = {(float(r.temp_left_c), float(r.temp_right_c)): r for _, r in gothic.iterrows()}

    def predict_row(row: pd.Series) -> float:
        temp_bin = row["temp_bin"]
        if pd.isna(temp_bin):
            return np.nan
        key = (float(temp_bin.left), float(temp_bin.right))
        if key not in thr_map:
            return np.nan

        r = thr_map[key]
        checks = [
            cond_for_feature(row["red"], r.red_threshold, r.red_direction),
            cond_for_feature(row["green"], r.green_threshold, r.green_direction),
            cond_for_feature(row["blue"], r.blue_threshold, r.blue_direction),
        ]
        n_true = sum(int(c) for c in checks)
        rule = str(r.rule).strip().lower()
        cloudy = n_true >= 1 if rule == "union" else n_true >= 2
        return int(cloudy)

    pred = train_df.apply(predict_row, axis=1)
    eval_df = train_df.loc[pred.notna()].copy()
    eval_df["gothic_transfer_pred"] = pred.loc[eval_df.index].astype(int)
    metrics = confusion_metrics(eval_df["cloud_binary"].to_numpy(), eval_df["gothic_transfer_pred"].to_numpy())
    return eval_df, metrics


def plot_confusion(cm: np.ndarray, title: str, out_path: Path, cmap: str = "Blues") -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    im = ax.imshow(cm, cmap=cmap, vmin=0, vmax=max(1, int(cm.max())))
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Clear", "Cloudy"])
    ax.set_yticklabels(["Clear", "Cloudy"])
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Observed")
    ax.set_title(title)
    for r in range(2):
        for c in range(2):
            value = int(cm[r, c])
            ax.text(c, r, f"{value}", ha="center", va="center", color="white" if value > cm.max() / 2 else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    station = add_clear_sky_longwave(load_station())
    rgb = load_rgb()
    station_lw = apply_lw_flag(station, DOC_A, DOC_B)

    station_merge = station_lw[
        ["time", "air_temp_c", "cloud_binary", "lw_csi", "lw_clear_sky", "obs_lw", "cos_sza"]
    ].dropna(subset=["cloud_binary"])

    train_df = pd.merge(
        sort_by_time(rgb),
        sort_by_time(station_merge),
        on="time",
        how="inner",
    )
    train_df = train_df.dropna(subset=["red", "green", "blue", "cloud_binary", "air_temp_c"]).copy()
    train_df["cloud_binary"] = train_df["cloud_binary"].astype(int)
    train_df = assign_temp_bins(train_df)

    sb_eval_df, sb_thresholds, sb_metrics = train_sb_per_bin(train_df)
    gothic_eval_df, gothic_metrics = apply_gothic_thresholds(train_df)
    plot_confusion(
        sb_metrics["cm"],
        (
            "SB per-bin DT vs LW cloud flag\n"
            f"n={len(sb_eval_df)} | acc={sb_metrics['accuracy']:.3f} | macro_f1={sb_metrics['macro_f1']:.3f}"
        ),
        OUT_DIR / "sb_per_bin_dt_confusion_lw_all_timesteps.png",
        cmap="Blues",
    )
    plot_confusion(
        gothic_metrics["cm"],
        (
            "Gothic thresholds at SB vs LW cloud flag\n"
            f"n={len(gothic_eval_df)} | acc={gothic_metrics['accuracy']:.3f} | macro_f1={gothic_metrics['macro_f1']:.3f}"
        ),
        OUT_DIR / "sb_gothic_transfer_confusion_lw_all_timesteps.png",
        cmap="Oranges",
    )

    summary_out = pd.DataFrame(
        [
            {
                "run": "SB per-bin DT vs LW flag",
                "n_eval": len(sb_eval_df),
                "accuracy": sb_metrics["accuracy"],
                "macro_f1": sb_metrics["macro_f1"],
                "clear_f1": sb_metrics["clear_f1"],
                "cloudy_f1": sb_metrics["cloudy_f1"],
            },
            {
                "run": "Gothic thresholds at SB vs LW flag",
                "n_eval": len(gothic_eval_df),
                "accuracy": gothic_metrics["accuracy"],
                "macro_f1": gothic_metrics["macro_f1"],
                "clear_f1": gothic_metrics["clear_f1"],
                "cloudy_f1": gothic_metrics["cloudy_f1"],
            },
        ]
    )

    train_df.to_csv(OUT_DIR / "senator_beck_rgb_lw_training_table_pyrgeom.csv", index=False)
    sb_thresholds.to_csv(OUT_DIR / "sb_temp_bin_rgb_thresholds_from_lw_dt.csv", index=False)
    summary_out.to_csv(OUT_DIR / "sb_lw_vs_rgb_summary_metrics.csv", index=False)

    metadata_lines = [
        f"Station file: {SBSP_CSV}",
        f"RGB file: {RGB_CSV}",
        f"Gothic threshold file: {GOTHIC_THRESH_CSV}",
        "LW variable: Pyrgeom_W",
        "Air temperature source: (LoAir_Min_C + LoAir_Max_C)/2",
        "RH source: Lo_RH",
        "LW clear-sky baseline: Brutsaert emissivity parameterization",
        f"Doc constant a: {DOC_A}",
        f"Doc constant b: {DOC_B}",
        f"LW training rows after RGB merge: {len(train_df)}",
        f"LW class counts: {train_df['cloud_binary'].value_counts().sort_index().to_dict()}",
    ]
    (OUT_DIR / "sb_lw_metadata.txt").write_text("\n".join(metadata_lines) + "\n")

    print("Saved output directory:", OUT_DIR)
    print("LW constants:", {"a": DOC_A, "b": DOC_B})
    print("LW training rows:", len(train_df))
    print("LW class counts:", train_df["cloud_binary"].value_counts().sort_index().to_dict())
    print("SB per-bin DT metrics:", {k: v for k, v in sb_metrics.items() if k != "cm"})
    print("Gothic transfer metrics:", {k: v for k, v in gothic_metrics.items() if k != "cm"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
