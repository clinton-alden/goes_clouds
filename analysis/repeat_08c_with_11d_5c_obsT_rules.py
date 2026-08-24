#!/usr/bin/env python3
"""Repeat notebook 08c TSI/RGB/ACM comparisons using 11d 5 C obs-T rules."""

from __future__ import annotations

import glob
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path("/glade/u/home/cdalden/goes_work/analysis")
OUT_11D = ROOT / "output_11d_gothic"
OUT_DIR = ROOT / "output_08c_11d_5c_obsT_rules"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RGB_DIR = Path("/glade/derecho/scratch/cdalden/gothic/goes16/rgb_composite")
RULES_CSV = OUT_11D / "gothic_rgb_5cbin_decision_tree_rules_obsT_kt050_090.csv"
OBS_TEMP_CSV = OUT_11D / "gothic_observed_air_temp_5min.csv"
ACM_RGB_08C_CSV = ROOT / "output_08c/colorado_domain_cloud_fraction_14z_00z.csv"
TSI_PATH = Path("/glade/u/home/cdalden/scratch/surface_obs/colorado/sail_tsi_cloud_frac.csv")

START_DATE = "2021-10-01"
END_DATE = "2023-06-30"
LOCAL_TZ = "America/Denver"
LOCAL_START = "09:00"
LOCAL_END = "14:00"
TEMP_MATCH_TOLERANCE = pd.Timedelta("3min")
TSI_MATCH_TOLERANCE = pd.Timedelta("3min")
ACM_MATCH_TOLERANCE = pd.Timedelta("3min")

DOMAIN_LON_MIN = -107.04666
DOMAIN_LON_MAX = -106.93268
DOMAIN_LAT_MIN = 38.91294
DOMAIN_LAT_MAX = 38.95834
COND_RE = re.compile(r"(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)")


def corr_or_nan(df: pd.DataFrame, x: str, y: str) -> float:
    valid = df[[x, y]].dropna()
    if len(valid) < 2:
        return np.nan
    return float(valid[x].corr(valid[y]))


def mae(df: pd.DataFrame, x: str, y: str) -> float:
    valid = df[[x, y]].dropna()
    if valid.empty:
        return np.nan
    return float(np.abs(valid[x] - valid[y]).mean())


def rmse(df: pd.DataFrame, x: str, y: str) -> float:
    valid = df[[x, y]].dropna()
    if valid.empty:
        return np.nan
    return float(np.sqrt(((valid[x] - valid[y]) ** 2).mean()))


def sort_by_time(df: pd.DataFrame, time_col: str = "t") -> pd.DataFrame:
    # Work around a pandas datetime argsort issue seen in the current py3.14 env.
    order = np.argsort(df[time_col].to_numpy(dtype="datetime64[ns]"))
    return df.iloc[order].reset_index(drop=True)


def latlon_to_goes_scan(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    from pyproj import CRS, Transformer

    h = 35786023.0
    crs_geo = CRS.from_epsg(4326)
    crs_goes = CRS.from_proj4(
        "+proj=geos +h=35786023 +a=6378137 +b=6356752.31414 "
        "+lon_0=-75 +sweep=x +units=m +no_defs"
    )
    transformer = Transformer.from_crs(crs_geo, crs_goes, always_xy=True)
    x_m, y_m = transformer.transform(lon, lat)
    return np.asarray(x_m) / h, np.asarray(y_m) / h


def domain_scan_bounds() -> tuple[float, float, float, float]:
    lons = np.array([DOMAIN_LON_MIN, DOMAIN_LON_MIN, DOMAIN_LON_MAX, DOMAIN_LON_MAX])
    lats = np.array([DOMAIN_LAT_MIN, DOMAIN_LAT_MAX, DOMAIN_LAT_MIN, DOMAIN_LAT_MAX])
    xs, ys = latlon_to_goes_scan(lons, lats)
    return float(xs.min()), float(xs.max()), float(ys.min()), float(ys.max())


def load_observed_temp() -> pd.DataFrame:
    temp = pd.read_csv(OBS_TEMP_CSV, parse_dates=["time_local"])
    temp = temp.dropna(subset=["observed_air_temp_c"])
    temp = sort_by_time(temp, "time_local")
    return temp[["time_local", "observed_air_temp_c"]]


def assign_temp_bin(temp_c: float) -> str | float:
    if pd.isna(temp_c):
        return np.nan
    edges = [-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0]
    clipped = min(max(float(temp_c), edges[0]), np.nextafter(edges[-1], -np.inf))
    for left, right in zip(edges[:-1], edges[1:]):
        if left <= clipped < right:
            return f"[{int(left)}, {int(right)})"
    return np.nan


def rule_mask(values: dict[str, np.ndarray], rule: str) -> np.ndarray:
    out = np.ones_like(values["red"], dtype=bool)
    for feature, op, threshold in COND_RE.findall(rule):
        value = values[feature]
        threshold_f = float(threshold)
        if op == "<=":
            out &= value <= threshold_f
        elif op == ">":
            out &= value > threshold_f
    return out


def list_rgb_files() -> list[Path]:
    start_key = START_DATE.replace("-", "")
    end_key = END_DATE.replace("-", "")
    files = []
    for path in sorted(glob.glob(str(RGB_DIR / "goes16_C02_C05_C13_rgb_gothic_*.nc"))):
        date_key = Path(path).stem.rsplit("_", 1)[-1]
        if start_key <= date_key <= end_key:
            files.append(Path(path))
    if not files:
        raise FileNotFoundError(f"No RGB files found in {RGB_DIR}")
    return files


def build_rgb_5c_pixel_cloud_fraction(rules: pd.DataFrame) -> pd.DataFrame:
    out_csv = OUT_DIR / "gothic_rgb_11d_5c_obsT_pixel_cloud_fraction.csv"
    if out_csv.exists():
        df = pd.read_csv(out_csv, parse_dates=["t", "time_local"])
        print(f"Reusing 5 C obs-T RGB cloud fraction: {out_csv} ({len(df)} rows)", flush=True)
        return df

    x_min, x_max, y_min, y_max = domain_scan_bounds()
    files = list_rgb_files()
    temp = load_observed_temp()
    cloudy_rules_by_bin = {
        temp_bin: grp.loc[grp["prediction"].eq(1), "rule"].astype(str).tolist()
        for temp_bin, grp in rules.groupby("temp_bin")
    }
    frames = []

    for i, path in enumerate(files, start=1):
        if i == 1 or i % 25 == 0 or i == len(files):
            print(f"Applying 5 C obs-T rules to RGB file {i}/{len(files)}: {path.name}", flush=True)

        with xr.open_dataset(path) as ds:
            x_slice = slice(x_min, x_max) if ds["x"][0] < ds["x"][-1] else slice(x_max, x_min)
            y_slice = slice(y_min, y_max) if ds["y"][0] < ds["y"][-1] else slice(y_max, y_min)
            sub = ds.sel(x=x_slice, y=y_slice)
            if sub.sizes.get("x", 0) == 0 or sub.sizes.get("y", 0) == 0:
                raise ValueError(f"Domain selected no pixels in {path}")

            times_utc = pd.DatetimeIndex(pd.to_datetime(sub["t"].values))
            times_local = times_utc.tz_localize("UTC").tz_convert(LOCAL_TZ).tz_localize(None)
            keep = pd.Series(True, index=times_local).between_time(LOCAL_START, LOCAL_END)
            keep = keep.reindex(times_local, fill_value=False).to_numpy()
            if not keep.any():
                continue

            sub = sub.isel(t=np.where(keep)[0])
            times_utc = times_utc[keep]
            times_local = pd.DatetimeIndex(times_local[keep])
            match = pd.merge_asof(
                sort_by_time(pd.DataFrame({"time_local": times_local}), "time_local"),
                temp,
                on="time_local",
                direction="nearest",
                tolerance=TEMP_MATCH_TOLERANCE,
            )
            valid_temp = match["observed_air_temp_c"].notna().to_numpy()
            if not valid_temp.any():
                continue

            sub = sub.isel(t=np.where(valid_temp)[0])
            times_utc = times_utc[valid_temp]
            times_local = times_local[valid_temp]
            obs_temp = match.loc[valid_temp, "observed_air_temp_c"].to_numpy(float)
            temp_bins = pd.Series(obs_temp).map(assign_temp_bin).to_numpy()

            red = sub["red"].values
            green = sub["green"].values
            blue = sub["blue"].values
            cloud_frac = np.full(red.shape[0], np.nan, dtype=float)

            for temp_bin in pd.unique(temp_bins):
                if pd.isna(temp_bin):
                    continue
                idx = np.where(temp_bins == temp_bin)[0]
                rules_for_bin = cloudy_rules_by_bin.get(temp_bin, [])
                if not rules_for_bin:
                    continue
                values = {"red": red[idx], "green": green[idx], "blue": blue[idx]}
                cloud = np.zeros_like(values["red"], dtype=bool)
                for rule in rules_for_bin:
                    cloud |= rule_mask(values, rule)
                cloud_frac[idx] = cloud.mean(axis=(1, 2))

            frames.append(
                pd.DataFrame(
                    {
                        "t": times_utc.astype("datetime64[ns]"),
                        "time_local": times_local.astype("datetime64[ns]"),
                        "observed_air_temp_c": obs_temp,
                        "temp_bin": temp_bins,
                        "rgb_11d_5c_cloud_frac": cloud_frac,
                        "rgb_11d_5c_n_pixels": sub.sizes["x"] * sub.sizes["y"],
                        "rgb_path": str(path),
                    }
                )
            )

    if not frames:
        raise RuntimeError("No RGB 5 C cloud-fraction rows produced")
    df = pd.concat(frames, ignore_index=True).dropna(subset=["rgb_11d_5c_cloud_frac"])
    df = sort_by_time(df).drop_duplicates("t").reset_index(drop=True)
    df.to_csv(out_csv, index=False)
    print(f"Saved 5 C obs-T RGB cloud fraction: {out_csv} ({len(df)} rows)", flush=True)
    return df


def season_label(month: int) -> str:
    if month in {11, 12, 1, 2, 3}:
        return "Nov-Mar"
    if month in {4, 5, 6}:
        return "Apr-Jun"
    if month == 10:
        return "Oct"
    return "Other"


def build_compare_df() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rules = pd.read_csv(RULES_CSV)
    rgb_pred = build_rgb_5c_pixel_cloud_fraction(rules)

    acm_df = pd.read_csv(ACM_RGB_08C_CSV, parse_dates=["date", "t"])
    acm_df["t"] = pd.to_datetime(acm_df["t"]).astype("datetime64[ns]")
    acm_df = acm_df[["date", "date_key", "t", "acm_cloud_frac", "acm_n_pixels", "acm_path"]]
    acm_df = sort_by_time(acm_df)

    tsi_df = pd.read_csv(TSI_PATH, parse_dates=["t"])
    tsi_df["t"] = pd.to_datetime(tsi_df["t"]).astype("datetime64[ns]")
    tsi_df.loc[tsi_df["tsi_frac"] < 0, "tsi_frac"] = np.nan
    tsi_df = sort_by_time(tsi_df)

    compare = pd.merge_asof(
        sort_by_time(rgb_pred),
        acm_df,
        on="t",
        direction="nearest",
        tolerance=ACM_MATCH_TOLERANCE,
    )
    compare = pd.merge_asof(
        sort_by_time(compare),
        tsi_df[["t", "tsi_frac"]],
        on="t",
        direction="nearest",
        tolerance=TSI_MATCH_TOLERANCE,
    )
    compare["date"] = compare["t"].dt.normalize()
    compare["season"] = compare["t"].dt.month.map(season_label)
    compare["acm_residual"] = compare["acm_cloud_frac"] - compare["tsi_frac"]
    compare["rgb_11d_5c_residual"] = compare["rgb_11d_5c_cloud_frac"] - compare["tsi_frac"]
    compare.to_csv(OUT_DIR / "gothic_08c_compare_with_11d_5c_obsT_rules.csv", index=False)

    season_rows = []
    for season, df in compare.groupby("season", sort=False):
        valid = df.loc[df["tsi_frac"].notna()]
        season_rows.append(
            {
                "season": season,
                "rows": len(df),
                "valid_tsi_rows": len(valid),
                "mean_tsi_cloud_fraction": valid["tsi_frac"].mean(),
                "mean_rgb_11d_5c_cloud_fraction": valid["rgb_11d_5c_cloud_frac"].mean(),
                "mean_acm_cloud_fraction": valid["acm_cloud_frac"].mean(),
                "rgb_11d_5c_mean_residual": valid["rgb_11d_5c_residual"].mean(),
                "acm_mean_residual": valid["acm_residual"].mean(),
                "rgb_11d_5c_mae": mae(valid, "rgb_11d_5c_cloud_frac", "tsi_frac"),
                "acm_mae": mae(valid, "acm_cloud_frac", "tsi_frac"),
                "rgb_11d_5c_rmse": rmse(valid, "rgb_11d_5c_cloud_frac", "tsi_frac"),
                "acm_rmse": rmse(valid, "acm_cloud_frac", "tsi_frac"),
                "corr_rgb_11d_5c_tsi": corr_or_nan(valid, "rgb_11d_5c_cloud_frac", "tsi_frac"),
                "corr_acm_tsi": corr_or_nan(valid, "acm_cloud_frac", "tsi_frac"),
            }
        )
    seasonal_summary = pd.DataFrame(season_rows)
    seasonal_summary.to_csv(OUT_DIR / "gothic_08c_11d_5c_obsT_seasonal_summary.csv", index=False)
    return compare, seasonal_summary, rgb_pred


def save_figures(compare_df: pd.DataFrame) -> None:
    hist_df = compare_df.loc[compare_df["tsi_frac"].notna()].copy()
    bins = np.linspace(-1, 1, 41)
    period_defs = [("Nov-Mar", [11, 12, 1, 2, 3]), ("Apr-Jun", [4, 5, 6])]

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True, sharey=True)
    for ax, (period_label, months) in zip(axes, period_defs):
        period_df = hist_df.loc[hist_df["t"].dt.month.isin(months)]
        rgb_tsi = period_df["rgb_11d_5c_residual"].dropna()
        acm_tsi = period_df["acm_residual"].dropna()
        ax.hist(
            rgb_tsi,
            bins=bins,
            alpha=0.6,
            density=True,
            label=f"RGB 5 C 11d - TSI (mean={rgb_tsi.mean():.2f})",
        )
        ax.hist(
            acm_tsi,
            bins=bins,
            alpha=0.6,
            density=True,
            label=f"ACM - TSI (mean={acm_tsi.mean():.2f})",
        )
        ax.axvline(0, color="k", ls="--", lw=1)
        ax.set_title(f"{period_label}: residuals (n={len(period_df)})", fontweight="bold")
        ax.set_ylabel("Density", fontweight="bold")
        ax.legend()
    axes[-1].set_xlabel("Residual cloud fraction", fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gothic_08c_11d_5c_obsT_seasonal_residual_hist.png", dpi=200)
    plt.close(fig)

    valid = hist_df.loc[
        hist_df[["tsi_frac", "rgb_11d_5c_cloud_frac", "acm_cloud_frac"]].notna().all(axis=1)
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    axes[0].scatter(valid["tsi_frac"], valid["rgb_11d_5c_cloud_frac"], s=12, alpha=0.5, label="RGB 5 C 11d vs TSI")
    axes[0].scatter(valid["tsi_frac"], valid["acm_cloud_frac"], s=12, alpha=0.5, label="ACM vs TSI")
    axes[0].plot([0, 1], [0, 1], "k--", lw=1)
    axes[0].set_title("GOES cloud fraction vs TSI using 5 C 11d RGB rules")
    axes[0].set_xlabel("TSI cloud fraction")
    axes[0].set_ylabel("GOES cloud fraction")
    axes[0].set_xlim(-0.05, 1.05)
    axes[0].set_ylim(-0.05, 1.05)
    axes[0].legend()

    axes[1].hist(valid["rgb_11d_5c_residual"], bins=bins, alpha=0.6, label="RGB 5 C 11d - TSI")
    axes[1].hist(valid["acm_residual"], bins=bins, alpha=0.6, label="ACM - TSI")
    axes[1].axvline(0, color="k", ls="--", lw=1)
    axes[1].set_title("Residual distributions")
    axes[1].set_xlabel("Residual cloud fraction")
    axes[1].set_ylabel("Count")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gothic_08c_11d_5c_obsT_scatter_residual_hist.png", dpi=200)
    plt.close(fig)


def main() -> int:
    compare, seasonal_summary, rgb_pred = build_compare_df()
    save_figures(compare)
    print(f"RGB 5 C obs-T predictions: {len(rgb_pred)}")
    print(f"Compare rows: {len(compare)}")
    print(f"Valid TSI rows: {compare['tsi_frac'].notna().sum()}")
    print(seasonal_summary.to_string(index=False))
    print(f"Outputs: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
