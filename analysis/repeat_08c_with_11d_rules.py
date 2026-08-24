#!/usr/bin/env python3
"""Repeat notebook 08c comparisons using 11d RGB decision-tree rules."""

from __future__ import annotations

import re
import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path("/glade/u/home/cdalden/goes_work/analysis")
OUT_11D = ROOT / "output_11d_gothic"
OUT_DIR = ROOT / "output_08c_11d_rules"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RGB_CSV = OUT_11D / "gothic_rgb_SW_domain.csv"
RGB_DIR = Path("/glade/derecho/scratch/cdalden/gothic/goes16/rgb_composite")
RULES_CSV = OUT_11D / "gothic_rgb_tempbin_decision_tree_rules.csv"
ACM_RGB_08C_CSV = ROOT / "output_08c/colorado_domain_cloud_fraction_14z_00z.csv"
TSI_PATH = Path("/glade/u/home/cdalden/scratch/surface_obs/colorado/sail_tsi_cloud_frac.csv")
ERA5_DIR = Path("/glade/derecho/scratch/cdalden/gothic/era5/t2m_hourly")

START_DATE = "2021-10-01"
END_DATE = "2023-06-30"
LOCAL_TZ = "America/Denver"
LOCAL_START = "09:00"
LOCAL_END = "14:00"
TSI_MATCH_TOLERANCE = pd.Timedelta("3min")
ACM_MATCH_TOLERANCE = pd.Timedelta("3min")

ERA5_PIXEL_LAT = 38.95781
ERA5_PIXEL_LON = -106.98795
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


def resolve_time_coord_name(da: xr.DataArray) -> str:
    if "time" in da.coords:
        return "time"
    if "valid_time" in da.coords:
        return "valid_time"
    for dim in da.dims:
        if "time" in dim:
            return dim
    raise KeyError("Could not find ERA5 time coordinate")


def load_era5_temp_series() -> pd.Series:
    months = pd.period_range(START_DATE, END_DATE, freq="M")
    frames = []
    for month in months:
        path = ERA5_DIR / f"era5_t2m_gothic_{month.year}{month.month:02d}.nc"
        if not path.exists():
            raise FileNotFoundError(path)
        with xr.open_dataset(path) as ds:
            da = (
                ds["t2m"]
                .sel(latitude=ERA5_PIXEL_LAT, longitude=ERA5_PIXEL_LON, method="nearest")
                - 273.15
            )
            time_name = resolve_time_coord_name(da)
            frames.append(
                pd.DataFrame(
                    {
                        "time_utc": pd.to_datetime(da[time_name].values),
                        "era5_land_t2m_c": da.values.astype(float),
                    }
                )
            )

    era5 = pd.concat(frames, ignore_index=True).drop_duplicates("time_utc").sort_values("time_utc")
    era5["time_local"] = (
        era5["time_utc"].dt.tz_localize("UTC").dt.tz_convert(LOCAL_TZ).dt.tz_localize(None)
    )
    temp_series = (
        era5.set_index("time_local")["era5_land_t2m_c"].sort_index().groupby(level=0).mean()
    )
    return temp_series


def interpolate_temp_at_times(temp_series: pd.Series, times_local: pd.Series) -> pd.DataFrame:
    times = pd.DatetimeIndex(pd.to_datetime(times_local)).sort_values().unique()
    temp_at_times = (
        temp_series.reindex(temp_series.index.union(times))
        .sort_index()
        .interpolate(method="time")
        .reindex(times)
        .ffill()
        .bfill()
    )
    return pd.DataFrame({"time_local": times, "era5_land_t2m_c": temp_at_times.to_numpy(float)})


def assign_temp_bin(temp_c: float) -> str | float:
    if pd.isna(temp_c):
        return np.nan
    bins = [(-20.0, -10.0), (-10.0, 0.0), (0.0, 10.0), (10.0, 20.0)]
    clipped = min(max(float(temp_c), bins[0][0]), np.nextafter(bins[-1][1], -np.inf))
    for left, right in bins:
        if left <= clipped < right:
            return f"[{int(left)}, {int(right)})"
    return np.nan


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


def build_rgb_11d_pixel_cloud_fraction(rules: pd.DataFrame) -> pd.DataFrame:
    out_csv = OUT_DIR / "gothic_rgb_11d_pixel_cloud_fraction.csv"
    if out_csv.exists():
        df = pd.read_csv(out_csv, parse_dates=["t", "time_local"])
        print(f"Reusing pixel-level 11d RGB cloud fraction: {out_csv} ({len(df)} rows)", flush=True)
        return df

    x_min, x_max, y_min, y_max = domain_scan_bounds()
    files = list_rgb_files()
    temp_series = load_era5_temp_series()
    cloudy_rules_by_bin = {
        temp_bin: grp.loc[grp["prediction"].eq(1), "rule"].astype(str).tolist()
        for temp_bin, grp in rules.groupby("temp_bin")
    }
    frames = []

    for i, path in enumerate(files, start=1):
        if i == 1 or i % 25 == 0 or i == len(files):
            print(f"Applying 11d pixel rules to RGB file {i}/{len(files)}: {path.name}", flush=True)

        with xr.open_dataset(path) as ds:
            x_slice = slice(x_min, x_max) if ds["x"][0] < ds["x"][-1] else slice(x_max, x_min)
            y_slice = slice(y_min, y_max) if ds["y"][0] < ds["y"][-1] else slice(y_max, y_min)
            sub = ds.sel(x=x_slice, y=y_slice)
            if sub.sizes.get("x", 0) == 0 or sub.sizes.get("y", 0) == 0:
                raise ValueError(f"Domain selected no pixels in {path}")

            times_utc = pd.DatetimeIndex(pd.to_datetime(sub["t"].values))
            times_local = times_utc.tz_localize("UTC").tz_convert(LOCAL_TZ).tz_localize(None)
            keep = pd.Series(True, index=times_local).between_time(
                LOCAL_START, LOCAL_END
            ).reindex(times_local, fill_value=False).to_numpy()
            if not keep.any():
                continue

            sub = sub.isel(t=np.where(keep)[0])
            times_utc = times_utc[keep]
            times_local = times_local[keep]
            temp_df = interpolate_temp_at_times(temp_series, pd.Series(times_local))
            temp_bins = temp_df["era5_land_t2m_c"].map(assign_temp_bin).to_numpy()

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
                        "time_local": pd.to_datetime(times_local).astype("datetime64[ns]"),
                        "era5_land_t2m_c": temp_df["era5_land_t2m_c"].to_numpy(float),
                        "temp_bin": temp_bins,
                        "rgb_11d_cloud_frac": cloud_frac,
                        "rgb_11d_n_pixels": sub.sizes["x"] * sub.sizes["y"],
                        "rgb_path": str(path),
                    }
                )
            )

    if not frames:
        raise RuntimeError("No RGB 11d cloud-fraction rows produced")
    df = pd.concat(frames, ignore_index=True).dropna(subset=["rgb_11d_cloud_frac"])
    df = df.sort_values("t").drop_duplicates("t").reset_index(drop=True)
    df.to_csv(out_csv, index=False)
    print(f"Saved pixel-level 11d RGB cloud fraction: {out_csv} ({len(df)} rows)", flush=True)
    return df


def build_compare_df() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rules = pd.read_csv(RULES_CSV)
    rgb_pred = build_rgb_11d_pixel_cloud_fraction(rules)

    acm_df = pd.read_csv(ACM_RGB_08C_CSV, parse_dates=["date", "t"])
    acm_df["t"] = pd.to_datetime(acm_df["t"]).astype("datetime64[ns]")
    acm_df = acm_df[["date", "date_key", "t", "acm_cloud_frac", "acm_n_pixels", "acm_path"]]

    tsi_df = pd.read_csv(TSI_PATH, parse_dates=["t"])
    tsi_df["t"] = pd.to_datetime(tsi_df["t"]).astype("datetime64[ns]")
    tsi_df.loc[tsi_df["tsi_frac"] < 0, "tsi_frac"] = np.nan

    compare = pd.merge_asof(
        rgb_pred.sort_values("t"),
        acm_df.sort_values("t"),
        on="t",
        direction="nearest",
        tolerance=ACM_MATCH_TOLERANCE,
    )
    compare = pd.merge_asof(
        compare.sort_values("t"),
        tsi_df[["t", "tsi_frac"]].sort_values("t"),
        on="t",
        direction="nearest",
        tolerance=TSI_MATCH_TOLERANCE,
    )
    compare["date"] = compare["t"].dt.normalize()
    compare["acm_residual"] = compare["acm_cloud_frac"] - compare["tsi_frac"]
    compare["rgb_11d_residual"] = compare["rgb_11d_cloud_frac"] - compare["tsi_frac"]
    compare.to_csv(OUT_DIR / "gothic_08c_compare_with_11d_rules.csv", index=False)

    summary = pd.DataFrame(
        {
            "metric": [
                "Rows",
                "Rows with valid TSI",
                "Mean ACM cloud fraction",
                "Mean RGB 11d cloud value",
                "Mean TSI cloud fraction",
                "Corr(ACM, TSI)",
                "Corr(RGB 11d, TSI)",
                "ACM mean residual",
                "RGB 11d mean residual",
                "ACM MAE",
                "RGB 11d MAE",
                "ACM RMSE",
                "RGB 11d RMSE",
            ],
            "value": [
                len(compare),
                compare["tsi_frac"].notna().sum(),
                compare["acm_cloud_frac"].mean(),
                compare["rgb_11d_cloud_frac"].mean(),
                compare["tsi_frac"].mean(),
                corr_or_nan(compare, "acm_cloud_frac", "tsi_frac"),
                corr_or_nan(compare, "rgb_11d_cloud_frac", "tsi_frac"),
                compare["acm_residual"].mean(),
                compare["rgb_11d_residual"].mean(),
                mae(compare, "acm_cloud_frac", "tsi_frac"),
                mae(compare, "rgb_11d_cloud_frac", "tsi_frac"),
                rmse(compare, "acm_cloud_frac", "tsi_frac"),
                rmse(compare, "rgb_11d_cloud_frac", "tsi_frac"),
            ],
        }
    )
    summary.to_csv(OUT_DIR / "gothic_08c_11d_summary.csv", index=False)
    return compare, summary, rgb_pred


def save_figures(compare_df: pd.DataFrame) -> None:
    valid_rgb = compare_df[["tsi_frac", "rgb_11d_cloud_frac"]].dropna()
    valid_acm = compare_df[["tsi_frac", "acm_cloud_frac"]].dropna()

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    axes[0].scatter(valid_rgb["tsi_frac"], valid_rgb["rgb_11d_cloud_frac"], s=12, alpha=0.5, label="RGB 11d vs TSI")
    axes[0].scatter(valid_acm["tsi_frac"], valid_acm["acm_cloud_frac"], s=12, alpha=0.5, label="ACM vs TSI")
    axes[0].plot([0, 1], [0, 1], "k--", lw=1)
    axes[0].set_title("GOES cloud value vs TSI using 11d RGB rules")
    axes[0].set_xlabel("TSI cloud fraction")
    axes[0].set_ylabel("GOES cloud value/fraction")
    axes[0].set_xlim(-0.05, 1.05)
    axes[0].set_ylim(-0.05, 1.05)
    axes[0].legend()

    axes[1].plot(compare_df["t"], compare_df["acm_cloud_frac"], label="ACM", lw=1.2, alpha=0.9)
    axes[1].plot(compare_df["t"], compare_df["rgb_11d_cloud_frac"], label="RGB 11d", lw=1.2, alpha=0.9)
    axes[1].plot(compare_df["t"], compare_df["tsi_frac"], label="TSI", lw=1.2, alpha=0.7)
    axes[1].set_title("Time series")
    axes[1].set_xlabel("Time")
    axes[1].set_ylabel("Cloud value/fraction")
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].legend()
    plt.tight_layout()
    fig.savefig(OUT_DIR / "gothic_08c_11d_scatter_timeseries.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    resid_rgb = compare_df[["tsi_frac", "rgb_11d_residual"]].dropna()
    resid_acm = compare_df[["tsi_frac", "acm_residual"]].dropna()
    axes[0].scatter(resid_rgb["tsi_frac"], resid_rgb["rgb_11d_residual"], s=12, alpha=0.5, label="RGB 11d residual")
    axes[0].scatter(resid_acm["tsi_frac"], resid_acm["acm_residual"], s=12, alpha=0.5, label="ACM residual")
    axes[0].axhline(0, color="k", ls="--", lw=1)
    axes[0].set_title("Residual plot: GOES minus TSI")
    axes[0].set_xlabel("TSI cloud fraction")
    axes[0].set_ylabel("Residual")
    axes[0].legend()
    axes[1].hist(resid_rgb["rgb_11d_residual"], bins=30, alpha=0.6, label="RGB 11d residual")
    axes[1].hist(resid_acm["acm_residual"], bins=30, alpha=0.6, label="ACM residual")
    axes[1].axvline(0, color="k", ls="--", lw=1)
    axes[1].set_title("Residual distributions")
    axes[1].set_xlabel("Residual")
    axes[1].set_ylabel("Count")
    axes[1].legend()
    plt.tight_layout()
    fig.savefig(OUT_DIR / "gothic_08c_11d_residual_scatter_hist.png", dpi=200)
    plt.close(fig)

    period_defs = [("Nov-Mar", [11, 12, 1, 2, 3]), ("Apr-Jun", [4, 5, 6])]
    hist_df = compare_df.loc[compare_df["tsi_frac"].notna()].copy()
    bins = np.linspace(-1, 1, 41)
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, sharey=True)
    for ax, (period_label, months) in zip(axes, period_defs):
        period_df = hist_df.loc[hist_df["t"].dt.month.isin(months)]
        rgb_tsi = period_df["rgb_11d_residual"].dropna()
        acm_tsi = period_df["acm_residual"].dropna()
        ax.hist(rgb_tsi, bins=bins, alpha=0.6, density=True, label=f"RGB 11d - TSI (mean={rgb_tsi.mean():.2f})")
        ax.hist(acm_tsi, bins=bins, alpha=0.6, density=True, label=f"ACM - TSI (mean={acm_tsi.mean():.2f})")
        ax.axvline(0, color="k", ls="--", lw=1)
        ax.set_title(f"{period_label}: residuals (n={len(period_df)})", fontweight="bold")
        ax.set_xlabel("Residual cloud fraction", fontweight="bold")
        ax.set_ylabel("Density", fontweight="bold")
        ax.legend()
    plt.tight_layout()
    fig.savefig(OUT_DIR / "gothic_08c_11d_seasonal_residual_hist.png", dpi=200)
    plt.close(fig)

    plot_df = compare_df.loc[
        compare_df[["tsi_frac", "rgb_11d_cloud_frac", "acm_cloud_frac"]].notna().all(axis=1)
    ].copy()
    plot_df["rgb_11d_abs_error"] = np.abs(plot_df["rgb_11d_cloud_frac"] - plot_df["tsi_frac"])
    plot_df["acm_abs_error"] = np.abs(plot_df["acm_cloud_frac"] - plot_df["tsi_frac"])
    bin_edges = np.linspace(0, 1, 11)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    plot_df["tsi_bin"] = pd.cut(plot_df["tsi_frac"], bins=bin_edges, include_lowest=True)

    def binned_summary(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
        grouped = df.groupby("tsi_bin", observed=False)[value_col]
        return pd.DataFrame(
            {
                "median": grouped.median(),
                "q25": grouped.quantile(0.25),
                "q75": grouped.quantile(0.75),
                "count": grouped.count(),
            }
        )

    rgb_summary = binned_summary(plot_df, "rgb_11d_cloud_frac")
    acm_summary = binned_summary(plot_df, "acm_cloud_frac")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    ax = axes[0]
    ax.plot([0, 1], [0, 1], color="k", ls="--", lw=1.25, label="1:1 line")
    ax.plot(bin_centers, rgb_summary["median"], color="tab:orange", lw=2, marker="o", label="RGB 11d median")
    ax.fill_between(bin_centers, rgb_summary["q25"], rgb_summary["q75"], color="tab:orange", alpha=0.25, label="RGB 11d IQR")
    ax.plot(bin_centers, acm_summary["median"], color="tab:blue", lw=2, marker="o", label="ACM median")
    ax.fill_between(bin_centers, acm_summary["q25"], acm_summary["q75"], color="tab:blue", alpha=0.25, label="ACM IQR")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("TSI cloud fraction", fontweight="bold")
    ax.set_ylabel("GOES cloud value/fraction", fontweight="bold")
    ax.set_title(f"Binned median vs TSI (n={len(plot_df)})", fontweight="bold")
    ax.legend(frameon=True)

    ax = axes[1]
    for err_col, label, color in [
        ("rgb_11d_abs_error", "RGB 11d absolute error", "tab:orange"),
        ("acm_abs_error", "ACM absolute error", "tab:blue"),
    ]:
        vals = np.sort(plot_df[err_col].dropna().to_numpy())
        ecdf = np.arange(1, len(vals) + 1) / len(vals)
        ax.plot(vals, ecdf, lw=2, color=color, label=label)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Absolute error vs TSI", fontweight="bold")
    ax.set_ylabel("Empirical CDF", fontweight="bold")
    ax.set_title("Agreement with TSI", fontweight="bold")
    ax.legend(frameon=True, loc="lower right")
    plt.tight_layout()
    fig.savefig(OUT_DIR / "gothic_08c_11d_binned_median_ecdf.png", dpi=200)
    plt.close(fig)


def main() -> int:
    compare, summary, rgb_pred = build_compare_df()
    save_figures(compare)
    print(f"RGB 11d predictions: {len(rgb_pred)}")
    print(f"Compare rows: {len(compare)}")
    print(f"Valid TSI rows: {compare['tsi_frac'].notna().sum()}")
    print(summary.to_string(index=False))
    print(f"Outputs: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
