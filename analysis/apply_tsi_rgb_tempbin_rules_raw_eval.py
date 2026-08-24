#!/usr/bin/env python
"""Apply even-day-trained RGB temperature-bin rules to raw odd-day RGB pixels."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from sklearn.metrics import classification_report


RGB_DIR = Path("/glade/u/home/cdalden/scratch/colorado/goes16/rgb_composite")
ACM_DIR = Path("/glade/u/home/cdalden/scratch/colorado_acm/goes16/daily_nc_latlon")
TSI_CSV = Path("/glade/u/home/cdalden/scratch/surface_obs/colorado/sail_tsi_cloud_frac.csv")
ERA5_PIXEL_5MIN_CSV = Path("analysis/era5_land_t2m_pixel_5min.csv")
RULES_CSV = Path("analysis/tsi_rgb_decision_tree_era5_t2m_10c_rules.csv")

OUT_MASK_NC = Path("analysis/tsi_rgb_tempbin_raw_odd_day_cloud_mask_10km.nc")
OUT_FRACTION_CSV = Path("analysis/tsi_rgb_tempbin_raw_odd_day_cloud_fraction_vs_tsi_acm.csv")
OUT_METRICS_CSV = Path("analysis/tsi_rgb_tempbin_raw_odd_day_cloud_fraction_metrics.csv")
OUT_MONTHLY_METRICS_CSV = Path(
    "analysis/tsi_rgb_tempbin_raw_odd_day_cloud_fraction_monthly_metrics.csv"
)
OUT_SEASONAL_METRICS_CSV = Path(
    "analysis/tsi_rgb_tempbin_raw_odd_day_cloud_fraction_seasonal_metrics.csv"
)

LAT_MIN = 38.91617092558588
LAT_MAX = 39.00624907441412
LON_MIN = -107.05032765812774
LON_MAX = -106.93495234187225

TSI_CLEAR_THRESHOLD = 0.1
TSI_CLOUDY_THRESHOLD = 0.9
TSI_MATCH_TOLERANCE = pd.Timedelta("15s")
ACM_MATCH_TOLERANCE = pd.Timedelta("3min")
COND_RE = re.compile(r"(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)")


def oriented_slice(coord: xr.DataArray, lower: float, upper: float) -> slice:
    first = float(coord.values[0])
    last = float(coord.values[-1])
    return slice(lower, upper) if first <= last else slice(upper, lower)


def parse_filename_date(path: Path) -> pd.Timestamp:
    token = path.stem.rsplit("_", 1)[-1]
    return pd.to_datetime(token, format="%Y%m%d")


def corrected_times(path: Path, ds: xr.Dataset) -> pd.DatetimeIndex:
    file_date = parse_filename_date(path)
    internal_times = pd.DatetimeIndex(pd.to_datetime(ds["t"].values))
    time_of_day = internal_times - internal_times.normalize()
    return pd.DatetimeIndex(file_date + time_of_day)


def load_cloudy_rules(path: Path) -> dict[str, list[str]]:
    rules = pd.read_csv(path)
    rules = rules.loc[rules["predicted_label"].eq("cloud")].copy()
    if rules.empty:
        raise ValueError(f"No cloudy rules found in {path}")
    return {
        str(temp_bin): group["rule"].astype(str).tolist()
        for temp_bin, group in rules.groupby("temp_bin")
    }


def temp_bin_midpoints(labels: list[str]) -> dict[str, float]:
    out = {}
    for label in labels:
        left, right = re.findall(r"-?\d+", label)
        out[label] = (float(left) + float(right)) / 2.0
    return out


def nearest_trained_bin(temp_c: float, midpoints: dict[str, float]) -> str:
    return min(midpoints, key=lambda label: abs(float(temp_c) - midpoints[label]))


def ten_c_bin_label(temp_c: float) -> str:
    left = np.floor(float(temp_c) / 10.0) * 10.0
    right = left + 10.0
    return f"[{int(left)}, {int(right)})"


def apply_rule(values: dict[str, np.ndarray], rule: str) -> np.ndarray:
    mask = np.ones_like(values["red"], dtype=bool)
    for feature, op, threshold_text in COND_RE.findall(rule):
        threshold = float(threshold_text)
        if op == "<=":
            mask &= values[feature] <= threshold
        elif op == ">":
            mask &= values[feature] > threshold
        else:
            raise ValueError(f"Unsupported operator in rule: {op}")
    return mask


def classify_pixels(
    red: np.ndarray,
    green: np.ndarray,
    blue: np.ndarray,
    bins: np.ndarray,
    cloudy_rules_by_bin: dict[str, list[str]],
) -> np.ndarray:
    cloud = np.zeros(red.shape, dtype=bool)
    for temp_bin in pd.unique(pd.Series(bins)):
        idx = np.where(bins == temp_bin)[0]
        rules = cloudy_rules_by_bin.get(str(temp_bin), [])
        if len(idx) == 0 or not rules:
            continue
        values = {"red": red[idx], "green": green[idx], "blue": blue[idx]}
        bin_cloud = np.zeros_like(values["red"], dtype=bool)
        for rule in rules:
            bin_cloud |= apply_rule(values, rule)
        cloud[idx] = bin_cloud
    invalid = np.isnan(red) | np.isnan(green) | np.isnan(blue)
    return np.where(invalid, 255, cloud.astype(np.uint8))


def acm_cloud_fraction_for_file(path: Path) -> pd.DataFrame:
    with xr.open_dataset(path) as ds:
        lat_slice = oriented_slice(ds["latitude"], LAT_MIN, LAT_MAX)
        lon_slice = oriented_slice(ds["longitude"], LON_MIN, LON_MAX)
        sub = ds.sel(latitude=lat_slice, longitude=lon_slice)
        if sub.sizes.get("latitude", 0) == 0 or sub.sizes.get("longitude", 0) == 0:
            raise ValueError(f"No ACM pixels selected in {path}")

        bcm = sub["BCM"].values
        valid_pixel_count = np.isfinite(bcm).sum(axis=(1, 2))
        cloud_pixel_count = (bcm == 1).sum(axis=(1, 2))
        cloud_fraction = np.divide(
            cloud_pixel_count,
            valid_pixel_count,
            out=np.full(cloud_pixel_count.shape, np.nan, dtype=float),
            where=valid_pixel_count > 0,
        )
        return pd.DataFrame(
            {
                "time": pd.to_datetime(sub["time"].values),
                "acm_cloud_fraction": cloud_fraction,
                "acm_cloud_pixel_count": cloud_pixel_count,
                "acm_valid_pixel_count": valid_pixel_count,
                "acm_file": str(path),
            }
        )


def pair_metrics(
    df: pd.DataFrame,
    estimate_col: str,
    obs_col: str = "tsi_frac",
    metric_set_prefix: str = "",
) -> dict[str, float | int | str]:
    valid = df.dropna(subset=[estimate_col, obs_col])
    if valid.empty:
        return {"metric_set": metric_set_prefix.rstrip("_"), "n": 0}
    err = valid[estimate_col] - valid[obs_col]
    return {
        "metric_set": metric_set_prefix.rstrip("_"),
        "n": len(valid),
        "correlation": float(valid[estimate_col].corr(valid[obs_col])),
        "bias": float(err.mean()),
        "mae": float(err.abs().mean()),
        "rmse": float(np.sqrt((err**2).mean())),
        "tsi_mean": float(valid[obs_col].mean()),
        "estimate_mean": float(valid[estimate_col].mean()),
    }


def regression_metrics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def add_metrics(metric_set: str, valid: pd.DataFrame, estimate_col: str) -> None:
        valid = valid.dropna(subset=[estimate_col, "tsi_frac"])
        if valid.empty:
            return
        err = valid[estimate_col] - valid["tsi_frac"]
        row = {
            "metric_set": metric_set,
            "estimate": estimate_col,
            "n": len(valid),
            "correlation": float(valid[estimate_col].corr(valid["tsi_frac"])),
            "mae": float(err.abs().mean()),
            "rmse": float(np.sqrt((err**2).mean())),
            "bias": float(err.mean()),
            "tsi_mean": float(valid["tsi_frac"].mean()),
            "estimate_mean": float(valid[estimate_col].mean()),
        }
        strict = valid.loc[
            (valid["tsi_frac"] < TSI_CLEAR_THRESHOLD)
            | (valid["tsi_frac"] > TSI_CLOUDY_THRESHOLD)
        ].copy()
        if not strict.empty:
            strict["tsi_cloud_binary"] = (strict["tsi_frac"] > TSI_CLOUDY_THRESHOLD).astype(int)
            strict["estimate_cloud_binary"] = (strict[estimate_col] >= 0.5).astype(int)
            report = classification_report(
                strict["tsi_cloud_binary"],
                strict["estimate_cloud_binary"],
                output_dict=True,
                zero_division=0,
            )
            row.update(
                {
                    "strict_binary_n": len(strict),
                    "clear_f1": report["0"]["f1-score"],
                    "cloud_f1": report["1"]["f1-score"],
                    "macro_f1": report["macro avg"]["f1-score"],
                    "weighted_f1": report["weighted avg"]["f1-score"],
                }
            )
        rows.append(row)

    for estimate_col in ["goes_rgb_cloud_fraction", "acm_cloud_fraction"]:
        valid_all = df.dropna(subset=[estimate_col, "tsi_frac"])
        add_metrics("all_matched_tsi", valid_all, estimate_col)
        if "eval_temp_bin_10c" in df:
            for temp_bin, group in valid_all.groupby("eval_temp_bin_10c"):
                add_metrics(f"eval_temp_bin={temp_bin}", group, estimate_col)
    return pd.DataFrame(rows)


def grouped_pair_metrics(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    rows = []
    for group_value, group in df.groupby(group_col, dropna=False):
        for estimate_col in ["goes_rgb_cloud_fraction", "acm_cloud_fraction"]:
            row = pair_metrics(
                group,
                estimate_col=estimate_col,
                metric_set_prefix=str(group_value),
            )
            row[group_col] = group_value
            row["estimate"] = estimate_col
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    cloudy_rules_by_bin = load_cloudy_rules(RULES_CSV)
    trained_midpoints = temp_bin_midpoints(list(cloudy_rules_by_bin))

    era5 = pd.read_csv(ERA5_PIXEL_5MIN_CSV, parse_dates=["time"])
    era5_series = era5.set_index("time")["era5_land_t2m_c"].sort_index()

    tsi = pd.read_csv(TSI_CSV).rename(columns={"t": "time"})
    tsi["time"] = pd.to_datetime(tsi["time"], errors="coerce")
    tsi = tsi.dropna(subset=["time"])
    tsi = tsi.loc[tsi["tsi_frac"].between(0, 1)]
    tsi = tsi.set_index("time").sort_index().reset_index()

    masks = []
    fraction_frames = []
    acm_frames = []
    files = sorted(RGB_DIR.glob("*.nc"))

    for i, path in enumerate(files, start=1):
        file_date = parse_filename_date(path)
        if file_date.day % 2 == 0:
            continue
        acm_path = ACM_DIR / f"goes16_acm_colorado_{file_date:%Y%m%d}.nc"
        if not acm_path.exists():
            print(f"[skip] missing ACM {acm_path}", flush=True)
            continue
        if i == 1 or i % 50 == 0 or i == len(files):
            print(f"[odd-day eval] {i}/{len(files)} {path.name}", flush=True)

        with xr.open_dataset(path) as ds:
            lat_slice = oriented_slice(ds["latitude"], LAT_MIN, LAT_MAX)
            lon_slice = oriented_slice(ds["longitude"], LON_MIN, LON_MAX)
            sub = ds.sel(latitude=lat_slice, longitude=lon_slice)
            if sub.sizes.get("latitude", 0) == 0 or sub.sizes.get("longitude", 0) == 0:
                raise ValueError(f"No pixels selected in {path}")

            times = corrected_times(path, sub)
            temp_at_times = (
                era5_series.reindex(era5_series.index.union(times))
                .sort_index()
                .interpolate(method="time")
                .reindex(times)
                .ffill()
                .bfill()
            )
            temp_bins = np.array(
                [nearest_trained_bin(temp, trained_midpoints) for temp in temp_at_times.values],
                dtype=object,
            )
            original_temp_bins = np.array(
                [ten_c_bin_label(temp) for temp in temp_at_times.values],
                dtype=object,
            )

            red = sub["red"].values
            green = sub["green"].values
            blue = sub["blue"].values
            cloud_mask = classify_pixels(red, green, blue, temp_bins, cloudy_rules_by_bin)
            valid_pixel_count = (cloud_mask != 255).sum(axis=(1, 2))
            cloud_pixel_count = (cloud_mask == 1).sum(axis=(1, 2))
            cloud_fraction = np.divide(
                cloud_pixel_count,
                valid_pixel_count,
                out=np.full(cloud_pixel_count.shape, np.nan, dtype=float),
                where=valid_pixel_count > 0,
            )

            masks.append(
                xr.DataArray(
                    cloud_mask.astype(np.uint8),
                    dims=("time", "latitude", "longitude"),
                    coords={
                        "time": times,
                        "latitude": sub["latitude"].values,
                        "longitude": sub["longitude"].values,
                    },
                )
            )
            fraction_frames.append(
                pd.DataFrame(
                    {
                        "time": times,
                        "goes_rgb_cloud_fraction": cloud_fraction,
                        "cloud_pixel_count": cloud_pixel_count,
                        "valid_pixel_count": valid_pixel_count,
                        "era5_land_t2m_c": temp_at_times.values,
                        "original_temp_bin_10c": original_temp_bins,
                        "eval_temp_bin_10c": temp_bins,
                        "routed_to_nearest_trained_bin": original_temp_bins != temp_bins,
                        "rgb_file": str(path),
                    }
                )
            )
            acm_frames.append(acm_cloud_fraction_for_file(acm_path))

    if not masks:
        raise RuntimeError("No odd-day RGB files processed")

    mask_da = xr.concat(masks, dim="time").sortby("time")
    out_ds = xr.Dataset(
        data_vars={
            "cloud_binary": mask_da,
        },
        attrs={
            "title": "Odd-day eval RGB cloud mask from even-day-trained TSI/ERA5 temperature-bin rules",
            "cloud_binary_values": "0=clear, 1=cloud, 255=invalid RGB pixel",
            "domain": f"lat {LAT_MIN:.6f} to {LAT_MAX:.6f}, lon {LON_MIN:.6f} to {LON_MAX:.6f}",
            "rules_source": str(RULES_CSV),
        },
    )
    out_ds["cloud_binary"].attrs["long_name"] = "RGB decision-tree cloud mask"
    out_ds.to_netcdf(
        OUT_MASK_NC,
        encoding={"cloud_binary": {"zlib": True, "complevel": 4, "dtype": "uint8"}},
    )

    fractions = pd.concat(fraction_frames, ignore_index=True).sort_values("time")
    acm_fractions = pd.concat(acm_frames, ignore_index=True)
    acm_fractions = acm_fractions.set_index("time").sort_index().reset_index()
    fractions = pd.merge_asof(
        fractions,
        tsi[["time", "tsi_frac"]],
        on="time",
        direction="nearest",
        tolerance=TSI_MATCH_TOLERANCE,
    )
    fractions = pd.merge_asof(
        fractions,
        acm_fractions,
        on="time",
        direction="nearest",
        tolerance=ACM_MATCH_TOLERANCE,
    )
    fractions.to_csv(OUT_FRACTION_CSV, index=False)

    metrics = regression_metrics(fractions)
    metrics.to_csv(OUT_METRICS_CSV, index=False)
    matched = fractions.dropna(subset=["tsi_frac"]).copy()
    matched["month"] = matched["time"].dt.to_period("M").astype(str)
    matched["season"] = matched["time"].dt.month.map(
        {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
    )
    grouped_pair_metrics(matched, "month").to_csv(OUT_MONTHLY_METRICS_CSV, index=False)
    grouped_pair_metrics(matched, "season").to_csv(OUT_SEASONAL_METRICS_CSV, index=False)

    print(f"Wrote {OUT_MASK_NC} ({mask_da.sizes['time']} timesteps)")
    print(f"Wrote {OUT_FRACTION_CSV} ({len(fractions)} rows)")
    print(f"Wrote {OUT_METRICS_CSV}")
    print(f"Wrote {OUT_MONTHLY_METRICS_CSV}")
    print(f"Wrote {OUT_SEASONAL_METRICS_CSV}")
    if not metrics.empty:
        print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
