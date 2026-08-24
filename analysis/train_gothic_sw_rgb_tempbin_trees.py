#!/usr/bin/env python3
"""Train Gothic RGB decision trees against SW-derived cloud labels.

Outputs:
  - gothic_rgb_SW_domain.csv
  - gothic_sw_cloud_binary_5min.csv
  - gothic_observed_air_temp_5min.csv
  - gothic_rgb_sw_training_table.csv
  - gothic_rgb_tempbin_decision_tree_rules.csv
  - gothic_rgb_tempbin_model_metrics.csv
"""

from __future__ import annotations

import glob
import os
import re
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier, _tree

try:
    from metsim import constants as metsim_constants
    from metsim.disaggregate import shortwave as metsim_shortwave
    from metsim.physics import solar_geom
except ModuleNotFoundError:
    # This machine has metsim unpacked in the conda package cache, but the active
    # kernel does not expose it or numba. Metsim only needs numba's jit decorator
    # for acceleration here, so a no-op shim is enough for this notebook workflow.
    metsim_pkg = Path("/glade/work/cdalden/.conda/pkgs/metsim-2.4.4-pyhd8ed1ab_0/site-packages")
    if metsim_pkg.exists():
        sys.path.append(str(metsim_pkg))
    if "numba" not in sys.modules:
        numba_stub = types.ModuleType("numba")

        def _jit(*args, **kwargs):
            if args and callable(args[0]) and len(args) == 1 and not kwargs:
                return args[0]
            return lambda func: func

        numba_stub.jit = _jit
        sys.modules["numba"] = numba_stub
    from metsim import constants as metsim_constants
    from metsim.disaggregate import shortwave as metsim_shortwave
    from metsim.physics import solar_geom


ROOT = Path("/glade/u/home/cdalden/goes_work/analysis")
OUT_DIR = ROOT / "output_11d_gothic"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RGB_DIR = Path("/glade/derecho/scratch/cdalden/gothic/goes16/rgb_composite")
SW_NC = Path("/glade/u/home/cdalden/scratch/surface_obs/colorado/sos_full_dataset_30min.nc")
SW_DAILY_DIR = Path("/glade/u/home/cdalden/scratch/surface_obs/colorado/gucsebsM1.b1")
ERA5_T2M_DIR = Path("/glade/derecho/scratch/cdalden/gothic/era5/t2m_hourly")

START_DATE = "2021-10-01"
END_DATE = "2023-06-30"
LOCAL_TZ = "America/Denver"
LOCAL_START = "09:00"
LOCAL_END = "14:00"
TIME_STEP_MIN = 5

SW_VAR = "Rsw_in_9m_d"
SW_DAILY_VAR = "down_short_hemisp"
OBS_TEMP_VAR = "T_2m_c"
ERA5_TEMP_VAR = "t2m"
TEMP_COL = "observed_air_temp_c"
CLOUDY_RATIO = 0.40
CLEAR_RATIO = 0.80

GOTHIC_LAT = 38.96
GOTHIC_LON = -106.99
GOTHIC_ELEV_M = 2915.0
METSIM_LAPSE_RATE = 0.0065
METSIM_TBASE = 0.95

DOMAIN_LON_MIN = -107.04666
DOMAIN_LON_MAX = -106.93268
DOMAIN_LAT_MIN = 38.91294
DOMAIN_LAT_MAX = 38.95834

TREE_FEATURES = ["red", "green", "blue"]
TEMP_BINS = [(-np.inf, -10.0), (-10.0, 0.0), (0.0, 10.0), (10.0, np.inf)]
TREE_MAX_DEPTH = 2
MIN_ROWS_PER_BIN = 40
MIN_CLASS_COUNT_PER_BIN = 10


def local_time_overlap(
    rgb: pd.DataFrame,
    sw_cloud: pd.DataFrame,
) -> tuple[pd.Timestamp, pd.Timestamp]:
    rgb_start = pd.to_datetime(rgb["time_local"]).min()
    rgb_end = pd.to_datetime(rgb["time_local"]).max()
    sw_start = pd.to_datetime(sw_cloud["time_local"]).min()
    sw_end = pd.to_datetime(sw_cloud["time_local"]).max()
    overlap_start = max(rgb_start, sw_start)
    overlap_end = min(rgb_end, sw_end)
    if pd.isna(overlap_start) or pd.isna(overlap_end) or overlap_start > overlap_end:
        raise ValueError(
            "No overlapping local times between RGB domain means and SW cloud labels: "
            f"RGB={rgb_start} to {rgb_end}, SW={sw_start} to {sw_end}"
        )
    return overlap_start, overlap_end


def clip_to_time_overlap(
    rgb: pd.DataFrame,
    sw_cloud: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    overlap_start, overlap_end = local_time_overlap(rgb, sw_cloud)
    rgb_overlap = rgb[
        pd.to_datetime(rgb["time_local"]).between(overlap_start, overlap_end, inclusive="both")
    ].copy()
    sw_overlap = sw_cloud[
        pd.to_datetime(sw_cloud["time_local"]).between(overlap_start, overlap_end, inclusive="both")
    ].copy()
    print(
        "Using RGB/SW local-time overlap: "
        f"{overlap_start} to {overlap_end} "
        f"({len(rgb_overlap)} RGB rows, {len(sw_overlap)} SW-label rows)",
        flush=True,
    )
    return rgb_overlap, sw_overlap


def filter_local_time_window(df: pd.DataFrame, time_col: str) -> pd.DataFrame:
    timestamps = pd.to_datetime(df[time_col])
    start_time = pd.Timestamp(LOCAL_START).time()
    end_time = pd.Timestamp(LOCAL_END).time()
    mask = timestamps.dt.time.between(start_time, end_time)
    return df.loc[mask].copy()


def daylight_saving_offset_hours(days: pd.DatetimeIndex, local_tz: str) -> np.ndarray:
    noons = pd.DatetimeIndex(days) + pd.Timedelta(hours=12)
    return np.array(
        [noon.tz_localize(local_tz).dst().total_seconds() / 3600 for noon in noons]
    )


def build_metsim_clear_sky(
    times_local: pd.Series,
    lat: float,
    lon: float,
    elev_m: float,
    lapse_rate: float,
    time_step_min: int,
    tbase: float,
    local_tz: str,
) -> pd.DataFrame:
    metsim_constants.TBASE = tbase
    try:
        solar_geom.recompile()
    except AttributeError:
        pass

    tiny_rad_fract, daylength, flat_potrad, tt_max0 = solar_geom(elev_m, lat, lapse_rate)
    all_days = pd.date_range(times_local.min().normalize(), times_local.max().normalize(), freq="D")
    yday_zero_based = all_days.dayofyear.to_numpy() - 1
    clear_daily = flat_potrad[yday_zero_based] * tt_max0[yday_zero_based]
    params = {"time_step": time_step_min, "method": "mtclim", "utc_offset": False, "lon": lon}
    clear_subdaily = metsim_shortwave(
        clear_daily,
        daylength[yday_zero_based],
        yday_zero_based + 1,
        tiny_rad_fract,
        params,
    )

    clear_times_solar = pd.date_range(
        all_days.min(), periods=len(clear_subdaily), freq=f"{time_step_min}min"
    )
    day_offsets = pd.DataFrame(
        {
            "date_solar": all_days,
            "dst_offset_hours": daylight_saving_offset_hours(all_days, local_tz),
        }
    )
    clear_df = pd.DataFrame({"time_solar": clear_times_solar, "sw_clear": clear_subdaily})
    clear_df["date_solar"] = clear_df["time_solar"].dt.normalize()
    clear_df = clear_df.merge(day_offsets, on="date_solar", how="left")
    clear_df["time_local"] = clear_df["time_solar"] + pd.to_timedelta(
        clear_df["dst_offset_hours"], unit="h"
    )
    return clear_df[["time_local", "sw_clear"]]


def latlon_to_goes_scan(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    try:
        from pyproj import CRS, Transformer
    except ImportError as exc:
        raise ImportError("pyproj is required to map the lat/lon box onto GOES x/y scan angles") from exc

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


def collect_rgb_domain_mean() -> pd.DataFrame:
    out = OUT_DIR / "gothic_rgb_SW_domain.csv"
    if out.exists():
        rgb = pd.read_csv(out, parse_dates=["time_utc", "time_local", "time_local_5min"])
        print(f"Reusing RGB domain means: {out} ({len(rgb)} rows)", flush=True)
        return rgb

    x_min, x_max, y_min, y_max = domain_scan_bounds()
    rows = []
    files = sorted(
        p
        for p in glob.glob(str(RGB_DIR / "goes16_C02_C05_C13_rgb_gothic_*.nc"))
        if START_DATE.replace("-", "") <= Path(p).stem.rsplit("_", 1)[-1] <= END_DATE.replace("-", "")
    )
    if not files:
        raise FileNotFoundError(f"No Gothic RGB files found in {RGB_DIR}")

    for idx, path in enumerate(files, start=1):
        if idx == 1 or idx % 25 == 0 or idx == len(files):
            print(f"Reading RGB file {idx}/{len(files)}: {Path(path).name}", flush=True)
        with xr.open_dataset(path) as ds:
            x_slice = slice(x_min, x_max) if ds["x"][0] < ds["x"][-1] else slice(x_max, x_min)
            y_slice = slice(y_min, y_max) if ds["y"][0] < ds["y"][-1] else slice(y_max, y_min)
            subset = ds.sel(x=x_slice, y=y_slice)
            if subset.sizes.get("x", 0) == 0 or subset.sizes.get("y", 0) == 0:
                raise ValueError(
                    f"Lat/lon domain selected no pixels in {path}; "
                    f"scan bounds x=({x_min}, {x_max}), y=({y_min}, {y_max})"
                )
            mean_ds = subset[TREE_FEATURES].mean(dim=("y", "x"), skipna=True)
            df = mean_ds.to_dataframe().reset_index()
            df["source_file"] = Path(path).name
            df["n_pixels"] = subset.sizes["x"] * subset.sizes["y"]
            rows.append(df)

    rgb = pd.concat(rows, ignore_index=True)
    rgb = rgb.rename(columns={"t": "time_utc"})
    rgb["time_utc"] = pd.to_datetime(rgb["time_utc"])
    rgb["time_local"] = (
        rgb["time_utc"].dt.tz_localize("UTC").dt.tz_convert(LOCAL_TZ).dt.tz_localize(None)
    )
    rgb["time_local_5min"] = rgb["time_local"].dt.floor(f"{TIME_STEP_MIN}min")
    rgb = rgb.set_index("time_local_5min").between_time(LOCAL_START, LOCAL_END).reset_index()
    rgb = rgb.sort_values("time_local_5min").drop_duplicates("time_local_5min")

    rgb.to_csv(out, index=False)
    print(f"Saved RGB domain means: {out} ({len(rgb)} rows)", flush=True)
    return rgb


def read_sw_daily_files() -> pd.DataFrame:
    files = sorted(SW_DAILY_DIR.glob("gucsebsM1.b1.*.custom.cdf"))
    rows = []
    start_utc = pd.Timestamp(START_DATE) - pd.Timedelta(days=1)
    end_utc = pd.Timestamp(END_DATE) + pd.Timedelta(days=1)

    for path in files:
        match = re.search(r"\.(\d{8})\.", path.name)
        if not match:
            continue
        file_day_utc = pd.Timestamp(match.group(1))
        if not start_utc <= file_day_utc <= end_utc:
            continue

        with xr.open_dataset(path, decode_times=True) as ds:
            if SW_DAILY_VAR not in ds:
                raise KeyError(f"Expected {SW_DAILY_VAR} in {path}")
            values = np.asarray(ds[SW_DAILY_VAR].to_numpy(), dtype=float)

        time_utc = pd.date_range(file_day_utc, periods=len(values), freq="30min")
        time_local = time_utc.tz_localize("UTC").tz_convert(LOCAL_TZ).tz_localize(None)
        rows.append(
            pd.DataFrame(
                {
                    "time_local": time_local,
                    "sw_obs": values,
                    "sw_obs_var": SW_DAILY_VAR,
                    "sw_obs_source": str(path),
                }
            )
        )

    if not rows:
        raise FileNotFoundError(f"No {SW_DAILY_VAR} daily files found in {SW_DAILY_DIR}")

    sw = pd.concat(rows, ignore_index=True)
    sw = sw.dropna(subset=["time_local", "sw_obs"]).sort_values("time_local")
    sw = sw[
        (sw["time_local"] >= pd.Timestamp(START_DATE))
        & (sw["time_local"] <= pd.Timestamp(END_DATE) + pd.Timedelta(days=1))
    ].copy()
    sw = sw.drop_duplicates("time_local", keep="last")
    print(
        f"Loaded SW daily files from {SW_DAILY_DIR}: "
        f"{len(sw)} valid {SW_DAILY_VAR} rows from "
        f"{sw['time_local'].min()} to {sw['time_local'].max()}",
        flush=True,
    )
    return sw


def build_sw_cloud_binary() -> pd.DataFrame:
    sw = read_sw_daily_files()
    sw_series = sw.set_index("time_local")["sw_obs"].sort_index()
    sw_5min_times = pd.date_range(
        sw_series.index.min().floor("D"),
        sw_series.index.max().ceil("D"),
        freq=f"{TIME_STEP_MIN}min",
    )
    sw_5min_times = sw_5min_times[
        (sw_5min_times >= pd.Timestamp(START_DATE))
        & (sw_5min_times <= pd.Timestamp(END_DATE) + pd.Timedelta(days=1))
    ]
    sw_at_rgb = (
        sw_series.reindex(sw_series.index.union(sw_5min_times))
        .sort_index()
        .interpolate(method="time")
        .reindex(sw_5min_times)
    )
    sw_5 = pd.DataFrame({"time_local": sw_5min_times, "sw_obs": sw_at_rgb.to_numpy(dtype=float)})
    sw_5 = sw_5.dropna(subset=["sw_obs"])
    sw_5 = sw_5.set_index("time_local").between_time(LOCAL_START, LOCAL_END).reset_index()

    clear_df = build_metsim_clear_sky(
        times_local=sw_5["time_local"],
        lat=GOTHIC_LAT,
        lon=GOTHIC_LON,
        elev_m=GOTHIC_ELEV_M,
        lapse_rate=METSIM_LAPSE_RATE,
        time_step_min=TIME_STEP_MIN,
        tbase=METSIM_TBASE,
        local_tz=LOCAL_TZ,
    )
    sw_cloud = sw_5.merge(clear_df, on="time_local", how="inner")
    sw_cloud.loc[sw_cloud["sw_clear"] <= 0, "sw_clear"] = np.nan
    sw_cloud["sw_ratio"] = sw_cloud["sw_obs"] / sw_cloud["sw_clear"]
    sw_cloud["cloud_binary"] = np.select(
        [sw_cloud["sw_ratio"] < CLOUDY_RATIO, sw_cloud["sw_ratio"] > CLEAR_RATIO],
        [1, 0],
        default=np.nan,
    )
    sw_cloud["sw_flag"] = np.select(
        [sw_cloud["cloud_binary"].eq(1), sw_cloud["cloud_binary"].eq(0)],
        ["cloudy", "clear"],
        default="ambiguous",
    )
    sw_cloud = sw_cloud.dropna(subset=["cloud_binary"]).copy()
    sw_cloud["cloud_binary"] = sw_cloud["cloud_binary"].astype(int)

    out = OUT_DIR / "gothic_sw_cloud_binary_5min.csv"
    sw_cloud.to_csv(out, index=False)
    print(f"Saved SW cloud binary: {out} ({len(sw_cloud)} rows)")
    return sw_cloud


def collect_era5_air_temp_at_times(times_local: pd.Series) -> pd.DataFrame:
    times = pd.DatetimeIndex(pd.to_datetime(times_local)).sort_values().unique()
    files = sorted(ERA5_T2M_DIR.glob("era5_t2m_gothic_*.nc"))
    if not files:
        raise FileNotFoundError(f"No ERA5 {ERA5_TEMP_VAR} files found in {ERA5_T2M_DIR}")

    rows = []
    for path in files:
        with xr.open_dataset(path, decode_times=True) as ds:
            if ERA5_TEMP_VAR not in ds:
                raise KeyError(f"Expected {ERA5_TEMP_VAR} in {path}")
            da = ds[ERA5_TEMP_VAR]
            time_name = "valid_time" if "valid_time" in da.coords else "time"
            spatial_dims = [dim for dim in da.dims if dim != time_name]
            if spatial_dims:
                da = da.mean(dim=spatial_dims, skipna=True)
            era5_times_local = (
                pd.DatetimeIndex(pd.to_datetime(da[time_name].values))
                .tz_localize("UTC")
                .tz_convert(LOCAL_TZ)
                .tz_localize(None)
            )
            rows.append(
                pd.DataFrame(
                    {
                        "time_local": era5_times_local,
                        TEMP_COL: np.asarray(da.to_numpy(), dtype=float) - 273.15,
                        "observed_air_temp_var": ERA5_TEMP_VAR,
                        "observed_air_temp_source": str(path),
                    }
                )
            )

    temp = pd.concat(rows, ignore_index=True)
    temp = temp.dropna(subset=["time_local", TEMP_COL]).sort_values("time_local")
    temp = temp[
        (temp["time_local"] >= times.min() - pd.Timedelta(hours=1))
        & (temp["time_local"] <= times.max() + pd.Timedelta(hours=1))
    ].copy()
    if temp.empty:
        raise ValueError(
            f"No ERA5 {ERA5_TEMP_VAR} temperatures overlap requested times "
            f"{times.min()} to {times.max()}"
        )

    temp_series = temp.set_index("time_local")[TEMP_COL].sort_index().groupby(level=0).mean()
    temp_at_times = (
        temp_series.reindex(temp_series.index.union(times))
        .sort_index()
        .interpolate(method="time")
        .reindex(times)
    )
    out_df = pd.DataFrame(
        {
            "time_local": times,
            TEMP_COL: temp_at_times.to_numpy(dtype=float),
            "observed_air_temp_var": ERA5_TEMP_VAR,
            "observed_air_temp_source": str(ERA5_T2M_DIR),
        }
    )

    out = OUT_DIR / "gothic_observed_air_temp_5min.csv"
    out_df.to_csv(out, index=False)
    print(f"Saved ERA5 Gothic air temperature: {out} ({len(out_df)} rows)")
    return out_df


def temp_bin_label(left: float, right: float) -> str:
    if np.isneginf(left):
        return f"<{int(right)}"
    if np.isposinf(right):
        return f">={int(left)}"
    return f"[{int(left)}, {int(right)})"


def assign_temp_bin(temp_c: float) -> str:
    if pd.isna(temp_c):
        return np.nan
    for left, right in TEMP_BINS:
        if left <= float(temp_c) < right:
            return temp_bin_label(left, right)
    raise RuntimeError(f"Could not bin temperature {temp_c}")


def export_tree_rules(clf: DecisionTreeClassifier, temp_bin: str, metrics_row: dict) -> list[dict]:
    tree = clf.tree_
    rows = []

    def walk(node: int, clauses: list[str]) -> None:
        if tree.feature[node] != _tree.TREE_UNDEFINED:
            feature = TREE_FEATURES[tree.feature[node]]
            threshold = float(tree.threshold[node])
            walk(node=tree.children_left[node], clauses=clauses + [f"{feature} <= {threshold:.8g}"])
            walk(node=tree.children_right[node], clauses=clauses + [f"{feature} > {threshold:.8g}"])
            return

        counts = tree.value[node][0]
        pred = int(np.argmax(counts))
        rows.append(
            {
                **metrics_row,
                "leaf_id": int(node),
                "rule": " and ".join(clauses) if clauses else "always",
                "prediction": pred,
                "prediction_label": "cloudy" if pred == 1 else "clear",
                "leaf_clear_weight": float(counts[0]),
                "leaf_cloudy_weight": float(counts[1]),
            }
        )

    walk(0, [])
    return rows


def train_tempbin_trees(train_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    train_df = train_df.dropna(subset=TREE_FEATURES + ["cloud_binary", TEMP_COL]).copy()
    train_df["temp_bin"] = train_df[TEMP_COL].map(assign_temp_bin)

    rules = []
    metrics = []
    for left, right in TEMP_BINS:
        temp_bin = temp_bin_label(left, right)
        grp = train_df[train_df["temp_bin"].eq(temp_bin)].copy()
        n_rows = len(grp)
        counts = grp["cloud_binary"].value_counts()
        n_clear = int(counts.get(0, 0))
        n_cloudy = int(counts.get(1, 0))
        metrics_row = {
            "temp_bin": temp_bin,
            "temp_left_c": left,
            "temp_right_c": right,
            "n_rows": n_rows,
            "n_clear": n_clear,
            "n_cloudy": n_cloudy,
        }

        if n_rows < MIN_ROWS_PER_BIN or min(n_clear, n_cloudy) < MIN_CLASS_COUNT_PER_BIN:
            metrics.append({**metrics_row, "status": "skipped_low_sample"})
            continue

        X = grp[TREE_FEATURES]
        y = grp["cloud_binary"].astype(int)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.25, random_state=42, stratify=y
        )
        clf = DecisionTreeClassifier(
            max_depth=TREE_MAX_DEPTH,
            min_samples_split=20,
            class_weight="balanced",
            random_state=42,
        )
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        rep = classification_report(
            y_test,
            y_pred,
            labels=[0, 1],
            target_names=["Clear", "Cloudy"],
            output_dict=True,
            zero_division=0,
        )
        full_metrics = {
            **metrics_row,
            "accuracy": float(rep["accuracy"]),
            "clear_precision": float(rep["Clear"]["precision"]),
            "clear_recall": float(rep["Clear"]["recall"]),
            "clear_f1": float(rep["Clear"]["f1-score"]),
            "cloudy_precision": float(rep["Cloudy"]["precision"]),
            "cloudy_recall": float(rep["Cloudy"]["recall"]),
            "cloudy_f1": float(rep["Cloudy"]["f1-score"]),
            "macro_f1": float(rep["macro avg"]["f1-score"]),
            "status": "trained",
        }
        metrics.append(full_metrics)
        rules.extend(export_tree_rules(clf, temp_bin, full_metrics))

    metrics_df = pd.DataFrame(metrics)
    rules_df = pd.DataFrame(rules)
    metrics_out = OUT_DIR / "gothic_rgb_tempbin_model_metrics.csv"
    rules_out = OUT_DIR / "gothic_rgb_tempbin_decision_tree_rules.csv"
    metrics_df.to_csv(metrics_out, index=False)
    rules_df.to_csv(rules_out, index=False)
    print(f"Saved temp-bin metrics: {metrics_out}")
    print(f"Saved temp-bin decision tree rules: {rules_out}")
    return metrics_df, rules_df


def main() -> int:
    rgb = collect_rgb_domain_mean()
    sw_cloud = build_sw_cloud_binary()
    rgb, sw_cloud = clip_to_time_overlap(rgb, sw_cloud)
    obs_temp = collect_era5_air_temp_at_times(sw_cloud["time_local"])

    rgb_for_merge = (
        rgb.sort_values("time_local")
        .drop_duplicates("time_local")
        .rename(columns={"time_local": "rgb_time_local", "time_utc": "rgb_time_utc"})
    )
    sw_cloud["time_local"] = pd.to_datetime(sw_cloud["time_local"]).astype("datetime64[ns]")
    rgb_for_merge["rgb_time_local"] = pd.to_datetime(rgb_for_merge["rgb_time_local"]).astype(
        "datetime64[ns]"
    )
    train_df = pd.merge_asof(
        sw_cloud.sort_values("time_local"),
        rgb_for_merge,
        left_on="time_local",
        right_on="rgb_time_local",
        direction="nearest",
        tolerance=pd.Timedelta(minutes=3),
    )
    train_df = train_df.dropna(subset=TREE_FEATURES).copy()
    train_df["rgb_time_delta_seconds"] = (
        train_df["rgb_time_local"] - train_df["time_local"]
    ).dt.total_seconds()
    train_df = train_df.merge(obs_temp, on="time_local", how="inner")
    train_df = filter_local_time_window(train_df, "time_local")
    train_out = OUT_DIR / "gothic_rgb_sw_training_table.csv"
    train_df.to_csv(train_out, index=False)
    print(
        f"Saved training table: {train_out} ({len(train_df)} rows, "
        f"{LOCAL_START}-{LOCAL_END} {LOCAL_TZ})"
    )

    metrics_df, rules_df = train_tempbin_trees(train_df)
    print("\nTraining rows by temperature bin:")
    print(metrics_df.to_string(index=False))
    print("\nRule count by temperature bin:")
    if rules_df.empty:
        print("No trees were trained.")
    else:
        print(rules_df.groupby("temp_bin").size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
