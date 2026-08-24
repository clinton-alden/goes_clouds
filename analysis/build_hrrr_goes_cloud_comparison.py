#!/usr/bin/env python3
"""Collocate June 2022 GOES cloud masks to HRRR and compute statistics."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import xarray as xr


HOURS = (0, 1, *range(14, 24))
DOMAINS = {
    "full": (-109.0, -104.0, 37.0, 41.0),
    "zoom": (-107.25, -106.75, 38.75, 39.25),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    base = Path("/glade/derecho/scratch/cdalden/hrrr_goes_june2022")
    parser.add_argument("--hrrr-dir", type=Path, default=base / "hrrr_daily")
    parser.add_argument("--goes-dir", type=Path, default=base / "goes_masks")
    parser.add_argument("--output-dir", type=Path, default=base / "comparison")
    parser.add_argument("--cloud-threshold", type=float, default=0.5)
    return parser.parse_args()


def nearest_hrrr_mapping(
    goes_lat: np.ndarray,
    goes_lon: np.ndarray,
    hrrr_lat: np.ndarray,
    hrrr_lon: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    glon, glat = np.meshgrid(goes_lon, goes_lat)
    mean_lat = np.deg2rad(np.nanmean(hrrr_lat))
    h_points = np.column_stack((hrrr_lon.ravel() * np.cos(mean_lat), hrrr_lat.ravel()))
    g_points = np.column_stack((glon.ravel() * np.cos(mean_lat), glat.ravel()))
    tree = cKDTree(h_points)
    distance, mapping = tree.query(g_points, k=1)
    weights = np.cos(np.deg2rad(glat.ravel()))
    # A generous 0.08-degree guard prevents accidental extrapolation while
    # retaining edge pixels whose nearest 3-km HRRR center lies just outside.
    valid = np.isfinite(distance) & (distance <= 0.08)
    return mapping.astype(np.int64), weights.astype(np.float64), valid


def upscale_goes_hourly(
    mask: xr.DataArray,
    valid_times: pd.DatetimeIndex,
    mapping: np.ndarray,
    weights: np.ndarray,
    mapping_valid: np.ndarray,
    hrrr_shape: tuple[int, int],
) -> np.ndarray:
    source_times = pd.DatetimeIndex(pd.to_datetime(mask["t"].values)).floor("h")
    values = np.asarray(mask.values, dtype=np.float32)
    output = np.full((len(valid_times), *hrrr_shape), np.nan, dtype=np.float32)
    n_cells = int(np.prod(hrrr_shape))
    for time_index, valid_time in enumerate(valid_times):
        selected = np.where(source_times == valid_time)[0]
        if not len(selected):
            continue
        fine_fraction = np.nanmean(values[selected], axis=0).ravel()
        valid = mapping_valid & np.isfinite(fine_fraction)
        numerator = np.bincount(
            mapping[valid], weights=fine_fraction[valid] * weights[valid], minlength=n_cells
        )
        denominator = np.bincount(mapping[valid], weights=weights[valid], minlength=n_cells)
        coarse = np.divide(
            numerator,
            denominator,
            out=np.full(n_cells, np.nan, dtype=np.float64),
            where=denominator > 0,
        )
        output[time_index] = coarse.reshape(hrrr_shape)
    return output


def metrics(observed: np.ndarray, modeled: np.ndarray, threshold: float) -> dict[str, float | int]:
    observed = np.asarray(np.ma.filled(observed, np.nan), dtype=np.float64)
    modeled = np.asarray(np.ma.filled(modeled, np.nan), dtype=np.float64)
    valid = np.isfinite(observed) & np.isfinite(modeled)
    observed = observed[valid]
    modeled = modeled[valid]
    if not len(observed):
        return {key: np.nan for key in ("bias", "mae", "rmse", "correlation", "accuracy", "pod", "far", "csi")} | {"n": 0}
    residual = modeled - observed
    correlation = np.corrcoef(observed, modeled)[0, 1] if len(observed) > 1 and np.std(observed) > 0 and np.std(modeled) > 0 else np.nan
    obs_cloud = observed >= threshold
    mod_cloud = modeled >= threshold
    # code: 0=correct negative, 1=false alarm, 2=miss, 3=hit
    contingency = np.bincount(
        2 * obs_cloud.astype(np.uint8) + mod_cloud.astype(np.uint8), minlength=4
    )
    correct_negatives, false_alarms, misses, hits = map(int, contingency)
    if hits + misses + false_alarms + correct_negatives != len(observed):
        raise RuntimeError(
            "Binary contingency counts do not cover every valid sample: "
            f"n={len(observed)}, counts={(hits, misses, false_alarms, correct_negatives)}, "
            f"types={(type(observed), type(obs_cloud), type(mod_cloud))}"
        )
    return {
        "n": int(len(observed)),
        "bias": float(np.mean(residual)),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "correlation": float(correlation),
        "accuracy": (hits + correct_negatives) / len(observed),
        "pod": hits / (hits + misses) if hits + misses else np.nan,
        "far": false_alarms / (hits + false_alarms) if hits + false_alarms else np.nan,
        "csi": hits / (hits + misses + false_alarms) if hits + misses + false_alarms else np.nan,
        "hits": hits,
        "misses": misses,
        "false_alarms": false_alarms,
        "correct_negatives": correct_negatives,
    }


def weighted_domain_mean(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> np.ndarray:
    selected = np.where(mask, weights, np.nan)
    valid = np.isfinite(values)
    numerator = np.nansum(values * selected, axis=(-2, -1))
    denominator = np.nansum(np.where(valid, selected, np.nan), axis=(-2, -1))
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    hrrr_files = sorted(args.hrrr_dir.glob("hrrr_tcc_colorado_202206*.nc"))
    goes_files = sorted(args.goes_dir.glob("*_202206*_cloud_binary_tempbin10c.nc"))
    if len(hrrr_files) != 30 or len(goes_files) != 30:
        raise RuntimeError(f"Expected 30 HRRR and 30 GOES files; found {len(hrrr_files)} and {len(goes_files)}")
    goes_by_date = {path.name.split("_colorado_")[1][:8]: path for path in goes_files}

    daily_hrrr = []
    daily_goes = []
    mapping = weights = mapping_valid = None
    hrrr_lat = hrrr_lon = None
    comparisons = None
    for index, hrrr_path in enumerate(hrrr_files, start=1):
        date_key = hrrr_path.stem.rsplit("_", 1)[-1]
        goes_path = goes_by_date[date_key]
        print(f"[{index:02d}/30] {date_key}", flush=True)
        with xr.open_dataset(hrrr_path) as hds, xr.open_dataset(goes_path) as gds:
            hrrr = np.asarray(hds["hrrr_tcc"].values, dtype=np.float32) / 100.0
            valid_times = pd.DatetimeIndex(pd.to_datetime(hds["valid_time"].values))
            if hrrr_lat is None:
                hrrr_lat = np.asarray(hds["latitude"].values, dtype=np.float32)
                hrrr_lon = np.asarray(hds["longitude"].values, dtype=np.float32)
                comparisons = np.asarray(hds["comparison"].values).astype(str)
                mapping, weights, mapping_valid = nearest_hrrr_mapping(
                    np.asarray(gds["latitude"].values),
                    np.asarray(gds["longitude"].values),
                    hrrr_lat,
                    hrrr_lon,
                )
            goes = upscale_goes_hourly(
                gds["cloud_binary"], valid_times, mapping, weights, mapping_valid, hrrr_lat.shape
            )
        daily_hrrr.append(hrrr)
        daily_goes.append(goes)

    times = pd.DatetimeIndex(
        [pd.Timestamp(2022, 6, day) + pd.Timedelta(hours=hour) for day in range(1, 31) for hour in HOURS]
    )
    hrrr_all = np.concatenate(daily_hrrr, axis=1)
    goes_all = np.concatenate(daily_goes, axis=0)
    output_ds = xr.Dataset(
        {
            "hrrr_cloud_fraction": (("comparison", "time", "y", "x"), hrrr_all),
            "goes_cloud_fraction": (("time", "y", "x"), goes_all),
        },
        coords={
            "comparison": comparisons,
            "time": times,
            "latitude": (("y", "x"), hrrr_lat),
            "longitude": (("y", "x"), hrrr_lon),
        },
        attrs={
            "title": "Hourly HRRR and upscaled GOES RGB cloud fractions for June 2022",
            "goes_upscaling": "GOES hourly pixel means assigned to nearest HRRR grid center with cosine-latitude weighting",
        },
    )
    output_ds["hrrr_cloud_fraction"].attrs["units"] = "1"
    output_ds["goes_cloud_fraction"].attrs["units"] = "1"
    collocated_path = args.output_dir / "hrrr_goes_cloud_collocated_june2022.nc"
    output_ds.to_netcdf(
        collocated_path,
        encoding={
            "hrrr_cloud_fraction": {"zlib": True, "complevel": 4, "dtype": "float32"},
            "goes_cloud_fraction": {"zlib": True, "complevel": 4, "dtype": "float32"},
        },
    )

    # Compute statistics from the persisted product that downstream users read.
    # The netCDF4 backend can apply encoding in-place to arrays owned by the
    # Dataset during serialization, so do not rely on the pre-write buffers.
    with xr.open_dataset(collocated_path) as persisted:
        persisted.load()
        hrrr_all = np.array(persisted["hrrr_cloud_fraction"].values, copy=True)
        goes_all = np.array(persisted["goes_cloud_fraction"].values, copy=True)
        hrrr_lat = np.array(persisted["latitude"].values, copy=True)
        hrrr_lon = np.array(persisted["longitude"].values, copy=True)
        times = pd.DatetimeIndex(np.array(persisted["time"].values, copy=True))
        time_hours = np.array(persisted["time"].dt.hour.values, copy=True)

    area_weights = np.cos(np.deg2rad(hrrr_lat)).astype(np.float64)
    summary_rows = []
    timeseries_rows = []
    for domain_name, (lon_min, lon_max, lat_min, lat_max) in DOMAINS.items():
        domain_mask = (
            (hrrr_lon >= lon_min) & (hrrr_lon <= lon_max) &
            (hrrr_lat >= lat_min) & (hrrr_lat <= lat_max)
        )
        goes_mean = weighted_domain_mean(goes_all, area_weights, domain_mask)
        for comparison_index, comparison in enumerate(comparisons):
            hrrr_values = hrrr_all[comparison_index]
            hrrr_mean = weighted_domain_mean(hrrr_values, area_weights, domain_mask)
            for time_index, timestamp in enumerate(times):
                timeseries_rows.append({
                    "time": timestamp,
                    "utc_hour": timestamp.hour,
                    "domain": domain_name,
                    "comparison": comparison,
                    "goes_cloud_fraction": goes_mean[time_index],
                    "hrrr_cloud_fraction": hrrr_mean[time_index],
                })
            for hour in HOURS:
                time_sel = time_hours == hour
                spatial_obs = goes_all[time_sel][:, domain_mask].ravel()
                spatial_mod = hrrr_values[time_sel][:, domain_mask].ravel()
                row = {"domain": domain_name, "comparison": comparison, "utc_hour": hour, "aggregation": "collocated_grid_cells"}
                row.update(metrics(spatial_obs, spatial_mod, args.cloud_threshold))
                summary_rows.append(row)

                row = {"domain": domain_name, "comparison": comparison, "utc_hour": hour, "aggregation": "domain_mean_time_series"}
                row.update(metrics(goes_mean[time_sel], hrrr_mean[time_sel], args.cloud_threshold))
                summary_rows.append(row)

    pd.DataFrame(timeseries_rows).to_csv(args.output_dir / "hourly_domain_timeseries.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(args.output_dir / "hourly_summary_statistics.csv", index=False)
    print(f"Wrote {collocated_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
