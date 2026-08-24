#!/usr/bin/env python3
"""Compute statistics from a persisted HRRR-GOES comparison file."""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from build_hrrr_goes_cloud_comparison import DOMAINS, HOURS, metrics, weighted_domain_mean

ROOT = Path("/glade/derecho/scratch/cdalden/hrrr_goes_june2022/comparison")


def main() -> int:
    with xr.open_dataset(ROOT / "hrrr_goes_cloud_collocated_june2022.nc") as ds:
        ds.load()
        hrrr = np.array(ds.hrrr_cloud_fraction.values, copy=True)
        goes = np.array(ds.goes_cloud_fraction.values, copy=True)
        lat = np.array(ds.latitude.values, copy=True)
        lon = np.array(ds.longitude.values, copy=True)
        times = pd.DatetimeIndex(np.array(ds.time.values, copy=True))
        comparisons = np.array(ds.comparison.values).astype(str)
    weights = np.cos(np.deg2rad(lat)).astype(np.float64)
    hours = times.hour.to_numpy(copy=True)
    summary, series = [], []
    for domain, (xmin, xmax, ymin, ymax) in DOMAINS.items():
        mask = (lon >= xmin) & (lon <= xmax) & (lat >= ymin) & (lat <= ymax)
        goes_mean = weighted_domain_mean(goes, weights, mask)
        for ci, comparison in enumerate(comparisons):
            modeled = hrrr[ci]
            modeled_mean = weighted_domain_mean(modeled, weights, mask)
            for ti, timestamp in enumerate(times):
                series.append({"time": timestamp, "utc_hour": timestamp.hour,
                               "domain": domain, "comparison": comparison,
                               "goes_cloud_fraction": goes_mean[ti],
                               "hrrr_cloud_fraction": modeled_mean[ti]})
            for hour in HOURS:
                selected = hours == hour
                row = {"domain": domain, "comparison": comparison,
                       "utc_hour": hour, "aggregation": "collocated_grid_cells"}
                row.update(metrics(goes[selected][:, mask].ravel(),
                                   modeled[selected][:, mask].ravel(), 0.5))
                summary.append(row)
                row = {"domain": domain, "comparison": comparison,
                       "utc_hour": hour, "aggregation": "domain_mean_time_series"}
                row.update(metrics(goes_mean[selected], modeled_mean[selected], 0.5))
                summary.append(row)
    pd.DataFrame(series).to_csv(ROOT / "hourly_domain_timeseries.csv", index=False)
    pd.DataFrame(summary).to_csv(ROOT / "hourly_summary_statistics.csv", index=False)
    print(f"Wrote summaries in {ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
