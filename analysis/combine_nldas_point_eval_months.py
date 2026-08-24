#!/usr/bin/env python
"""Combine monthly Gothic NLDAS point-eval CSVs and remake annual plots."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from compare_2022_nldas_point_sw import plot_density_and_bias


DEFAULT_OUT_DIR = Path("analysis/output_17_nldas_point_eval")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2022)
    parser.add_argument("--hours", default="14-23,0")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def monthly_path(out_dir: Path, year: int, month: int, hours: str) -> Path:
    import calendar

    start = f"{year:04d}{month:02d}01"
    end = f"{year:04d}{month:02d}{calendar.monthrange(year, month)[1]:02d}"
    tag = f"{start}_{end}_h{hours.replace(',', '-')}"
    return out_dir / f"gothic_nldas_point_sw_{tag}.csv"


def compute_metrics(compare: pd.DataFrame) -> pd.DataFrame:
    compare = compare.copy()
    compare["month"] = compare["time_utc"].dt.to_period("M").astype(str)
    rows = []
    for source in ("nldas2", "nldas3"):
        if source not in compare:
            continue
        for month, group in compare.groupby("month"):
            valid = group[[source, "gothic_obs_sw"]].dropna()
            if valid.empty:
                continue
            err = valid[source] - valid["gothic_obs_sw"]
            rows.append(
                {
                    "source": source,
                    "month": month,
                    "n": len(valid),
                    "bias_wm2": err.mean(),
                    "mae_wm2": err.abs().mean(),
                    "rmse_wm2": (err.pow(2).mean()) ** 0.5,
                    "corr": valid[source].corr(valid["gothic_obs_sw"]) if len(valid) > 1 else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    frames = []
    missing = []
    for month in range(1, 13):
        path = monthly_path(args.out_dir, args.year, month, args.hours)
        if path.exists():
            frames.append(pd.read_csv(path, parse_dates=["time_utc"]))
        else:
            missing.append(path)

    if not frames:
        raise FileNotFoundError(f"No monthly point-eval CSVs found in {args.out_dir}")

    compare = (
        pd.concat(frames, ignore_index=True)
        .drop_duplicates("time_utc", keep="last")
        .sort_values("time_utc")
    )
    metrics = compute_metrics(compare)

    tag = f"{args.year}_h{args.hours.replace(',', '-')}"
    compare_csv = args.out_dir / f"gothic_nldas_point_sw_compare_{tag}.csv"
    metrics_csv = args.out_dir / f"gothic_nldas_point_sw_monthly_bias_{tag}.csv"
    figure_png = args.out_dir / f"gothic_nldas_point_sw_density_bias_{tag}.png"

    compare.to_csv(compare_csv, index=False)
    metrics.to_csv(metrics_csv, index=False)
    plot_density_and_bias(compare, metrics, figure_png, args.hours)

    print(f"Wrote {compare_csv}")
    print(f"Wrote {metrics_csv}")
    print(f"Wrote {figure_png}")
    if missing:
        print("Missing monthly CSVs:")
        for path in missing:
            print(path)
    print(metrics)


if __name__ == "__main__":
    main()
