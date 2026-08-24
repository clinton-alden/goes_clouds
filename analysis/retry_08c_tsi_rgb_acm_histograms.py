#!/usr/bin/env python3
"""Retry notebook 08c TSI/RGB/ACM histogram comparisons.

This uses the existing domain cloud-fraction table from notebook 08c:
TSI cloud fraction, RGB cloud fraction over the domain, and GOES ACM cloud
fraction. It does not use ERA5 or decision-tree temperature bins.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("/glade/u/home/cdalden/goes_work/analysis")
GOES_CF_CSV = ROOT / "output_08c/colorado_domain_cloud_fraction_14z_00z.csv"
TSI_CSV = Path("/glade/u/home/cdalden/scratch/surface_obs/colorado/sail_tsi_cloud_frac.csv")
OUT_DIR = ROOT / "output_08c_retry_tsi_rgb_acm"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TSI_MATCH_TOLERANCE = pd.Timedelta("3min")


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


def build_compare() -> pd.DataFrame:
    goes = pd.read_csv(GOES_CF_CSV, parse_dates=["date", "t"])
    goes["t"] = pd.to_datetime(goes["t"]).astype("datetime64[ns]")
    goes = goes.sort_values("t")

    tsi = pd.read_csv(TSI_CSV, parse_dates=["t"])
    tsi["t"] = pd.to_datetime(tsi["t"]).astype("datetime64[ns]")
    tsi.loc[tsi["tsi_frac"] < 0, "tsi_frac"] = np.nan
    tsi = tsi.sort_values("t")

    compare = pd.merge_asof(
        goes,
        tsi[["t", "tsi_frac"]],
        on="t",
        direction="nearest",
        tolerance=TSI_MATCH_TOLERANCE,
    )
    compare["rgb_minus_tsi"] = compare["rgb_cloud_frac"] - compare["tsi_frac"]
    compare["acm_minus_tsi"] = compare["acm_cloud_frac"] - compare["tsi_frac"]
    compare.to_csv(OUT_DIR / "colorado_08c_tsi_rgb_acm_compare.csv", index=False)
    return compare


def save_summary(compare: pd.DataFrame) -> pd.DataFrame:
    summary = pd.DataFrame(
        {
            "metric": [
                "Rows",
                "Rows with valid TSI",
                "Mean TSI cloud fraction",
                "Mean RGB cloud fraction",
                "Mean ACM cloud fraction",
                "Corr(RGB, TSI)",
                "Corr(ACM, TSI)",
                "RGB mean residual",
                "ACM mean residual",
                "RGB MAE",
                "ACM MAE",
                "RGB RMSE",
                "ACM RMSE",
            ],
            "value": [
                len(compare),
                compare["tsi_frac"].notna().sum(),
                compare["tsi_frac"].mean(),
                compare["rgb_cloud_frac"].mean(),
                compare["acm_cloud_frac"].mean(),
                corr_or_nan(compare, "rgb_cloud_frac", "tsi_frac"),
                corr_or_nan(compare, "acm_cloud_frac", "tsi_frac"),
                compare["rgb_minus_tsi"].mean(),
                compare["acm_minus_tsi"].mean(),
                mae(compare, "rgb_cloud_frac", "tsi_frac"),
                mae(compare, "acm_cloud_frac", "tsi_frac"),
                rmse(compare, "rgb_cloud_frac", "tsi_frac"),
                rmse(compare, "acm_cloud_frac", "tsi_frac"),
            ],
        }
    )
    summary.to_csv(OUT_DIR / "colorado_08c_tsi_rgb_acm_summary.csv", index=False)
    return summary


def save_histograms(compare: pd.DataFrame) -> None:
    valid = compare.loc[
        compare[["tsi_frac", "rgb_cloud_frac", "acm_cloud_frac"]].notna().all(axis=1)
    ].copy()
    bins_cf = np.linspace(0, 1, 41)
    bins_resid = np.linspace(-1, 1, 41)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(valid["tsi_frac"], bins=bins_cf, alpha=0.55, density=True, label="TSI")
    ax.hist(valid["rgb_cloud_frac"], bins=bins_cf, alpha=0.55, density=True, label="RGB domain fraction")
    ax.hist(valid["acm_cloud_frac"], bins=bins_cf, alpha=0.55, density=True, label="GOES ACM fraction")
    ax.set_xlabel("Cloud fraction")
    ax.set_ylabel("Density")
    ax.set_title(f"Cloud-fraction distributions (n={len(valid)})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "tsi_rgb_acm_cloud_fraction_hist.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(valid["rgb_minus_tsi"], bins=bins_resid, alpha=0.6, density=True, label="RGB - TSI")
    ax.hist(valid["acm_minus_tsi"], bins=bins_resid, alpha=0.6, density=True, label="ACM - TSI")
    ax.axvline(0, color="k", linestyle="--", linewidth=1)
    ax.set_xlabel("Residual cloud fraction")
    ax.set_ylabel("Density")
    ax.set_title("Residual distributions")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "rgb_acm_minus_tsi_residual_hist.png", dpi=200)
    plt.close(fig)

    period_defs = [("Nov-Mar", [11, 12, 1, 2, 3]), ("Apr-Jun", [4, 5, 6])]
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True, sharey=True)
    for ax, (label, months) in zip(axes, period_defs):
        period = valid.loc[valid["t"].dt.month.isin(months)]
        ax.hist(period["rgb_minus_tsi"], bins=bins_resid, alpha=0.6, density=True, label="RGB - TSI")
        ax.hist(period["acm_minus_tsi"], bins=bins_resid, alpha=0.6, density=True, label="ACM - TSI")
        ax.axvline(0, color="k", linestyle="--", linewidth=1)
        ax.set_title(f"{label} residuals (n={len(period)})")
        ax.set_ylabel("Density")
        ax.legend()
    axes[-1].set_xlabel("Residual cloud fraction")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "seasonal_rgb_acm_minus_tsi_residual_hist.png", dpi=200)
    plt.close(fig)


def main() -> int:
    compare = build_compare()
    summary = save_summary(compare)
    save_histograms(compare)
    print(summary.to_string(index=False))
    print(f"Outputs: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
