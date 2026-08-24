#!/usr/bin/env python3
"""Build independent RGB-rule and CUES-shortwave binary time series outside Notebook 13."""

from pathlib import Path
import argparse
import pandas as pd

from mammoth_cues_binary_analysis import (
    apply_rgb_rules,
    attach_shortwave,
    load_clean_cues_sw,
    load_era5_cues,
    rebuild_rgb_means,
)

ROOT = Path(__file__).resolve().parent
RGB_DIR = Path("/glade/derecho/scratch/cdalden/mammoth/goes18/rgb_composite")
ERA5_DIR = Path("/glade/derecho/scratch/cdalden/mammoth/era5_land/t2m_hourly")
RULES = ROOT / "tsi_rgb_decision_tree_era5_t2m_10c_rules.csv"
CUES = ROOT / "CUES_1min_data_atmos_radiation_soiltemp_precip_2015to2025.nc"
OUT = ROOT / "output_13_mammoth_domain_cloud_fraction"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-rgb", action="store_true", help="Reread raw cubes even when yearly means exist")
    parser.add_argument("--cloudy-kt", type=float, default=0.55, help="Cloud if k_t is below this threshold")
    parser.add_argument("--clear-kt", type=float, default=0.85, help="Clear if k_t is above this threshold")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    yearly = []
    metadata_rows = []
    for year in (2022, 2023, 2024):
        yearly_path = OUT / f"mammoth_cues_rgb_means_raw_rebuild_{year}.csv"
        if yearly_path.exists() and not args.rebuild_rgb:
            table = pd.read_csv(yearly_path, parse_dates=["time"])
            metadata = {"yearly_rgb_cache": str(yearly_path), "source_file_count": "previous raw rebuild"}
            print(f"Reusing completed raw RGB rebuild: {yearly_path}")
        else:
            table, metadata = rebuild_rgb_means(RGB_DIR, yearly_path, year=year)
        yearly.append(table)
        metadata_rows.append({"year": year, **metadata})

    # All finite RGB scans are retained, including nighttime and partial days.
    rgb = pd.concat(yearly, ignore_index=True).drop_duplicates("time")
    rgb["era5_temp_c"], _ = load_era5_cues(ERA5_DIR, rgb.time)
    missing_era5 = int(rgb.era5_temp_c.isna().sum())
    if missing_era5:
        print(f"Dropping {missing_era5:,} RGB rows without matched ERA5-Land temperature")
        rgb = rgb.dropna(subset=["era5_temp_c"]).reset_index(drop=True)
    rgb_binary = apply_rgb_rules(rgb, RULES)
    rgb_path = OUT / "mammoth_cues_rgb_binary_all_available.csv"
    rgb_binary.to_csv(rgb_path, index=False)

    start = rgb_binary.time.min() - pd.Timedelta(minutes=5)
    end = rgb_binary.time.max() + pd.Timedelta(minutes=5)
    sw, cues_attrs = load_clean_cues_sw(CUES, start=start, end=end)
    sw_on_rgb_clock = attach_shortwave(
        rgb_binary[["time"]], sw, cloudy_kt=args.cloudy_kt, clear_kt=args.clear_kt
    )
    sw_binary = sw_on_rgb_clock.loc[
        sw_on_rgb_clock.sw_cloud_binary.notna(),
        ["time", "sw_obs_5min", "sw_clear", "cos_sza", "k_t", "sw_cloud_binary"],
    ].copy()
    sw_binary["sw_cloud_binary"] = sw_binary.sw_cloud_binary.astype(int)
    sw_path = OUT / "mammoth_cues_sw_binary_all_available.csv"
    sw_binary.to_csv(sw_path, index=False)

    pd.DataFrame(metadata_rows).to_csv(OUT / "mammoth_cues_rgb_rebuild_metadata.csv", index=False)
    print(f"RGB binary: {rgb_path} ({len(rgb_binary):,} rows)")
    print(f"SW binary:  {sw_path} ({len(sw_binary):,} unambiguous rows)")
    print(f"Shortwave rule: cloud if k_t < {args.cloudy_kt}; clear if k_t > {args.clear_kt}")
    print("CUES source time zone:", cues_attrs.get("time_zone"), cues_attrs.get("time_zone_offset"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
    load_era5_cues,
