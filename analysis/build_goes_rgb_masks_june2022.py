#!/usr/bin/env python3
"""Regenerate June 2022 Colorado RGB cloud masks without rendering GIFs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COLORADO_PROCESSING = ROOT / "processing" / "colorado"
sys.path.insert(0, str(COLORADO_PROCESSING))

from apply_tempbin_thresholds_colorado import build_cloud_mask, infer_output_paths  # noqa: E402


DEFAULT_HOURS = {0, 1, *range(14, 24)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2022)
    parser.add_argument("--month", type=int, default=6)
    parser.add_argument("--day", type=int, help="Optional single day of month")
    parser.add_argument(
        "--rgb-dir",
        type=Path,
        default=Path("/glade/u/home/cdalden/scratch/colorado/goes16/rgb_composite"),
    )
    parser.add_argument(
        "--era5-file",
        type=Path,
        default=Path(
            "/glade/derecho/scratch/cdalden/tmp/colorado/era5_land/t2m_hourly/"
            "era5land_t2m_colorado_202206.nc"
        ),
    )
    parser.add_argument(
        "--threshold-csv",
        type=Path,
        default=ROOT / "analysis/output_12_rgb_threshold_transfer/gothic_temp_bin_rgb_thresholds_10c.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/glade/derecho/scratch/cdalden/hrrr_goes_june2022/goes_masks"),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.era5_file.exists():
        raise FileNotFoundError(args.era5_file)
    if not args.threshold_csv.exists():
        raise FileNotFoundError(args.threshold_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rgb_files = sorted(args.rgb_dir.glob(f"goes16_C02_C05_C13_rgb_colorado_{args.year}{args.month:02d}*.nc"))
    if args.day is not None:
        wanted = f"{args.year}{args.month:02d}{args.day:02d}"
        rgb_files = [path for path in rgb_files if wanted in path.name]
        if len(rgb_files) != 1:
            raise RuntimeError(f"Expected one RGB file for {wanted}, found {len(rgb_files)}")
    elif len(rgb_files) != 30:
        raise RuntimeError(f"Expected 30 June RGB files, found {len(rgb_files)}")

    for index, rgb_file in enumerate(rgb_files, start=1):
        mask_path, _ = infer_output_paths(rgb_file, args.output_dir, args.output_dir)
        if mask_path.exists() and not args.overwrite:
            print(f"[{index:02d}/{len(rgb_files):02d} skip] {mask_path}", flush=True)
            continue
        print(f"[{index:02d}/{len(rgb_files):02d} build] {rgb_file.name}", flush=True)
        build_cloud_mask(
            rgb_path=rgb_file,
            era5_path=args.era5_file,
            threshold_csv=args.threshold_csv,
            mask_path=mask_path,
            target_hours=DEFAULT_HOURS,
        )
        print(f"[{index:02d}/{len(rgb_files):02d} wrote] {mask_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
