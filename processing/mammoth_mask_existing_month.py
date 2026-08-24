#!/usr/bin/env python3
"""Mask every available Mammoth RGB day in one month with the exact tree."""

from __future__ import annotations

import argparse
import calendar
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("/glade/derecho/scratch/cdalden/mammoth")
PYTHON = "/glade/work/cdalden/conda-envs/goes_downloading/bin/python"
DOMAIN = "mammoth_1deg"
ERA_DIR = BASE / "era5_land/t2m_hourly_mammoth_1deg"
RULES = ROOT / "analysis/output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv"


def segments(year: int, month: int) -> list[tuple[str, int, int]]:
    end = calendar.monthrange(year, month)[1]
    if (year, month) < (2023, 1):
        return [("goes17", 1, end)]
    if (year, month) == (2023, 1):
        return [("goes17", 1, 2), ("goes18", 3, end)]
    return [("goes18", 1, end)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("year", type=int)
    parser.add_argument("month", type=int)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--python", default=PYTHON)
    args = parser.parse_args()

    if not ((2019, 3) <= (args.year, args.month) <= (2025, 12)):
        raise ValueError("Mammoth range is March 2019 through December 2025")

    total_rgb = 0
    total_masks = 0
    for satellite, first, last in segments(args.year, args.month):
        rgb_dir = BASE / satellite / "rgb_composite"
        mask_dir = BASE / satellite / "cloud_mask_tempbin_tree"
        expected_rgb = [
            rgb_dir / (
                f"{satellite}_C02_C05_C13_rgb_{DOMAIN}_"
                f"{args.year}{args.month:02d}{day:02d}.nc"
            )
            for day in range(first, last + 1)
        ]
        available = [path for path in expected_rgb if path.is_file() and path.stat().st_size]
        if not available:
            print(f"NO RGB {satellite} {args.year}-{args.month:02d} days {first}-{last}", flush=True)
            continue

        command = [
            args.python,
            str(ROOT / "processing/run_tree_mask_month.py"),
            "--year", str(args.year),
            "--month", str(args.month),
            "--day-start", str(first),
            "--day-end", str(last),
            "--satellite", satellite,
            "--domain", DOMAIN,
            "--base", str(BASE),
            "--era-dir", str(ERA_DIR),
            "--rules", str(RULES),
            "--workers", str(args.workers),
            "--python", args.python,
            "--allow-missing-rgb",
        ]
        subprocess.run(command, check=True)
        total_rgb += len(available)
        total_masks += sum(
            (mask_dir / f"{path.stem}_cloud_binary_tempbin_tree.nc").is_file()
            for path in available
        )

    if total_rgb == 0:
        raise RuntimeError(f"No Mammoth RGB files found for {args.year}-{args.month:02d}")
    if total_masks != total_rgb:
        raise RuntimeError(f"Only {total_masks}/{total_rgb} available RGB days were masked")

    expected_days = calendar.monthrange(args.year, args.month)[1]
    if total_masks == expected_days:
        era = ERA_DIR / f"era5land_t2m_{DOMAIN}_{args.year}{args.month:02d}.nc"
        if era.exists():
            era.unlink()
            print(f"DELETED disposable ERA5 month: {era}", flush=True)

    print(
        f"COMPLETE AVAILABLE Mammoth {args.year}-{args.month:02d}: "
        f"{total_masks}/{expected_days} days masked",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
