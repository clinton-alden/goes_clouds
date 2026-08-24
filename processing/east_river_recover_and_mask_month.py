#!/usr/bin/env python3
"""Recover any missing East River RGB days, then build exact-tree masks."""

from __future__ import annotations

import argparse
import calendar
import subprocess
from pathlib import Path

import pandas as pd
import xarray as xr


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("/glade/derecho/scratch/cdalden/east_river_goes")
PYTHON = "/glade/work/cdalden/conda-envs/goes_downloading/bin/python"
BOUNDS = (-107.25, 38.70, -106.75, 39.25)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("year", type=int)
    p.add_argument("month", type=int)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--python", default=PYTHON)
    p.add_argument(
        "--skip-unrecoverable-days",
        action="store_true",
        help="Continue the month when a day still fails after all recovery attempts",
    )
    args = p.parse_args()
    end = calendar.monthrange(args.year, args.month)[1]
    rgb_dir = BASE / "goes16/rgb_composite"
    validator = ROOT / "processing/goes_rgb_domain/validate_rgb_dates.py"
    first = f"{args.year:04d}-{args.month:02d}-01"
    last = f"{args.year:04d}-{args.month:02d}-{end:02d}"
    validate_month = [
        args.python, str(validator), str(rgb_dir), "--start", first, "--end", last,
        "--goes", "goes16", "--domain", "east_river", "--bounds", *map(str, BOUNDS),
    ]

    skipped_days: list[int] = []
    if subprocess.run(validate_month).returncode != 0:
        print(f"RGB recovery required for {args.year}-{args.month:02d}", flush=True)
        for day in range(1, end + 1):
            command = [
                    args.python, str(ROOT / "processing/recover_goes_rgb_day.py"),
                    "--base", str(BASE), "--goes", "goes16", "--domain", "east_river",
                    "--year", str(args.year), "--month", str(args.month), "--day", str(day),
                    "--bounds", *map(str, BOUNDS),
                    "--dem", str(BASE / "static/SRTMGL3_east_river.tif"),
                    "--python", args.python, "--attempts", "3",
                ]
            try:
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError:
                if not args.skip_unrecoverable_days:
                    raise
                skipped_days.append(day)
                print(
                    f"SKIPPING unrecoverable RGB day {args.year}-{args.month:02d}-{day:02d}; "
                    "continuing with the rest of the month",
                    flush=True,
                )
        if not skipped_days:
            subprocess.run(validate_month, check=True)
    else:
        print(f"RGB month already complete: {args.year}-{args.month:02d}", flush=True)

    mask_command = [
        args.python, str(ROOT / "processing/run_tree_mask_month.py"),
        "--year", str(args.year), "--month", str(args.month),
        "--satellite", "goes16", "--domain", "east_river", "--base", str(BASE),
        "--era-dir", str(BASE / "era5_land/t2m_hourly"),
        "--rules", str(ROOT / "analysis/output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv"),
        "--workers", str(args.workers), "--python", args.python,
    ]
    if args.skip_unrecoverable_days:
        mask_command.append("--allow-missing-rgb")
    subprocess.run(mask_command, check=True)

    # ERA5 is disposable for this workflow, but delete it only after all daily
    # products are demonstrably complete and carry the exact-tree provenance.
    mask_dir = BASE / "goes16/cloud_mask_tempbin_tree"
    expected_hours = {0, 1, *range(14, 24)}
    for day in range(1, end + 1):
        token = f"{args.year}{args.month:02d}{day:02d}"
        rgb = rgb_dir / f"goes16_C02_C05_C13_rgb_east_river_{token}.nc"
        if not rgb.is_file() or not rgb.stat().st_size:
            if args.skip_unrecoverable_days:
                continue
            raise RuntimeError(f"Cannot validate masks: missing RGB {rgb}")
        mask = mask_dir / f"goes16_C02_C05_C13_rgb_east_river_{token}_cloud_binary_tempbin_tree.nc"
        if not mask.is_file() or not mask.stat().st_size:
            raise RuntimeError(f"Cannot delete ERA5: missing mask {mask}")
        with xr.open_dataset(mask) as ds:
            times = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
            if len(times) != 144 or set(times.hour) != expected_hours:
                raise RuntimeError(f"Cannot delete ERA5: invalid timestamps in {mask}")
            if ds.attrs.get("tree_logic") != "AND within each leaf; OR across cloudy leaves":
                raise RuntimeError(f"Cannot delete ERA5: wrong rule provenance in {mask}")

    era = BASE / "era5_land/t2m_hourly" / f"era5land_t2m_east_river_{args.year}{args.month:02d}.nc"
    if era.exists():
        era.unlink()
        print(f"DELETED disposable ERA5 month after complete mask validation: {era}", flush=True)
    if skipped_days:
        labels = ", ".join(f"{args.year}-{args.month:02d}-{day:02d}" for day in skipped_days)
        print(f"COMPLETE AVAILABLE DAYS; UNRECOVERABLE DAYS SKIPPED: {labels}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
