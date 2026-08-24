#!/usr/bin/env python3
"""Resume a Colorado month, skip failed days, and record an audit manifest."""

from __future__ import annotations

import argparse
import calendar
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "processing"))
from stream_goes_rgb_mask_day import outputs_are_valid  # noqa: E402

BASE = Path("/glade/derecho/scratch/cdalden/colorado_goes")
BOUNDS = (-109.0, 37.0, -104.0, 41.0)
DOMAIN = "colorado_5x4"
RULES = ROOT / "analysis/output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("year", type=int)
    parser.add_argument("month", type=int)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--satellite", default="goes16")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--day-start", type=int, default=1)
    parser.add_argument("--day-end", type=int)
    args = parser.parse_args()

    if not ((2017, 1) <= (args.year, args.month) <= (2024, 12)):
        raise ValueError("Supported Colorado range is January 2017–December 2024")

    last = calendar.monthrange(args.year, args.month)[1]
    first = max(1, args.day_start)
    end = min(args.day_end or last, last)
    day_script = ROOT / "processing/stream_goes_rgb_mask_day_strict.py"
    rgb_dir = args.base / args.satellite / "rgb_composite_packed"
    mask_dir = args.base / args.satellite / "cloud_mask_tempbin_tree_packed"
    status_dir = args.base / "status" / "monthly"
    status_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, dict[str, str]] = {}
    for day in range(first, end + 1):
        token = f"{args.year}{args.month:02d}{day:02d}"
        date = f"{args.year:04d}-{args.month:02d}-{day:02d}"
        # The public GOES-16 ABI archive has no usable data for this workflow
        # in January-February 2017. Record these dates without repeated requests.
        if (args.year, args.month) < (2017, 3):
            results[token] = {"status": "unavailable", "reason": "pre-archive GOES-16 ABI data"}
            print(f"UNAVAILABLE {token}: pre-archive GOES-16 ABI data", flush=True)
            continue

        command = [
            args.python, str(day_script), "--date", date, "--base", str(args.base),
            "--satellite", args.satellite, "--domain", DOMAIN,
            "--bounds", *map(str, BOUNDS),
            "--dem", str(args.base / "static/SRTMGL3_colorado_5x4.tif"),
            "--era-dir", str(args.base / "era5_land/t2m_hourly"),
            "--rules", str(RULES), "--python", args.python, "--workers", str(args.workers),
        ]
        completed = subprocess.run(command, check=False)
        results[token] = {
            "status": "complete" if completed.returncode == 0 else "failed",
            "reason": "" if completed.returncode == 0 else f"day process exit {completed.returncode}",
        }
        if completed.returncode:
            print(f"SKIP FAILED DAY {token}; continuing month", flush=True)

    # Re-audit products independently of child exit codes.
    valid: list[str] = []
    invalid: list[str] = []
    for day in range(first, end + 1):
        token = f"{args.year}{args.month:02d}{day:02d}"
        rgb = rgb_dir / f"{args.satellite}_C02_C05_C13_rgb_{DOMAIN}_{token}.nc"
        mask = mask_dir / f"{args.satellite}_C02_C05_C13_rgb_{DOMAIN}_{token}_cloud_binary_tempbin_tree.nc"
        if outputs_are_valid(rgb, mask, BOUNDS):
            valid.append(token)
            results[token] = {"status": "complete", "reason": "strict output validation passed"}
        else:
            invalid.append(token)

    manifest = status_dir / f"colorado_{args.year}{args.month:02d}_strict.json"
    temporary = manifest.with_suffix(".json.part")
    temporary.write_text(json.dumps({
        "year": args.year, "month": args.month, "domain": DOMAIN,
        "required_scans_per_day": 144,
        "tree_logic": "AND within each leaf; OR across cloudy leaves",
        "valid_days": valid, "invalid_or_unavailable_days": invalid,
        "days": results,
    }, indent=2) + "\n")
    temporary.replace(manifest)

    if not invalid and first == 1 and end == last:
        era = args.base / "era5_land/t2m_hourly" / f"era5land_t2m_{DOMAIN}_{args.year}{args.month:02d}.nc"
        era.unlink(missing_ok=True)
    print(
        f"MONTH FINISHED Colorado {args.year}-{args.month:02d}: "
        f"strict_valid={len(valid)}/{end-first+1}, skipped={len(invalid)}, manifest={manifest}",
        flush=True,
    )
    # A month with gaps is a successful audited pass; the manifest carries gaps.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
