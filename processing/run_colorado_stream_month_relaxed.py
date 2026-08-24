#!/usr/bin/env python3
"""Resume a Colorado month using common timestamps and >=132 scans/day."""

from __future__ import annotations

import argparse
import calendar
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "processing"))
from stream_goes_rgb_mask_day_relaxed import relaxed_output_count  # noqa: E402

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
    parser.add_argument("--minimum-scans", type=int, default=132)
    parser.add_argument("--day-start", type=int, default=1)
    parser.add_argument("--day-end", type=int)
    args = parser.parse_args()

    if not ((2017, 1) <= (args.year, args.month) <= (2024, 12)):
        raise ValueError("Supported Colorado range is January 2017–December 2024")
    if not 1 <= args.minimum_scans <= 144:
        raise ValueError("minimum-scans must be between 1 and 144")

    last = calendar.monthrange(args.year, args.month)[1]
    first = max(1, args.day_start)
    end = min(args.day_end or last, last)
    day_script = ROOT / "processing/stream_goes_rgb_mask_day_relaxed.py"
    rgb_dir = args.base / args.satellite / "rgb_composite_packed"
    mask_dir = args.base / args.satellite / "cloud_mask_tempbin_tree_packed"
    status_dir = args.base / "status" / "monthly_relaxed"
    status_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, dict[str, object]] = {}
    for day in range(first, end + 1):
        token = f"{args.year}{args.month:02d}{day:02d}"
        if (args.year, args.month) < (2017, 3):
            results[token] = {"status": "unavailable", "valid_scans": 0}
            continue
        command = [
            args.python, str(day_script), "--date", f"{args.year:04d}-{args.month:02d}-{day:02d}",
            "--base", str(args.base), "--satellite", args.satellite, "--domain", DOMAIN,
            "--bounds", *map(str, BOUNDS),
            "--dem", str(args.base / "static/SRTMGL3_colorado_5x4.tif"),
            "--era-dir", str(args.base / "era5_land/t2m_hourly"),
            "--rules", str(RULES), "--python", args.python, "--workers", str(args.workers),
            "--minimum-scans", str(args.minimum_scans),
        ]
        completed = subprocess.run(command, check=False)
        rgb = rgb_dir / f"{args.satellite}_C02_C05_C13_rgb_{DOMAIN}_{token}.nc"
        mask = mask_dir / f"{args.satellite}_C02_C05_C13_rgb_{DOMAIN}_{token}_cloud_binary_tempbin_tree.nc"
        n = relaxed_output_count(rgb, mask, list(BOUNDS), args.minimum_scans)
        results[token] = {
            "status": "complete" if completed.returncode == 0 and n else "failed",
            "valid_scans": n,
        }
        if not n:
            print(f"SKIP DAY {token}: fewer than {args.minimum_scans} usable aligned scans", flush=True)

    valid = [token for token, result in results.items() if result["status"] == "complete"]
    skipped = [token for token, result in results.items() if result["status"] != "complete"]
    manifest = status_dir / f"colorado_{args.year}{args.month:02d}_min{args.minimum_scans}.json"
    temporary = manifest.with_suffix(".json.part")
    temporary.write_text(json.dumps({
        "year": args.year, "month": args.month, "domain": DOMAIN,
        "minimum_common_aligned_scans": args.minimum_scans,
        "nominal_scans_per_day": 144,
        "tree_logic": "AND within each leaf; OR across cloudy leaves",
        "valid_days": valid, "skipped_days": skipped, "days": results,
    }, indent=2) + "\n")
    temporary.replace(manifest)

    if not skipped and first == 1 and end == last:
        era = args.base / "era5_land/t2m_hourly" / f"era5land_t2m_{DOMAIN}_{args.year}{args.month:02d}.nc"
        era.unlink(missing_ok=True)
    print(
        f"MONTH FINISHED Colorado {args.year}-{args.month:02d}: "
        f">={args.minimum_scans} scans={len(valid)}/{end-first+1}, skipped={len(skipped)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
