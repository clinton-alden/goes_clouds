#!/usr/bin/env python3
"""Run one Colorado day transactionally and require 144 aligned ABI scans."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "processing"))
from stream_goes_rgb_mask_day import channel_files, outputs_are_valid  # noqa: E402


def scan_keys(paths: list[Path]) -> set[str]:
    """Return nominal YYYYJJJHHMM scan starts, shared by all ABI bands."""
    keys: set[str] = set()
    for path in paths:
        match = re.search(r"_s(\d{11})", path.name)
        if not match:
            raise RuntimeError(f"Cannot parse ABI scan start from {path.name}")
        if match.group(1) in keys:
            raise RuntimeError(f"Duplicate nominal ABI scan time {match.group(1)}")
        keys.add(match.group(1))
    return keys


def raw_day_is_complete(raw_day: Path) -> tuple[bool, str]:
    keys = {channel: scan_keys(channel_files(raw_day, channel)) for channel in ("C02", "C05", "C13")}
    counts = {channel: len(value) for channel, value in keys.items()}
    aligned = keys["C02"] == keys["C05"] == keys["C13"]
    return aligned and counts == {"C02": 144, "C05": 144, "C13": 144}, f"counts={counts}, aligned={aligned}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--satellite", default="goes16")
    parser.add_argument("--domain", default="colorado_5x4")
    parser.add_argument("--bounds", type=float, nargs=4, default=(-109, 37, -104, 41))
    parser.add_argument("--dem", type=Path, required=True)
    parser.add_argument("--era-dir", type=Path, required=True)
    parser.add_argument("--rules", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    date = pd.Timestamp(args.date)
    token = date.strftime("%Y%m%d")
    raw_day = args.base / args.satellite / str(date.year) / str(date.month) / str(date.day)
    rgb = args.base / args.satellite / "rgb_composite_packed" / (
        f"{args.satellite}_C02_C05_C13_rgb_{args.domain}_{token}.nc"
    )
    mask = args.base / args.satellite / "cloud_mask_tempbin_tree_packed" / (
        f"{args.satellite}_C02_C05_C13_rgb_{args.domain}_{token}_cloud_binary_tempbin_tree.nc"
    )
    temp_rgb = args.base / "tmp" / f"{args.satellite}_C02_C05_C13_rgb_{args.domain}_{token}.nc"

    if outputs_are_valid(rgb, mask, args.bounds):
        print(f"STRICT VALID existing day: {token}", flush=True)
        return 0

    # Never let an earlier partial product or float-RGB cache be reused.
    rgb.unlink(missing_ok=True)
    mask.unlink(missing_ok=True)
    temp_rgb.unlink(missing_ok=True)

    command = [
        args.python, str(ROOT / "processing/stream_goes_rgb_mask_day.py"),
        "--date", args.date, "--base", str(args.base),
        "--satellite", args.satellite, "--domain", args.domain,
        "--bounds", *map(str, args.bounds), "--dem", str(args.dem),
        "--era-dir", str(args.era_dir), "--rules", str(args.rules),
        "--python", args.python, "--workers", str(args.workers), "--keep-work",
    ]
    try:
        subprocess.run(command, check=True)
        raw_valid, detail = raw_day_is_complete(raw_day)
        if not raw_valid:
            raise RuntimeError(f"Day does not contain 144 aligned scans: {detail}")
        if not outputs_are_valid(rgb, mask, args.bounds):
            raise RuntimeError("Packed RGB/mask failed strict size, orientation, packing, or tree validation")
    except Exception:
        # Derived partial days are worse than an explicit gap. Preserve raw data
        # for a later retry, but remove products that could be mistaken as valid.
        rgb.unlink(missing_ok=True)
        mask.unlink(missing_ok=True)
        temp_rgb.unlink(missing_ok=True)
        raise

    shutil.rmtree(raw_day, ignore_errors=True)
    temp_rgb.unlink(missing_ok=True)
    print(f"STRICT COMPLETE rgb={rgb} mask={mask}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
