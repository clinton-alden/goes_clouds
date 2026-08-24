#!/usr/bin/env python3
"""Process one Colorado day using at least 132 common aligned ABI scans."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "processing"))
def relaxed_output_count(rgb: Path, mask: Path, bounds: list[float], minimum: int) -> int:
    if not rgb.is_file() or not mask.is_file() or not rgb.stat().st_size or not mask.stat().st_size:
        return 0
    try:
        with xr.open_dataset(rgb) as r, xr.open_dataset(mask) as m:
            n = int(r.sizes.get("t", 0))
            alignment_ok = n == 144 or r.attrs.get("scan_alignment") == "nominal_start_intersection"
            valid = (
                minimum <= n <= 144
                and alignment_ok
                and dict(m.sizes) == dict(r.sizes)
                and r.sizes.get("latitude") == 961
                and r.sizes.get("longitude") == 1601
                and np.array_equal(r.t.values, m.t.values)
                and np.all(np.diff(r.latitude) > 0)
                and np.all(np.diff(r.longitude) > 0)
                and abs(float(r.longitude.min()) - bounds[0]) < 0.01
                and abs(float(r.longitude.max()) - bounds[2]) < 0.01
                and abs(float(r.latitude.min()) - bounds[1]) < 0.01
                and abs(float(r.latitude.max()) - bounds[3]) < 0.01
                and r.red.encoding.get("dtype") == np.dtype("uint8")
                and m.attrs.get("tree_logic") == "AND within each leaf; OR across cloudy leaves"
            )
            return n if valid else 0
    except Exception:
        return 0


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
    parser.add_argument("--minimum-scans", type=int, default=132)
    args = parser.parse_args()

    date = pd.Timestamp(args.date)
    token = date.strftime("%Y%m%d")
    raw_day = args.base / args.satellite / str(date.year) / str(date.month) / str(date.day)
    rgb = args.base / args.satellite / "rgb_composite_packed" / f"{args.satellite}_C02_C05_C13_rgb_{args.domain}_{token}.nc"
    mask = args.base / args.satellite / "cloud_mask_tempbin_tree_packed" / f"{args.satellite}_C02_C05_C13_rgb_{args.domain}_{token}_cloud_binary_tempbin_tree.nc"
    temp_rgb = args.base / "tmp" / f"{args.satellite}_C02_C05_C13_rgb_{args.domain}_{token}.nc"

    existing = relaxed_output_count(rgb, mask, args.bounds, args.minimum_scans)
    if existing:
        print(f"RELAXED VALID existing day {token}: {existing}/144 aligned scans", flush=True)
        return 0

    rgb.unlink(missing_ok=True)
    mask.unlink(missing_ok=True)
    temp_rgb.unlink(missing_ok=True)
    env = os.environ.copy()
    env["MIN_ALIGNED_SCANS"] = str(args.minimum_scans)
    command = [
        args.python, str(ROOT / "processing/stream_goes_rgb_mask_day.py"),
        "--date", args.date, "--base", str(args.base), "--satellite", args.satellite,
        "--domain", args.domain, "--bounds", *map(str, args.bounds),
        "--dem", str(args.dem), "--era-dir", str(args.era_dir), "--rules", str(args.rules),
        "--python", args.python, "--workers", str(args.workers), "--keep-work",
    ]
    try:
        subprocess.run(command, check=True, env=env)
        n = relaxed_output_count(rgb, mask, args.bounds, args.minimum_scans)
        if not n:
            raise RuntimeError("RGB/mask failed relaxed aligned-scan validation")
    except Exception:
        rgb.unlink(missing_ok=True)
        mask.unlink(missing_ok=True)
        temp_rgb.unlink(missing_ok=True)
        raise


    shutil.rmtree(raw_day, ignore_errors=True)
    temp_rgb.unlink(missing_ok=True)
    print(f"RELAXED COMPLETE {token}: {n}/144 aligned scans", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
