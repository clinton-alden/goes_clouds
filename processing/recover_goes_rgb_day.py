#!/usr/bin/env python3
"""Rebuild one GOES RGB day with targeted cleanup and bounded retries."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], env: dict[str, str]) -> None:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=env)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--goes", required=True)
    p.add_argument("--domain", required=True)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, required=True)
    p.add_argument("--day", type=int, required=True)
    p.add_argument("--bounds", type=float, nargs=4, required=True)
    p.add_argument("--dem", type=Path, required=True)
    p.add_argument("--python", default="python")
    p.add_argument("--attempts", type=int, default=3)
    args = p.parse_args()

    token = f"{args.year}{args.month:02d}{args.day:02d}"
    rgb = args.base / args.goes / "rgb_composite" / f"{args.goes}_C02_C05_C13_rgb_{args.domain}_{token}.nc"
    validator = ROOT / "processing/goes_rgb_domain/validate_rgb_dates.py"
    start = f"{args.year:04d}-{args.month:02d}-{args.day:02d}"
    validate = [
        args.python, str(validator), str(rgb.parent), "--start", start, "--end", start,
        "--goes", args.goes, "--domain", args.domain, "--bounds", *map(str, args.bounds),
    ]
    if subprocess.run(validate).returncode == 0:
        print(f"VALID existing {token}", flush=True)
        return 0

    env = os.environ.copy()
    env.update({
        "GOES_HOURS": "00,01,14,15,16,17,18,19,20,21,22,23",
        "LON_MIN": str(args.bounds[0]), "LAT_MIN": str(args.bounds[1]),
        "LON_MAX": str(args.bounds[2]), "LAT_MAX": str(args.bounds[3]),
        "SHARED_DEM_PATH": str(args.dem), "HDF5_USE_FILE_LOCKING": "FALSE",
        "ORTHO_MAX_WORKERS": "8", "RGB_MAX_WORKERS": "1",
        "PYTHONPATH": str(ROOT / "processing/gothic") + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""),
    })
    raw_day = args.base / args.goes / str(args.year) / str(args.month) / str(args.day)

    for attempt in range(1, args.attempts + 1):
        print(f"RECOVERY {token}: attempt {attempt}/{args.attempts}", flush=True)
        if raw_day.exists():
            shutil.rmtree(raw_day)
        for channel in ("C02", "C05", "C13"):
            zarr = args.base / args.goes / channel / f"{args.goes}_{channel}_{args.domain}_{token}.zarr"
            if zarr.exists():
                shutil.rmtree(zarr)
        rgb.unlink(missing_ok=True)
        try:
            for channel in ("C02", "C05", "C13"):
                run([
                    args.python, str(ROOT / "data_download/download-goes.py"),
                    "-B", f"noaa-{args.goes}", "-Y", str(args.year), "-M", str(args.month),
                    "-D", str(args.day), str(args.day), "-p", "ABI-L1b-RadC", "-c", channel,
                    "-b", *map(str, args.bounds), "-d", str(args.base),
                ], env)
            run([args.python, str(ROOT / "processing/east_river_goes/batch_ortho.py"), str(raw_day)], env)
            run([
                args.python, str(ROOT / "processing/gothic/zarr_v2.py"), f"{args.base}/",
                str(args.year), str(args.month), args.goes, args.domain,
                str(args.day), str(args.day),
            ], env)
            run([
                args.python, str(ROOT / "processing/gothic/rgb_v2.py"), str(args.base / args.goes),
                str(args.year), str(args.month), args.domain, args.goes,
                str(args.day), str(args.day),
            ], env)
            run(validate, env)
            print(f"RECOVERED {token}", flush=True)
            return 0
        except subprocess.CalledProcessError as exc:
            print(f"FAILED {token} attempt {attempt}: {exc}", flush=True)
    raise RuntimeError(f"Could not recover {token} after {args.attempts} clean attempts")


if __name__ == "__main__":
    raise SystemExit(main())
