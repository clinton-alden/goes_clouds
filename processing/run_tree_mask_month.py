#!/usr/bin/env python3
"""Build exact-tree masks for one (possibly partial) satellite month."""

from __future__ import annotations

import argparse
import calendar
import fcntl
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import xarray as xr


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, required=True)
    p.add_argument("--day-start", type=int, default=1)
    p.add_argument("--day-end", type=int)
    p.add_argument("--satellite", required=True)
    p.add_argument("--domain", required=True)
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--era-dir", type=Path, required=True)
    p.add_argument("--rules", type=Path, required=True)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--python", default="python")
    p.add_argument(
        "--allow-missing-rgb",
        action="store_true",
        help="Mask every available RGB day instead of skipping the entire month",
    )
    args = p.parse_args()

    last = calendar.monthrange(args.year, args.month)[1]
    end = min(args.day_end or last, last)
    days = range(args.day_start, end + 1)
    rgb_dir = args.base / args.satellite / "rgb_composite"
    out_dir = args.base / args.satellite / "cloud_mask_tempbin_tree"
    expected = [
        rgb_dir / f"{args.satellite}_C02_C05_C13_rgb_{args.domain}_{args.year}{args.month:02d}{day:02d}.nc"
        for day in days
    ]
    missing = [path for path in expected if not path.is_file() or path.stat().st_size == 0]
    if missing:
        if args.allow_missing_rgb:
            print(
                f"MISSING RGB: processing {len(expected) - len(missing)}/{len(expected)} "
                f"available days for {args.satellite} {args.year}-{args.month:02d}",
                flush=True,
            )
            expected = [path for path in expected if path not in missing]
            if not expected:
                raise RuntimeError("No RGB files are available to mask")
        else:
            print(
                f"SKIPPED {args.satellite} {args.year}-{args.month:02d} days "
                f"{args.day_start}-{end}: {len(missing)}/{len(expected)} RGB files missing",
                flush=True,
            )
            return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name("build_rgb_tempbin_tree_mask.py")

    def output_for(rgb: Path) -> Path:
        return out_dir / f"{rgb.stem}_cloud_binary_tempbin_tree.nc"

    def valid_output(rgb: Path) -> bool:
        output = output_for(rgb)
        if not output.is_file() or not output.stat().st_size:
            return False
        try:
            with xr.open_dataset(output) as ds:
                times = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
                return (
                    len(times) == 144
                    and set(times.hour) == {0, 1, *range(14, 24)}
                    and ds.attrs.get("tree_logic") == "AND within each leaf; OR across cloudy leaves"
                )
        except Exception:
            return False

    def command(rgb: Path) -> list[str]:
        return [
            args.python, str(script), str(rgb), str(output_for(rgb)),
            "--era-dir", str(args.era_dir), "--rules", str(args.rules),
        ]

    # Serialize the first build across array tasks sharing an ERA month. This
    # prevents duplicate CDS requests for split satellite months such as Jan 2023.
    args.era_dir.mkdir(parents=True, exist_ok=True)
    lock_path = args.era_dir / f".{args.year}{args.month:02d}.download.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        first = expected[0]
        if not valid_output(first):
            subprocess.run(command(first), check=True)

    remaining = [rgb for rgb in expected[1:] if not valid_output(rgb)]
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(subprocess.run, command(rgb), check=True): rgb for rgb in remaining}
        for future in as_completed(futures):
            future.result()

    present = sum(valid_output(rgb) for rgb in expected)
    if present != len(expected):
        raise RuntimeError(f"Only {present}/{len(expected)} masks exist after processing")
    print(
        f"COMPLETE {args.satellite} {args.year}-{args.month:02d} days "
        f"{args.day_start}-{end}: {present} exact-tree masks",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
