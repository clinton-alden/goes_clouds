#!/usr/bin/env python3
"""Create monthly exact-tree mask manifests for East River and Mammoth."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
HEADER = ["domain", "satellite", "year", "month", "day_start", "day_end", "base", "era_dir"]


def add_period(rows, domain, satellite, start, end, base, era):
    for month in pd.period_range(start, end, freq="M"):
        rows.append([domain, satellite, month.year, month.month, 1, month.days_in_month, base, era])


def write(path: Path, rows):
    with path.open("w", newline="") as stream:
        csv.writer(stream).writerows([HEADER, *rows])
    print(f"Wrote {len(rows)} tasks: {path}")


def ready_rows(rows):
    ready = []
    for row in rows:
        domain, satellite, year, month, day_start, day_end, base, _ = row
        rgb_dir = Path(base) / satellite / "rgb_composite"
        paths = [
            rgb_dir / f"{satellite}_C02_C05_C13_rgb_{domain}_{year}{month:02d}{day:02d}.nc"
            for day in range(day_start, day_end + 1)
        ]
        if all(path.is_file() and path.stat().st_size > 0 for path in paths):
            ready.append(row)
    return ready


def main() -> int:
    east_base = "/glade/derecho/scratch/cdalden/east_river_goes"
    east_era = f"{east_base}/era5_land/t2m_hourly"
    east = []
    add_period(east, "east_river", "goes16", "2017-04", "2024-12", east_base, east_era)
    write(ROOT / "processing/east_river_tree_mask_months.csv", east)
    write(ROOT / "processing/east_river_tree_mask_ready.csv", ready_rows(east))

    mammoth_base = "/glade/derecho/scratch/cdalden/mammoth"
    mammoth_era = f"{mammoth_base}/era5_land/t2m_hourly_mammoth_1deg"
    mammoth = []
    add_period(mammoth, "mammoth_1deg", "goes17", "2019-03", "2022-12", mammoth_base, mammoth_era)
    mammoth.append(["mammoth_1deg", "goes17", 2023, 1, 1, 2, mammoth_base, mammoth_era])
    mammoth.append(["mammoth_1deg", "goes18", 2023, 1, 3, 31, mammoth_base, mammoth_era])
    add_period(mammoth, "mammoth_1deg", "goes18", "2023-02", "2025-12", mammoth_base, mammoth_era)
    write(ROOT / "processing/mammoth_tree_mask_months.csv", mammoth)
    write(ROOT / "processing/mammoth_tree_mask_ready.csv", ready_rows(mammoth))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
