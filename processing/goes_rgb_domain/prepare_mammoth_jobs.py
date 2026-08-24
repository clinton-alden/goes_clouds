#!/usr/bin/env python3
"""Create and preview the Mammoth GOES-17/18 task manifest; optionally submit."""

import argparse
import calendar
import csv
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "mammoth_201903_202512_tasks.csv"


def tasks():
    rows = []
    for absolute in range(2019 * 12 + 2, 2023 * 12):
        year, month0 = divmod(absolute, 12)
        month = month0 + 1
        rows.append(("goes17", year, month, 1, calendar.monthrange(year, month)[1]))
    rows.extend((("goes17", 2023, 1, 1, 2), ("goes18", 2023, 1, 3, 31)))
    for absolute in range(2023 * 12 + 1, 2026 * 12):
        year, month0 = divmod(absolute, 12)
        month = month0 + 1
        rows.append(("goes18", year, month, 1, calendar.monthrange(year, month)[1]))
    return rows


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--submit", action="store_true")
    p.add_argument("--concurrency", type=int, default=12)
    a = p.parse_args()
    rows = tasks()
    with MANIFEST.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("goes", "year", "month", "day_start", "day_end"))
        w.writerows(rows)
    command = ["qsub", "-N", "mammoth_rgb", "-J", f"1-{len(rows)}%{a.concurrency}",
               "-v", f"MANIFEST={MANIFEST}", str(HERE / "task_array.pbs")]
    print(f"Wrote {len(rows)} tasks: {MANIFEST}")
    print("Bounds: -119.53390,37.13740,-118.53390,38.13740")
    print("GOES-17: 2019-03-01 through 2023-01-02")
    print("GOES-18: 2023-01-03 through 2025-12-31")
    print("Hours: 00,01,14-23 UTC")
    print("Command:", " ".join(command))
    if not a.submit:
        print("PREVIEW ONLY: no jobs submitted")
        return 0
    print("Submitted:", subprocess.run(command, check=True, text=True, capture_output=True).stdout.strip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
