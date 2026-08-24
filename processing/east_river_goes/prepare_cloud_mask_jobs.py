#!/usr/bin/env python3
"""Preview or submit the Apr 2017-Dec 2024 East River cloud-mask array."""

import argparse
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--concurrency", type=int, default=12)
    args = parser.parse_args()
    if not 1 <= args.concurrency <= 12:
        parser.error("--concurrency must be between 1 and 12")
    pbs = Path(__file__).with_name("cloud_mask_month.pbs").resolve()
    command = ["qsub", "-N", "east_river_mask", "-J", f"1-93%{args.concurrency}", str(pbs)]
    print("Period: 2017-04 through 2024-12 (93 months)")
    print("Hours: 00, 01, 14-23 UTC")
    print(f"Concurrency: {args.concurrency}")
    print("Command:", " ".join(command))
    if not args.submit:
        print("PREVIEW ONLY: no jobs submitted")
        return 0
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    print("Submitted:", result.stdout.strip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
