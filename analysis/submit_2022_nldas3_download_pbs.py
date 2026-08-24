#!/usr/bin/env python
"""Submit one PBS job per month to download 2022 NLDAS-3 Gothic point subsets."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


REPO = Path("/glade/u/home/cdalden/goes_work")
LOG_DIR = REPO / "analysis" / "logs_nldas3_download"
MONTH_RUNNER = REPO / "analysis" / "run_month_nldas3_download.sh"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2022)
    parser.add_argument("--hours", default="14-23,0")
    parser.add_argument("--months", default="1-12", help="Months, e.g. 1-12 or 5,6,7")
    parser.add_argument("--walltime", default="06:00:00")
    parser.add_argument("--select", default="select=1:ncpus=1:mem=8GB")
    parser.add_argument("--queue", default="casper")
    parser.add_argument("--account", default="P48500028")
    parser.add_argument("--force-nldas3", action="store_true", help="Overwrite existing cached NLDAS3 hourly files")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def parse_months(spec: str) -> list[int]:
    months = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            start, end = [int(x) for x in part.split("-", 1)]
            months.extend(range(start, end + 1))
        elif part:
            months.append(int(part))
    return months


def submit_month(args: argparse.Namespace, month: int) -> str:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    job_name = f"nldas3dl{args.year}{month:02d}"
    cmd = [
        "qsub",
        "-q",
        args.queue,
        "-A",
        args.account,
        "-N",
        job_name,
        "-l",
        args.select,
        "-l",
        f"walltime={args.walltime}",
        "-o",
        str(LOG_DIR / f"{job_name}.out"),
        "-e",
        str(LOG_DIR / f"{job_name}.err"),
        "--",
        str(MONTH_RUNNER),
        str(args.year),
        f"{month:02d}",
        args.hours,
        "1" if args.force_nldas3 else "0",
    ]
    if args.dry_run:
        print(" ".join(cmd))
        return "DRYRUN"
    result = subprocess.run(cmd, text=True, capture_output=True, check=True)
    job_id = result.stdout.strip()
    print(f"{month:02d}: {job_id}")
    return job_id


def main() -> None:
    args = parse_args()
    job_ids = [submit_month(args, month) for month in parse_months(args.months)]
    print("Submitted jobs:")
    for job_id in job_ids:
        print(job_id)


if __name__ == "__main__":
    main()
