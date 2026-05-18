#!/usr/bin/env python
"""Write daily ACM zarr files from orthorectified NetCDF inputs."""

import argparse
import os

import utils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert orthorectified ACM NetCDF files to daily zarr files."
    )
    parser.add_argument("--base-dir", required=True, help="GOES base directory.")
    parser.add_argument("--year", required=True, type=int)
    parser.add_argument("--month", required=True, type=int)
    parser.add_argument("--start-day", required=True, type=int)
    parser.add_argument("--end-day", required=True, type=int)
    parser.add_argument("--goes", required=True, help="GOES model (e.g., goes16).")
    parser.add_argument("--domain", required=True, help="Domain/site token for output names.")
    parser.add_argument(
        "--channel",
        default="ACMC",
        help="Channel token expected in ACM filenames and zarr output names.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.start_day > args.end_day:
        raise ValueError(
            f"start-day ({args.start_day}) must be <= end-day ({args.end_day})"
        )

    base_dir = args.base_dir.rstrip("/") + "/"
    os.makedirs(os.path.join(base_dir, args.goes, args.channel), exist_ok=True)

    print(
        f"Writing zarr for {args.domain} {args.goes} {args.year}-{args.month:02d} "
        f"days {args.start_day:02d}-{args.end_day:02d} channel={args.channel}"
    )
    utils.goes_nc_to_zarr(
        base_dir,
        [args.channel],
        args.start_day,
        args.end_day,
        args.month,
        args.year,
        args.domain,
        args.goes,
        surprise=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
