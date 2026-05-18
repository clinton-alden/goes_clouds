#!/usr/bin/env python3
"""
Create GOES Zarr outputs for an explicit day range.

Usage:
    python zarr_v2_days.py <directory> <year> <month> <start_day> <end_day> <goes_model> <domain>
"""

import os
import sys

import utils


def main():
    if len(sys.argv) != 8:
        print(
            "Usage: python zarr_v2_days.py <directory> <year> <month> <start_day> <end_day> <goes_model> <domain>"
        )
        sys.exit(1)

    in_dir = sys.argv[1]
    year = int(sys.argv[2])
    month = int(sys.argv[3])
    start_day = int(sys.argv[4])
    end_day = int(sys.argv[5])
    goes_model = sys.argv[6]
    domain = sys.argv[7]

    channel_list = ["C02", "C05", "C13"]

    print(
        f"Processing {goes_model} for {year}-{month:02d} (days {start_day}-{end_day})"
    )

    utils.goes_nc_to_zarr(
        in_dir,
        channel_list,
        start_day,
        end_day,
        month,
        year,
        domain,
        goes_model,
        surprise=True,
    )

    for day in range(start_day, end_day + 1):
        nc_dir = os.path.join(in_dir, goes_model, str(year), str(month), str(day))

        for root, _, files in os.walk(nc_dir):
            print(f"Inspecting directory: {root}")
            print(f"Files found: {files}")
            for file in files:
                if file.endswith(".nc"):
                    file_path = os.path.join(root, file)
                    try:
                        os.remove(file_path)
                        print(f"Deleted: {file_path}")
                    except Exception as exc:
                        print(f"Failed to delete {file_path}: {exc}")
                        raise


if __name__ == "__main__":
    main()
