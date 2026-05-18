#!/usr/bin/env python3
"""
RGB Processing Script v2
Source: https://github.com/clinton-alden/goes_clouds/blob/master/processing/rgb_v2.py

Processes GOES Zarr files to create RGB composites for a full month.
Automatically cleans up Zarr files after processing.

Usage:
    python rgb_v2.py <in_dir> <year> <month> [domain] [goes]

Example:
    python rgb_v2.py /storage/cdalden/goes/colorado/ 2022 7 colorado goes16
"""

import sys
import os
import shutil
import calendar
import utils


def main():
    if len(sys.argv) < 4:
        print("Usage: python rgb_v2.py <in_dir> <year> <month> [domain] [goes]")
        sys.exit(1)

    in_dir = sys.argv[1].rstrip("/") + "/"
    year = int(sys.argv[2])
    month = int(sys.argv[3])
    domain = sys.argv[4] if len(sys.argv) > 4 else "scripps"
    goes = sys.argv[5] if len(sys.argv) > 5 else "goes18"

    # Auto-compute end of month
    start_day = 1
    end_day = calendar.monthrange(year, month)[1]

    print(f"Processing {year}-{month:02d} ({start_day} → {end_day})")
    print(f"Input dir: {in_dir}")
    print(f"Domain: {domain}, GOES: {goes}")

    channel_list = ["C02", "C05", "C13"]

    for day in range(start_day, end_day + 1):
        date = f"{year}{month:02d}{day:02d}"
        print(f"\n=== Starting {date} ===")

        try:
            # Convert raw → RGB
            utils.goes_rad_to_rgb(in_dir, date, goes, location=domain)
            print(f"Finished RGB for {date}")

        except Exception as e:
            print(f'RGB for {date} failed with error: {e}')

        # Cleanup Zarr files
        for channel in channel_list:
            zarr_path = f"{in_dir}{channel}/{goes}_{channel}_{domain}_{date}.zarr"

            if os.path.exists(zarr_path):
                try:
                    shutil.rmtree(zarr_path)
                    print(f"Deleted: {zarr_path}")
                except Exception as e:
                    print(f"Failed to delete {zarr_path}: {e}")
            else:
                print(f"Not found: {zarr_path}")


if __name__ == "__main__":
    main()