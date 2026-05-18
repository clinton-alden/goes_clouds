#!/usr/bin/env python3
import sys
import os
import shutil
import calendar
import utils
from multiprocessing import Pool, cpu_count

# Mammoth Mountain test version

# Module-level worker function (must be at top level for pickling)
def process_day(args):
    """Process a single day - worker function for multiprocessing"""
    day, year, month, in_dir, goes, domain = args
    date = f"{year}{month:02d}{day:02d}"
    try:
        utils.goes_rad_to_rgb(in_dir, date, goes, location=domain)
        return f"SUCCESS: {date}"
    except Exception as e:
        return f"ERROR {date}: {e}"

def main():
    if len(sys.argv) < 4:
        print("Usage: python rgb_v2.py <in_dir> <year> <month> [domain] [goes]")
        sys.exit(1)

    in_dir = sys.argv[1].rstrip("/") + "/"
    year = int(sys.argv[2])
    month = int(sys.argv[3])
    domain = sys.argv[4] if len(sys.argv) > 4 else "mammoth"
    goes = sys.argv[5] if len(sys.argv) > 5 else "goes18"

    # Auto-compute end of month
    start_day = 1
    end_day = calendar.monthrange(year, month)[1]

    print(f"Processing {year}-{month:02d} ({start_day} → {end_day})")
    print(f"Input dir: {in_dir}")
    print(f"Domain: {domain}, GOES: {goes}")

    channel_list = ["C02", "C05", "C13"]

    # Process days in parallel
    days = list(range(start_day, end_day + 1))
    max_workers = min(cpu_count(), 4)
    
    print(f"Processing {len(days)} days with {max_workers} workers")
    # Prepare arguments for each day
    args_list = [(day, year, month, in_dir, goes, domain) for day in days]
    
    with Pool(processes=max_workers) as pool:
        results = pool.map(process_day, args_list)
    
    for result in results:
        print(result)


if __name__ == "__main__":
    main()
