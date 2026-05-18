#!/usr/bin/env python3
import argparse
import calendar
import datetime as dt
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import utils


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Backfill missing daily RGB files for one month by checking existing "
            "rgb_composite outputs and required C02/C05/C13 zarr inputs."
        )
    )
    parser.add_argument("--base-dir", required=True, help="Base directory, e.g. /glade/derecho/scratch/cdalden/colorado")
    parser.add_argument("--goes", default="goes16", help="GOES platform (default: goes16)")
    parser.add_argument("--domain", default="colorado", help="Domain name used in filenames (default: colorado)")
    parser.add_argument("--year", type=int, required=True, help="Target year")
    parser.add_argument("--month", type=int, required=True, help="Target month (1-12)")
    parser.add_argument("--start-date", default="2021-10-01", help="Inclusive global start date (YYYY-MM-DD)")
    parser.add_argument("--end-date", default="2023-06-15", help="Inclusive global end date (YYYY-MM-DD)")
    parser.add_argument("--workers", type=int, default=8, help="Parallel workers for missing days in this month")
    parser.add_argument("--dry-run", action="store_true", help="Only print missing dates; do not run RGB generation")
    return parser.parse_args()


def rgb_file_path(goes_root: str, goes: str, domain: str, yyyymmdd: str) -> str:
    name = f"{goes}_C02_C05_C13_rgb_{domain}_{yyyymmdd}.nc"
    return os.path.join(goes_root, "rgb_composite", name)


def zarr_file_path(goes_root: str, goes: str, domain: str, channel: str, yyyymmdd: str) -> str:
    name = f"{goes}_{channel}_{domain}_{yyyymmdd}.zarr"
    return os.path.join(goes_root, channel, name)


def required_zarr_exists(goes_root: str, goes: str, domain: str, yyyymmdd: str) -> bool:
    for channel in ("C02", "C05", "C13"):
        if not os.path.exists(zarr_file_path(goes_root, goes, domain, channel, yyyymmdd)):
            return False
    return True


def process_one_date(goes_root: str, goes: str, domain: str, yyyymmdd: str) -> tuple[str, str, str]:
    out_path = rgb_file_path(goes_root, goes, domain, yyyymmdd)
    if os.path.exists(out_path):
        return (yyyymmdd, "skip", "rgb already exists")
    try:
        utils.goes_rad_to_rgb(goes_root + "/", yyyymmdd, goes, location=domain)
        return (yyyymmdd, "ok", "created")
    except Exception as exc:  # noqa: BLE001
        return (yyyymmdd, "error", str(exc))


def main() -> int:
    args = parse_args()

    goes_root = os.path.join(args.base_dir, args.goes)
    start_date = dt.date.fromisoformat(args.start_date)
    end_date = dt.date.fromisoformat(args.end_date)

    if start_date > end_date:
        raise ValueError("start-date must be <= end-date")
    if not (1 <= args.month <= 12):
        raise ValueError("month must be 1-12")
    if not os.path.isdir(goes_root):
        raise FileNotFoundError(f"GOES root not found: {goes_root}")

    month_days = calendar.monthrange(args.year, args.month)[1]
    candidate_dates: list[str] = []
    for day in range(1, month_days + 1):
        d = dt.date(args.year, args.month, day)
        if start_date <= d <= end_date:
            candidate_dates.append(d.strftime("%Y%m%d"))

    existing_rgb = 0
    missing_rgb_with_inputs: list[str] = []
    missing_rgb_without_inputs: list[str] = []

    for yyyymmdd in candidate_dates:
        if os.path.exists(rgb_file_path(goes_root, args.goes, args.domain, yyyymmdd)):
            existing_rgb += 1
            continue
        if required_zarr_exists(goes_root, args.goes, args.domain, yyyymmdd):
            missing_rgb_with_inputs.append(yyyymmdd)
        else:
            missing_rgb_without_inputs.append(yyyymmdd)

    print(
        f"[{args.year}-{args.month:02d}] in-range dates={len(candidate_dates)} "
        f"existing_rgb={existing_rgb} "
        f"missing_with_inputs={len(missing_rgb_with_inputs)} "
        f"missing_without_inputs={len(missing_rgb_without_inputs)}"
    )

    if missing_rgb_without_inputs:
        print("missing_without_inputs_dates=" + ",".join(missing_rgb_without_inputs))

    if not missing_rgb_with_inputs:
        print(f"[{args.year}-{args.month:02d}] nothing to process")
        return 0

    print("process_dates=" + ",".join(missing_rgb_with_inputs))

    if args.dry_run:
        return 0

    workers = max(1, args.workers)
    print(f"[{args.year}-{args.month:02d}] launching workers={workers}")

    ok_count = 0
    skip_count = 0
    err_count = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [
            ex.submit(process_one_date, goes_root, args.goes, args.domain, d)
            for d in missing_rgb_with_inputs
        ]
        for fut in as_completed(futures):
            yyyymmdd, status, msg = fut.result()
            print(f"{yyyymmdd} status={status} msg={msg}")
            if status == "ok":
                ok_count += 1
            elif status == "skip":
                skip_count += 1
            else:
                err_count += 1

    print(
        f"[{args.year}-{args.month:02d}] summary created={ok_count} "
        f"skipped={skip_count} errors={err_count}"
    )
    return 1 if err_count > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())