#!/usr/bin/env python
"""Download hourly NLDAS-2 primary-forcing LWdown subsets around Gothic."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from urllib.parse import quote

import pandas as pd


DEFAULT_OUT_DIR = Path("/glade/derecho/scratch/cdalden/tmp/nldas2_lw")
BASE_URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi"
PRODUCT = "NLDAS_FORA0125_H"
VARIABLE = "LWdown"
GOTHIC_BBOX = (-107.1, 38.8, -106.8, 39.1)


def parse_hours(spec: str) -> list[int]:
    hours: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = [int(x) for x in part.split("-", 1)]
            if start <= end:
                hours.extend(range(start, end + 1))
            else:
                hours.extend(list(range(start, 24)) + list(range(0, end + 1)))
        else:
            hours.append(int(part))
    out = []
    for hour in hours:
        if hour < 0 or hour > 23:
            raise ValueError(f"Invalid hour {hour}")
        if hour not in out:
            out.append(hour)
    return out


def build_times(start: str, end: str, hours: list[int]) -> pd.DatetimeIndex:
    times = []
    for day in pd.date_range(pd.Timestamp(start), pd.Timestamp(end), freq="D"):
        for hour in hours:
            offset_day = day + pd.Timedelta(days=1) if hour == 0 and 0 in hours and hours.index(0) > 0 else day
            times.append(offset_day + pd.Timedelta(hours=hour))
    return pd.DatetimeIndex(times).drop_duplicates().sort_values()


def output_name(dt: pd.Timestamp) -> str:
    return f"{PRODUCT}.A{dt:%Y%m%d}.{dt:%H}00.020.{VARIABLE}.SUB.nc4"


def existing_match(out_dir: Path, dt: pd.Timestamp) -> Path | None:
    matches = sorted(out_dir.glob(f"*A{dt:%Y%m%d}.{dt:%H}00*{VARIABLE}*"))
    return matches[0] if matches else None


def nldas2_url(dt: pd.Timestamp, bbox: tuple[float, float, float, float]) -> str:
    year = dt.strftime("%Y")
    doy = f"{dt.dayofyear:03d}"
    base_file = f"{PRODUCT}.A{dt:%Y%m%d}.{dt:%H}00.020.nc"
    filename = f"/data/NLDAS/{PRODUCT}.2.0/{year}/{doy}/{base_file}"
    label = output_name(dt)
    west, south, east, north = bbox
    params = {
        "FILENAME": filename,
        "SERVICE": "L34RS_LDAS",
        "VARIABLES": VARIABLE,
        "VERSION": "1.02",
        "FORMAT": "bmM0Lw",
        "DATASET_VERSION": "2.0",
        "BBOX": f"{west},{south},{east},{north}",
        "LABEL": label,
    }
    query = "&".join(f"{key}={quote(value, safe='')}" for key, value in params.items())
    return f"{BASE_URL}?{query}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", required=True, help="First date, YYYY-MM-DD")
    parser.add_argument("--end", required=True, help="Last date, YYYY-MM-DD")
    parser.add_argument("--hours", default="14-23,0")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--bbox",
        default=",".join(str(value) for value in GOTHIC_BBOX),
        help="west,south,east,north bounding box",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def parse_bbox(spec: str) -> tuple[float, float, float, float]:
    parts = [float(part.strip()) for part in spec.split(",")]
    if len(parts) != 4:
        raise ValueError("--bbox must be west,south,east,north")
    return tuple(parts)


def download_one(dt: pd.Timestamp, args: argparse.Namespace) -> Path:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / output_name(dt)
    existing = existing_match(args.out_dir, dt)
    if existing and not args.force:
        print(f"Exists {existing}", flush=True)
        return existing

    url = nldas2_url(dt, parse_bbox(args.bbox))
    if args.dry_run:
        print(f"Would download {dt} -> {out}", flush=True)
        print(url, flush=True)
        return out

    tmp = out.with_suffix(out.suffix + ".part")
    auth_args = []
    netrc = Path.home() / ".netrc"
    if netrc.exists():
        auth_args = ["--netrc-file", str(netrc)]
    else:
        print(
            f"WARNING: {netrc} does not exist. GES DISC usually requires Earthdata "
            "credentials in ~/.netrc; cookie-only auth may fail with HTTP 401.",
            flush=True,
        )

    cmd = [
        "curl",
        "-f",
        "-L",
        "--retry",
        "4",
        "--retry-delay",
        "10",
        *auth_args,
        "-c",
        str(Path.home() / ".urs_cookies"),
        "-b",
        str(Path.home() / ".urs_cookies"),
        "-o",
        str(tmp),
        url,
    ]
    print(f"Downloading {dt} -> {out}", flush=True)
    subprocess.run(cmd, check=True)
    tmp.replace(out)
    return out


def main() -> None:
    args = parse_args()
    times = build_times(args.start, args.end, parse_hours(args.hours))
    print(f"Requested {len(times)} NLDAS2 LW hours from {times.min()} to {times.max()}", flush=True)
    for dt in times:
        download_one(pd.Timestamp(dt), args)


if __name__ == "__main__":
    main()
