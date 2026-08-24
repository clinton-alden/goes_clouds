#!/usr/bin/env python3
"""Download Terra/Aqua MODIS C6.1 cloud-mask swaths over Gothic, Colorado.

The CMR catalog and Earthdata Cloud file URLs are public.  An optional Earthdata
bearer token may be supplied with --token-file or EARTHDATA_TOKEN if NASA's
access policy changes. Existing nonempty files are skipped, so interrupted runs
can be restarted safely.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path


CMR_URL = "https://cmr.earthdata.nasa.gov/search/granules.json"
COLLECTIONS = {
    "MOD35_L2": "C1443561895-LAADS",  # Terra, Collection 6.1
    "MYD35_L2": "C1443715587-LAADS",  # Aqua, Collection 6.1
}
GOTHIC_BBOX = (-107.08, 38.89, -106.90, 39.03)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", default="2021-10-01")
    parser.add_argument("--end", default="2023-06-30", help="Exclusive end date")
    parser.add_argument(
        "--output-dir", default="/glade/u/home/cdalden/scratch/colorado_modis"
    )
    parser.add_argument("--products", nargs="+", choices=COLLECTIONS, default=list(COLLECTIONS))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--token-file", type=Path)
    parser.add_argument("--catalog-only", action="store_true")
    return parser.parse_args()


def request_json(url: str, params: dict[str, str]) -> tuple[dict, dict]:
    req = urllib.request.Request(url + "?" + urllib.parse.urlencode(params))
    req.add_header("Accept", "application/json")
    with urllib.request.urlopen(req, timeout=120) as response:
        return json.load(response), dict(response.headers)


def catalog(product: str, start: str, end: str) -> list[dict]:
    entries: list[dict] = []
    page_num = 1
    while True:
        payload, _ = request_json(
            CMR_URL,
            {
                "collection_concept_id": COLLECTIONS[product],
                "bounding_box": ",".join(map(str, GOTHIC_BBOX)),
                "temporal": f"{start}T00:00:00Z,{end}T00:00:00Z",
                "page_size": "2000",
                "page_num": str(page_num),
            },
        )
        page = payload["feed"]["entry"]
        entries.extend(page)
        if len(page) < 2000:
            return entries
        page_num += 1


def data_url(entry: dict) -> str:
    for link in entry.get("links", []):
        if (
            link.get("rel", "").endswith("/data#")
            and link.get("type") == "application/x-hdfeos"
            and link.get("href", "").startswith("https://")
        ):
            return link["href"]
    raise ValueError(f"No HDF-EOS data URL for {entry['producer_granule_id']}")


def load_token(args: argparse.Namespace) -> str:
    if os.environ.get("EARTHDATA_TOKEN"):
        return os.environ["EARTHDATA_TOKEN"].strip()
    if args.token_file:
        return args.token_file.expanduser().read_text().strip()
    default = Path.home() / ".earthdata_token"
    if default.exists():
        return default.read_text().strip()
    return ""


def destination(root: Path, entry: dict) -> Path:
    name = entry["producer_granule_id"]
    product, acquisition = name.split(".")[:2]
    stamp = datetime.strptime(acquisition, "A%Y%j")
    return root / product / f"{stamp.year:04d}" / f"{stamp.timetuple().tm_yday:03d}" / name


def download_one(entry: dict, root: Path, token: str, retries: int) -> tuple[str, str]:
    target = destination(root, entry)
    if target.exists() and target.stat().st_size > 0:
        return "skipped", str(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    req = urllib.request.Request(data_url(entry))
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    req.add_header("User-Agent", "goes-work-modis-downloader/1.0")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=300) as response, partial.open("wb") as out:
                while chunk := response.read(1024 * 1024):
                    out.write(chunk)
            if partial.stat().st_size == 0:
                raise OSError("server returned an empty file")
            partial.replace(target)
            return "downloaded", str(target)
        except (OSError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            partial.unlink(missing_ok=True)
            if attempt + 1 == retries:
                return "failed", f"{target}: {exc}"
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def main() -> int:
    args = parse_args()
    root = Path(args.output_dir)
    all_entries: list[dict] = []
    for product in args.products:
        found = catalog(product, args.start, args.end)
        print(f"{product}: {len(found)} intersecting granules", flush=True)
        all_entries.extend(found)

    root.mkdir(parents=True, exist_ok=True)
    manifest = root / "cmr_manifest.jsonl"
    with manifest.open("w") as stream:
        for entry in all_entries:
            stream.write(json.dumps(entry, separators=(",", ":")) + "\n")
    print(f"Manifest: {manifest} ({len(all_entries)} total granules)", flush=True)
    if args.catalog_only:
        return 0

    token = load_token(args)
    counts = {"downloaded": 0, "skipped": 0, "failed": 0}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(download_one, e, root, token, args.retries) for e in all_entries]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            status, detail = future.result()
            counts[status] += 1
            if status == "failed":
                print(f"FAILED: {detail}", file=sys.stderr, flush=True)
            if index % 25 == 0 or index == len(futures):
                print(f"Progress {index}/{len(futures)}: {counts}", flush=True)
    print(f"Complete: {counts}", flush=True)
    return 1 if counts["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
