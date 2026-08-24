#!/usr/bin/env python
"""Parallel downloader for NASA GES DISC NLDAS subset URL lists.

The GES DISC subset text file contains one URL per hour plus usually a README.
This script filters to the data URLs, parses the LABEL parameter for a clean
local filename, and downloads in parallel with an Earthdata bearer token.

Token sources, in order:
  1. --token-file
  2. NASA_TOKEN environment variable

Using this script avoids putting the token in each wget command where it is
visible in process listings.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url_list", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("."))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--token-file", type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--min-existing-bytes", type=int, default=1024)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=180)
    return parser.parse_args()


def read_token(args: argparse.Namespace) -> str:
    if args.token_file:
        token = args.token_file.read_text().strip()
    else:
        token = os.environ.get("NASA_TOKEN", "").strip()
    if not token:
        raise RuntimeError("Set NASA_TOKEN or pass --token-file with an Earthdata bearer token.")
    if token.lower().startswith("bearer "):
        token = token.split(None, 1)[1]
    return token


def iter_data_urls(path: Path) -> list[str]:
    urls = []
    for raw in path.read_text().splitlines():
        url = raw.strip()
        if not url or "HTTP_services.cgi" not in url:
            continue
        urls.append(url)
    return urls


def label_from_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    qs = urllib.parse.parse_qs(parsed.query)
    label = qs.get("LABEL", [None])[0]
    if label:
        return Path(label).name

    filename = qs.get("FILENAME", [None])[0]
    if filename:
        return Path(filename).name + ".SUB.nc4"

    raise ValueError(f"Could not infer output label from URL: {url}")


def download_one(
    url: str,
    out_dir: Path,
    token: str,
    force: bool,
    min_existing_bytes: int,
    retries: int,
    timeout: int,
    dry_run: bool,
) -> tuple[str, str]:
    name = label_from_url(url)
    out = out_dir / name
    part = out.with_suffix(out.suffix + ".part")

    if out.exists() and not force and out.stat().st_size >= min_existing_bytes:
        return "exists", name
    if dry_run:
        return "would_download", name

    headers = {"Authorization": f"Bearer {token}"}
    request = urllib.request.Request(url, headers=headers)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                if response.status >= 400:
                    raise urllib.error.HTTPError(
                        url, response.status, response.reason, response.headers, None
                    )
                with part.open("wb") as fh:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        fh.write(chunk)
            if part.stat().st_size < min_existing_bytes:
                raise RuntimeError(f"Downloaded file is suspiciously small: {part.stat().st_size} bytes")
            part.replace(out)
            return "downloaded", name
        except Exception as exc:  # noqa: BLE001 - report final exception below
            last_error = exc
            try:
                part.unlink()
            except FileNotFoundError:
                pass
            if attempt < retries:
                time.sleep(min(60, 2**attempt))
    return "failed", f"{name}: {last_error}"


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    urls = iter_data_urls(args.url_list)
    print(f"Data URLs: {len(urls)}")
    print(f"Output dir: {args.out_dir}")
    print(f"Workers: {args.workers}")
    token = "DRYRUN" if args.dry_run else read_token(args)

    counts: dict[str, int] = {}
    failures: list[str] = []
    with cf.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                download_one,
                url,
                args.out_dir,
                token,
                args.force,
                args.min_existing_bytes,
                args.retries,
                args.timeout,
                args.dry_run,
            )
            for url in urls
        ]
        for idx, future in enumerate(cf.as_completed(futures), start=1):
            status, detail = future.result()
            counts[status] = counts.get(status, 0) + 1
            if status == "failed":
                failures.append(detail)
            if idx == 1 or idx % 50 == 0 or idx == len(futures):
                summary = " ".join(f"{key}={value}" for key, value in sorted(counts.items()))
                print(f"{idx}/{len(futures)} {summary}", flush=True)

    if failures:
        print("Failures:")
        for failure in failures[:50]:
            print(failure)
        if len(failures) > 50:
            print(f"... {len(failures) - 50} more failures")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
