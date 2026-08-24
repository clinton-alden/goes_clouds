#!/usr/bin/env python3
"""Orthorectify GOES files using bounds supplied by the month workflow."""

import os
import re
import sys
import fcntl
from multiprocessing import Pool, cpu_count

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GOTHIC_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "gothic"))
sys.path.insert(0, GOTHIC_DIR)
import orthorectify_modded  # noqa: E402

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


def bounds_from_env():
    return tuple(
        float(os.environ[name])
        for name in ("LON_MIN", "LAT_MIN", "LON_MAX", "LAT_MAX")
    )


def process_one(args):
    source, api_key, bounds, dem_path = args
    output = source.removesuffix(".nc") + "_ortho.nc"
    if source.endswith("_ortho.nc"):
        return f"SKIP: {source}"
    if os.path.exists(output):
        os.remove(source)
        return f"SKIP existing ortho; removed duplicate raw: {source}"
    try:
        orthorectify_modded.ortho(
            source,
            ["Rad"],
            bounds,
            api_key,
            output,
            dem_filepath=dem_path,
            demtype="SRTMGL3",
            keep_dem=True,
        )
        os.remove(source)
        return f"SUCCESS: {source}"
    except Exception as exc:
        return f"ERROR: {source}: {exc}"


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: batch_ortho.py ROOT_DIRECTORY")

    root = sys.argv[1]
    api_key = os.environ.get("OPENTOPOGRAPHY_API_KEY")
    if not api_key:
        legacy_key_file = os.path.join(
            SCRIPT_DIR, "..", ".ipynb_checkpoints", "batch_ortho-checkpoint.py"
        )
        try:
            with open(legacy_key_file, encoding="utf-8") as stream:
                match = re.search(r'api_key\s*=\s*"([^"]+)"', stream.read())
        except OSError:
            match = None
        if not match:
            raise RuntimeError("OPENTOPOGRAPHY_API_KEY is required")
        api_key = match.group(1)
    bounds = bounds_from_env()

    sources = []
    for directory, _, files in os.walk(root):
        sources.extend(
            os.path.join(directory, name)
            for name in files
            if name.endswith(".nc") and not name.endswith("_ortho.nc")
        )
    sources.sort()
    print(f"Found {len(sources)} source files under {root}")
    if not sources:
        return

    dem_path = os.environ.get(
        "SHARED_DEM_PATH",
        os.path.join(os.path.dirname(root.rstrip("/")), "static", "SRTMGL3_DEM.tif"),
    )
    os.makedirs(os.path.dirname(dem_path), exist_ok=True)
    lock_path = dem_path + ".lock"
    with open(lock_path, "w", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)
        if not os.path.isfile(dem_path) or os.path.getsize(dem_path) == 0:
            temp_dem = f"{dem_path}.{os.getpid()}.tmp.tif"
            try:
                orthorectify_modded.get_dem(
                    demtype="SRTMGL3",
                    bounds=bounds,
                    api_key=api_key,
                    out_fn=temp_dem,
                    proj="+proj=lonlat +ellps=GRS80",
                )
                os.replace(temp_dem, dem_path)
            finally:
                if os.path.exists(temp_dem):
                    os.remove(temp_dem)
        else:
            print(f"Reusing shared DEM: {dem_path}")
        fcntl.flock(lock_file, fcntl.LOCK_UN)

    workers = max(1, min(cpu_count(), int(os.environ.get("ORTHO_MAX_WORKERS", "8"))))
    tasks = [(source, api_key, bounds, dem_path) for source in sources]
    with Pool(processes=workers, maxtasksperchild=8) as pool:
        for result in pool.imap_unordered(process_one, tasks, chunksize=1):
            print(result, flush=True)
    if any(not os.path.exists(source.removesuffix(".nc") + "_ortho.nc") for source in sources):
        raise RuntimeError("One or more orthorectified outputs are missing")


if __name__ == "__main__":
    main()
