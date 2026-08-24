#!/usr/bin/env python3
"""Download, orthorectify, assemble, mask, pack, and clean one GOES day."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
GOTHIC = ROOT / "processing/gothic"
sys.path.insert(0, str(GOTHIC))
from utils import _abi_scan_midpoint_from_filename  # noqa: E402


def channel_files(raw_day: Path, channel: str) -> list[Path]:
    return sorted(raw_day.glob(f"ABI-L1b-RadC/*/{channel}/*_ortho.nc"))


def nominal_scan_key(path: Path) -> str:
    """Return the ABI nominal YYYYJJJHHMM scan start shared by all bands."""
    match = re.search(r"_s(\d{11})", path.name)
    if not match:
        raise RuntimeError(f"Cannot parse ABI scan start from {path.name}")
    return match.group(1)


def normalize(value: np.ndarray, low: float, high: float, invert: bool = False) -> np.ndarray:
    result = np.clip((value - low) / (high - low), 0.0, 1.0)
    return 1.0 - result if invert else result


def assemble_float_rgb(raw_day: Path, output: Path) -> None:
    files = {channel: channel_files(raw_day, channel) for channel in ("C02", "C05", "C13")}
    counts = {channel: len(paths) for channel, paths in files.items()}
    maps = {
        channel: {nominal_scan_key(path): path for path in paths}
        for channel, paths in files.items()
    }
    if any(len(maps[channel]) != counts[channel] for channel in maps):
        raise RuntimeError(f"Duplicate nominal scan times in orthorectified inputs: {counts}")
    common = sorted(set(maps["C02"]) & set(maps["C05"]) & set(maps["C13"]))
    minimum = int(os.environ.get("MIN_ALIGNED_SCANS", "1"))
    if len(common) < minimum:
        raise RuntimeError(
            f"Only {len(common)} common aligned scans; require {minimum}; channel counts={counts}"
        )

    red, green, blue, times = [], [], [], []
    latitude = longitude = None
    for key in common:
        c02, c05, c13 = (maps[channel][key] for channel in ("C02", "C05", "C13"))
        with xr.open_dataset(c02) as d2, xr.open_dataset(c05) as d5, xr.open_dataset(c13) as d13:
            if latitude is None:
                latitude = np.asarray(d2.latitude, dtype=np.float64)
                longitude = np.asarray(d2.longitude, dtype=np.float64)
            green.append(normalize(np.asarray(d2.ref, dtype=np.float32), 0.0, 0.78))
            blue.append(normalize(np.asarray(d5.ref, dtype=np.float32), 0.01, 0.59))
            red.append(normalize(np.asarray(d13.tb, dtype=np.float32), 219.65, 280.65, invert=True))
            times.append(_abi_scan_midpoint_from_filename(str(c02)))

    ds = xr.Dataset(
        {name: (("t", "latitude", "longitude"), np.stack(values).astype(np.float32))
         for name, values in (("red", red), ("green", green), ("blue", blue))},
        coords={"t": np.asarray(times), "latitude": latitude, "longitude": longitude},
        attrs={
            "rgb_recipe": "NOAA Day Cloud Phase Distinction",
            "processing": "direct orthorectified channels; no intermediate Zarr",
            "scan_alignment": "nominal_start_intersection",
            "source_channel_counts": str(counts),
            "common_aligned_scans": len(common),
        },
    ).sortby("t")
    output.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(output, encoding={name: {"zlib": True, "complevel": 1, "dtype": "float32"} for name in ("red", "green", "blue")})


def pack_rgb(source: Path, output: Path) -> None:
    with xr.open_dataset(source) as ds:
        packed = xr.Dataset(
            {name: (("t", "latitude", "longitude"), np.clip(ds[name].values, 0, 1).astype(np.float32))
             for name in ("red", "green", "blue")},
            coords={"t": ds.t, "latitude": ds.latitude, "longitude": ds.longitude},
            attrs={**ds.attrs, "packing": "uint8; physical_value=stored_value/255"},
        )
        for name in ("red", "green", "blue"):
            packed[name].attrs.update(valid_min=np.float32(0), valid_max=np.float32(1))
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_suffix(".nc.part")
        packed.to_netcdf(temporary, encoding={name: {
            "zlib": True, "complevel": 4, "dtype": "uint8",
            "scale_factor": np.float64(1 / 255), "add_offset": np.float64(0),
        } for name in ("red", "green", "blue")})
        temporary.replace(output)


def outputs_are_valid(rgb: Path, mask: Path, bounds: list[float] | tuple[float, ...]) -> bool:
    if not rgb.is_file() or not mask.is_file() or not rgb.stat().st_size or not mask.stat().st_size:
        return False
    try:
        with xr.open_dataset(rgb) as r, xr.open_dataset(mask) as m:
            expected_dims = {"t": 144, "latitude": 961, "longitude": 1601}
            return (
                dict(r.sizes) == expected_dims
                and dict(m.sizes) == expected_dims
                and np.array_equal(r.t.values, m.t.values)
                and np.all(np.diff(r.latitude) > 0)
                and np.all(np.diff(r.longitude) > 0)
                and abs(float(r.longitude.min()) - bounds[0]) < 0.01
                and abs(float(r.longitude.max()) - bounds[2]) < 0.01
                and abs(float(r.latitude.min()) - bounds[1]) < 0.01
                and abs(float(r.latitude.max()) - bounds[3]) < 0.01
                and r.red.encoding.get("dtype") == np.dtype("uint8")
                and m.attrs.get("tree_logic") == "AND within each leaf; OR across cloudy leaves"
            )
    except Exception:
        return False


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--date", required=True, help="YYYY-MM-DD")
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--satellite", default="goes16")
    p.add_argument("--domain", default="colorado_5x4")
    p.add_argument("--bounds", type=float, nargs=4, default=(-109, 37, -104, 41))
    p.add_argument("--dem", type=Path, required=True)
    p.add_argument("--era-dir", type=Path, required=True)
    p.add_argument("--rules", type=Path, required=True)
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--keep-work", action="store_true")
    a = p.parse_args()
    date = pd.Timestamp(a.date)
    token = date.strftime("%Y%m%d")
    raw_day = a.base / a.satellite / str(date.year) / str(date.month) / str(date.day)
    rgb_dir = a.base / a.satellite / "rgb_composite_packed"
    mask_dir = a.base / a.satellite / "cloud_mask_tempbin_tree_packed"
    temp_rgb = a.base / "tmp" / f"{a.satellite}_C02_C05_C13_rgb_{a.domain}_{token}.nc"
    rgb = rgb_dir / f"{a.satellite}_C02_C05_C13_rgb_{a.domain}_{token}.nc"
    mask = mask_dir / f"{a.satellite}_C02_C05_C13_rgb_{a.domain}_{token}_cloud_binary_tempbin_tree.nc"

    if outputs_are_valid(rgb, mask, a.bounds):
        print(f"VALID existing packed RGB and exact-tree mask: {token}")
        return 0

    env = os.environ.copy()
    env.update({
        "GOES_HOURS": "00,01,14,15,16,17,18,19,20,21,22,23",
        "LON_MIN": str(a.bounds[0]), "LAT_MIN": str(a.bounds[1]),
        "LON_MAX": str(a.bounds[2]), "LAT_MAX": str(a.bounds[3]),
        "SHARED_DEM_PATH": str(a.dem), "ORTHO_MAX_WORKERS": str(a.workers),
        "ORTHO_GRID_DX_DEG": "0.003125", "ORTHO_GRID_DY_DEG": "0.004166666667",
        "HDF5_USE_FILE_LOCKING": "FALSE", "OMP_NUM_THREADS": "1",
        "ORTHO_COMPACT_OUTPUT": "1",
    })
    completed = False
    try:
        def download_channel(channel: str) -> None:
            subprocess.run([
                a.python, str(ROOT / "data_download/download-goes.py"),
                "-B", f"noaa-{a.satellite}", "-Y", str(date.year), "-M", str(date.month),
                "-D", str(date.day), str(date.day), "-p", "ABI-L1b-RadC", "-c", channel,
                "-b", *map(str, a.bounds), "-d", str(a.base),
            ], check=True, env=env)
        existing_paths = {channel: channel_files(raw_day, channel) for channel in ("C02", "C05", "C13")}
        existing_counts = {channel: len(paths) for channel, paths in existing_paths.items()}
        existing_keys = {
            channel: {nominal_scan_key(path) for path in paths}
            for channel, paths in existing_paths.items()
        }
        existing_common = set.intersection(*existing_keys.values()) if existing_keys else set()
        minimum = int(env.get("MIN_ALIGNED_SCANS", "144"))
        if len(existing_common) >= minimum:
            print(
                f"Reusing orthorectified day: channel counts={existing_counts}, "
                f"common aligned={len(existing_common)}, required={minimum}",
                flush=True,
            )
        else:
            with ThreadPoolExecutor(max_workers=3) as pool:
                list(pool.map(download_channel, ("C02", "C05", "C13")))
        for attempt in range(1, 4):
            try:
                subprocess.run([a.python, str(ROOT / "processing/east_river_goes/batch_ortho.py"), str(raw_day)], check=True, env=env)
                break
            except subprocess.CalledProcessError:
                if attempt == 3:
                    raise
                bad_raw = [path for path in raw_day.rglob("*.nc") if not path.name.endswith("_ortho.nc") and not path.with_name(path.stem + "_ortho.nc").exists()]
                print(f"Orthorectification attempt {attempt} left {len(bad_raw)} bad/missing scans; redownloading", flush=True)
                for path in bad_raw:
                    path.unlink(missing_ok=True)
                with ThreadPoolExecutor(max_workers=3) as pool:
                    list(pool.map(download_channel, ("C02", "C05", "C13")))
        if temp_rgb.is_file() and temp_rgb.stat().st_size:
            print(f"Reusing assembled float RGB: {temp_rgb}", flush=True)
        else:
            assemble_float_rgb(raw_day, temp_rgb)
        subprocess.run([
            a.python, str(ROOT / "processing/build_rgb_tempbin_tree_mask.py"), str(temp_rgb), str(mask),
            "--era-dir", str(a.era_dir), "--rules", str(a.rules), "--binary-only",
        ], check=True, env=env)
        pack_rgb(temp_rgb, rgb)
        completed = True
    finally:
        if completed and not a.keep_work:
            shutil.rmtree(raw_day, ignore_errors=True)
            temp_rgb.unlink(missing_ok=True)

    with xr.open_dataset(rgb) as r, xr.open_dataset(mask) as m:
        if dict(r.sizes) != dict(m.sizes):
            raise RuntimeError(f"RGB/mask dimensions differ: {dict(r.sizes)} vs {dict(m.sizes)}")
        if not np.all(np.diff(r.latitude) > 0) or not np.all(np.diff(r.longitude) > 0):
            raise RuntimeError("Output coordinates are not south-to-north/west-to-east")
        if m.attrs.get("tree_logic") != "AND within each leaf; OR across cloudy leaves":
            raise RuntimeError("Mask provenance does not contain exact decision-tree logic")
    print(f"COMPLETE rgb={rgb} mask={mask}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
