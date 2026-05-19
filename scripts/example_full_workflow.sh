#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# End-to-end GOES VINTAGE example.
#
# Edit the variables below, or override them from the shell:
#
#   DOMAIN=mammoth \
#   LON_MIN=-119.5 LAT_MIN=36.5 LON_MAX=-118.5 LAT_MAX=37.9 \
#   START_DATE=2022-03-25 END_DATE=2022-03-31 \
#   BASE_DIR=/path/to/scratch/mammoth \
#   ./scripts/example_full_workflow.sh
#
# DOMAIN is only an output filename label. The lon/lat bounds define the
# actual processing area.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Load repo-local private credentials when available. This file is git-ignored.
if [[ -f "${REPO_DIR}/private_env.sh" ]]; then
  # shellcheck source=/dev/null
  source "${REPO_DIR}/private_env.sh"
fi

if [[ -z "${PYTHON_BIN:-}" && -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
  PYTHON_BIN="${CONDA_PREFIX}/bin/python"
else
  PYTHON_BIN="${PYTHON_BIN:-python}"
fi

DOMAIN="${DOMAIN:-colorado}"
GOES="${GOES:-goes16}"
START_DATE="${START_DATE:-2020-06-30}"
END_DATE="${END_DATE:-2020-06-30}"

LON_MIN="${LON_MIN:--109.0}"
LAT_MIN="${LAT_MIN:-37.0}"
LON_MAX="${LON_MAX:--104.0}"
LAT_MAX="${LAT_MAX:-41.0}"

# GOES_HOURS accepts a single hour ("20"), a range ("18-22"), or a comma list.
# GOES_TIMESTEPS_PER_HOUR=1 keeps the example small. Set to "None" for all scans.
GOES_HOURS="${GOES_HOURS:-20}"
GOES_TIMESTEPS_PER_HOUR="${GOES_TIMESTEPS_PER_HOUR:-1}"

BASE_DIR="${BASE_DIR:-${REPO_DIR}/demo_output/${DOMAIN}}"
THRESHOLD_CSV="${THRESHOLD_CSV:-${REPO_DIR}/thresholds/gothic_temp_bin_rgb_thresholds_10c.csv}"

# Set RUN_WORKFLOW=0 to print the workflow config and skip downloads/processing.
RUN_WORKFLOW="${RUN_WORKFLOW:-1}"
OVERWRITE_MASK="${OVERWRITE_MASK:-1}"
KEEP_MASK_DIAGNOSTICS="${KEEP_MASK_DIAGNOSTICS:-0}"

# Set CREATE_PNG=0 to skip the final RGB/VINTAGE-mask side-by-side PNG.
CREATE_PNG="${CREATE_PNG:-1}"
PLOT_DATE="${PLOT_DATE:-${START_DATE}}"
PLOT_TIME_UTC="${PLOT_TIME_UTC:-}"
PLOT_PATH="${PLOT_PATH:-}"

export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"
export DOMAIN GOES START_DATE END_DATE
export LON_MIN LAT_MIN LON_MAX LAT_MAX
export GOES_HOURS GOES_TIMESTEPS_PER_HOUR
export BASE_DIR THRESHOLD_CSV
export RUN_WORKFLOW OVERWRITE_MASK
export KEEP_MASK_DIAGNOSTICS
export CREATE_PNG PLOT_DATE PLOT_TIME_UTC PLOT_PATH

"${PYTHON_BIN}" - <<'PY'
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import workflow as wf


def env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() not in {"0", "false", "no", "off", ""}


def env_optional_int(name: str) -> int | None:
    raw = os.environ.get(name, "").strip()
    if raw == "" or raw.lower() == "none":
        return None
    return int(raw)


config = wf.WorkflowConfig(
    domain=os.environ["DOMAIN"],
    goes=os.environ["GOES"],
    start_date=os.environ["START_DATE"],
    end_date=os.environ["END_DATE"],
    lon_min=float(os.environ["LON_MIN"]),
    lat_min=float(os.environ["LAT_MIN"]),
    lon_max=float(os.environ["LON_MAX"]),
    lat_max=float(os.environ["LAT_MAX"]),
    base_dir=Path(os.environ["BASE_DIR"]),
    goes_hours=os.environ["GOES_HOURS"],
    goes_timesteps_per_hour=env_optional_int("GOES_TIMESTEPS_PER_HOUR"),
    threshold_csv=Path(os.environ["THRESHOLD_CSV"]),
    overwrite_mask=env_bool("OVERWRITE_MASK", True),
    keep_mask_diagnostics=env_bool("KEEP_MASK_DIAGNOSTICS", False),
    run=env_bool("RUN_WORKFLOW", True),
)

print("=== GOES WORKFLOW CONFIG ===")
print(f"output label: {config.domain}")
print(f"GOES satellite: {config.goes}")
print(f"dates: {config.start_date} to {config.end_date}")
print(f"bounds: {config.lon_min}, {config.lat_min}, {config.lon_max}, {config.lat_max}")
print(f"GOES hours: {config.goes_hours}")
print(f"timesteps per hour: {config.goes_timesteps_per_hour}")
print(f"mask UTC window: {config.start_hour_utc:g}-{config.end_hour_utc:g}")
print(f"base dir: {config.base_dir.resolve()}")
print(f"keep mask diagnostics: {config.keep_mask_diagnostics}")

wf.validate_credentials(config)
wf.download_goes(config)
wf.orthorectify(config)
wf.build_zarr(config)
rgb_paths = wf.build_rgb(config)
mask_paths = wf.apply_mask(config)

print("=== WORKFLOW OUTPUTS ===")
for path in rgb_paths + mask_paths:
    print(path)

if not env_bool("CREATE_PNG", True):
    raise SystemExit(0)

if not config.run:
    print("Dry run: skipping PNG creation")
    raise SystemExit(0)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plot_date = pd.Timestamp(os.environ.get("PLOT_DATE") or config.start_date)
rgb_path = config.rgb_path(plot_date)
mask_path = config.mask_path(plot_date)

plot_path_raw = os.environ.get("PLOT_PATH", "").strip()
if plot_path_raw:
    plot_path = Path(plot_path_raw)
else:
    plot_path = (
        config.base_dir
        / config.goes
        / "plots"
        / f"{config.goes}_rgb_mask_{config.domain}_{config.date_ymd(plot_date)}.png"
    )
plot_path.parent.mkdir(parents=True, exist_ok=True)

with xr.open_dataset(rgb_path) as rgb_ds, xr.open_dataset(mask_path) as mask_ds:
    requested_time = os.environ.get("PLOT_TIME_UTC", "").strip()
    plot_time = (
        pd.Timestamp(requested_time)
        if requested_time
        else pd.Timestamp(mask_ds["t"].values[0])
    )
    rgb_frame = rgb_ds.sel(t=plot_time, method="nearest")
    mask_frame = mask_ds.sel(t=plot_time, method="nearest")
    actual_time = pd.Timestamp(mask_frame["t"].values)

    rgb_image = np.stack(
        [
            rgb_frame["red"].values,
            rgb_frame["green"].values,
            rgb_frame["blue"].values,
        ],
        axis=-1,
    )
    rgb_image = np.clip(np.nan_to_num(rgb_image, nan=0.0), 0.0, 1.0)
    vintage_mask = mask_frame["vintage_mask"].values
    lon = rgb_ds["longitude"].values
    lat = rgb_ds["latitude"].values
    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
axes[0].imshow(rgb_image, origin="lower", extent=extent, aspect="auto")
axes[0].set_title(f"RGB {actual_time:%Y-%m-%d %H:%M UTC}")
axes[0].set_xlabel("Longitude")
axes[0].set_ylabel("Latitude")

masked = np.ma.masked_invalid(vintage_mask)
axes[1].imshow(masked, origin="lower", extent=extent, aspect="auto", cmap="gray_r", vmin=0, vmax=1)
axes[1].set_title("GOES VINTAGE mask")
axes[1].set_xlabel("Longitude")
axes[1].set_ylabel("Latitude")

fig.savefig(plot_path, dpi=160)
plt.close(fig)
print(f"Wrote PNG: {plot_path}")
PY
