#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Repair already-processed Colorado RGB/GIF outputs by writing new RGB files
# with the latitude coordinate values reversed, then regenerating GIFs from
# those corrected files using the corrected plotting orientation.
#
# What this script does:
# 1. Finds existing Colorado RGB composite files.
# 2. Writes a new tagged RGB file with latitude coordinate values flipped.
# 3. Regenerates the GIF from that tagged RGB using the corrected plotting code.
#
# New corrected outputs are written with a _latcorrected suffix.
# This script does not create backup copies.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado}"
GOES="${GOES:-goes16}"
DOMAIN="${DOMAIN:-colorado}"
GIF_START="${GIF_START:-1400}"
GIF_END="${GIF_END:-2355}"
DATE_GLOB="${DATE_GLOB:-*}"
EXCLUDE_DATES="${EXCLUDE_DATES:-20230105}"
OUTPUT_TAG="${OUTPUT_TAG:-_latcorrected}"

RGB_DIR="${BASE_DIR}/${GOES}/rgb_composite"
GIF_DIR="${BASE_DIR}/${GOES}/gif_loops"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "ERROR: PYTHON_BIN not executable: ${PYTHON_BIN}" >&2
  exit 2
fi

if [[ ! -d "${RGB_DIR}" ]]; then
  echo "ERROR: RGB directory not found: ${RGB_DIR}" >&2
  exit 2
fi

mkdir -p "${GIF_DIR}" "${SCRIPT_DIR}/plots" "${SCRIPT_DIR}/gifs"

shopt -s nullglob
pattern="${RGB_DIR}/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${DATE_GLOB}.nc"
rgb_files=(${pattern})
shopt -u nullglob

if [[ "${#rgb_files[@]}" -eq 0 ]]; then
  echo "No RGB files matched: ${pattern}"
  exit 0
fi

echo "Found ${#rgb_files[@]} RGB files to repair/regenerate"
echo "Excluded dates: ${EXCLUDE_DATES}"

for rgb_path in "${rgb_files[@]}"; do
  rgb_name="$(basename "${rgb_path}")"
  date="${rgb_name##*_}"
  date="${date%.nc}"

  skip_date=0
  for excluded in ${EXCLUDE_DATES}; do
    if [[ "${date}" == "${excluded}" ]]; then
      skip_date=1
      break
    fi
  done
  if [[ "${skip_date}" -eq 1 ]]; then
    echo "${date}: skipping because it is excluded"
    continue
  fi

  rgb_out="${rgb_path%.nc}${OUTPUT_TAG}.nc"
  gif_path="${GIF_DIR}/goes_rgb_${DOMAIN}_${date}_${GIF_START}_${GIF_END}${OUTPUT_TAG}.gif"
  echo "${date}: writing ${rgb_out} and regenerating ${gif_path}"

  if ! "${PYTHON_BIN}" - <<PY
import os
import shutil
import sys
import xarray as xr

root_dir = "${ROOT_DIR}"
proc_dir = "${SCRIPT_DIR}"
rgb_path = "${rgb_path}"
rgb_out = "${rgb_out}"
gif_path = "${gif_path}"
date = "${date}"
domain = "${DOMAIN}"
gif_start = "${GIF_START}"
gif_end = "${GIF_END}"

sys.path.insert(0, os.path.join(root_dir, "analysis"))
from analysis_utils import make_gif

os.makedirs(os.path.join(proc_dir, "plots"), exist_ok=True)
os.makedirs(os.path.join(proc_dir, "gifs"), exist_ok=True)
os.chdir(proc_dir)

with xr.open_dataset(rgb_path) as ds:
    fixed = ds.copy(deep=True)
    fixed = fixed.assign_coords(latitude=ds["latitude"].values[::-1])
    tmp_rgb = rgb_out + ".tmp"
    if os.path.exists(tmp_rgb):
        os.remove(tmp_rgb)
    fixed.to_netcdf(tmp_rgb, mode="w", format="NETCDF4")

os.replace(tmp_rgb, rgb_out)

with xr.open_dataset(rgb_out) as ds_fixed:
    make_gif(ds_fixed, date, gif_start, gif_end, subset=domain, mask=False)

src = os.path.join(proc_dir, "gifs", f"goes_rgb_{domain}_{date}_{gif_start}_{gif_end}.gif")
if not os.path.exists(src):
    raise FileNotFoundError(src)
shutil.move(src, gif_path)
print(f"Wrote {rgb_out}")
print(f"Wrote {gif_path}")
PY
  then
    echo "${date}: GIF regeneration failed" >&2
    continue
  fi
  echo "${date}: DONE"
done

echo "Repair/regeneration pass complete."
