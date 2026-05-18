#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Generate missing GIF loops for existing Colorado RGB composites.
# Only creates GIFs for RGB files that do not already have a matching GIF.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado}"
GOES="${GOES:-goes16}"
DOMAIN="${DOMAIN:-colorado}"
GIF_START="${GIF_START:-1400}"
GIF_END="${GIF_END:-2355}"
DATE_GLOB="${DATE_GLOB:-*}"

RGB_DIR="${BASE_DIR}/${GOES}/rgb_composite"
GIF_DIR="${BASE_DIR}/${GOES}/gif_loops"
TMP_ROOT="${BASE_DIR}/${GOES}/gif_tmp"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "ERROR: PYTHON_BIN not executable: ${PYTHON_BIN}" >&2
  exit 2
fi

if [[ ! -d "${RGB_DIR}" ]]; then
  echo "ERROR: RGB directory not found: ${RGB_DIR}" >&2
  exit 2
fi

mkdir -p "${GIF_DIR}" "${TMP_ROOT}"

shopt -s nullglob
rgb_files=(${RGB_DIR}/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${DATE_GLOB}.nc)
shopt -u nullglob

if [[ "${#rgb_files[@]}" -eq 0 ]]; then
  echo "No RGB files matched: ${RGB_DIR}/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${DATE_GLOB}.nc"
  exit 0
fi

echo "Found ${#rgb_files[@]} RGB files to inspect"

for rgb_path in "${rgb_files[@]}"; do
  rgb_name="$(basename "${rgb_path}" .nc)"
  suffix="${rgb_name#${GOES}_C02_C05_C13_rgb_${DOMAIN}_}"

  date="${suffix}"
  tag=""
  if [[ "${suffix}" == *_latcorrected ]]; then
    date="${suffix%_latcorrected}"
    tag="_latcorrected"
  fi

  gif_path="${GIF_DIR}/goes_rgb_${DOMAIN}_${date}_${GIF_START}_${GIF_END}${tag}.gif"
  if [[ -f "${gif_path}" ]]; then
    echo "${date}${tag}: GIF already exists, skipping"
    continue
  fi

  work_dir="$(mktemp -d "${TMP_ROOT}/${date}${tag}.XXXXXX")"
  echo "${date}${tag}: generating ${gif_path}"

  if ! "${PYTHON_BIN}" - <<PY
import os
import shutil
import sys
import xarray as xr

root_dir = "${ROOT_DIR}"
rgb_path = "${rgb_path}"
gif_path = "${gif_path}"
date = "${date}"
domain = "${DOMAIN}"
gif_start = "${GIF_START}"
gif_end = "${GIF_END}"
work_dir = "${work_dir}"

sys.path.insert(0, os.path.join(root_dir, "analysis"))
from analysis_utils import make_gif

plots_dir = os.path.join(work_dir, "plots")
gifs_dir = os.path.join(work_dir, "gifs")
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(gifs_dir, exist_ok=True)
os.chdir(work_dir)

with xr.open_dataset(rgb_path) as ds:
    make_gif(ds, date, gif_start, gif_end, subset=domain, mask=False)

src = os.path.join(gifs_dir, f"goes_rgb_{domain}_{date}_{gif_start}_{gif_end}.gif")
if not os.path.exists(src):
    raise FileNotFoundError(src)
shutil.move(src, gif_path)
print(f"Wrote {gif_path}")
PY
  then
    echo "${date}${tag}: GIF generation failed" >&2
    rm -rf "${work_dir}"
    continue
  fi

  rm -rf "${work_dir}"
  echo "${date}${tag}: DONE"
done

echo "Missing-GIF generation pass complete."
