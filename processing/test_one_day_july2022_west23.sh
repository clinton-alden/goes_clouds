#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# One-day end-to-end GOES workflow test:
# download -> keep 16Z-23Z only -> ortho -> zarr -> rgb -> gif loop
# Defaults target July 2022 Colorado test domain.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
DOWNLOAD_SCRIPT="${DOWNLOAD_SCRIPT:-${ROOT_DIR}/data_download/download-goes.py}"

GOES="${GOES:-goes16}"
YEAR="${YEAR:-2022}"
MONTH="${MONTH:-7}"
DAY="${DAY:-15}"
DOMAIN="${DOMAIN:-colorado}"

# Colorado bounds: lon [-109, -104], lat [37, 41]
LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

# Base dir should be a high-capacity scratch location.
BASE_DIR="${BASE_DIR:-/glade/derecho/scratch/cdalden/goes_test_west23}"
KEEP_HOURS_START="${KEEP_HOURS_START:-16}"
KEEP_HOURS_END="${KEEP_HOURS_END:-23}"
GOES_HOURS="${GOES_HOURS:-16-23}"

DAY_PAD=$(printf "%02d" "$((10#${DAY}))")
MONTH_INT=$((10#${MONTH}))
DAY_INT=$((10#${DAY}))

channels=(C02 C05 C13)

echo "=== CONFIG ==="
echo "PYTHON_BIN=${PYTHON_BIN}"
echo "DOWNLOAD_SCRIPT=${DOWNLOAD_SCRIPT}"
echo "BASE_DIR=${BASE_DIR}"
echo "GOES=${GOES}"
echo "DATE=${YEAR}-${MONTH_INT}-${DAY_INT}"
echo "DOMAIN=${DOMAIN}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}"
echo "HOURS=${KEEP_HOURS_START}Z-${KEEP_HOURS_END}Z"
echo "GOES_HOURS=${GOES_HOURS}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "ERROR: PYTHON_BIN is not executable: ${PYTHON_BIN}" >&2
  exit 2
fi

if [[ ! -f "${DOWNLOAD_SCRIPT}" ]]; then
  echo "ERROR: download script not found: ${DOWNLOAD_SCRIPT}" >&2
  exit 2
fi

mkdir -p "${BASE_DIR}"
mkdir -p "${SCRIPT_DIR}/plots" "${SCRIPT_DIR}/gifs"

pushd "${SCRIPT_DIR}" >/dev/null

echo "=== STEP 0: RESET DAY DIRECTORY ==="
rm -rf "${BASE_DIR}/${GOES}/${YEAR}/${MONTH_INT}/${DAY_INT}"

echo "=== STEP 1: DOWNLOAD ==="
export GOES_HOURS
for channel in "${channels[@]}"; do
  echo "Downloading ${channel}..."
  "${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" \
    -B "noaa-${GOES}" \
    -Y "${YEAR}" \
    -M "${MONTH_INT}" \
    -D "${DAY_INT}" "${DAY_INT}" \
    -p ABI-L1b-RadC \
    -c "${channel}" \
    -b "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}" \
    -d "${BASE_DIR}"
done

echo "=== STEP 1b: KEEP ONLY 16Z-23Z ==="
DAY_DIR="${BASE_DIR}/${GOES}/${YEAR}/${MONTH_INT}/${DAY_INT}/ABI-L1b-RadC"
if [[ ! -d "${DAY_DIR}" ]]; then
  echo "ERROR: expected day directory missing: ${DAY_DIR}" >&2
  exit 3
fi

for hour in $(seq -w 0 23); do
  h=$((10#${hour}))
  if [[ "${h}" -lt "${KEEP_HOURS_START}" || "${h}" -gt "${KEEP_HOURS_END}" ]]; then
    dir_path="${DAY_DIR}/${hour}"
    if [[ -d "${dir_path}" ]]; then
      chmod -R u+w "${dir_path}" 2>/dev/null || true
      rm -rf "${dir_path}"/* 2>/dev/null || true
    fi
  fi
done

echo "=== STEP 2: ORTHO ==="
"${PYTHON_BIN}" "${SCRIPT_DIR}/batch_ortho.py" "${BASE_DIR}/${GOES}/${YEAR}/${MONTH_INT}/" "${DOMAIN}"

echo "=== STEP 2b: CLEAN EXISTING DAY OUTPUTS ==="
for channel in "${channels[@]}"; do
  zarr_day_path="${BASE_DIR}/${GOES}/${channel}/${GOES}_${channel}_${DOMAIN}_${YEAR}$(printf '%02d' ${MONTH_INT})${DAY_PAD}.zarr"
  if [[ -d "${zarr_day_path}" ]]; then
    rm -rf "${zarr_day_path}"
  fi
done
rgb_day_file="${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${YEAR}$(printf '%02d' ${MONTH_INT})${DAY_PAD}.nc"
if [[ -f "${rgb_day_file}" ]]; then
  rm -f "${rgb_day_file}"
fi

echo "=== STEP 3: ZARR (single month, one downloaded day) ==="
"${PYTHON_BIN}" "${SCRIPT_DIR}/zarr_v2.py" "${BASE_DIR}/" "${YEAR}" "${MONTH_INT}" "${GOES}" "${DOMAIN}"

echo "=== STEP 4: RGB (single day) ==="
echo "=== STEP 3b: DEDUPE ZARR TIMESTAMPS ==="
"${PYTHON_BIN}" - <<PY
import os
import shutil
import xarray as xr

base = "${BASE_DIR}/${GOES}"
date = "${YEAR}$(printf '%02d' ${MONTH_INT})${DAY_PAD}"
domain = "${DOMAIN}"
goes = "${GOES}"

for channel in ("C02", "C05", "C13"):
  zarr_path = f"{base}/{channel}/{goes}_{channel}_{domain}_{date}.zarr"
  if not os.path.isdir(zarr_path):
    print(f"Missing zarr, skipping dedupe: {zarr_path}")
    continue

  ds = xr.open_dataset(zarr_path)
  if "t" not in ds.coords:
    ds.close()
    print(f"No t coord, skipping dedupe: {zarr_path}")
    continue

  t_index = ds.indexes["t"]
  if t_index.is_unique:
    ds.close()
    print(f"t already unique: {zarr_path}")
    continue

  keep = ~t_index.duplicated(keep="first")
  ds_u = ds.isel(t=keep)
  tmp = zarr_path + ".tmp"
  if os.path.isdir(tmp):
    shutil.rmtree(tmp)
  ds_u.to_zarr(tmp, mode="w")
  ds.close()
  shutil.rmtree(zarr_path)
  os.rename(tmp, zarr_path)
  print(f"Deduped t in {zarr_path}")
PY

"${PYTHON_BIN}" - <<PY
import utils

base = "${BASE_DIR}/${GOES}/"
date = "${YEAR}$(printf '%02d' ${MONTH_INT})${DAY_PAD}"
print(f"Generating RGB for {date} from {base}")
utils.goes_rad_to_rgb(base, date, "${GOES}", location="${DOMAIN}")
print("RGB done")
PY

echo "=== STEP 5: PLOT LOOP (GIF) ==="
"${PYTHON_BIN}" - <<PY
import os
import xarray as xr
import sys

sys.path.insert(0, os.path.join("${ROOT_DIR}", "analysis"))
from analysis_utils import make_gif

date = "${YEAR}$(printf '%02d' ${MONTH_INT})${DAY_PAD}"
rgb_path = os.path.join("${BASE_DIR}", "${GOES}", "rgb_composite", f"${GOES}_C02_C05_C13_rgb_${DOMAIN}_{date}.nc")
print(f"Opening {rgb_path}")

if not os.path.exists(rgb_path):
    raise FileNotFoundError(f"RGB output missing: {rgb_path}")

os.makedirs("${SCRIPT_DIR}/plots", exist_ok=True)
os.makedirs("${SCRIPT_DIR}/gifs", exist_ok=True)
os.chdir("${SCRIPT_DIR}")

ds = xr.open_dataset(rgb_path)
make_gif(ds, date, "1600", "2355", subset="${DOMAIN}", mask=False)
ds.close()

print("GIF creation done")
PY

echo "=== DONE ==="
echo "RGB file: ${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${YEAR}$(printf '%02d' ${MONTH_INT})${DAY_PAD}.nc"
echo "GIF dir: ${SCRIPT_DIR}/gifs"

popd >/dev/null
