#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
DOWNLOAD_SCRIPT="${DOWNLOAD_SCRIPT:-${ROOT_DIR}/data_download/download-goes.py}"
ORTHO_SCRIPT="${ORTHO_SCRIPT:-${SCRIPT_DIR}/batch_ortho.py}"
ZARR_SCRIPT="${ZARR_SCRIPT:-${SCRIPT_DIR}/zarr_v2_days.py}"

GOES="${GOES:-goes17}"
DOMAIN="${DOMAIN:-sierra}"
BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/sierra}"
YEAR="${YEAR:-2020}"
MONTH="${MONTH:-1}"
GOES_HOURS="${GOES_HOURS:-14-23}"
GIF_START="${GIF_START:-1400}"
GIF_END="${GIF_END:-2355}"

LON_MIN="${LON_MIN:--121}"
LAT_MIN="${LAT_MIN:-35}"
LON_MAX="${LON_MAX:--118}"
LAT_MAX="${LAT_MAX:-40}"

LOG_DIR="${BASE_DIR}/${GOES}/logs"
GIF_DIR="${BASE_DIR}/${GOES}/gif_loops"
mkdir -p "${LOG_DIR}"
mkdir -p "${GIF_DIR}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "ERROR: PYTHON_BIN is not executable: ${PYTHON_BIN}" >&2
  exit 2
fi

if [[ ! -f "${DOWNLOAD_SCRIPT}" ]]; then
  echo "ERROR: download script not found: ${DOWNLOAD_SCRIPT}" >&2
  exit 2
fi

if [[ ! -f "${ORTHO_SCRIPT}" ]]; then
  echo "ERROR: ortho script not found: ${ORTHO_SCRIPT}" >&2
  exit 2
fi

if [[ ! -f "${ZARR_SCRIPT}" ]]; then
  echo "ERROR: zarr script not found: ${ZARR_SCRIPT}" >&2
  exit 2
fi

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${MONTH})[1])
PY
)"

month_pad="$(printf '%02d' "${MONTH}")"
month_dir="${BASE_DIR}/${GOES}/${YEAR}/${MONTH}"
month_log="${LOG_DIR}/sierra_${GOES}_${YEAR}${month_pad}.log"
channels=(C02 C05 C13)

echo "=== MONTH CONFIG ===" | tee "${month_log}"
echo "BASE_DIR=${BASE_DIR}" | tee -a "${month_log}"
echo "GOES=${GOES} DOMAIN=${DOMAIN}" | tee -a "${month_log}"
echo "YEAR=${YEAR} MONTH=${MONTH}" | tee -a "${month_log}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}" | tee -a "${month_log}"
echo "UTC HOURS=${GOES_HOURS}" | tee -a "${month_log}"
echo "GIF WINDOW=${GIF_START}-${GIF_END}" | tee -a "${month_log}"
echo "DAYS=${days_in_month}" | tee -a "${month_log}"

mkdir -p "${BASE_DIR}"
mkdir -p "${SCRIPT_DIR}/plots" "${SCRIPT_DIR}/gifs"

for day in $(seq 1 "${days_in_month}"); do
  day_pad="$(printf '%02d' "${day}")"
  day_dir="${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${day}/ABI-L1b-RadC"
  date_ymd="${YEAR}${month_pad}${day_pad}"
  day_log="${LOG_DIR}/sierra_${GOES}_${date_ymd}.log"
  rgb_day_file="${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${date_ymd}.nc"
  gif_day_file="${GIF_DIR}/goes_rgb_${DOMAIN}_${date_ymd}_${GIF_START}_${GIF_END}.gif"

  echo | tee -a "${month_log}"
  echo "=== DAY ${YEAR}-${month_pad}-${day_pad} ===" | tee -a "${month_log}" "${day_log}"

  rm -rf "${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${day}" >>"${day_log}" 2>&1 || true
  for channel in "${channels[@]}"; do
    rm -rf "${BASE_DIR}/${GOES}/${channel}/${GOES}_${channel}_${DOMAIN}_${date_ymd}.zarr" >>"${day_log}" 2>&1 || true
  done
  rm -f "${rgb_day_file}" "${gif_day_file}" >>"${day_log}" 2>&1 || true

  echo "Downloading channels for ${date_ymd}" | tee -a "${month_log}" "${day_log}"
  for channel in "${channels[@]}"; do
    GOES_HOURS="${GOES_HOURS}" \
      "${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" \
        -B "noaa-${GOES}" \
        -Y "${YEAR}" \
        -M "${MONTH}" \
        -D "${day}" "${day}" \
        -p ABI-L1b-RadC \
        -c "${channel}" \
        -b "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}" \
        -d "${BASE_DIR}" >>"${day_log}" 2>&1
  done

  day_nc_count=$(find "${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${day}" -type f -name '*.nc' 2>/dev/null | wc -l || true)
  if [[ "${day_nc_count}" -eq 0 ]]; then
    echo "ERROR: no downloaded NetCDF files found for ${date_ymd}" | tee -a "${month_log}" "${day_log}"
    exit 1
  fi

  echo "Orthorectifying ${date_ymd}" | tee -a "${month_log}" "${day_log}"
  "${PYTHON_BIN}" "${ORTHO_SCRIPT}" "${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/" "${DOMAIN}" >>"${day_log}" 2>&1

  echo "Building zarr for ${date_ymd}" | tee -a "${month_log}" "${day_log}"
  "${PYTHON_BIN}" "${ZARR_SCRIPT}" "${BASE_DIR}/" "${YEAR}" "${MONTH}" "${day}" "${day}" "${GOES}" "${DOMAIN}" >>"${day_log}" 2>&1

  echo "Deduping zarr timestamps for ${date_ymd}" | tee -a "${month_log}" "${day_log}"
  "${PYTHON_BIN}" - <<PY >>"${day_log}" 2>&1
import os
import shutil
import xarray as xr

base = "${BASE_DIR}/${GOES}"
date = "${date_ymd}"
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

  echo "Building RGB for ${date_ymd}" | tee -a "${month_log}" "${day_log}"
  "${PYTHON_BIN}" - <<PY >>"${day_log}" 2>&1
import utils

base = "${BASE_DIR}/${GOES}/"
date = "${date_ymd}"
print(f"Generating RGB for {date} from {base}")
utils.goes_rad_to_rgb(base, date, "${GOES}", location="${DOMAIN}")
print("RGB done")
PY

  echo "Building GIF for ${date_ymd}" | tee -a "${month_log}" "${day_log}"
  "${PYTHON_BIN}" - <<PY >>"${day_log}" 2>&1
import os
import shutil
import sys
import xarray as xr

sys.path.insert(0, os.path.join("${ROOT_DIR}", "analysis"))
from analysis_utils import make_gif

date = "${date_ymd}"
rgb_path = "${rgb_day_file}"
gif_out = "${gif_day_file}"

if not os.path.exists(rgb_path):
    raise FileNotFoundError(f"RGB output missing: {rgb_path}")

os.makedirs("${SCRIPT_DIR}/plots", exist_ok=True)
os.makedirs("${SCRIPT_DIR}/gifs", exist_ok=True)
os.makedirs("${GIF_DIR}", exist_ok=True)
os.chdir("${SCRIPT_DIR}")

with xr.open_dataset(rgb_path) as ds:
    make_gif(ds, date, "${GIF_START}", "${GIF_END}", subset="${DOMAIN}", mask=False)

src = os.path.join("${SCRIPT_DIR}", "gifs", f"goes_rgb_${DOMAIN}_{date}_${GIF_START}_${GIF_END}.gif")
shutil.copy2(src, gif_out)
print(f"GIF written to {gif_out}")
PY

  echo "DONE ${date_ymd}" | tee -a "${month_log}" "${day_log}"
done

echo "=== DONE ${YEAR}-${month_pad} ===" | tee -a "${month_log}"
