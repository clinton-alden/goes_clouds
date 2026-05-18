#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python}"
DOWNLOAD_SCRIPT="${DOWNLOAD_SCRIPT:-${SCRIPT_DIR}/download-goes.py}"
ORTHO_SCRIPT="${ORTHO_SCRIPT:-${SCRIPT_DIR}/batch_ortho.py}"
ZARR_SCRIPT="${ZARR_SCRIPT:-${SCRIPT_DIR}/zarr_v2_days.py}"

GOES="${GOES:-goes16}"
DOMAIN="${DOMAIN:-colorado}"
BASE_DIR="${BASE_DIR:?BASE_DIR is required, e.g. /path/to/scratch/colorado}"
YEAR="${YEAR:?YEAR is required}"
MONTH="${MONTH:?MONTH is required}"
GOES_HOURS="${GOES_HOURS:-14-23}"

LON_MIN="${LON_MIN:?LON_MIN is required}"
LAT_MIN="${LAT_MIN:?LAT_MIN is required}"
LON_MAX="${LON_MAX:?LON_MAX is required}"
LAT_MAX="${LAT_MAX:?LAT_MAX is required}"

if [[ -z "${OPENTOPOGRAPHY_API_KEY:-}" ]]; then
  echo "ERROR: OPENTOPOGRAPHY_API_KEY is required for orthorectification DEM downloads." >&2
  echo "Create your own OpenTopography key and export OPENTOPOGRAPHY_API_KEY before running." >&2
  exit 2
fi

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${MONTH})[1])
PY
)"

month_pad="$(printf '%02d' "${MONTH}")"
log_dir="${BASE_DIR}/${GOES}/logs"
mkdir -p "${log_dir}" "${BASE_DIR}"
month_log="${log_dir}/${DOMAIN}_${GOES}_${YEAR}${month_pad}.log"
channels=(C02 C05 C13)

echo "=== MONTH CONFIG ===" | tee "${month_log}"
echo "BASE_DIR=${BASE_DIR}" | tee -a "${month_log}"
echo "GOES=${GOES} DOMAIN=${DOMAIN}" | tee -a "${month_log}"
echo "YEAR=${YEAR} MONTH=${MONTH}" | tee -a "${month_log}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}" | tee -a "${month_log}"
echo "UTC HOURS=${GOES_HOURS}" | tee -a "${month_log}"

for day in $(seq 1 "${days_in_month}"); do
  day_pad="$(printf '%02d' "${day}")"
  date_ymd="${YEAR}${month_pad}${day_pad}"
  day_log="${log_dir}/${DOMAIN}_${GOES}_${date_ymd}.log"
  rgb_day_file="${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${date_ymd}.nc"

  echo | tee -a "${month_log}"
  echo "=== DAY ${YEAR}-${month_pad}-${day_pad} ===" | tee -a "${month_log}" "${day_log}"

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
    if "t" not in ds.coords or ds.indexes["t"].is_unique:
        ds.close()
        continue
    keep = ~ds.indexes["t"].duplicated(keep="first")
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
  PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}" "${PYTHON_BIN}" - <<PY >>"${day_log}" 2>&1
import utils
utils.goes_rad_to_rgb("${BASE_DIR}/${GOES}/", "${date_ymd}", "${GOES}", location="${DOMAIN}")
PY

  if [[ ! -f "${rgb_day_file}" ]]; then
    echo "ERROR: expected RGB output missing: ${rgb_day_file}" | tee -a "${month_log}" "${day_log}"
    exit 1
  fi

  echo "DONE ${date_ymd}" | tee -a "${month_log}" "${day_log}"
done

echo "=== DONE ${YEAR}-${month_pad} ===" | tee -a "${month_log}"
