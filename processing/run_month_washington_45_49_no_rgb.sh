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
DOMAIN="${DOMAIN:-washington}"
BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/washington}"
YEAR="${YEAR:-2022}"
MONTH="${MONTH:-5}"
GOES_HOURS="${GOES_HOURS:-14-23}"

LON_MIN="${LON_MIN:--125}"
LAT_MIN="${LAT_MIN:-45}"
LON_MAX="${LON_MAX:--120}"
LAT_MAX="${LAT_MAX:-49}"

LOG_DIR="${BASE_DIR}/${GOES}/logs"
mkdir -p "${LOG_DIR}"

if [[ ! -f "${DOWNLOAD_SCRIPT}" ]]; then
  echo "Download script not found: ${DOWNLOAD_SCRIPT}" >&2
  exit 2
fi

if [[ ! -f "${ORTHO_SCRIPT}" ]]; then
  echo "Ortho script not found: ${ORTHO_SCRIPT}" >&2
  exit 2
fi

if [[ ! -f "${ZARR_SCRIPT}" ]]; then
  echo "Zarr script not found: ${ZARR_SCRIPT}" >&2
  exit 2
fi

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${MONTH})[1])
PY
)"

echo "=== MONTH CONFIG ==="
echo "BASE_DIR=${BASE_DIR}"
echo "GOES=${GOES} DOMAIN=${DOMAIN}"
echo "YEAR=${YEAR} MONTH=${MONTH}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}"
echo "UTC HOURS=${GOES_HOURS}"
echo "DAYS=${days_in_month}"

for day in $(seq 1 "${days_in_month}"); do
  date_str="$(printf "%04d%02d%02d" "${YEAR}" "${MONTH}" "${day}")"
  day_dir="${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${day}"
  log_path="${LOG_DIR}/${date_str}.log"

  echo "=== ${date_str}: START ===" | tee -a "${log_path}"
  echo "DAY_DIR=${day_dir}" | tee -a "${log_path}"

  rm -rf "${day_dir}"
  mkdir -p "${day_dir}"

  download_failed=0
  for channel in C02 C05 C13; do
    echo "${date_str}: download ${channel}" | tee -a "${log_path}"
    if ! GOES_HOURS="${GOES_HOURS}" \
      "${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" \
        -B "noaa-${GOES}" \
        -Y "${YEAR}" \
        -M "${MONTH}" \
        -D "${day}" "${day}" \
        -p ABI-L1b-RadC \
        -c "${channel}" \
        -b "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}" \
        -d "${BASE_DIR}" >>"${log_path}" 2>&1; then
      download_failed=1
      echo "${date_str}: download failed for ${channel}" | tee -a "${log_path}"
      break
    fi
  done

  if [[ "${download_failed}" -ne 0 ]]; then
    echo "${date_str}: FAILED during download, continuing to next day" | tee -a "${log_path}"
    continue
  fi

  raw_count=$(find "${day_dir}" -type f -name '*.nc' ! -name '*_ortho.nc' | wc -l || true)
  if [[ "${raw_count}" -eq 0 ]]; then
    echo "${date_str}: no raw NetCDF files downloaded" | tee -a "${log_path}"
    continue
  fi

  echo "${date_str}: orthorectify ${raw_count} files" | tee -a "${log_path}"
  if ! "${PYTHON_BIN}" "${ORTHO_SCRIPT}" "${day_dir}/" "${DOMAIN}" >>"${log_path}" 2>&1; then
    echo "${date_str}: ortho failed, continuing to next day" | tee -a "${log_path}"
    continue
  fi

  ortho_count=$(find "${day_dir}" -type f -name '*_ortho.nc' | wc -l || true)
  if [[ "${ortho_count}" -eq 0 ]]; then
    echo "${date_str}: no orthorectified files found" | tee -a "${log_path}"
    continue
  fi

  echo "${date_str}: create zarr outputs" | tee -a "${log_path}"
  if ! "${PYTHON_BIN}" "${ZARR_SCRIPT}" "${BASE_DIR}/" "${YEAR}" "${MONTH}" "${day}" "${day}" "${GOES}" "${DOMAIN}" >>"${log_path}" 2>&1; then
    echo "${date_str}: zarr failed, continuing to next day" | tee -a "${log_path}"
    continue
  fi

  echo "${date_str}: DONE" | tee -a "${log_path}"
done
