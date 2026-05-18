#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
DOWNLOAD_SCRIPT="${DOWNLOAD_SCRIPT:-${ROOT_DIR}/data_download/download-goes.py}"
ORTHO_SCRIPT="${ORTHO_SCRIPT:-${SCRIPT_DIR}/batch_ortho_acm.py}"
DAILY_NC_SCRIPT="${DAILY_NC_SCRIPT:-${SCRIPT_DIR}/acm_daily_nc.py}"

GOES="${GOES:-goes16}"
DOMAIN="${DOMAIN:-colorado}"
BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado_acm}"
YEAR="${YEAR:-2021}"
MONTH="${MONTH:-10}"
GOES_HOURS="${GOES_HOURS:-00-23}"
ORTHO_ACM_MAX_WORKERS="${ORTHO_ACM_MAX_WORKERS:-4}"
GOES_DATA_VARS="${GOES_DATA_VARS:-ACM}"

LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

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

if [[ ! -f "${DAILY_NC_SCRIPT}" ]]; then
  echo "Daily NetCDF script not found: ${DAILY_NC_SCRIPT}" >&2
  exit 2
fi

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${MONTH})[1])
PY
)"

month_log="${LOG_DIR}/acm_${YEAR}$(printf '%02d' "${MONTH}").log"

echo "=== MONTH CONFIG ===" | tee "${month_log}"
echo "BASE_DIR=${BASE_DIR}" | tee -a "${month_log}"
echo "GOES=${GOES} DOMAIN=${DOMAIN}" | tee -a "${month_log}"
echo "YEAR=${YEAR} MONTH=${MONTH}" | tee -a "${month_log}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}" | tee -a "${month_log}"
echo "UTC HOURS=${GOES_HOURS}" | tee -a "${month_log}"
echo "DAYS=${days_in_month}" | tee -a "${month_log}"

for day in $(seq 1 "${days_in_month}"); do
  date_str="$(printf "%04d%02d%02d" "${YEAR}" "${MONTH}" "${day}")"
  day_dir="${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${day}"
  raw_dir="${day_dir}/ABI-L2-ACMC"

  echo "=== ${date_str}: START ===" | tee -a "${month_log}"
  echo "DAY_DIR=${day_dir}" | tee -a "${month_log}"

  mkdir -p "${day_dir}"

  echo "${date_str}: download ABI-L2-ACMC" | tee -a "${month_log}"
  if ! GOES_HOURS="${GOES_HOURS}" \
    "${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" \
      -B "noaa-${GOES}" \
      -Y "${YEAR}" \
      -M "${MONTH}" \
      -D "${day}" "${day}" \
      -p ABI-L2-ACMC \
      -c C02 \
      -b "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}" \
      -d "${BASE_DIR}" >>"${month_log}" 2>&1; then
    echo "${date_str}: download failed, continuing to next day" | tee -a "${month_log}"
    continue
  fi

  raw_count=$(find "${raw_dir}" -type f -name '*.nc' ! -name '*_ortho.nc' | wc -l || true)
  if [[ "${raw_count}" -eq 0 ]]; then
    echo "${date_str}: no raw ACM files downloaded" | tee -a "${month_log}"
    continue
  fi

  echo "${date_str}: orthorectify ${raw_count} files" | tee -a "${month_log}"
  if ! GOES_DATA_VARS="${GOES_DATA_VARS}" ORTHO_ACM_MAX_WORKERS="${ORTHO_ACM_MAX_WORKERS}" \
    "${PYTHON_BIN}" "${ORTHO_SCRIPT}" "${day_dir}/" "${DOMAIN}" >>"${month_log}" 2>&1; then
    echo "${date_str}: ortho failed, continuing to next day" | tee -a "${month_log}"
    continue
  fi

  ortho_count=$(find "${raw_dir}" -type f -name '*_ortho.nc' | wc -l || true)
  if [[ "${ortho_count}" -eq 0 ]]; then
    echo "${date_str}: no orthorectified files found" | tee -a "${month_log}"
    continue
  fi

  echo "${date_str}: DONE" | tee -a "${month_log}"
done

echo "=== MONTH PACKAGING START ===" | tee -a "${month_log}"
"${PYTHON_BIN}" "${DAILY_NC_SCRIPT}" \
  --base-dir "${BASE_DIR}" \
  --year "${YEAR}" \
  --month "${MONTH}" \
  --goes "${GOES}" \
  --domain "${DOMAIN}" >>"${month_log}" 2>&1
echo "=== MONTH PACKAGING DONE ===" | tee -a "${month_log}"
