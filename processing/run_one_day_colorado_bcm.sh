#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
DOWNLOAD_SCRIPT="${DOWNLOAD_SCRIPT:-${ROOT_DIR}/data_download/download-goes.py}"
ORTHO_SCRIPT="${ORTHO_SCRIPT:-${SCRIPT_DIR}/batch_ortho_acm.py}"
PLOT_SCRIPT="${PLOT_SCRIPT:-${SCRIPT_DIR}/plot_acm_gif.py}"

GOES="${GOES:-goes16}"
DOMAIN="${DOMAIN:-colorado}"
BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado}"
YEAR="${YEAR:-2022}"
MONTH="${MONTH:-7}"
DAY="${DAY:-2}"
START_HOUR="${START_HOUR:-14}"
END_HOUR="${END_HOUR:-23}"
MASK_VAR="${MASK_VAR:-BCM}"

LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

DATE=$(printf "%04d%02d%02d" "${YEAR}" "${MONTH}" "${DAY}")
DAY_DIR="${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${DAY}"
RAW_DIR="${DAY_DIR}/ABI-L2-ACMC"
LOG_DIR="${BASE_DIR}/${GOES}/logs"
GIF_DIR="${BASE_DIR}/${GOES}/gif_loops"
TMP_DIR="${BASE_DIR}/${GOES}/gif_tmp/${DATE}_bcm"
LOG_PATH="${LOG_DIR}/${DATE}_bcm.log"
GIF_PATH="${GIF_DIR}/goes_${MASK_VAR,,}_${DOMAIN}_${DATE}_1400_2355.gif"

mkdir -p "${LOG_DIR}" "${GIF_DIR}" "${TMP_DIR}"
: > "${LOG_PATH}"

echo "=== BCM ONE-DAY CONFIG ===" | tee -a "${LOG_PATH}"
echo "BASE_DIR=${BASE_DIR}" | tee -a "${LOG_PATH}"
echo "DATE=${DATE}" | tee -a "${LOG_PATH}"
echo "GOES=${GOES} DOMAIN=${DOMAIN}" | tee -a "${LOG_PATH}"
echo "MASK_VAR=${MASK_VAR}" | tee -a "${LOG_PATH}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}" | tee -a "${LOG_PATH}"
echo "UTC HOURS=${START_HOUR}-${END_HOUR}" | tee -a "${LOG_PATH}"

rm -rf "${DAY_DIR}"
mkdir -p "${DAY_DIR}"

echo "${DATE}: download ABI-L2-ACMC" | tee -a "${LOG_PATH}"
GOES_HOURS="${START_HOUR}-${END_HOUR}" \
"${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" \
  -B "noaa-${GOES}" \
  -Y "${YEAR}" \
  -M "${MONTH}" \
  -D "${DAY}" "${DAY}" \
  -p ABI-L2-ACMC \
  -c C02 \
  -b "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}" \
  -d "${BASE_DIR}" >>"${LOG_PATH}" 2>&1

raw_count=$(find "${RAW_DIR}" -type f -name '*.nc' | wc -l || true)
if [[ "${raw_count}" -eq 0 ]]; then
  echo "${DATE}: no BCM/ACM NetCDF files downloaded" | tee -a "${LOG_PATH}"
  exit 1
fi

echo "${DATE}: orthorectify ${raw_count} files for ${MASK_VAR}" | tee -a "${LOG_PATH}"
GOES_DATA_VARS="${MASK_VAR}" ORTHO_ACM_MAX_WORKERS="${ORTHO_ACM_MAX_WORKERS:-4}" \
  "${PYTHON_BIN}" "${ORTHO_SCRIPT}" "${DAY_DIR}/" "${DOMAIN}" >>"${LOG_PATH}" 2>&1

ortho_count=$(find "${RAW_DIR}" -type f -name '*_ortho.nc' | wc -l || true)
if [[ "${ortho_count}" -eq 0 ]]; then
  echo "${DATE}: no orthorectified files found" | tee -a "${LOG_PATH}"
  exit 1
fi

echo "${DATE}: render ${MASK_VAR} gif" | tee -a "${LOG_PATH}"
"${PYTHON_BIN}" "${PLOT_SCRIPT}" \
  --input-dir "${DAY_DIR}" \
  --output "${GIF_PATH}" \
  --date "${DATE}" \
  --start-hour "${START_HOUR}" \
  --end-hour "${END_HOUR}" \
  --domain "${DOMAIN}" \
  --tmp-dir "${TMP_DIR}" >>"${LOG_PATH}" 2>&1

rmdir "${TMP_DIR}" 2>/dev/null || true

echo "${DATE}: DONE" | tee -a "${LOG_PATH}"
echo "GIF: ${GIF_PATH}" | tee -a "${LOG_PATH}"
