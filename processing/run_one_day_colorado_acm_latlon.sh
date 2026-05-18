#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
DOWNLOAD_SCRIPT="${DOWNLOAD_SCRIPT:-${ROOT_DIR}/data_download/download-goes.py}"
CLIP_SCRIPT="${CLIP_SCRIPT:-${SCRIPT_DIR}/clip_acm_to_latlon.py}"
DAILY_SCRIPT="${DAILY_SCRIPT:-${SCRIPT_DIR}/acm_daily_latlon_nc.py}"
PLOT_SCRIPT="${PLOT_SCRIPT:-${SCRIPT_DIR}/plot_acm_daily_gif.py}"

GOES="${GOES:-goes16}"
DOMAIN="${DOMAIN:-colorado}"
BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado_acm}"
YEAR="${YEAR:-2021}"
MONTH="${MONTH:-10}"
DAY="${DAY:-1}"
GOES_HOURS="${GOES_HOURS:-14-23}"
MASK_VAR="${MASK_VAR:-AUTO}"

LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

DATE="$(printf "%04d%02d%02d" "${YEAR}" "${MONTH}" "${DAY}")"
DAY_DIR="${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${DAY}"
RAW_DIR="${DAY_DIR}/ABI-L2-ACMC"
CLIP_DIR="${RAW_DIR}/clipped"
LOG_DIR="${BASE_DIR}/${GOES}/logs_latlon"
GIF_DIR="${BASE_DIR}/${GOES}/gif_loops_latlon"
TMP_DIR="${BASE_DIR}/${GOES}/gif_tmp/${DATE}_acm"
DAILY_PATH="${BASE_DIR}/${GOES}/daily_nc_latlon/${GOES}_acm_${DOMAIN}_${DATE}.nc"
GIF_PATH="${GIF_DIR}/${GOES}_acm_${DOMAIN}_${DATE}.gif"
LOG_PATH="${LOG_DIR}/${DATE}_acm_latlon.log"

mkdir -p "${LOG_DIR}" "${GIF_DIR}" "${TMP_DIR}"
: > "${LOG_PATH}"

echo "=== ONE DAY ACM LATLON CONFIG ===" | tee -a "${LOG_PATH}"
echo "DATE=${DATE}" | tee -a "${LOG_PATH}"
echo "BASE_DIR=${BASE_DIR}" | tee -a "${LOG_PATH}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}" | tee -a "${LOG_PATH}"
echo "UTC HOURS=${GOES_HOURS}" | tee -a "${LOG_PATH}"

mkdir -p "${DAY_DIR}"

echo "${DATE}: download ABI-L2-ACMC" | tee -a "${LOG_PATH}"
GOES_HOURS="${GOES_HOURS}" \
"${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" \
  -B "noaa-${GOES}" \
  -Y "${YEAR}" \
  -M "${MONTH}" \
  -D "${DAY}" "${DAY}" \
  -p ABI-L2-ACMC \
  -c C02 \
  -b "${LAT_MIN}" "${LAT_MAX}" "${LON_MIN}" "${LON_MAX}" \
  -d "${BASE_DIR}" >>"${LOG_PATH}" 2>&1

mapfile -t raw_files < <(find "${RAW_DIR}" -type f -name '*.nc' ! -path '*/clipped/*' | sort)
if [[ "${#raw_files[@]}" -eq 0 ]]; then
  echo "${DATE}: no raw ACM files downloaded" | tee -a "${LOG_PATH}"
  exit 1
fi

mkdir -p "${CLIP_DIR}"
echo "${DATE}: clip ${#raw_files[@]} files to lat/lon domain" | tee -a "${LOG_PATH}"
for raw_file in "${raw_files[@]}"; do
  out_file="${CLIP_DIR}/$(basename "${raw_file%.nc}")_clip.nc"
  if [[ -f "${out_file}" ]]; then
    continue
  fi
  "${PYTHON_BIN}" "${CLIP_SCRIPT}" \
    "${raw_file}" \
    "${out_file}" \
    --lon-min "${LON_MIN}" \
    --lon-max "${LON_MAX}" \
    --lat-min "${LAT_MIN}" \
    --lat-max "${LAT_MAX}" \
    --mask-var "${MASK_VAR}" >>"${LOG_PATH}" 2>&1
done

echo "${DATE}: build daily file" | tee -a "${LOG_PATH}"
"${PYTHON_BIN}" "${DAILY_SCRIPT}" \
  --base-dir "${BASE_DIR}" \
  --year "${YEAR}" \
  --month "${MONTH}" \
  --goes "${GOES}" \
  --domain "${DOMAIN}" \
  --mask-var "${MASK_VAR}" >>"${LOG_PATH}" 2>&1

if [[ ! -f "${DAILY_PATH}" ]]; then
  echo "${DATE}: expected daily file missing: ${DAILY_PATH}" | tee -a "${LOG_PATH}"
  exit 1
fi

echo "${DATE}: render gif" | tee -a "${LOG_PATH}"
"${PYTHON_BIN}" "${PLOT_SCRIPT}" \
  --input "${DAILY_PATH}" \
  --output "${GIF_PATH}" \
  --domain "${DOMAIN}" \
  --mask-var "${MASK_VAR}" \
  --tmp-dir "${TMP_DIR}" >>"${LOG_PATH}" 2>&1

rmdir "${TMP_DIR}" 2>/dev/null || true

echo "${DATE}: DONE" | tee -a "${LOG_PATH}"
echo "DAILY=${DAILY_PATH}" | tee -a "${LOG_PATH}"
echo "GIF=${GIF_PATH}" | tee -a "${LOG_PATH}"
