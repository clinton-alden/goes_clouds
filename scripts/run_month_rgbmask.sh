#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python}"
APPLY_SCRIPT="${APPLY_SCRIPT:-${SCRIPT_DIR}/apply_tempbin_thresholds.py}"

YEAR="${YEAR:?YEAR is required}"
MONTH="${MONTH:?MONTH is required}"
GOES="${GOES:-goes16}"
DOMAIN="${DOMAIN:-colorado}"
RGB_DIR="${RGB_DIR:?RGB_DIR is required}"
THRESHOLD_CSV="${THRESHOLD_CSV:-${REPO_DIR}/thresholds/gothic_temp_bin_rgb_thresholds_10c.csv}"
ERA5_DIR="${ERA5_DIR:?ERA5_DIR is required}"
OUTPUT_BASE="${OUTPUT_BASE:?OUTPUT_BASE is required}"
MASK_DIR="${MASK_DIR:-${OUTPUT_BASE}/cloud_mask_tempbin_10c}"
GIF_DIR="${GIF_DIR:-${OUTPUT_BASE}/gif_loops_tempbin_10c}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
FRAME_DURATION="${FRAME_DURATION:-0.15}"
SKIP_DOWNLOAD="${SKIP_DOWNLOAD:-0}"

mkdir -p "${MASK_DIR}" "${GIF_DIR}" "${LOG_DIR}"

month_pad="$(printf '%02d' "${MONTH}")"
month_num="$((10#${MONTH}))"
month_log="${LOG_DIR}/rgbmask_${DOMAIN}_${GOES}_${YEAR}${month_pad}.log"

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${month_num})[1])
PY
)"

echo "=== MONTH RGB MASK CONFIG ===" | tee "${month_log}"
echo "YEAR=${YEAR} MONTH=${MONTH}" | tee -a "${month_log}"
echo "RGB_DIR=${RGB_DIR}" | tee -a "${month_log}"
echo "ERA5_DIR=${ERA5_DIR}" | tee -a "${month_log}"
echo "MASK_DIR=${MASK_DIR}" | tee -a "${month_log}"
echo "GIF_DIR=${GIF_DIR}" | tee -a "${month_log}"

for day in $(seq 1 "${days_in_month}"); do
  date_str="$(printf '%04d%02d%02d' "${YEAR}" "${month_num}" "${day}")"
  rgb_file="${RGB_DIR}/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${date_str}.nc"

  if [[ ! -f "${rgb_file}" ]]; then
    echo "[skip] ${date_str}: missing RGB file ${rgb_file}" | tee -a "${month_log}"
    continue
  fi

  args=(
    "${APPLY_SCRIPT}"
    --rgb-file "${rgb_file}"
    --threshold-csv "${THRESHOLD_CSV}"
    --era5-dir "${ERA5_DIR}"
    --mask-dir "${MASK_DIR}"
    --gif-dir "${GIF_DIR}"
    --frame-duration "${FRAME_DURATION}"
    --domain "${DOMAIN}"
    --overwrite
  )
  if [[ "${SKIP_DOWNLOAD}" == "1" ]]; then
    args+=(--skip-download)
  fi

  echo "[build] ${date_str}" | tee -a "${month_log}"
  if ! "${PYTHON_BIN}" "${args[@]}" >>"${month_log}" 2>&1; then
    echo "[fail] ${date_str}" | tee -a "${month_log}"
    continue
  fi
done

echo "=== MONTH RGB MASK DONE ===" | tee -a "${month_log}"
