#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
APPLY_SCRIPT="${APPLY_SCRIPT:-${SCRIPT_DIR}/apply_tempbin_thresholds_colorado.py}"

YEAR="${YEAR:?YEAR is required}"
MONTH="${MONTH:?MONTH is required}"

RGB_DIR="${RGB_DIR:-/glade/u/home/cdalden/scratch/colorado/goes16/rgb_composite}"
THRESHOLD_CSV="${THRESHOLD_CSV:-/glade/u/home/cdalden/goes_work/analysis/output_12_rgb_threshold_transfer/gothic_temp_bin_rgb_thresholds_10c.csv}"
ERA5_DIR="${ERA5_DIR:-/glade/derecho/scratch/cdalden/tmp/colorado/era5_land/t2m_hourly}"
OUTPUT_BASE="${OUTPUT_BASE:-/glade/u/home/cdalden/scratch/colorado_rgbmasks}"
MASK_DIR="${MASK_DIR:-${OUTPUT_BASE}/cloud_mask_tempbin_10c}"
GIF_DIR="${GIF_DIR:-${OUTPUT_BASE}/gif_loops_tempbin_10c}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
FRAME_DURATION="${FRAME_DURATION:-0.15}"

mkdir -p "${MASK_DIR}" "${GIF_DIR}" "${LOG_DIR}"

month_pad="$(printf '%02d' "${MONTH}")"
month_num="$((10#${MONTH}))"
MONTH_LOG="${LOG_DIR}/rgbmask_${YEAR}${month_pad}.log"

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${month_num})[1])
PY
)"

echo "=== MONTH RGB MASK CONFIG ===" | tee "${MONTH_LOG}"
echo "YEAR=${YEAR} MONTH=${MONTH}" | tee -a "${MONTH_LOG}"
echo "RGB_DIR=${RGB_DIR}" | tee -a "${MONTH_LOG}"
echo "ERA5_DIR=${ERA5_DIR}" | tee -a "${MONTH_LOG}"
echo "MASK_DIR=${MASK_DIR}" | tee -a "${MONTH_LOG}"
echo "GIF_DIR=${GIF_DIR}" | tee -a "${MONTH_LOG}"

for day in $(seq 1 "${days_in_month}"); do
  date_str="$(printf '%04d%02d%02d' "${YEAR}" "${month_num}" "${day}")"
  rgb_file="${RGB_DIR}/goes16_C02_C05_C13_rgb_colorado_${date_str}.nc"

  if [[ ! -f "${rgb_file}" ]]; then
    echo "[skip] ${date_str}: missing RGB file ${rgb_file}" | tee -a "${MONTH_LOG}"
    continue
  fi

  echo "[build] ${date_str}" | tee -a "${MONTH_LOG}"
  if ! "${PYTHON_BIN}" "${APPLY_SCRIPT}" \
    --rgb-file "${rgb_file}" \
    --threshold-csv "${THRESHOLD_CSV}" \
    --era5-dir "${ERA5_DIR}" \
    --mask-dir "${MASK_DIR}" \
    --gif-dir "${GIF_DIR}" \
    --frame-duration "${FRAME_DURATION}" \
    --overwrite \
    --skip-download >>"${MONTH_LOG}" 2>&1; then
    echo "[fail] ${date_str}" | tee -a "${MONTH_LOG}"
    continue
  fi
done

echo "=== MONTH RGB MASK DONE ===" | tee -a "${MONTH_LOG}"
