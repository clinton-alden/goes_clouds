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
THRESHOLD_CSV="${THRESHOLD_CSV:-${REPO_DIR}/thresholds/gothic_vintage_rgb_tree_rules_5c_sw_kt050_090.csv}"
ERA5_DIR="${ERA5_DIR:?ERA5_DIR is required}"
OUTPUT_BASE="${OUTPUT_BASE:?OUTPUT_BASE is required}"
MASK_DIR="${MASK_DIR:-${OUTPUT_BASE}/vintage_mask}"
GIF_DIR="${GIF_DIR:-${OUTPUT_BASE}/vintage_gif_loops}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
FRAME_DURATION="${FRAME_DURATION:-0.15}"
SKIP_DOWNLOAD="${SKIP_DOWNLOAD:-0}"
KEEP_MASK_DIAGNOSTICS="${KEEP_MASK_DIAGNOSTICS:-0}"
GOES_HOURS="${GOES_HOURS:-14-23}"
LON_MIN="${LON_MIN:?LON_MIN is required}"
LAT_MIN="${LAT_MIN:?LAT_MIN is required}"
LON_MAX="${LON_MAX:?LON_MAX is required}"
LAT_MAX="${LAT_MAX:?LAT_MAX is required}"

mkdir -p "${MASK_DIR}" "${GIF_DIR}" "${LOG_DIR}"

month_pad="$(printf '%02d' "${MONTH}")"
month_num="$((10#${MONTH}))"
month_log="${LOG_DIR}/rgbmask_${DOMAIN}_${GOES}_${YEAR}${month_pad}.log"

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${month_num})[1])
PY
)"

read -r START_HOUR_UTC END_HOUR_UTC < <("${PYTHON_BIN}" - <<PY
value = "${GOES_HOURS}".strip()
if not value:
    hours = list(range(24))
elif "-" in value and "," not in value:
    start, end = (int(part) for part in value.split("-", 1))
    hours = list(range(start, end + 1))
else:
    hours = [int(part) for part in value.split(",") if part.strip()]
if any(hour < 0 or hour > 23 for hour in hours):
    raise SystemExit(f"GOES_HOURS values must be between 0 and 23: {value!r}")
print(min(hours), min(max(hours) + 1, 24))
PY
)

echo "=== MONTH VINTAGE MASK CONFIG ===" | tee "${month_log}"
echo "YEAR=${YEAR} MONTH=${MONTH}" | tee -a "${month_log}"
echo "RGB_DIR=${RGB_DIR}" | tee -a "${month_log}"
echo "ERA5_DIR=${ERA5_DIR}" | tee -a "${month_log}"
echo "MASK_DIR=${MASK_DIR}" | tee -a "${month_log}"
echo "GIF_DIR=${GIF_DIR}" | tee -a "${month_log}"
echo "GOES_HOURS=${GOES_HOURS}" | tee -a "${month_log}"
echo "MASK WINDOW UTC=${START_HOUR_UTC}-${END_HOUR_UTC}" | tee -a "${month_log}"
echo "KEEP_MASK_DIAGNOSTICS=${KEEP_MASK_DIAGNOSTICS}" | tee -a "${month_log}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}" | tee -a "${month_log}"

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
    --lon-min "${LON_MIN}"
    --lat-min "${LAT_MIN}"
    --lon-max "${LON_MAX}"
    --lat-max "${LAT_MAX}"
    --start-hour-utc "${START_HOUR_UTC}"
    --end-hour-utc "${END_HOUR_UTC}"
    --overwrite
  )
  if [[ "${SKIP_DOWNLOAD}" == "1" ]]; then
    args+=(--skip-download)
  fi
  if [[ "${KEEP_MASK_DIAGNOSTICS}" == "1" ]]; then
    args+=(--keep-diagnostics)
  fi

  echo "[build] ${date_str}" | tee -a "${month_log}"
  if ! "${PYTHON_BIN}" "${args[@]}" >>"${month_log}" 2>&1; then
    echo "[fail] ${date_str}" | tee -a "${month_log}"
    continue
  fi
done

echo "=== MONTH VINTAGE MASK DONE ===" | tee -a "${month_log}"
