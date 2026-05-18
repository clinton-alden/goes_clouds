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
GOES_HOURS="${GOES_HOURS:-00-23}"
MASK_VAR="${MASK_VAR:-AUTO}"

LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

LOG_DIR="${BASE_DIR}/${GOES}/logs_latlon"
GIF_DIR="${BASE_DIR}/${GOES}/gif_loops_latlon"
TMP_ROOT="${BASE_DIR}/${GOES}/gif_tmp_monthly"
mkdir -p "${LOG_DIR}"
mkdir -p "${GIF_DIR}" "${TMP_ROOT}"
MONTH_LOG="${LOG_DIR}/acm_latlon_${YEAR}$(printf '%02d' "${MONTH}").log"

days_in_month="$("${PYTHON_BIN}" - <<PY
import calendar
print(calendar.monthrange(${YEAR}, ${MONTH})[1])
PY
)"

echo "=== MONTH ACM LATLON CONFIG ===" | tee "${MONTH_LOG}"
echo "YEAR=${YEAR} MONTH=${MONTH}" | tee -a "${MONTH_LOG}"
echo "BASE_DIR=${BASE_DIR}" | tee -a "${MONTH_LOG}"
echo "BOUNDS=${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}" | tee -a "${MONTH_LOG}"
echo "UTC HOURS=${GOES_HOURS}" | tee -a "${MONTH_LOG}"

for day in $(seq 1 "${days_in_month}"); do
  date_str="$(printf "%04d%02d%02d" "${YEAR}" "${MONTH}" "${day}")"
  day_dir="${BASE_DIR}/${GOES}/${YEAR}/${MONTH}/${day}"
  raw_dir="${day_dir}/ABI-L2-ACMC"
  clip_dir="${raw_dir}/clipped"

  echo "${date_str}: download" | tee -a "${MONTH_LOG}"
  mkdir -p "${day_dir}"
  if ! GOES_HOURS="${GOES_HOURS}" \
    "${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" \
      -B "noaa-${GOES}" \
      -Y "${YEAR}" \
      -M "${MONTH}" \
      -D "${day}" "${day}" \
      -p ABI-L2-ACMC \
      -c C02 \
      -b "${LAT_MIN}" "${LAT_MAX}" "${LON_MIN}" "${LON_MAX}" \
      -d "${BASE_DIR}" >>"${MONTH_LOG}" 2>&1; then
    echo "${date_str}: download failed" | tee -a "${MONTH_LOG}"
    continue
  fi

  mapfile -t raw_files < <(find "${raw_dir}" -type f -name '*.nc' ! -path '*/clipped/*' | sort)
  if [[ "${#raw_files[@]}" -eq 0 ]]; then
    echo "${date_str}: no raw files downloaded" | tee -a "${MONTH_LOG}"
    continue
  fi

  mkdir -p "${clip_dir}"
  echo "${date_str}: clip ${#raw_files[@]} files" | tee -a "${MONTH_LOG}"
  for raw_file in "${raw_files[@]}"; do
    out_file="${clip_dir}/$(basename "${raw_file%.nc}")_clip.nc"
    if [[ -f "${out_file}" ]]; then
      continue
    fi
    if ! "${PYTHON_BIN}" "${CLIP_SCRIPT}" \
      "${raw_file}" \
      "${out_file}" \
      --lon-min "${LON_MIN}" \
      --lon-max "${LON_MAX}" \
      --lat-min "${LAT_MIN}" \
      --lat-max "${LAT_MAX}" \
      --mask-var "${MASK_VAR}" >>"${MONTH_LOG}" 2>&1; then
      echo "${date_str}: clip failed for $(basename "${raw_file}")" | tee -a "${MONTH_LOG}"
    fi
  done
done

echo "=== MONTH DAILY PACKAGING START ===" | tee -a "${MONTH_LOG}"
"${PYTHON_BIN}" "${DAILY_SCRIPT}" \
  --base-dir "${BASE_DIR}" \
  --year "${YEAR}" \
  --month "${MONTH}" \
  --goes "${GOES}" \
  --domain "${DOMAIN}" \
  --mask-var "${MASK_VAR}" >>"${MONTH_LOG}" 2>&1
echo "=== MONTH DAILY PACKAGING DONE ===" | tee -a "${MONTH_LOG}"

echo "=== MONTH GIF RENDER START ===" | tee -a "${MONTH_LOG}"
for day in $(seq 1 "${days_in_month}"); do
  date_str="$(printf "%04d%02d%02d" "${YEAR}" "${MONTH}" "${day}")"
  daily_path="${BASE_DIR}/${GOES}/daily_nc_latlon/${GOES}_acm_${DOMAIN}_${date_str}.nc"
  gif_path="${GIF_DIR}/${GOES}_acm_${DOMAIN}_${date_str}.gif"
  tmp_dir="${TMP_ROOT}/${date_str}"

  if [[ ! -f "${daily_path}" ]]; then
    echo "${date_str}: daily file missing, skipping gif" | tee -a "${MONTH_LOG}"
    continue
  fi

  mkdir -p "${tmp_dir}"
  echo "${date_str}: render gif" | tee -a "${MONTH_LOG}"
  if ! "${PYTHON_BIN}" "${PLOT_SCRIPT}" \
    --input "${daily_path}" \
    --output "${gif_path}" \
    --domain "${DOMAIN}" \
    --mask-var "${MASK_VAR}" \
    --tmp-dir "${tmp_dir}" >>"${MONTH_LOG}" 2>&1; then
    echo "${date_str}: gif render failed" | tee -a "${MONTH_LOG}"
  fi
  rmdir "${tmp_dir}" 2>/dev/null || true
done
echo "=== MONTH GIF RENDER DONE ===" | tee -a "${MONTH_LOG}"
