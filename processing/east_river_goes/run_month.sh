#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=/glade/u/home/cdalden/goes_work/processing/east_river_goes
ROOT_DIR=/glade/u/home/cdalden/goes_work
GOTHIC_DIR=${ROOT_DIR}/processing/gothic
PYTHON_BIN=${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}
DOWNLOAD_SCRIPT=${ROOT_DIR}/data_download/download-goes.py
BASE_DIR=${BASE_DIR:-/glade/derecho/scratch/cdalden/east_river_goes}
GOES=${GOES:-goes16}
DOMAIN=${DOMAIN:-east_river}
YEAR=${YEAR:?YEAR is required}
MONTH=${MONTH:?MONTH is required}
GOES_HOURS=${GOES_HOURS:-00,01,14,15,16,17,18,19,20,21,22,23}
LON_MIN=${LON_MIN:--107.25}; LAT_MIN=${LAT_MIN:-38.70}
LON_MAX=${LON_MAX:--106.75}; LAT_MAX=${LAT_MAX:-39.25}
SHARED_DEM_PATH=${SHARED_DEM_PATH:-${BASE_DIR}/static/SRTMGL3_east_river.tif}
export GOES_HOURS LON_MIN LAT_MIN LON_MAX LAT_MAX SHARED_DEM_PATH

month_num=$((10#${MONTH})); month_pad=$(printf '%02d' "${month_num}")
days_in_month=$(date -d "${YEAR}-${month_pad}-01 +1 month -1 day" +%d)
month_dir=${BASE_DIR}/${GOES}/${YEAR}/${month_num}
rgb_dir=${BASE_DIR}/${GOES}/rgb_composite
log_dir=${BASE_DIR}/logs
mkdir -p "${log_dir}" "${rgb_dir}"
exec > >(tee -a "${log_dir}/east_river_${YEAR}${month_pad}.log") 2>&1

rgb_count=$(find "${rgb_dir}" -maxdepth 1 -type f -name "${GOES}_C02_C05_C13_rgb_${DOMAIN}_${YEAR}${month_pad}*.nc" | wc -l)
if [[ "${rgb_count}" -eq "${days_in_month}" ]]; then
  if "${PYTHON_BIN}" "${SCRIPT_DIR}/validate_rgb_month.py" "${rgb_dir}" "${YEAR}" "${month_num}" \
      --goes "${GOES}" --domain "${DOMAIN}" --bounds "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}"; then
    echo "DONE $(date --iso-8601=seconds): ${YEAR}-${month_pad}; already complete and oriented; RGB days=${rgb_count}"
    exit 0
  fi
  echo "Existing RGB month failed geographic validation; rebuilding"
fi

echo "START $(date --iso-8601=seconds): ${YEAR}-${month_pad}; existing RGB=${rgb_count}/${days_in_month}"

# A previous attempt may already have complete daily channel Zarr stores even
# when its final RGB write failed. Rebuild from those first; this is cheap and
# avoids redownloading or re-orthorectifying valid data.
PYTHONPATH="${GOTHIC_DIR}${PYTHONPATH:+:${PYTHONPATH}}" RGB_MAX_WORKERS=${RGB_MAX_WORKERS:-8} \
  "${PYTHON_BIN}" "${GOTHIC_DIR}/rgb_v2.py" "${BASE_DIR}/${GOES}" "${YEAR}" "${month_num}" "${DOMAIN}" "${GOES}"
rgb_count=$(find "${rgb_dir}" -maxdepth 1 -type f -name "${GOES}_C02_C05_C13_rgb_${DOMAIN}_${YEAR}${month_pad}*.nc" | wc -l)
if [[ "${rgb_count}" -eq "${days_in_month}" ]]; then
  if "${PYTHON_BIN}" "${SCRIPT_DIR}/validate_rgb_month.py" "${rgb_dir}" "${YEAR}" "${month_num}" \
      --goes "${GOES}" --domain "${DOMAIN}" --bounds "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}"; then
    echo "DONE $(date --iso-8601=seconds): ${YEAR}-${month_pad}; recovered from existing Zarr; RGB days=${rgb_count}"
    exit 0
  fi
fi

for channel in C02 C05 C13; do
  "${PYTHON_BIN}" "${DOWNLOAD_SCRIPT}" -B "noaa-${GOES}" -Y "${YEAR}" -M "${month_num}" \
    -D 1 "${days_in_month}" -p ABI-L1b-RadC -c "${channel}" \
    -b "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}" -d "${BASE_DIR}"
done

raw_count=$(find "${month_dir}" -type f -name '*.nc' ! -name '*_ortho.nc' | wc -l)
[[ "${raw_count}" -gt 0 ]] || { echo "ERROR: no raw files in ${month_dir}" >&2; exit 1; }
"${PYTHON_BIN}" "${SCRIPT_DIR}/batch_ortho.py" "${month_dir}"
PYTHONPATH="${GOTHIC_DIR}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" \
  "${GOTHIC_DIR}/zarr_v2.py" "${BASE_DIR}/" "${YEAR}" "${month_num}" "${GOES}" "${DOMAIN}"
PYTHONPATH="${GOTHIC_DIR}${PYTHONPATH:+:${PYTHONPATH}}" RGB_MAX_WORKERS=${RGB_MAX_WORKERS:-8} \
  "${PYTHON_BIN}" "${GOTHIC_DIR}/rgb_v2.py" "${BASE_DIR}/${GOES}" "${YEAR}" "${month_num}" "${DOMAIN}" "${GOES}"

missing=()
for day in $(seq -w 1 "${days_in_month}"); do
  file=${rgb_dir}/${GOES}_C02_C05_C13_rgb_${DOMAIN}_${YEAR}${month_pad}${day}.nc
  [[ -s "${file}" ]] || missing+=("${YEAR}${month_pad}${day}")
done
if (( ${#missing[@]} )); then
  echo "ERROR: missing/nonempty RGB days: ${missing[*]}" >&2
  exit 1
fi
"${PYTHON_BIN}" "${SCRIPT_DIR}/validate_rgb_month.py" "${rgb_dir}" "${YEAR}" "${month_num}" \
  --goes "${GOES}" --domain "${DOMAIN}" --bounds "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}"
echo "DONE $(date --iso-8601=seconds): ${YEAR}-${month_pad}; RGB days=${days_in_month}"
