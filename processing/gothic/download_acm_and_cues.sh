#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

: "${TRACE:=0}"
export HDF5_USE_FILE_LOCKING=FALSE
if [ "${TRACE}" -ne 0 ]; then
  set -x
fi

host_name="${HOSTNAME:-$(hostname)}"
if [[ "${host_name}" =~ ^crlogin[0-9]+$ ]] && [ "${ALLOW_LOGIN_RUN:-0}" != "1" ]; then
  echo "ERROR: refusing to run on interactive login node (${host_name})."
  echo "Submit this workflow with qsub from a PBS script."
  echo "If you must override for debugging, set ALLOW_LOGIN_RUN=1 (not recommended)."
  exit 97
fi

usage() {
  cat <<USAGE
Usage: $0 [DOMAIN] [GOES] [GOES_BASE_DIR]

DOMAIN             default: gothic
GOES               default: goes16
GOES_BASE_DIR      default: /glade/derecho/scratch/cdalden/<domain>

Env overrides:
  RGB_COMPOSITE_DIR   directory containing completed RGB files
  CUES_SOURCE_PATH    input CUES file (default CUES file in analysis/)
  CUES_OUTPUT_DIR     where month cue chunks are written (default <base>/cues/<domain>/)
  DOWNLOAD_SCRIPT     download-goes.py (default ../../data_download/download-goes.py)
  ORTHO_SCRIPT        batch_ortho.py (default ./batch_ortho.py)
  ZARR_SCRIPT         ACM daily zarr writer (default ./zarr_acm_daily.py)
  CUES_EXTRACT_SCRIPT python cue extractor (default ./extract_cues_month.py)
  ORTHO_MAX_WORKERS   passed to batch_ortho.py (default 28)
  GOES_DATA_VARS      comma separated data vars for orthorectification (default ACM)
  ZARR_CHANNEL        channel token for ACM zarr files (default ACMC)
  TEST_DAY            optional one-day test mode (1-31)
  DRY_RUN             if 1, print commands only (default 0)
USAGE
  exit 1
}

DOMAIN=${1:-${DOMAIN:-gothic}}
GOES=${2:-${GOES:-goes16}}
GOES_BASE_DIR=${3:-${GOES_BASE_DIR:-/glade/derecho/scratch/cdalden/${DOMAIN}}}
RGB_COMPOSITE_DIR=${RGB_COMPOSITE_DIR:-${GOES_BASE_DIR}/${GOES}/rgb_composite}
CUES_SOURCE_PATH=${CUES_SOURCE_PATH:-/glade/u/home/cdalden/goes_work/analysis/CUES_1min_data_atmos_radiation_soiltemp_precip_2015to2025.nc}
CUES_OUTPUT_DIR=${CUES_OUTPUT_DIR:-${GOES_BASE_DIR}/cues/${DOMAIN}}
DOWNLOAD_SCRIPT=${DOWNLOAD_SCRIPT:-../../data_download/download-goes.py}
ORTHO_SCRIPT=${ORTHO_SCRIPT:-./batch_ortho.py}
ZARR_SCRIPT=${ZARR_SCRIPT:-./zarr_acm_daily.py}
CUES_EXTRACT_SCRIPT=${CUES_EXTRACT_SCRIPT:-./extract_cues_month.py}
ORTHO_MAX_WORKERS=${ORTHO_MAX_WORKERS:-28}
GOES_DATA_VARS=${GOES_DATA_VARS:-ACM}
ZARR_CHANNEL=${ZARR_CHANNEL:-ACMC}
PBS_ARRAY_INDEX=${PBS_ARRAY_INDEX:-1}
TEST_DAY=${TEST_DAY:-}
DRY_RUN=${DRY_RUN:-0}
PYTHON=${PYTHON:-python}
SKIP_DOWNLOAD_IF_PRESENT=${SKIP_DOWNLOAD_IF_PRESENT:-1}

if [ "${DOMAIN}" = "cues" ] && [ ! -d "${RGB_COMPOSITE_DIR}" ]; then
  alt_rgb_dir=/glade/u/home/cdalden/scratch/mammoth/goes16/rgb_composite
  if [ -d "${alt_rgb_dir}" ]; then
    echo "INFO: DOMAIN=cues using mammoth RGB dir: ${alt_rgb_dir}"
    RGB_COMPOSITE_DIR="${alt_rgb_dir}"
  fi
fi

if [ ! -d "${RGB_COMPOSITE_DIR}" ]; then
  echo "ERROR: RGB_COMPOSITE_DIR does not exist: ${RGB_COMPOSITE_DIR}"
  usage
fi

RGB_TOKEN_DOMAIN="${DOMAIN}"
if [ "${DOMAIN}" = "cues" ]; then
  RGB_TOKEN_DOMAIN="mammoth"
fi

first_file=$(find "${RGB_COMPOSITE_DIR}" -maxdepth 1 -type f -name "${GOES}_*_rgb_${RGB_TOKEN_DOMAIN}_*.nc" | sort | head -n 1 || true)
last_file=$(find "${RGB_COMPOSITE_DIR}" -maxdepth 1 -type f -name "${GOES}_*_rgb_${RGB_TOKEN_DOMAIN}_*.nc" | sort | tail -n 1 || true)
if [ -z "${first_file}" ] || [ -z "${last_file}" ]; then
  echo "ERROR: no RGB files found under ${RGB_COMPOSITE_DIR} for domain ${DOMAIN}"
  exit 2
fi

first_token="$(basename "${first_file}" | sed -E 's/.*_([0-9]{8})\.nc/\1/')"
last_token="$(basename "${last_file}" | sed -E 's/.*_([0-9]{8})\.nc/\1/')"
if [ "${first_token}" == "${first_file}" ] || [ "${last_token}" == "${last_file}" ]; then
  echo "ERROR: could not parse YYYYMMDD token from ${first_file} or ${last_file}"
  exit 3
fi

start_year="${first_token:0:4}"
start_month="${first_token:4:2}"
end_year="${last_token:0:4}"
end_month="${last_token:4:2}"
start_ym=$((10#${start_year} * 100 + 10#${start_month}))
end_ym=$((10#${end_year} * 100 + 10#${end_month}))

target_offset=$((PBS_ARRAY_INDEX - 1))
target_month_num=$((10#${start_month} + target_offset))
target_year=$((10#${start_year} + (target_month_num - 1) / 12))
target_month=$(((target_month_num - 1) % 12 + 1))
target_month_padded=$(printf "%02d" "${target_month}")
target_ym=$((target_year * 100 + target_month))

if [ "${target_ym}" -gt "${end_ym}" ]; then
  echo "INFO: array index ${PBS_ARRAY_INDEX} maps to ${target_year}-${target_month_padded} > ${end_year}-${end_month}; done."
  exit 0
fi

case "${DOMAIN}" in
  gothic)
    LON_MIN=-107.08
    LAT_MIN=38.89
    LON_MAX=-106.90
    LAT_MAX=39.03
    ;;
  mammoth)
    LON_MIN=-119.32
    LAT_MIN=37.41
    LON_MAX=-118.75
    LAT_MAX=37.86
    ;;
  cues)
    LON_MIN=-119.32
    LAT_MIN=37.41
    LON_MAX=-118.75
    LAT_MAX=37.86
    ;;
  colorado)
    LON_MIN=-109
    LAT_MIN=37
    LON_MAX=-104
    LAT_MAX=41
    ;;
  scripps)
    LON_MIN=-118
    LAT_MIN=32.5
    LON_MAX=-117
    LAT_MAX=33.5
    ;;
  *)
    echo "ERROR: unsupported DOMAIN=${DOMAIN} for built-in bounds"
    usage
    ;;
esac

if ! num_days=$(date -d "${target_year}-${target_month_padded}-01 +1 month -1 day" +%d 2>/dev/null); then
  num_days=$(cal "${target_month}" "${target_year}" | awk 'NF {d=$NF} END {print d}')
fi

if [ -n "${TEST_DAY}" ]; then
  if ! [[ "${TEST_DAY}" =~ ^[0-9]{1,2}$ ]] || [ "${TEST_DAY}" -lt 1 ] || [ "${TEST_DAY}" -gt 31 ]; then
    echo "ERROR: TEST_DAY must be integer day 1-31"
    exit 5
  fi
  if [ "${TEST_DAY}" -gt "${num_days}" ]; then
    echo "ERROR: TEST_DAY ${TEST_DAY} exceeds days in month ${target_year}-${target_month_padded}"
    exit 6
  fi
  start_day="${TEST_DAY}"
  end_day="${TEST_DAY}"
else
  start_day=1
  end_day="${num_days}"
fi

if [ "${target_month}" -lt 1 ] || [ "${target_month}" -gt 12 ]; then
  echo "ERROR: invalid target month ${target_month}"
  exit 4
fi

echo "Processing ${DOMAIN} ${GOES} ${target_year}-${target_month_padded} (${start_day}-${end_day})"

echo "Downloading ACM product: noaa-${GOES} ${target_year}-${target_month_padded}"
skip_download=0
if [ "${SKIP_DOWNLOAD_IF_PRESENT}" = "1" ] && [ "${start_day}" = "${end_day}" ]; then
  day_dir="${GOES_BASE_DIR}/${GOES}/${target_year}/${target_month}/${start_day}/ABI-L2-ACMC"
  existing_nc_count=$(find "${day_dir}" -type f -name '*.nc' 2>/dev/null | wc -l || true)
  if [ "${existing_nc_count}" -ge 250 ]; then
    skip_download=1
    echo "Skipping download: found ${existing_nc_count} existing ACM files in ${day_dir}"
  fi
fi

if [ "${skip_download}" -eq 0 ]; then
  if [ "${DRY_RUN}" = "1" ]; then
    echo "[DRY RUN] ${PYTHON} ${DOWNLOAD_SCRIPT} -B noaa-${GOES} -Y ${target_year} -M ${target_month} -D ${start_day} ${end_day} -p ABI-L2-ACMC -c C02 -b ${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX} -d ${GOES_BASE_DIR}"
  else
    "${PYTHON}" "${DOWNLOAD_SCRIPT}" \
      -B "noaa-${GOES}" \
      -Y "${target_year}" \
      -M "${target_month}" \
      -D "${start_day}" "${end_day}" \
      -p ABI-L2-ACMC \
      -c C02 \
      -b "${LON_MIN}" "${LAT_MIN}" "${LON_MAX}" "${LAT_MAX}" \
      -d "${GOES_BASE_DIR}"
  fi
fi

echo "Running orthorectification for ${target_year}-${target_month_padded}"
if [ "${DRY_RUN}" = "1" ]; then
  echo "[DRY RUN] mkdir -p ${GOES_BASE_DIR}/cues/${DOMAIN}"
  echo "[DRY RUN] mkdir -p ${GOES_BASE_DIR}/${GOES}/${target_year}/${target_month}/"
else
  mkdir -p "${GOES_BASE_DIR}/cues/${DOMAIN}"
  mkdir -p "${GOES_BASE_DIR}/${GOES}/${target_year}/${target_month}/"
fi

export GOES_DATA_VARS
export ORTHO_MAX_WORKERS
if [ "${DRY_RUN}" = "1" ]; then
  echo "[DRY RUN] ${PYTHON} ${ORTHO_SCRIPT} ${GOES_BASE_DIR}/${GOES}/${target_year}/${target_month}/ ${DOMAIN}"
else
  "${PYTHON}" "${ORTHO_SCRIPT}" \
    "${GOES_BASE_DIR}/${GOES}/${target_year}/${target_month}/" \
    "${DOMAIN}"
fi

echo "Running ACM daily zarr for ${target_year}-${target_month_padded} (${start_day}-${end_day})"
if [ "${DRY_RUN}" = "1" ]; then
  echo "[DRY RUN] ${PYTHON} ${ZARR_SCRIPT} --base-dir ${GOES_BASE_DIR} --year ${target_year} --month ${target_month} --start-day ${start_day} --end-day ${end_day} --goes ${GOES} --domain ${DOMAIN} --channel ${ZARR_CHANNEL}"
else
  "${PYTHON}" "${ZARR_SCRIPT}" \
    --base-dir "${GOES_BASE_DIR}" \
    --year "${target_year}" \
    --month "${target_month}" \
    --start-day "${start_day}" \
    --end-day "${end_day}" \
    --goes "${GOES}" \
    --domain "${DOMAIN}" \
    --channel "${ZARR_CHANNEL}"
fi

if [ -f "${CUES_SOURCE_PATH}" ]; then
  echo "Subsetting CUES month ${target_year}-${target_month_padded}"
  if [ "${DRY_RUN}" = "1" ]; then
    echo "[DRY RUN] ${PYTHON} ${CUES_EXTRACT_SCRIPT} --source ${CUES_SOURCE_PATH} --year ${target_year} --month ${target_month} --output-dir ${CUES_OUTPUT_DIR} --prefix cues_${DOMAIN}"
  else
    "${PYTHON}" "${CUES_EXTRACT_SCRIPT}" \
      --source "${CUES_SOURCE_PATH}" \
      --year "${target_year}" \
      --month "${target_month}" \
      --output-dir "${CUES_OUTPUT_DIR}" \
      --prefix "cues_${DOMAIN}"
  fi
else
  echo "WARNING: CUES source not found at ${CUES_SOURCE_PATH}; skipping CUES extraction."
fi
