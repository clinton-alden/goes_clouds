#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Gothic ACM (Aerosol Cloud Mask) version - 15km x 15km domain
# Centered at: 38.9600°N, 106.9900°W
# Bounds: -107.08, 38.89, -106.90, 39.03

: "${TRACE:=0}"
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
  cat <<EOF
Usage: $0 [DOMAIN] [GOES] [BASE_DIR] [YEAR] [MONTH]
  DOMAIN   default: gothic_acm
  GOES     default: goes16
  BASE_DIR default: /glade/derecho/scratch/cdalden/gothic_acm
  YEAR     optional, e.g. 2021
  MONTH    optional, e.g. 03 or 3
Env:
  DOWNLOAD_SCRIPT (default ../../data_download/download-goes.py)
  START_YEAR (default 2021, used when YEAR/MONTH not provided)
  END_YEAR   (default 2023, used when YEAR/MONTH not provided)
  START_DAY_OVERRIDE (optional day-of-month override)
  END_DAY_OVERRIDE   (optional day-of-month override)
EOF
  exit 1
}

# args / defaults
DOMAIN=${1:-gothic_acm}
GOES=${2:-goes16}
DOWNLOAD_SCRIPT=${DOWNLOAD_SCRIPT:-../../data_download/download-goes.py}
BASE_DIR=${3:-/glade/derecho/scratch/cdalden/gothic_acm}
YEAR_ARG=${4:-}
MONTH_ARG=${5:-}
START_YEAR=${START_YEAR:-2021}
END_YEAR=${END_YEAR:-2023}
START_DAY_OVERRIDE=${START_DAY_OVERRIDE:-}
END_DAY_OVERRIDE=${END_DAY_OVERRIDE:-}

# sanity
if [ ! -f "${DOWNLOAD_SCRIPT}" ]; then
  echo "ERROR: download script not found at: ${DOWNLOAD_SCRIPT}"
  echo "Set DOWNLOAD_SCRIPT or put the file at that path. Exiting."
  exit 2
fi

# months / years to process
months_years=()
if [ -n "${YEAR_ARG}" ] || [ -n "${MONTH_ARG}" ]; then
  if [ -z "${YEAR_ARG}" ] || [ -z "${MONTH_ARG}" ]; then
    echo "ERROR: YEAR and MONTH must both be set if either is provided"
    usage
  fi
  months_years=("${YEAR_ARG}-$(printf "%02d" "$((10#${MONTH_ARG}))")")
else
  for year in $(seq "${START_YEAR}" "${END_YEAR}"); do
    for month in $(seq 1 12); do
      months_years+=("${year}-$(printf "%02d" "${month}")")
    done
  done
fi

# ACM is a single product; this downloader still requires -c.
# Use a placeholder ABI channel token to satisfy CLI parsing.
product=ABI-L2-ACMC
channel=C02

# helper: zero-pad number (force base-10 to avoid octal interpretation of 08, 09)
pad2() { printf "%02d" "$((10#$1))"; }

# process a single month
process_month() {
    local year="$1"; local month="$2"
    echo "=== START ${year}-${month} (pid $$) ==="
    local month_padded
    month_padded=$(pad2 "${month}")
    local month_num
    month_num=$((10#${month_padded}))

    # compute number of days in month
    local num_days
    if num_days=$(date -d "${year}-${month_padded}-01 +1 month -1 day" +%d 2>/dev/null); then
        :
    else
        echo "Falling back to calendar for ${year}-${month_padded}"
        num_days=$(cal "${month}" "${year}" | awk 'NF {D=$NF} END {print D}')
    fi

    local start_day=1
    local end_day=${num_days}
    if [ -n "${START_DAY_OVERRIDE}" ]; then
      start_day=$((10#${START_DAY_OVERRIDE}))
    fi
    if [ -n "${END_DAY_OVERRIDE}" ]; then
      end_day=$((10#${END_DAY_OVERRIDE}))
    fi

    if [ "${start_day}" -lt 1 ] || [ "${end_day}" -lt "${start_day}" ] || [ "${end_day}" -gt "${num_days}" ]; then
      echo "ERROR: invalid day override for ${year}-${month_padded}: start_day=${start_day} end_day=${end_day} num_days=${num_days}"
      return 1
    fi

    echo "Processing ${year}-${month_padded} days ${start_day}-${end_day}"
    echo "BASE_DIR: ${BASE_DIR}"
    echo "GOES: ${GOES} DOMAIN: ${DOMAIN}"

    umask 022

    # Gothic bounds for ACM
    # Center ~38.960N, -106.99W over ~15km x 15km domain.
    echo "-> Download: GOES=${GOES} YEAR=${year} MONTH=${month} PRODUCT=${product}"
    python "${DOWNLOAD_SCRIPT}" \
        -B noaa-${GOES} \
        -Y $year \
        -M $month_num \
        -D $start_day $end_day \
        -p ${product} \
      -c ${channel} \
        -b -107.08 38.89 -106.90 39.03 \
        -d "${BASE_DIR}" \
        || {
            echo "WARNING: download failed for ${year}-${month} ${product}"
        }

    local month_dir
    month_dir="${BASE_DIR}/${GOES}/${year}/${month_num}"
    local month_nc_count
    month_nc_count=$(find "${month_dir}" -type f -name '*.nc' 2>/dev/null | wc -l || true)
    if [ "${month_nc_count}" -eq 0 ]; then
      echo "ERROR: no downloaded NetCDF files found in ${month_dir}; failing ${year}-${month_padded}."
      return 1
    fi

    echo "Found ${month_nc_count} NetCDF files for ${year}-${month_padded}"

    # STEP 2 - Ortho
    echo "-> Running ortho for ${year}-${month_padded}"
    if ! python ../batch_ortho_acm.py "${BASE_DIR}/${GOES}/${year}/$((10#$month))/" "${DOMAIN}"; then
      echo "batch_ortho_acm failed for ${year}-${month_padded}; skipping remaining steps for this month"
      return 1
    fi

    # STEP 3 - Zarr
    echo "-> Running zarr for ${year}-${month_padded}"
    if ! python ../zarr_acm_v2.py "${BASE_DIR}/" "${year}" "${month_num}" "${GOES}" "${DOMAIN}"; then
      echo "zarr_acm_v2 failed for ${year}-${month_padded}; skipping remaining steps for this month"
      return 1
    fi

    echo "=== DONE ${year}-${month_padded} ==="
}

# Run single month (for array jobs, no internal parallelism)
for ym in "${months_years[@]}"; do
    year="${ym%%-*}"
    month="${ym##*-}"
    process_month "${year}" "${month}"
    break  # Only process first month when called from PBS array
done
