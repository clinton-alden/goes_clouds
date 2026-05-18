#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Gothic test version - 15km x 15km domain
# Centered at: 38.9600°N, 106.9900°W
# Bounds: -107.08, 38.89, -106.90, 39.03

: "${TRACE:=0}"
if [ "${TRACE}" -ne 0 ]; then
  set -x
fi

usage() {
  cat <<EOF
Usage: $0 [DOMAIN] [GOES] [BASE_DIR] [YEAR] [MONTH]
  DOMAIN   default: gothic
  GOES     default: goes16
  BASE_DIR default: /glade/derecho/scratch/cdalden/gothic
  YEAR     optional, e.g. 2022
  MONTH    optional, e.g. 03 or 3
Env:
  DOWNLOAD_SCRIPT (default ../../data_download/download-goes.py)
  START_YEAR (default 2021, used when YEAR/MONTH not provided)
  END_YEAR   (default 2023, used when YEAR/MONTH not provided)
EOF
  exit 1
}

# args / defaults
DOMAIN=${1:-gothic}
GOES=${2:-goes16}
DOWNLOAD_SCRIPT=${DOWNLOAD_SCRIPT:-../../data_download/download-goes.py}
BASE_DIR=${3:-/glade/derecho/scratch/cdalden/gothic}
YEAR_ARG=${4:-}
MONTH_ARG=${5:-}
START_YEAR=${START_YEAR:-2021}
END_YEAR=${END_YEAR:-2023}

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

channels=(C02 C05 C13)

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

    echo "Processing ${year}-${month_padded} days ${start_day}-${end_day}"
    echo "BASE_DIR: ${BASE_DIR}"
    echo "GOES: ${GOES} DOMAIN: ${DOMAIN}"

    umask 022

    # Gothic bounds using same syntax/order as colorado script.
    # Center ~38.960N, -106.99W over ~15km x 15km domain.
    for channel in "${channels[@]}"; do
        echo "-> Download: GOES=${GOES} YEAR=${year} MONTH=${month} CHANNEL=${channel}"
        python "${DOWNLOAD_SCRIPT}" \
            -B noaa-${GOES} \
            -Y $year \
        -M $month_num \
            -D $start_day $end_day \
            -p ABI-L1b-RadC \
            -c $channel \
            -b -107.08 38.89 -106.90 39.03 \
            -d "${BASE_DIR}" \
            || {
                echo "WARNING: download failed for ${year}-${month} ${channel}"
            }
    done

    local month_dir
    month_dir="${BASE_DIR}/${GOES}/${year}/${month_num}"
    local month_nc_count
    month_nc_count=$(find "${month_dir}" -type f -name '*.nc' 2>/dev/null | wc -l || true)
    if [ "${month_nc_count}" -eq 0 ]; then
      echo "ERROR: no downloaded NetCDF files found in ${month_dir}; failing ${year}-${month_padded}."
      return 1
    fi

    # Keep hour 16 and later (16-23); delete 00-15 before ortho
    echo "-> Deleting hours 00-15 (keeping 16-23) for ${year}-${month_padded}"
    for day in $(seq -w "${start_day}" "${end_day}"); do
      for hour in $(seq -w 0 15); do
        dir_path="${BASE_DIR}/${GOES}/${year}/${month_num}/${day}/ABI-L1b-RadC/${hour}/"
        if [ -d "${dir_path}" ]; then
          chmod -R u+w "${dir_path}" 2>/dev/null || true
          rm -rf "${dir_path}"* 2>/dev/null || true
        fi
      done
    done

    # STEP 2 - Ortho
    echo "-> Running ortho for ${year}-${month_padded}"
    if ! python ./batch_ortho.py "${BASE_DIR}/${GOES}/${year}/$((10#$month))/" "${DOMAIN}"; then
      echo "batch_ortho failed for ${year}-${month_padded}; skipping remaining steps for this month"
      return 1
    fi

    # STEP 3 - Zarr
    echo "-> Running zarr for ${year}-${month_padded}"
    if ! python ./zarr_v2.py "${BASE_DIR}/" "${year}" "${month_num}" "${GOES}" "${DOMAIN}"; then
      echo "zarr_v2 failed for ${year}-${month_padded}; skipping remaining steps for this month"
      return 1
    fi
    
    # STEP 4 - RGB
    echo "-> Running RGB creation for ${year}-${month_padded}"
    python ./rgb_v2.py "${BASE_DIR}/${GOES}" "${year}" "${month_num}" "${DOMAIN}" "${GOES}" || echo "batch_rgb failed for ${year}-${month_padded}"

    echo "=== DONE ${year}-${month_padded} ==="
}

# Run single month (for array jobs, no internal parallelism)
for ym in "${months_years[@]}"; do
    year="${ym%%-*}"
    month="${ym##*-}"
    process_month "${year}" "${month}"
    break  # Only process first month when called from PBS array
done
