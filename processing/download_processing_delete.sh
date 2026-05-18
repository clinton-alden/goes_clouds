#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Enable optional tracing for debugging
: "${TRACE:=0}"
if [ "${TRACE}" -ne 0 ]; then
  set -x
fi

usage() {
  cat <<EOF
Usage: $0 [DOMAIN] [GOES] [BASE_DIR] [YEAR] [MONTH]
  DOMAIN   default: colorado
  GOES     default: goes16
  BASE_DIR default: /scratch/cdalden/goes/<domain>/<goes>/
  YEAR     optional, e.g. 2022
  MONTH    optional, e.g. 03 or 3
Env:
  MAX_PROCS (default 4)
  DOWNLOAD_SCRIPT (default ./download-goes.py)
  DEBUG_SINGLE (if set to 1, run only the first month non-parallel for debugging)
  START_YEAR (default 2021, used when YEAR/MONTH not provided)
  END_YEAR   (default 2023, used when YEAR/MONTH not provided)
EOF
  exit 1
}

# args / defaults
DOMAIN=${1:-colorado}
GOES=${2:-goes16}
DOWNLOAD_SCRIPT=${DOWNLOAD_SCRIPT:-../data_download/download-goes.py}
BASE_DIR=${3:-../../scratch/colorado}
MAX_PROCS=${MAX_PROCS:-4}
DEBUG_SINGLE=${DEBUG_SINGLE:-0}
RGB_CLOUD_LOCK=${RGB_CLOUD_LOCK:-/tmp/rgb_cloud_serial.lock}
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

  # Normalize month to 2 digits for path consistency.
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

# process a single month (no bg here)
process_month() {
    local year="$1"; local month="$2"
    echo "=== START ${year}-${month} (pid $$) ==="
    local month_padded
    month_padded=$(pad2 "${month}")

    # compute number of days in month (GNU date)
    local num_days
    if num_days=$(date -d "${year}-${month_padded}-01 +1 month -1 day" +%d 2>/dev/null); then
        :
    else
        echo "Falling back to calendar for ${year}-${month_padded}"
        # fallback using cal (BSD compat)
        num_days=$(cal "${month}" "${year}" | awk 'NF {D=$NF} END {print D}')
    fi

    local start_day=1
    local end_day=${num_days}

    echo "Processing ${year}-${month_padded} days ${start_day}-${end_day}"
    echo "BASE_DIR: ${BASE_DIR}"
    echo "GOES: ${GOES} DOMAIN: ${DOMAIN}"

    # Set umask to ensure downloaded files are writable
    umask 022

    for channel in "${channels[@]}"; do
        echo "-> Download: GOES=${GOES} YEAR=${year} MONTH=${month} CHANNEL=${channel}"
        # call the download script - pass the same args style as your working example
        python "${DOWNLOAD_SCRIPT}" \
            -B noaa-${GOES} \
            -Y $year \
            -M $month \
            -D $start_day $end_day \
            -p ABI-L1b-RadC \
            -c $channel \
            -b -109 37 -104 41 \
            -d "${BASE_DIR}" \
            || {
                echo "WARNING: download failed for ${year}-${month} ${channel}"
            }
    done
            # -b -107.5 38.3 -106 39.4 \
    # delete all hours except hour 18 - be careful about path structure the downloader uses
    echo "-> Keeping only hour 18 for ${year}-${month_padded}"
    for day in $(seq -w "${start_day}" "${end_day}"); do
      for hour in $(seq -w 0 17); do
        dir_path="${BASE_DIR}/${GOES}/${year}/${month}/${day}/ABI-L1b-RadC/${hour}/"
        if [ -d "${dir_path}" ]; then
          chmod -R u+w "${dir_path}" 2>/dev/null || true
          rm -rf "${dir_path}"* 2>/dev/null || true
        fi
      done
      for hour in $(seq -w 19 23); do
        dir_path="${BASE_DIR}/${GOES}/${year}/${month}/${day}/ABI-L1b-RadC/${hour}/"
        if [ -d "${dir_path}" ]; then
          chmod -R u+w "${dir_path}" 2>/dev/null || true
          rm -rf "${dir_path}"* 2>/dev/null || true
        fi
      done
    done

    # STEP 2 - Ortho
    echo "-> Running ortho for ${year}-${month_padded}"
    echo "${BASE_DIR}/${GOES}/${year}/${month}/"
    if ! python ./batch_ortho.py "${BASE_DIR}/${GOES}/${year}/$((10#$month))/" "${DOMAIN}"; then
      echo "batch_ortho failed for ${year}-${month_padded}; skipping remaining steps for this month"
      return 1
    fi

    # STEP 3 - Zarr
    echo "-> Running zarr for ${year}-${month_padded}"
    if ! python ./zarr_v2.py "${BASE_DIR}/" "${year}" "${month}" "${GOES}" "${DOMAIN}"; then
      echo "zarr_v2 failed for ${year}-${month_padded}; skipping remaining steps for this month"
      return 1
    fi
    
    # STEP 4/5 - serialize RGB and cloud frequency to reduce memory pressure
    exec {rgb_lock_fd}>"${RGB_CLOUD_LOCK}"
    flock -x "${rgb_lock_fd}"
    echo "-> Running RGB creation for ${year}-${month_padded} (serialized)"
    python ./rgb_v2.py "${BASE_DIR}/${GOES}" "${year}" "${month}" "${DOMAIN}" "${GOES}" || echo "batch_rgb failed for ${year}-${month_padded}"

    # echo "-> Calculating hourly cloud frequency and deleting RGB files for ${year}-${month_padded} (serialized)"
    # python daily_cloud_frequency.py "${BASE_DIR}" "${year}" "${month}" "${DOMAIN}" "${GOES}" || echo "cloud frequency calculation failed for ${year}-${month_padded}"
    # flock -u "${rgb_lock_fd}"
    # exec {rgb_lock_fd}>&-

    echo "=== DONE ${year}-${month_padded} ==="
}

# run month jobs in parallel, limited by MAX_PROCS
pids=()
for ym in "${months_years[@]}"; do
    year="${ym%%-*}"
    month="${ym##*-}"

    if [ "${DEBUG_SINGLE}" = "1" ]; then
        # Run only first month (non-parallel) for debugging
        process_month "${year}" "${month}"
        break
    fi

    # wait for a free slot
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_PROCS}" ]; do
        sleep 1
    done

    # run in background
    process_month "${year}" "${month}" &
    pids+=($!)
done

if [ "${DEBUG_SINGLE}" -ne 1 ]; then
  # wait for all background jobs
  echo "Waiting for $(jobs -rp | wc -l) background job(s)..."
  wait
  echo "All month jobs finished."
fi