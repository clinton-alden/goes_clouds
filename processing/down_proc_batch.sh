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
Usage: $0 [DOMAIN] [GOES] [BASE_DIR]
  DOMAIN   default: colorado
  GOES     default: goes16
  BASE_DIR default: /scratch/cdalden/goes/<domain>/<goes>/
Env:
  MAX_PROCS (default 4)
  DOWNLOAD_SCRIPT (default ./download-goes.py)
  DEBUG_SINGLE (if set to 1, run only the first month non-parallel for debugging)
EOF
  exit 1
}

# args / defaults
DOMAIN=${1:-colorado}
GOES=${2:-goes16}
DOWNLOAD_SCRIPT=${DOWNLOAD_SCRIPT:-../data_download/download-goes.py}
BASE_DIR=${3:-../../scratch/colorado}
YEAR=${4:-2022}
MONTH=${5:-04}
MAX_PROCS=${MAX_PROCS:-4}
DEBUG_SINGLE=${DEBUG_SINGLE:-0}
RGB_CLOUD_LOCK=${RGB_CLOUD_LOCK:-/tmp/rgb_cloud_serial.lock}

# sanity
if [ ! -f "${DOWNLOAD_SCRIPT}" ]; then
  echo "ERROR: download script not found at: ${DOWNLOAD_SCRIPT}"
  echo "Set DOWNLOAD_SCRIPT or put the file at that path. Exiting."
  exit 2
fi

# Process single month
months_years=("${YEAR}-${MONTH}")

channels=(C02 C05 C13)

# helper: zero-pad number
pad2() { printf "%02d" "$1"; }

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
    # local end_day=1  # for testing, limit to first 1 days

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
    # delete hours 00-13 - be careful about path structure the downloader uses
    echo "-> Deleting hours 00-13 for ${year}-${month_padded}"
    for day in $(seq -w "${start_day}" "${end_day}"); do
      for hour in $(seq -w 0 13); do
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
    python ./batch_ortho.py "${BASE_DIR}/${GOES}/${year}/$((10#$month))/" "${DOMAIN}" || echo "batch_ortho failed for ${year}-${month_padded}"

    # STEP 3 - Zarr
    echo "-> Running zarr for ${year}-${month_padded}"
    python ./zarr_v2.py "${BASE_DIR}/" "${year}" "${month}" "${GOES}" "${DOMAIN}" || echo "zarr_v2 failed for ${year}-${month_padded}"
    
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