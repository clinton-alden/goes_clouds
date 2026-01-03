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
BASE_DIR=${3:-/scratch/cdalden/goes/${DOMAIN}/${GOES}/}
MAX_PROCS=${MAX_PROCS:-4}
DEBUG_SINGLE=${DEBUG_SINGLE:-0}

# sanity
if [ ! -f "${DOWNLOAD_SCRIPT}" ]; then
  echo "ERROR: download script not found at: ${DOWNLOAD_SCRIPT}"
  echo "Set DOWNLOAD_SCRIPT or put the file at that path. Exiting."
  exit 2
fi

# months / years to process (you can edit this list or pass one via env)
months_years=(
    "2022-1"
    "2022-2"
    "2023-3"
)

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

    echo "Processing ${year}-${month_padded} days ${start_day}-${end_day}"
    echo "BASE_DIR: ${BASE_DIR}"
    echo "GOES: ${GOES} DOMAIN: ${DOMAIN}"

    for channel in "${channels[@]}"; do
        echo "-> Download: GOES=${GOES} YEAR=${year} MONTH=${month} CHANNEL=${channel}"
        # call the download script - pass the same args style as your working example
        python "${DOWNLOAD_SCRIPT}" \
            -B noaa-${GOES} \
            -Y $year \
            -M $month \
            -D $start_day $end_day \
            -p ABI-L1-RadC \
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
        # Note: adjust this path if your downloader writes files to a different structure.
        rm -rf "${BASE_DIR}/${year}/${month_padded}/${day}/ABI-L1-RadC/${hour}/"* 2>/dev/null || true
      done
    done

    # STEP 2 - Ortho
    echo "-> Running ortho for ${year}-${month_padded}"
    python batch_ortho.py "${BASE_DIR}/${year}/${month_padded}" "${DOMAIN}" || echo "batch_ortho failed for ${year}-${month_padded}"

    # STEP 3 - Zarr
    echo "-> Running zarr for ${year}-${month_padded}"
    python zarr_v2.py "${BASE_DIR}" "${year}" "${month_padded}" "${GOES}" || echo "zarr_v2 failed for ${year}-${month_padded}"

    # STEP 4 - RGB file creation
    echo "-> Running RGB creation for ${year}-${month_padded}"
    python batch_rgb.py "${BASE_DIR}" "${year}" "${month_padded}" "${DOMAIN}" "${GOES}" || echo "batch_rgb failed for ${year}-${month_padded}"

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