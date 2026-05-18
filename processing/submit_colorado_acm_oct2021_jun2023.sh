#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${PBS_SCRIPT:-${SCRIPT_DIR}/submit_month_colorado_acm.pbs}"
QSUB_BIN="${QSUB_BIN:-qsub}"

BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado_acm}"
DOMAIN="${DOMAIN:-colorado}"
GOES="${GOES:-goes16}"
GOES_HOURS="${GOES_HOURS:-00-23}"
LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

if [[ ! -f "${PBS_SCRIPT}" ]]; then
  echo "ERROR: PBS script not found: ${PBS_SCRIPT}" >&2
  exit 2
fi

echo "Queueing monthly Colorado ACM jobs from 2021-10 through 2023-06"
echo "PBS script: ${PBS_SCRIPT}"

for year in 2021 2022 2023; do
  for month in $(seq 1 12); do
    if [[ "${year}" -eq 2021 && "${month}" -lt 10 ]]; then
      continue
    fi
    if [[ "${year}" -eq 2023 && "${month}" -gt 6 ]]; then
      continue
    fi

    month_pad="$(printf "%02d" "${month}")"
    job_name="coacm${year}${month_pad}"
    qsub_cmd=(
      "${QSUB_BIN}"
      -N "${job_name}"
      -v "BASE_DIR=${BASE_DIR},YEAR=${year},MONTH=${month},DOMAIN=${DOMAIN},GOES=${GOES},GOES_HOURS=${GOES_HOURS},LON_MIN=${LON_MIN},LAT_MIN=${LAT_MIN},LON_MAX=${LON_MAX},LAT_MAX=${LAT_MAX}"
      "${PBS_SCRIPT}"
    )

    echo "Submitting ${year}-${month_pad}:"
    printf '  %q' "${qsub_cmd[@]}"
    printf '\n'

    "${qsub_cmd[@]}"
  done
done

echo
echo "All monthly Colorado ACM jobs have been queued."
