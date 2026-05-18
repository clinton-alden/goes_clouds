#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${PBS_SCRIPT:-${SCRIPT_DIR}/submit_month_colorado_acm_latlon.pbs}"
QSUB_BIN="${QSUB_BIN:-qsub}"
BATCH_SIZE="${BATCH_SIZE:-6}"

BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado_acm}"
DOMAIN="${DOMAIN:-colorado}"
GOES="${GOES:-goes16}"
GOES_HOURS="${GOES_HOURS:-00-23}"
MASK_VAR="${MASK_VAR:-AUTO}"
LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

MONTH_SPECS=(
  "2021:10"
  "2021:11"
  "2021:12"
  "2022:01"
  "2022:02"
  "2022:03"
  "2022:04"
  "2022:05"
  "2022:06"
  "2022:07"
  "2022:08"
  "2022:09"
  "2022:10"
  "2022:11"
  "2022:12"
  "2023:01"
  "2023:02"
  "2023:03"
  "2023:04"
  "2023:05"
  "2023:06"
)

join_by_colon() {
  local IFS=":"
  echo "$*"
}

if [[ ! -f "${PBS_SCRIPT}" ]]; then
  echo "ERROR: PBS script not found: ${PBS_SCRIPT}" >&2
  exit 2
fi

prev_batch_ids=()
total_months="${#MONTH_SPECS[@]}"
echo "Queueing ${total_months} monthly jobs in batches of ${BATCH_SIZE}"

for ((batch_start=0; batch_start<total_months; batch_start+=BATCH_SIZE)); do
  batch_end=$((batch_start + BATCH_SIZE - 1))
  if [[ "${batch_end}" -ge "${total_months}" ]]; then
    batch_end=$((total_months - 1))
  fi

  depend_args=()
  submitted_ids=()
  if [[ "${#prev_batch_ids[@]}" -gt 0 ]]; then
    dep_string="$(join_by_colon "${prev_batch_ids[@]}")"
    depend_args=(-W "depend=afterany:${dep_string}")
  fi

  for ((idx=batch_start; idx<=batch_end; idx++)); do
    spec="${MONTH_SPECS[idx]}"
    year="${spec%%:*}"
    month_pad="${spec##*:}"
    month_int=$((10#${month_pad}))
    job_name="coacmll${year}${month_pad}"

    qsub_cmd=(
      "${QSUB_BIN}"
      -N "${job_name}"
      "${depend_args[@]}"
      -v "BASE_DIR=${BASE_DIR},YEAR=${year},MONTH=${month_int},DOMAIN=${DOMAIN},GOES=${GOES},GOES_HOURS=${GOES_HOURS},MASK_VAR=${MASK_VAR},LON_MIN=${LON_MIN},LAT_MIN=${LAT_MIN},LON_MAX=${LON_MAX},LAT_MAX=${LAT_MAX}"
      "${PBS_SCRIPT}"
    )

    printf 'Submitting %s-%s: ' "${year}" "${month_pad}"
    printf '%q ' "${qsub_cmd[@]}"
    printf '\n'

    job_id="$("${qsub_cmd[@]}")"
    submitted_ids+=("${job_id}")
    echo "  job id: ${job_id}"
  done

  prev_batch_ids=("${submitted_ids[@]}")
done

echo "All monthly lat/lon ACM jobs queued."
