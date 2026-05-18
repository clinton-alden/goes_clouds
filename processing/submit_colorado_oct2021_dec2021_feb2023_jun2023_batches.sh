#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Queue these Colorado monthly jobs in batches of 4:
#   2021-10, 2021-11, 2021-12
#   2023-02, 2023-03, 2023-04, 2023-05, 2023-06
#
# Batch behavior matches the 2022 submitter:
# - up to 4 monthly jobs can run concurrently in a batch
# - the next batch is held until the previous batch exits

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${PBS_SCRIPT:-${SCRIPT_DIR}/submit_month_colorado_14z_00z.pbs}"

BATCH_SIZE="${BATCH_SIZE:-4}"

BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/colorado}"
DOMAIN="${DOMAIN:-colorado}"
GOES="${GOES:-goes16}"
GOES_HOURS="${GOES_HOURS:-14-23}"
GIF_START="${GIF_START:-1400}"
GIF_END="${GIF_END:-2355}"
LON_MIN="${LON_MIN:--109}"
LAT_MIN="${LAT_MIN:-37}"
LON_MAX="${LON_MAX:--104}"
LAT_MAX="${LAT_MAX:-41}"

QSUB_BIN="${QSUB_BIN:-qsub}"

MONTH_SPECS=(
  "2021:10"
  "2021:11"
  "2021:12"
  "2023:02"
  "2023:03"
  "2023:04"
  "2023:05"
  "2023:06"
)

if [[ ! -f "${PBS_SCRIPT}" ]]; then
  echo "ERROR: PBS script not found: ${PBS_SCRIPT}" >&2
  exit 2
fi

join_by_colon() {
  local IFS=":"
  echo "$*"
}

prev_batch_ids=()
total_months="${#MONTH_SPECS[@]}"

echo "Queueing ${total_months} targeted monthly jobs in batches of ${BATCH_SIZE}"
echo "PBS script: ${PBS_SCRIPT}"
echo "Months:"
for spec in "${MONTH_SPECS[@]}"; do
  echo "  ${spec/:/-}"
done

for ((batch_start=0; batch_start<total_months; batch_start+=BATCH_SIZE)); do
  batch_end=$((batch_start + BATCH_SIZE - 1))
  if [[ "${batch_end}" -ge "${total_months}" ]]; then
    batch_end=$((total_months - 1))
  fi

  submitted_ids=()
  depend_args=()

  if [[ "${#prev_batch_ids[@]}" -gt 0 ]]; then
    dep_string="$(join_by_colon "${prev_batch_ids[@]}")"
    depend_args=(-W "depend=afterany:${dep_string}")
    echo
    echo "Batch $((batch_start + 1))-$((batch_end + 1)) depends on: ${dep_string}"
  else
    echo
    echo "Submitting first batch: entries $((batch_start + 1))-$((batch_end + 1))"
  fi

  for ((idx=batch_start; idx<=batch_end; idx++)); do
    spec="${MONTH_SPECS[idx]}"
    year="${spec%%:*}"
    month_pad="${spec##*:}"
    month_int=$((10#${month_pad}))
    job_name="co${year}${month_pad}_g16_14z00z"

    qsub_cmd=(
      "${QSUB_BIN}"
      -N "${job_name}"
      "${depend_args[@]}"
      -v "BASE_DIR=${BASE_DIR},YEAR=${year},MONTH=${month_int},DOMAIN=${DOMAIN},GOES=${GOES},GOES_HOURS=${GOES_HOURS},GIF_START=${GIF_START},GIF_END=${GIF_END},LON_MIN=${LON_MIN},LAT_MIN=${LAT_MIN},LON_MAX=${LON_MAX},LAT_MAX=${LAT_MAX}"
      "${PBS_SCRIPT}"
    )

    echo "Submitting ${year}-${month_pad}:"
    printf '  %q' "${qsub_cmd[@]}"
    printf '\n'

    job_id="$("${qsub_cmd[@]}")"
    submitted_ids+=("${job_id}")
    echo "  job id: ${job_id}"
  done

  prev_batch_ids=("${submitted_ids[@]}")
done

echo
echo "All targeted months have been queued."
echo "Each batch can run up to ${BATCH_SIZE} months concurrently."
echo "The next batch starts only after the previous batch exits."
