#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Queue 2022 monthly jobs in batches of 4.
# Each batch of 4 months can run concurrently.
# Later batches are submitted with PBS dependencies on the previous batch.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${PBS_SCRIPT:-${SCRIPT_DIR}/submit_month_colorado_14z_00z.pbs}"

YEAR="${YEAR:-2022}"
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

if [[ ! -f "${PBS_SCRIPT}" ]]; then
  echo "ERROR: PBS script not found: ${PBS_SCRIPT}" >&2
  exit 2
fi

join_by_colon() {
  local IFS=":"
  echo "$*"
}

prev_batch_ids=()

echo "Queueing ${YEAR} monthly jobs in batches of ${BATCH_SIZE}"
echo "PBS script: ${PBS_SCRIPT}"

for batch_start in $(seq 1 "${BATCH_SIZE}" 12); do
  batch_end=$((batch_start + BATCH_SIZE - 1))
  if [[ "${batch_end}" -gt 12 ]]; then
    batch_end=12
  fi

  submitted_ids=()

  depend_args=()
  if [[ "${#prev_batch_ids[@]}" -gt 0 ]]; then
    dep_string="$(join_by_colon "${prev_batch_ids[@]}")"
    depend_args=(-W "depend=afterany:${dep_string}")
    echo
    echo "Batch ${batch_start}-${batch_end} depends on: ${dep_string}"
  else
    echo
    echo "Submitting first batch: months ${batch_start}-${batch_end}"
  fi

  for month in $(seq "${batch_start}" "${batch_end}"); do
    month_pad="$(printf "%02d" "${month}")"
    job_name="co${YEAR}${month_pad}_g16_14z00z"

    qsub_cmd=(
      "${QSUB_BIN}"
      -N "${job_name}"
      "${depend_args[@]}"
      -v "BASE_DIR=${BASE_DIR},YEAR=${YEAR},MONTH=${month},DOMAIN=${DOMAIN},GOES=${GOES},GOES_HOURS=${GOES_HOURS},GIF_START=${GIF_START},GIF_END=${GIF_END},LON_MIN=${LON_MIN},LAT_MIN=${LAT_MIN},LON_MAX=${LON_MAX},LAT_MAX=${LAT_MAX}"
      "${PBS_SCRIPT}"
    )

    echo "Submitting month ${month_pad}:"
    printf '  %q' "${qsub_cmd[@]}"
    printf '\n'

    job_id="$("${qsub_cmd[@]}")"
    submitted_ids+=("${job_id}")
    echo "  job id: ${job_id}"
  done

  prev_batch_ids=("${submitted_ids[@]}")
done

echo
echo "All 12 months have been queued."
echo "Each batch of 4 months can run concurrently."
echo "The next batch starts only after the previous batch exits."
