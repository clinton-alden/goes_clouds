#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${PBS_SCRIPT:-${SCRIPT_DIR}/submit_month_washington_45_49_no_rgb.pbs}"
QSUB_BIN="${QSUB_BIN:-qsub}"

YEARS_STR="${YEARS:-2021 2022 2023}"
MONTHS_STR="${MONTHS:-5 6 7 8 9}"
BATCH_SIZE="${BATCH_SIZE:-5}"

BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/washington}"
DOMAIN="${DOMAIN:-washington}"
GOES="${GOES:-goes17}"
GOES_HOURS="${GOES_HOURS:-14-23}"
LON_MIN="${LON_MIN:--125}"
LAT_MIN="${LAT_MIN:-45}"
LON_MAX="${LON_MAX:--120}"
LAT_MAX="${LAT_MAX:-49}"

if [[ ! -f "${PBS_SCRIPT}" ]]; then
  echo "ERROR: PBS script not found: ${PBS_SCRIPT}" >&2
  exit 2
fi

IFS=' ' read -r -a YEARS <<< "${YEARS_STR}"
IFS=' ' read -r -a MONTHS <<< "${MONTHS_STR}"

join_by_colon() {
  local IFS=":"
  echo "$*"
}

prev_batch_ids=()

echo "Queueing Washington May-Sep monthly jobs in yearly batches of ${BATCH_SIZE}"
echo "PBS script: ${PBS_SCRIPT}"
echo "Years: ${YEARS[*]}"
echo "Months: ${MONTHS[*]}"

for year in "${YEARS[@]}"; do
  submitted_ids=()
  depend_args=()

  if [[ "${#prev_batch_ids[@]}" -gt 0 ]]; then
    dep_string="$(join_by_colon "${prev_batch_ids[@]}")"
    depend_args=(-W "depend=afterany:${dep_string}")
    echo
    echo "Submitting ${year} batch after: ${dep_string}"
  else
    echo
    echo "Submitting first batch for ${year}"
  fi

  for month in "${MONTHS[@]}"; do
    month_pad="$(printf "%02d" "${month}")"
    job_name="wa${year}${month_pad}_g17_14z00z"

    qsub_cmd=(
      "${QSUB_BIN}"
      -N "${job_name}"
      "${depend_args[@]}"
      -v "BASE_DIR=${BASE_DIR},YEAR=${year},MONTH=${month},DOMAIN=${DOMAIN},GOES=${GOES},GOES_HOURS=${GOES_HOURS},LON_MIN=${LON_MIN},LAT_MIN=${LAT_MIN},LON_MAX=${LON_MAX},LAT_MAX=${LAT_MAX}"
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
echo "All Washington May-Sep jobs have been queued."
echo "Only one year (5 month jobs) is eligible to run at a time."
