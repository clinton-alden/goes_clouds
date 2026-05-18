#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${PBS_SCRIPT:-${SCRIPT_DIR}/submit_month_sierra_goes17.pbs}"

BASE_DIR="${BASE_DIR:-/glade/u/home/cdalden/scratch/sierra}"
DOMAIN="${DOMAIN:-sierra}"
GOES="${GOES:-goes17}"
GOES_HOURS="${GOES_HOURS:-14-23}"
GIF_START="${GIF_START:-1400}"
GIF_END="${GIF_END:-2355}"
START_YEAR="${START_YEAR:-2020}"
START_MONTH="${START_MONTH:-1}"
END_YEAR="${END_YEAR:-2022}"
END_MONTH="${END_MONTH:-12}"
MAX_CONCURRENT="${MAX_CONCURRENT:-6}"
LON_MIN="${LON_MIN:--121}"
LAT_MIN="${LAT_MIN:-35}"
LON_MAX="${LON_MAX:--118}"
LAT_MAX="${LAT_MAX:-40}"

QSUB_BIN="${QSUB_BIN:-qsub}"

if [[ ! -f "${PBS_SCRIPT}" ]]; then
  echo "ERROR: PBS script not found: ${PBS_SCRIPT}" >&2
  exit 2
fi

start_total_months=$((10#${START_YEAR} * 12 + 10#${START_MONTH}))
end_total_months=$((10#${END_YEAR} * 12 + 10#${END_MONTH}))
if (( end_total_months < start_total_months )); then
  echo "ERROR: end window is earlier than start window" >&2
  exit 2
fi

total_months=$((end_total_months - start_total_months + 1))

month_specs=()
for ((offset=0; offset<total_months; offset++)); do
  absolute_month=$((start_total_months + offset - 1))
  year=$((absolute_month / 12))
  month=$((absolute_month % 12 + 1))
  month_specs+=("${year}:$(printf '%02d' "${month}")")
done

echo "Queueing ${total_months} Sierra monthly jobs as individual rolling jobs"
echo "PBS script: ${PBS_SCRIPT}"
echo "Range: ${START_YEAR}-$(printf '%02d' "${START_MONTH}") through ${END_YEAR}-$(printf '%02d' "${END_MONTH}")"
echo "Rolling concurrency limit: ${MAX_CONCURRENT}"
echo "Bounds: ${LON_MIN} ${LAT_MIN} ${LON_MAX} ${LAT_MAX}"
echo "GIF window: ${GIF_START}-${GIF_END}"

submitted_ids=()
for ((idx=0; idx<total_months; idx++)); do
  spec="${month_specs[idx]}"
  year="${spec%%:*}"
  month_pad="${spec##*:}"
  month_int=$((10#${month_pad}))
  job_name="si${year}${month_pad}_g17"

  depend_args=()
  if (( idx >= MAX_CONCURRENT )); then
    upstream_job_id="${submitted_ids[idx - MAX_CONCURRENT]}"
    depend_args=(-W "depend=afterany:${upstream_job_id}")
  fi

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

echo
echo "All Sierra monthly jobs have been queued as individual jobs."
echo "PBS will keep up to ${MAX_CONCURRENT} months active at a time via per-job dependencies."
