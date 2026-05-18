#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_TEMPLATE="${PBS_TEMPLATE:-${SCRIPT_DIR}/submit_mammoth_acm_cues_overlap.pbs}"
CONCURRENCY="${CONCURRENCY:-8}"

START_YEAR="${START_YEAR:-2022}"
START_MONTH="${START_MONTH:-05}"
END_YEAR="${END_YEAR:-2023}"
END_MONTH="${END_MONTH:-01}"

start_total=$((10#${START_YEAR} * 12 + (10#${START_MONTH} - 1)))
end_total=$((10#${END_YEAR} * 12 + (10#${END_MONTH} - 1)))
months=$((end_total - start_total + 1))

if [ "${months}" -le 0 ]; then
  echo "ERROR: invalid start/end date for overlap ACM+CUES run."
  exit 4
fi

qsub_cmd=(qsub -J "1-${months}%${CONCURRENCY}" "${PBS_TEMPLATE}")
echo "Submitting Mammoth ACM+CUES overlap array: ${months} month(s) (${START_YEAR}-${START_MONTH} -> ${END_YEAR}-${END_MONTH})"
printf '  %q' "${qsub_cmd[@]}"
echo
"${qsub_cmd[@]}"
