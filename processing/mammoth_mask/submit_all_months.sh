#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for year in 2022 2023 2024; do
  first_month=1
  if [[ "${year}" -eq 2022 ]]; then
    first_month=5
  fi
  for month in $(seq "${first_month}" 12); do
    job_id="$(qsub -N "mmask_${year}$(printf '%02d' "${month}")" -v "YEAR=${year},MONTH=${month}" "${SCRIPT_DIR}/submit_month_mammoth_mask.pbs")"
    echo "${year}-$(printf '%02d' "${month}") ${job_id}"
  done
done
