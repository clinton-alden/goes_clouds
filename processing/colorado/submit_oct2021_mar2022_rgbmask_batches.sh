#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${SCRIPT_DIR}/submit_month_colorado_rgbmask.pbs"

for year_month in 2021-10 2021-11 2021-12 2022-01 2022-02 2022-03; do
  year="${year_month%-*}"
  month="${year_month#*-}"
  qsub -v YEAR="${year}",MONTH="${month}" "${PBS_SCRIPT}"
done
