#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBS_SCRIPT="${SCRIPT_DIR}/submit_month_colorado_rgbmask.pbs"

for year_month in \
  2022-04 2022-05 2022-06 2022-07 2022-08 2022-09 \
  2022-10 2022-11 2022-12 \
  2023-01 2023-02 2023-03 2023-04 2023-05 2023-06
do
  year="${year_month%-*}"
  month="${year_month#*-}"
  qsub -v YEAR="${year}",MONTH="${month}" "${PBS_SCRIPT}"
done
