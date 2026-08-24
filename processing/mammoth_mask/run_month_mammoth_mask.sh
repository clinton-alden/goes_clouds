#!/bin/bash
set -euo pipefail

YEAR="${YEAR:?YEAR is required}"
MONTH="${MONTH:?MONTH is required}"
PYTHON_BIN="${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ERA5_DIR="${ERA5_DIR:-/glade/derecho/scratch/cdalden/mammoth/era5_land/t2m_hourly}"

next_year="${YEAR}"
next_month="$((10#${MONTH} + 1))"
if [[ "${next_month}" -eq 13 ]]; then
  next_month=1
  next_year="$((YEAR + 1))"
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/download_era5_land_mammoth.py" --year "${YEAR}" --month "${MONTH}" --output-dir "${ERA5_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/download_era5_land_mammoth.py" --year "${next_year}" --month "${next_month}" --output-dir "${ERA5_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/apply_tsi_rules_mammoth.py" --year "${YEAR}" --month "${MONTH}" --era5-dir "${ERA5_DIR}"
