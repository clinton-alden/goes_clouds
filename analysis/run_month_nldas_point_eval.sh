#!/usr/bin/env bash
# Run one monthly Gothic point SW comparison for NLDAS-2, NLDAS-3, and obs.
set -euo pipefail

YEAR="${1:?YEAR required}"
MONTH="${2:?MONTH required, 01-12}"
HOURS="${3:-14-23,0}"

REPO="/glade/u/home/cdalden/goes_work"
PY="/glade/work/cdalden/conda-envs/goes_analysis/bin/python"
START="${YEAR}-${MONTH}-01"
END="$("${PY}" - <<PY
import calendar
year = int("${YEAR}")
month = int("${MONTH}")
print(f"{year:04d}-{month:02d}-{calendar.monthrange(year, month)[1]:02d}")
PY
)"

cd "${REPO}"

echo "Running ${YEAR}-${MONTH}, hours ${HOURS}"
echo "Start=${START} End=${END}"

"${PY}" analysis/download_nldas2_hours.py \
  --start "${START}" \
  --end "${END}" \
  --hours "${HOURS}"

"${PY}" analysis/build_nldas_point_eval.py \
  --start "${START}" \
  --end "${END}" \
  --hours "${HOURS}" \
  --download-nldas3

echo "Finished ${YEAR}-${MONTH}"
