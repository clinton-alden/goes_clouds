#!/usr/bin/env bash
# Download one month of requested Gothic-point NLDAS-2 and NLDAS-3 LW subsets.
set -euo pipefail

YEAR="${1:?YEAR required}"
MONTH="${2:?MONTH required, 01-12}"
HOURS="${3:-14-23,0}"
FORCE="${4:-0}"
SOURCE="${5:-both}"

REPO="/glade/u/home/cdalden/goes_work"
PY_GOES="/glade/work/cdalden/conda-envs/goes_analysis/bin/python"
PY_DL="/glade/work/cdalden/conda-envs/goes_downloading/bin/python"
START="${YEAR}-${MONTH}-01"
END="$("${PY_GOES}" - <<PY
import calendar
year = int("${YEAR}")
month = int("${MONTH}")
print(f"{year:04d}-{month:02d}-{calendar.monthrange(year, month)[1]:02d}")
PY
)"

cd "${REPO}"

FORCE_ARGS=()
if [[ "${FORCE}" == "1" || "${FORCE}" == "true" || "${FORCE}" == "TRUE" ]]; then
  FORCE_ARGS+=(--force)
fi

echo "Downloading NLDAS LW ${YEAR}-${MONTH}, hours ${HOURS}, source ${SOURCE}"
echo "Start=${START} End=${END} Force=${FORCE}"

if [[ "${SOURCE}" == "both" || "${SOURCE}" == "nldas2" ]]; then
  "${PY_GOES}" analysis/download_nldas2_lw_hours.py \
    --start "${START}" \
    --end "${END}" \
    --hours "${HOURS}" \
    --out-dir /glade/derecho/scratch/cdalden/tmp/nldas2_lw \
    "${FORCE_ARGS[@]}"
fi

if [[ "${SOURCE}" == "both" || "${SOURCE}" == "nldas3" ]]; then
  "${PY_GOES}" - <<PY | while read -r block_start block_end; do
from analysis.download_nldas2_lw_hours import build_times, parse_hours
import pandas as pd

times = build_times("${START}", "${END}", parse_hours("${HOURS}"))
if len(times) == 0:
    raise SystemExit
block_start = times[0]
prev = times[0]
for dt in times[1:]:
    if dt - prev != pd.Timedelta(hours=1):
        print(block_start.strftime("%Y-%m-%dT%H:%M"), prev.strftime("%Y-%m-%dT%H:%M"))
        block_start = dt
    prev = dt
print(block_start.strftime("%Y-%m-%dT%H:%M"), prev.strftime("%Y-%m-%dT%H:%M"))
PY
    "${PY_DL}" analysis/download_nldas3_variable_subset.py \
      --start "${block_start}" \
      --end "${block_end}" \
      --out-dir /glade/derecho/scratch/cdalden/tmp/nldas3_lw \
      --variables LWdown \
      --point-lat 38.96 \
      --point-lon -106.99 \
      "${FORCE_ARGS[@]}"
  done
fi

echo "Finished NLDAS LW ${YEAR}-${MONTH}"
