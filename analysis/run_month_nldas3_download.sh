#!/usr/bin/env bash
# Download one month of requested Gothic point NLDAS-3 SW subsets.
set -euo pipefail

YEAR="${1:?YEAR required}"
MONTH="${2:?MONTH required, 01-12}"
HOURS="${3:-14-23,0}"
FORCE_NLDAS3="${4:-0}"

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

echo "Downloading NLDAS3 ${YEAR}-${MONTH}, hours ${HOURS}"
echo "Start=${START} End=${END}"
echo "Force NLDAS3=${FORCE_NLDAS3}"

FORCE_ARGS=()
if [[ "${FORCE_NLDAS3}" == "1" || "${FORCE_NLDAS3}" == "true" || "${FORCE_NLDAS3}" == "TRUE" ]]; then
  FORCE_ARGS+=(--force-nldas3)
fi

"${PY}" analysis/build_nldas_point_eval.py \
  --start "${START}" \
  --end "${END}" \
  --hours "${HOURS}" \
  --download-nldas3 \
  --nldas2-dir /glade/u/home/cdalden/scratch/nldas \
  --nldas3-dir /glade/derecho/scratch/cdalden/tmp/nldas3_sw \
  --obs-dir /glade/u/home/cdalden/scratch/surface_obs/colorado/gucsebsM1.b1 \
  --out-dir analysis/output_17_nldas_point_eval \
  "${FORCE_ARGS[@]}"

echo "Finished NLDAS3 ${YEAR}-${MONTH}"
