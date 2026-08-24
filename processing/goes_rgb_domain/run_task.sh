#!/usr/bin/env bash
set -euo pipefail

ROOT=/glade/u/home/cdalden/goes_work
PY=${PYTHON_BIN:-/glade/work/cdalden/conda-envs/goes_downloading/bin/python}
BASE_DIR=${BASE_DIR:?}
GOES=${GOES:?}; DOMAIN=${DOMAIN:?}; YEAR=${YEAR:?}; MONTH=${MONTH:?}
DAY_START=${DAY_START:?}; DAY_END=${DAY_END:?}
LON_MIN=${LON_MIN:?}; LAT_MIN=${LAT_MIN:?}; LON_MAX=${LON_MAX:?}; LAT_MAX=${LAT_MAX:?}
GOES_HOURS=${GOES_HOURS:-00,01,14,15,16,17,18,19,20,21,22,23}
SHARED_DEM_PATH=${SHARED_DEM_PATH:-${BASE_DIR}/static/SRTMGL3_${DOMAIN}.tif}
export GOES_HOURS LON_MIN LAT_MIN LON_MAX LAT_MAX SHARED_DEM_PATH

month=$((10#$MONTH)); first=$(printf '%04d-%02d-%02d' "$YEAR" "$month" "$DAY_START")
last=$(printf '%04d-%02d-%02d' "$YEAR" "$month" "$DAY_END")
rgb_dir=${BASE_DIR}/${GOES}/rgb_composite
month_dir=${BASE_DIR}/${GOES}/${YEAR}/${month}
mkdir -p "$rgb_dir"

validate() {
  "$PY" "$ROOT/processing/goes_rgb_domain/validate_rgb_dates.py" "$rgb_dir" \
    --start "$first" --end "$last" --goes "$GOES" --domain "$DOMAIN" \
    --bounds "$LON_MIN" "$LAT_MIN" "$LON_MAX" "$LAT_MAX"
}
if validate; then echo "DONE: existing outputs are complete"; exit 0; fi

for channel in C02 C05 C13; do
  "$PY" "$ROOT/data_download/download-goes.py" -B "noaa-${GOES}" -Y "$YEAR" -M "$month" \
    -D "$DAY_START" "$DAY_END" -p ABI-L1b-RadC -c "$channel" \
    -b "$LON_MIN" "$LAT_MIN" "$LON_MAX" "$LAT_MAX" -d "$BASE_DIR"
done

"$PY" "$ROOT/processing/east_river_goes/batch_ortho.py" "$month_dir"
PYTHONPATH="$ROOT/processing/gothic" "$PY" "$ROOT/processing/gothic/zarr_v2.py" \
  "$BASE_DIR/" "$YEAR" "$month" "$GOES" "$DOMAIN"
PYTHONPATH="$ROOT/processing/gothic" RGB_MAX_WORKERS=${RGB_MAX_WORKERS:-8} \
  "$PY" "$ROOT/processing/gothic/rgb_v2.py" "$BASE_DIR/$GOES" "$YEAR" "$month" "$DOMAIN" "$GOES"
validate
