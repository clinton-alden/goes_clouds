#!/usr/bin/env bash
set -euo pipefail

BASE_DIR=${BASE_DIR:-/glade/derecho/scratch/cdalden/colorado}
GOES=${GOES:-goes16}
DOMAIN=${DOMAIN:-colorado}

cd /glade/u/home/cdalden/goes_work/processing

start="2021-10"
end="2023-06"

cur="$start"
while :; do
  y="${cur%-*}"
  m="${cur#*-}"
  echo "===== ${cur} ====="
  python ./rgb_backfill_month.py \
    --base-dir "${BASE_DIR}" \
    --goes "${GOES}" \
    --domain "${DOMAIN}" \
    --year "${y}" \
    --month "${m}" \
    --start-date 2021-10-01 \
    --end-date 2023-06-15 \
    --workers 1 \
    --dry-run

  if [ "$cur" = "$end" ]; then
    break
  fi

  next=$(date -d "${cur}-01 +1 month" +%Y-%m)
  cur="$next"
done
