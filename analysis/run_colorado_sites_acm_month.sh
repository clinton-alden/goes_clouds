#!/usr/bin/env bash
set -euo pipefail
year="$1"; month=$((10#$2)); ndays=$(date -d "${year}-$(printf '%02d' "$month")-01 +1 month -1 day" +%d)
archive=/glade/u/home/cdalden/scratch/colorado_acm/goes16
stage=/glade/derecho/scratch/cdalden/tmp/colorado_sites_acm_stage
out=/glade/u/home/cdalden/goes_work/analysis/output_20_colorado_sites_hourly/acm_daily
py=/glade/work/cdalden/conda-envs/goes_downloading/bin/python
mkdir -p "$stage" "$out"
for day in $(seq 1 "$((10#$ndays))"); do
  key=$(printf '%04d%02d%02d' "$year" "$month" "$day")
  [ -s "$out/gothic_acm_hourly_${key}.csv" ] && [ -s "$out/table_mountain_acm_hourly_${key}.csv" ] && [ -s "$out/senator_beck_acm_hourly_${key}.csv" ] && continue
  existing="$archive/$year/$month/$day"
  input="$existing"; made_stage=0
  if ! find "$existing" -type f -path '*/ABI-L2-ACMC/*/*.nc' -print -quit 2>/dev/null | grep -q .; then
    input="$stage/$key"; mkdir -p "$input"; made_stage=1
    "$py" analysis/download_goes18_acm_day.py --satellite goes16 --date "${year}-$(printf '%02d' "$month")-$(printf '%02d' "$day")" --output-dir "$input"
  fi
  "$py" analysis/build_colorado_sites_acm_day.py --input-dir "$input" --output-dir "$out" --date-key "$key"
  if [ "$made_stage" -eq 1 ]; then find "$input" -type f -name 'OR_ABI-L2-ACMC-*.nc' -delete; find "$input" -depth -type d -empty -delete; fi
done
