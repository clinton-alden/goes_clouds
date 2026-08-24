#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 YEAR MONTH [goes17|goes18]" >&2
    exit 2
fi

year="$1"
month=$((10#$2))
satellite="${3:-goes18}"
if [ "${satellite}" != "goes17" ] && [ "${satellite}" != "goes18" ]; then
    echo "Satellite must be goes17 or goes18" >&2
    exit 2
fi
base="/glade/derecho/scratch/cdalden/mammoth"
month_root="${base}/${satellite}/${year}/${month}"
month_out="/glade/u/home/cdalden/goes_work/analysis/output_19_sw_temporal_cf_eval/mammoth_acm_hourly_${year}$(printf '%02d' "${month}").csv"
download_script="/glade/u/home/cdalden/goes_work/analysis/download_goes18_acm_day.py"
extract_script="/glade/u/home/cdalden/goes_work/analysis/build_mammoth_acm_hourly_cloud_fraction.py"
python_bin="/glade/work/cdalden/conda-envs/goes_downloading/bin/python"
ortho_script="/glade/u/home/cdalden/goes_work/processing/gothic/batch_ortho.py"

if [ -s "${month_out}" ]; then
    echo "Reusing completed monthly table ${month_out}"
    exit 0
fi
mkdir -p "${month_root}"
num_days=$(date -d "${year}-$(printf '%02d' "${month}")-01 +1 month -1 day" +%d)
for day in $(seq 1 "$((10#${num_days}))"); do
    date_key=$(printf '%04d%02d%02d' "${year}" "${month}" "${day}")
    day_stage="${month_root}/$(printf '%02d' "${day}")/ABI-L2-ACMC"
    ortho_count=$(find "${day_stage}" -type f -name '*_ortho.nc' 2>/dev/null | wc -l || true)
    if [ "${ortho_count}" -ge 12 ]; then
        echo "Reusing ${ortho_count} orthorectified scans for ${date_key}"
        continue
    fi
    mkdir -p "${day_stage}"
    "${python_bin}" "${download_script}" \
        --date "${year}-$(printf '%02d' "${month}")-$(printf '%02d' "${day}")" \
        --output-dir "${day_stage}" --satellite "${satellite}" \
        --workers "${ACM_DOWNLOAD_WORKERS:-8}"
    if ! find "${day_stage}" -type f -name 'OR_ABI-L2-ACMC-*.nc' ! -name '*_ortho.nc' -print -quit | grep -q .; then
        echo "No ${satellite} ACM available for ${date_key}; skipping"
        continue
    fi
done

export GOES_DATA_VARS=BCM
export ORTHO_MAX_WORKERS="${ORTHO_MAX_WORKERS:-16}"
export ORTHO_ACM_MAX_WORKERS="${ORTHO_ACM_MAX_WORKERS:-16}"
export ORTHO_DEM_FILE="${ORTHO_DEM_FILE:-${base}/static/SRTMGL3_mammoth_50km.tif}"
"${python_bin}" "${ortho_script}" "${month_root}" mammoth
"${python_bin}" "${extract_script}" --ortho-dir "${month_root}" --output "${month_out}"

# The monthly table is the retained analysis product; remove scan staging only
# after it was written successfully. A failed job therefore remains resumable.
find "${month_root}" -type f -name '*_ortho.nc' -delete
find "${month_root}" -depth -type d -empty -delete
