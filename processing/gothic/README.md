# Gothic Test Domain

This directory contains GOES download-to-RGB processing scripts for a small 15km x 15km domain centered on Gothic, Colorado.

## Domain Configuration

- **Center**: 38.9600°N, 106.9900°W
- **Bounds**: -107.08, 38.89, -106.90, 39.03 (lon_min, lat_min, lon_max, lat_max)
- **Size**: ~15km x ~15km

## Files

- `batch_ortho.py` - Orthorectification with Gothic domain
- `zarr_v2.py` - Convert to zarr format
- `rgb_v2.py` - Generate RGB composites
- `download_processing.sh` - Main pipeline script
- `submit_gothic_202110_202306.pbs` - PBS submission script for Oct 2021 to Jun 2023
- `batch_ortho.py` now also supports `GOES_DATA_VARS` (default `Rad`, set to `ACM` for cloud mask products).
- `download_acm_and_cues.sh` - Month-indexed ACM downloader + orthorectifier + CUES monthly extractor
- `extract_cues_month.py` - Subset CUES NetCDF by month
- `submit_gothic_acm_cues.pbs` - PBS job template for ACM+CUES month tasks
- `submit_gothic_acm_cues_array.sh` - Auto-submits array sized from existing Gothic RGB date range

## Usage

### Submit full backfill window (Oct 2021 to Jun 2023)
```bash
qsub submit_gothic_202110_202306.pbs
```

### Submit ACM+CUES work for all months already in Gothic RGB (auto array)
```bash
cd /glade/u/home/cdalden/goes_work/processing/gothic
./submit_gothic_acm_cues_array.sh
```

### Submit single month manually
```bash
cd /glade/u/home/cdalden/goes_work/processing/gothic
export PBS_ARRAY_INDEX=1
./download_acm_and_cues.sh
```

### Run custom month
```bash
cd /glade/u/home/cdalden/goes_work/processing/gothic
./download_processing.sh gothic goes16 /glade/derecho/scratch/cdalden/gothic 2022 06
```

## Output Location

Data will be saved to: `/glade/derecho/scratch/cdalden/gothic/goes16/`

## Notes

- Default GOES platform is GOES-East (`goes16`).
- The PBS array maps 21 indices to months 2021-10 through 2023-06.
- Parallelization:
	- Month-level via PBS array (`#PBS -J 1-21%8`)
	- Per-file orthorectification via `ORTHO_MAX_WORKERS` in `batch_ortho.py`
	- Per-day RGB generation via `RGB_MAX_WORKERS` in `rgb_v2.py`
