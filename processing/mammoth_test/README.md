# Mammoth Mountain Test Domain

This directory contains test versions of the GOES processing pipeline configured for a small 50km x 50km domain centered on Mammoth Mountain ski area in California.

## Domain Configuration

- **Center**: 37.6308°N, 119.0326°W (Mammoth Mountain)
- **Bounds**: -119.32, 37.41, -118.75, 37.86 (lon_min, lat_min, lon_max, lat_max)
- **Size**: ~50km x ~50km

## Files

- `batch_ortho.py` - Orthorectification with mammoth domain
- `zarr_v2.py` - Convert to zarr format
- `rgb_v2.py` - Generate RGB composites
- `download_processing.sh` - Main pipeline script
- `submit_mammoth_test.pbs` - PBS submission script for single test month

## Usage

### Test single month (June 2022)
```bash
qsub submit_mammoth_test.pbs
```

### Run custom month
```bash
cd /glade/u/home/cdalden/goes_work/processing/mammoth_test
./download_processing.sh mammoth goes16 ../../../scratch/mammoth 2022 06
```

## Output Location

Data will be saved to: `/glade/derecho/scratch/cdalden/mammoth/goes16/`

## Notes

- Default test configuration: June 2022
- Memory: 16GB (smaller domain than Colorado)
- Walltime: 2 hours (shorter for testing)
- Only processes hour 18 UTC to minimize data volume
