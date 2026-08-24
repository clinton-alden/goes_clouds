# East River GOES-16 RGB backfill

This workflow downloads and orthorectifies GOES-16 ABI CONUS radiances and
builds the repository's C02/C05/C13 Day Cloud Phase Distinction-style RGB files
for the East River domain.

## Default configuration

- Period: April 2017 through December 2024 (93 monthly jobs)
- Bounds: 38.70 to 39.25 N, 106.75 to 107.25 W
- Satellite/product: GOES-16 `ABI-L1b-RadC`
- Channels: C02, C05, C13
- Hours: 00, 01, and 14 through 23 UTC (the established 14Z-01Z window)
- Output: `/glade/derecho/scratch/cdalden/east_river_goes`
- Concurrency: one 93-month PBS array, throttled to 12 running subjobs
- DEM: one shared `${BASE_DIR}/static/SRTMGL3_east_river.tif`; concurrent
  monthly jobs serialize its initial download and reuse it thereafter

The submitter queues all months at once. PBS starts no more than 12 monthly
array subjobs concurrently.

## Temperature-bin cloud masks

`build_cloud_mask_month.py` creates daily binary masks from the existing RGB
files and the trained Gothic 10 C temperature-bin thresholds. Each monthly job
downloads and validates one dedicated East River ERA5-Land 2 m temperature
file, then reuses it for every day in that month. Months with incomplete daily
RGB inputs are logged as `SKIPPED` before any ERA5 request; rerunning the array
will pick them up after their RGB inputs are complete. Outputs are restart-safe and
validated for completeness, binary values, requested hours, and geographic
orientation. By default the daily files contain only `cloud_binary`; pass
`--include-diagnostics` to retain interpolated temperature and bin fields too.

Preview the 93-month array (no submission):

```bash
python prepare_cloud_mask_jobs.py
```

Submission is deliberately explicit and supports a 1-12 concurrency limit:

```bash
python prepare_cloud_mask_jobs.py --submit --concurrency 12
```

After every array subjob succeeds, a dependent cleanup job removes the two
one-time shell drivers (`run_month.sh` and `submit_east_river_batches.sh`). If
any month fails, the cleanup job does not run and the drivers remain for retry.

## Preview and submit

Preview commands without submitting:

```bash
./submit_east_river_batches.sh
```

After reviewing the configuration, submit with:

```bash
CONFIRM_SUBMIT=YES ./submit_east_river_batches.sh
```

Defaults can be overridden through environment variables documented by the
scripts. The submitter uses an exported `OPENTOPOGRAPHY_API_KEY` when present;
otherwise it reads the existing key from the legacy processing checkpoint and
exports it into the PBS job environment without printing it.
