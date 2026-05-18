# GOES RGB + ERA5 Cloud-Mask Workflow

This branch is a standalone, public workflow for generating GOES Day Cloud Phase RGB files and applying an ERA5-Land-temperature-dependent cloud mask. It intentionally excludes the private analysis notebooks, experiments, and site-specific test code from the research repository.

The workflow does six things:

1. Create the `goes_downloading` conda environment.
2. Download GOES ABI L1b radiance files for `C02`, `C05`, and `C13`.
3. Orthorectify GOES files to a latitude/longitude study domain.
4. Convert daily NetCDF files to Zarr.
5. Build daily RGB NetCDF files.
6. Download ERA5-Land 2 m temperature and apply temperature-bin RGB cloud thresholds.

## Repository Layout

```text
.
├── README.md
├── environment.yml
├── example_cdsapirc
├── example_private_env.sh
├── docs/
│   └── pbs_notes.md
├── notebooks/
│   └── demo_goes_rgb_mask.ipynb
├── scripts/
│   ├── apply_tempbin_thresholds.py
│   ├── batch_ortho.py
│   ├── download-goes.py
│   ├── orthorectify_modded.py
│   ├── run_month_goes_rgb.sh
│   ├── run_month_rgbmask.sh
│   ├── submit_month_goes_rgb.pbs
│   ├── submit_month_rgbmask.pbs
│   ├── utils.py
│   └── zarr_v2_days.py
└── thresholds/
    └── gothic_temp_bin_rgb_thresholds_10c.csv
```

## 1. Create The Environment

```bash
git clone https://github.com/<your-user>/<your-repo>.git
cd <your-repo>
git checkout deployed

conda env create -f environment.yml
conda activate goes_downloading
```

If `goes_ortho` has trouble installing on your machine, install the geospatial/compiler dependencies first, then retry:

```bash
conda install -c conda-forge gdal gcc_linux-64 gxx_linux-64
pip install goes_ortho
```

## Demo Notebook

For a guided one-day example, open:

```text
notebooks/demo_goes_rgb_mask.ipynb
```

The notebook lets a user set lat/lon bounds, dates, and UTC hours, then runs the full workflow and ends with a side-by-side RGB and cloud-mask plot. The default demo is Colorado on `2020-06-30` using one 5-minute daylight scan so it can run as a practical smoke test.

`GOES_HOURS` controls downloaded imagery. A single value like `20` selects the 20 UTC hour; a range like `18-22` is inclusive and selects five hours. `GOES_TIMESTEPS_PER_HOUR=1` keeps only one 5-minute scan per selected hour for demos. Set it to `None` to download every scan in each selected hour. GOES CONUS imagery is usually every 5 minutes, so one full hour is about 12 files per channel. Because the workflow downloads three channels (`C02`, `C05`, `C13`), `18-22` for one day is roughly `5 * 12 * 3 = 180` files before orthorectification.

The notebook uses the small wrapper API in `scripts/workflow.py`, so each major step is a one-line function:

```python
wf.download_goes(config)
wf.orthorectify(config)
wf.build_zarr(config)
wf.build_rgb(config)
wf.apply_mask(config)
```

These functions show compact progress bars and suppress noisy command output unless something fails.

## 2. Add Your Own API Keys

GOES data come from public NOAA AWS buckets and do not require credentials.

Orthorectification downloads DEM data from OpenTopography. Create your own OpenTopography API key, then export it before running any ortho step:

```bash
export OPENTOPOGRAPHY_API_KEY="<your-opentopography-api-key>"
```

ERA5-Land data come from the Copernicus Climate Data Store. Create your own Copernicus/CDS account and write your token to `~/.cdsapirc`:

```bash
cp example_cdsapirc ~/.cdsapirc
chmod 600 ~/.cdsapirc
```

Then edit `~/.cdsapirc` with your real key. Do not commit real API keys to GitHub. `example_private_env.sh` is only a template.

## 3. Configure Your Domain

The examples below use a Colorado domain, but the same workflow works for any region where GOES ABI CONUS imagery and DEM coverage are appropriate.

```bash
export DOMAIN=colorado
export GOES=goes16
export BASE_DIR=/path/to/scratch/colorado
export LON_MIN=-109
export LAT_MIN=37
export LON_MAX=-104
export LAT_MAX=41
export GOES_HOURS=14-23
export OPENTOPOGRAPHY_API_KEY="<your-opentopography-api-key>"
```

For another region, change `DOMAIN`, `GOES`, `BASE_DIR`, and the bounding box. `scripts/batch_ortho.py` reads `LON_MIN`, `LAT_MIN`, `LON_MAX`, and `LAT_MAX`, so custom domains do not require editing the script.

## 4. Run GOES Download, Ortho, Zarr, RGB

Run one month locally or inside an interactive HPC job:

```bash
YEAR=2022 MONTH=3 ./scripts/run_month_goes_rgb.sh
```

Expected RGB output:

```text
${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_YYYYMMDD.nc
```

To submit with PBS, pass the same variables through `qsub`. Adapt account/queue directives in the PBS file for your institution:

```bash
qsub \
  -v YEAR=2022,MONTH=3,BASE_DIR="${BASE_DIR}",DOMAIN="${DOMAIN}",GOES="${GOES}",GOES_HOURS="${GOES_HOURS}",LON_MIN="${LON_MIN}",LAT_MIN="${LAT_MIN}",LON_MAX="${LON_MAX}",LAT_MAX="${LAT_MAX}",OPENTOPOGRAPHY_API_KEY="${OPENTOPOGRAPHY_API_KEY}" \
  scripts/submit_month_goes_rgb.pbs
```

## 5. Download ERA5-Land And Apply The Cloud Mask

For one RGB file, let the cloud-mask script download ERA5-Land automatically:

```bash
python scripts/apply_tempbin_thresholds.py \
  --rgb-file "${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_20220325.nc" \
  --threshold-csv thresholds/gothic_temp_bin_rgb_thresholds_10c.csv \
  --era5-dir "${BASE_DIR}/era5_land/t2m_hourly" \
  --mask-dir "${BASE_DIR}/${GOES}/cloud_mask_tempbin_10c" \
  --gif-dir "${BASE_DIR}/${GOES}/gif_loops_tempbin_10c" \
  --domain "${DOMAIN}" \
  --overwrite
```

Expected ERA5-Land output:

```text
${BASE_DIR}/era5_land/t2m_hourly/era5land_t2m_${DOMAIN}_YYYYMM.nc
```

Expected cloud-mask output:

```text
${BASE_DIR}/${GOES}/cloud_mask_tempbin_10c/${GOES}_C02_C05_C13_rgb_${DOMAIN}_YYYYMMDD_cloud_binary_tempbin10c.nc
```

To process a whole month:

```bash
export RGB_DIR="${BASE_DIR}/${GOES}/rgb_composite"
export ERA5_DIR="${BASE_DIR}/era5_land/t2m_hourly"
export OUTPUT_BASE="${BASE_DIR}/${GOES}"

YEAR=2022 MONTH=3 ./scripts/run_month_rgbmask.sh
```

To submit with PBS:

```bash
qsub \
  -v YEAR=2022,MONTH=3,GOES="${GOES}",DOMAIN="${DOMAIN}",RGB_DIR="${RGB_DIR}",ERA5_DIR="${ERA5_DIR}",OUTPUT_BASE="${OUTPUT_BASE}" \
  scripts/submit_month_rgbmask.pbs
```

## 6. Validate The Setup

```bash
conda activate goes_downloading
python - <<'PY'
import cdsapi
import goes_ortho
import goespy
import xarray
print("GOES workflow imports OK")
PY

test -n "${OPENTOPOGRAPHY_API_KEY}" || echo "Set OPENTOPOGRAPHY_API_KEY before ortho"
```

## Notes For New Research Sites

- `C02`, `C05`, and `C13` are required for this RGB/cloud-mask method.
- `GOES_HOURS=14-23` keeps daytime imagery for western U.S. examples. Change this for other longitudes/seasons; the RGB mask is intended for daylight imagery only.
- `thresholds/gothic_temp_bin_rgb_thresholds_10c.csv` contains trained RGB thresholds from the original research workflow. Treat it as a starting point and validate/retrain for different regions, seasons, land covers, or snow conditions.
- Large monthly jobs can use substantial storage. Keep `BASE_DIR` on scratch or project storage rather than your home directory.
