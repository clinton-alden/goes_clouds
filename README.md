# GOES Visible-Infrared-Near infrared Threshold AlGorithm for cloud Evolution (VINTAGE) Workflow

This branch is a standalone workflow for generating GOES Day Cloud Phase RGB files and applying the ERA5-Land-temperature-dependent GOES VINTAGE mask. It intentionally excludes the private analysis notebooks, experiments, and site-specific test code from the research repository.

The workflow does six things:

1. Create the `goes_vintage_workflow` conda environment.
2. Download GOES ABI L1b radiance files for `C02`, `C05`, and `C13`.
3. Orthorectify GOES files to a latitude/longitude study domain.
4. Convert daily NetCDF files to Zarr.
5. Build daily RGB NetCDF files.
6. Download ERA5-Land 2 m temperature and apply the GOES VINTAGE mask.

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
    └── gothic_vintage_rgb_tree_rules_5c_sw_kt050_090.csv
```

## 1. Create The Environment

```bash
git clone https://github.com/clinton-alden/goes_clouds.git
cd goes_clouds
git checkout deployed

conda env create -f environment.yml
conda activate goes_vintage_workflow
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

The notebook lets a user set lat/lon bounds, dates, and UTC hours, then runs the full workflow and ends with a side-by-side RGB and VINTAGE-mask plot. The default demo is Colorado on `2020-06-30` using one 5-minute daylight scan so it can run as a practical smoke test.

`GOES_HOURS` controls both downloaded imagery and the later VINTAGE-mask time window. A single value like `20` selects the 20 UTC hour and masks `20:00-20:59`; a range like `18-22` is inclusive for download and masks `18:00-22:59`. `GOES_TIMESTEPS_PER_HOUR=1` keeps only one 5-minute scan per selected hour for demos. Set it to `None` to download every scan in each selected hour. GOES CONUS imagery is every 5 minutes, so one full hour is about 12 files per channel. Because the workflow downloads three channels (`C02`, `C05`, `C13`), `18-22` for one day is roughly `5 * 12 * 3 = 180` files before orthorectification.

The notebook uses the small wrapper API in `scripts/workflow.py`, so each major step is a one-line function:

```python
wf.download_goes(config)
wf.orthorectify(config)
wf.build_zarr(config)
wf.build_rgb(config)
wf.apply_mask(config)
```

These functions show compact progress bars and suppress noisy command output unless something fails.

The same full workflow can be run from a configurable shell script:

```bash
source private_env.sh
./scripts/example_full_workflow.sh
```

Override variables to batch a different label, date range, or bounding box:

```bash
DOMAIN=mammoth \
BASE_DIR=/path/to/scratch/mammoth \
START_DATE=2022-03-25 END_DATE=2022-03-31 \
LON_MIN=-119.5 LAT_MIN=36.5 LON_MAX=-118.5 LAT_MAX=37.9 \
GOES_HOURS=18-22 GOES_TIMESTEPS_PER_HOUR=None \
./scripts/example_full_workflow.sh
```

Set `CREATE_PNG=0` to skip the final RGB/VINTAGE-mask PNG.
By default, VINTAGE NetCDF outputs store the compact binary mask plus small
summary variables. Set `KEEP_MASK_DIAGNOSTICS=1` if you also want to store the
per-pixel ERA5 temperature and threshold-bin diagnostic arrays.

## 2. Add Your Own API Keys

GOES data come from public NOAA AWS buckets and do not require credentials.

Copy the private environment template, add your own keys, and source it before
running the workflow:

```bash
cp example_private_env.sh private_env.sh
chmod 600 private_env.sh
# edit private_env.sh with your real keys
source private_env.sh
```

`private_env.sh` is git-ignored. It should contain:

```bash
export OPENTOPOGRAPHY_API_KEY="<your-opentopography-api-key>"
export CDSAPI_URL="https://cds.climate.copernicus.eu/api"
export CDSAPI_KEY="<your-cds-api-key>"
```

## 3. Configure Your Output Label And Bounds

The examples below use a Colorado output label and bounding box, but the same workflow works for any region where GOES ABI CONUS imagery and DEM coverage are appropriate. `DOMAIN` is only used as an output filename suffix; `LON_MIN`, `LAT_MIN`, `LON_MAX`, and `LAT_MAX` define the actual processing area.

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

For another region, change `DOMAIN`, `GOES`, `BASE_DIR`, and the bounding box. `scripts/batch_ortho.py` requires `LON_MIN`, `LAT_MIN`, `LON_MAX`, and `LAT_MAX`; it does not infer bounds from the domain label.

## 4. Run GOES Download, Ortho, Zarr, RGB

Run one month locally or inside an interactive HPC job:

```bash
YEAR=2022 MONTH=3 ./scripts/run_month_goes_rgb.sh
```

Expected RGB output:

```text
${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_YYYYMMDD.nc
```

To submit with PBS, pass the same variables through `qsub`.  account/queue directives in the PBS file for your institution:

```bash
qsub \
  -v YEAR=2022,MONTH=3,BASE_DIR="${BASE_DIR}",DOMAIN="${DOMAIN}",GOES="${GOES}",GOES_HOURS="${GOES_HOURS}",LON_MIN="${LON_MIN}",LAT_MIN="${LAT_MIN}",LON_MAX="${LON_MAX}",LAT_MAX="${LAT_MAX}",OPENTOPOGRAPHY_API_KEY="${OPENTOPOGRAPHY_API_KEY}" \
  scripts/submit_month_goes_rgb.pbs
```

## 5. Download ERA5-Land And Apply The VINTAGE Mask

For one RGB file, let the VINTAGE-mask script download ERA5-Land automatically:

```bash
python scripts/apply_tempbin_thresholds.py \
  --rgb-file "${BASE_DIR}/${GOES}/rgb_composite/${GOES}_C02_C05_C13_rgb_${DOMAIN}_20220325.nc" \
  --threshold-csv thresholds/gothic_vintage_rgb_tree_rules_5c_sw_kt050_090.csv \
  --era5-dir "${BASE_DIR}/era5_land/t2m_hourly" \
  --mask-dir "${BASE_DIR}/${GOES}/vintage_mask" \
  --gif-dir "${BASE_DIR}/${GOES}/vintage_gif_loops" \
  --domain "${DOMAIN}" \
  --lon-min "${LON_MIN}" \
  --lat-min "${LAT_MIN}" \
  --lon-max "${LON_MAX}" \
  --lat-max "${LAT_MAX}" \
  --overwrite
```

Expected ERA5-Land output:

```text
${BASE_DIR}/era5_land/t2m_hourly/era5land_t2m_${DOMAIN}_YYYYMM.nc
```

Expected VINTAGE-mask output:

```text
${BASE_DIR}/${GOES}/vintage_mask/${GOES}_C02_C05_C13_rgb_${DOMAIN}_YYYYMMDD_vintage_mask.nc
```

The default NetCDF stores the compact binary mask. Add `--keep-diagnostics` to store
per-pixel ERA5 temperature and threshold-bin arrays for debugging or method
development.

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
  -v YEAR=2022,MONTH=3,GOES="${GOES}",DOMAIN="${DOMAIN}",GOES_HOURS="${GOES_HOURS}",RGB_DIR="${RGB_DIR}",ERA5_DIR="${ERA5_DIR}",OUTPUT_BASE="${OUTPUT_BASE}",LON_MIN="${LON_MIN}",LAT_MIN="${LAT_MIN}",LON_MAX="${LON_MAX}",LAT_MAX="${LAT_MAX}" \
  scripts/submit_month_rgbmask.pbs
```

## 6. Validate The Setup

```bash
conda activate goes_vintage_workflow
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

- `C02`, `C05`, and `C13` are required for the GOES VINTAGE method.
- `GOES_HOURS=14-23` keeps daytime imagery for western U.S. examples. Change this for other longitudes/seasons; the RGB mask is intended for daylight imagery only.
- `thresholds/gothic_vintage_rgb_tree_rules_5c_sw_kt050_090.csv` contains the current default Gothic VINTAGE RGB decision-tree rules trained from SW-derived cloud labels using `k_t < 0.50` for cloudy and `k_t > 0.90` for clear in 5 C temperature bins. Temperatures colder than `-15 C` use the `[-15, -10)` rules, and temperatures warmer than `20 C` use the `[15, 20)` rules. Treat these rules as a starting point and validate/retrain for different regions, seasons, land covers, or snow conditions.
- Large monthly jobs can use substantial storage. Keep `BASE_DIR` on scratch or project storage rather than your home directory.
