# PBS Notes

This branch includes two generic PBS templates:

- `scripts/submit_month_goes_rgb.pbs` runs GOES download, orthorectification, Zarr conversion, and RGB generation for one month.
- `scripts/submit_month_rgbmask.pbs` applies the GOES VINTAGE mask for one month.

Before submitting, adapt the PBS account, queue, walltime, CPU, and memory directives to your institution. The templates intentionally avoid hard-coded account names and private filesystem paths.

Orthorectification requires an OpenTopography API key. Export `OPENTOPOGRAPHY_API_KEY` before `qsub`, pass it with `qsub -v`, or use your cluster's preferred secure secret handling. Do not commit a real key in a PBS script.

ERA5-Land downloads use the Copernicus/CDS token in `~/.cdsapirc`. If your compute nodes cannot access your home directory, copy the CDS config to an approved private location and set `HOME` or `CDSAPI_RC` according to your cluster's policy.
