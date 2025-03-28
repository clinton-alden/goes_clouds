# Suppress the RuntimeWarning for the cfgrib engine
# export PYTHONWARNINGS="ignore:Engine 'cfgrib' loading failed"

python ./download-goes.py --bucket noaa-goes16 --year 2023 --month 6 --days 1 20 --product ABI-L1b-RadC --channel C13 --bounds -109 37 -104 41 --dir /storage/cdalden/goes

# python ./goes_nc_to_zarr.py