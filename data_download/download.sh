# Suppress the RuntimeWarning for the cfgrib engine
# export PYTHONWARNINGS="ignore:Engine 'cfgrib' loading failed"

# IMPORTANT - use goes_old environment
# I don't fully remember but I think this env manually imported some goes-ortho code
python ./download-goes.py --bucket noaa-goes16 --year 2022 --month 7 --days 13 13 --product ABI-L1b-RadC --channel C13 --bounds -124 48 -121 49 --dir /storage/cdalden/goes/washington

# python ./goes_nc_to_zarr.py