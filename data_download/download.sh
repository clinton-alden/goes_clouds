# Suppress the RuntimeWarning for the cfgrib engine
# export PYTHONWARNINGS="ignore:Engine 'cfgrib' loading failed"

# IMPORTANT - use goes_old environment
# I don't fully remember but I think this env manually imported some goes-ortho code
for channel in C2 C5 C13; do
    python ./download-goes.py --bucket noaa-goes16 --year 2022 --month 7 --days 13 31 --product ABI-L1b-RadC --channel $channel --bounds -124 48 -121 49 --dir /storage/cdalden/goes/washington
done

# python ./goes_nc_to_zarr.py