# Suppress the RuntimeWarning for the cfgrib engine
# export PYTHONWARNINGS="ignore:Engine 'cfgrib' loading failed"

# IMPORTANT - use goes_old environment
# I don't fully remember but I think this env manually imported some goes-ortho code

# IMPORTANT - channels need to be 3 digits (ie C02, C05, or C13)

for day in {13..31}; do
    for channel in C02 C05 C13; do
        python ./download-goes.py --bucket noaa-goes17 --year 2022 --month 7 --days $day $day --product ABI-L1b-RadC --channel $channel --bounds -124 47 -121 49 --dir /storage/cdalden/goes/washington
    done
done

# python ./goes_nc_to_zarr.py