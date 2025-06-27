# Suppress the RuntimeWarning for the cfgrib engine
# export PYTHONWARNINGS="ignore:Engine 'cfgrib' loading failed"


python ./download-goes.py --bucket noaa-goes16 --year 2022 --month 7 --days 13 14 --product ABI-L1b-RadC --channel C2 --bounds -124 48 -121 49 --dir /storage/cdalden/goes/washington
