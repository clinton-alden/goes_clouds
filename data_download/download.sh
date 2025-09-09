# Suppress the RuntimeWarning for the cfgrib engine
# export PYTHONWARNINGS="ignore:Engine 'cfgrib' loading failed"


# IMPORTANT - channels need to be 3 digits (ie C02, C05, or C13)


# for day in {1..31}; do
#     for channel in C02 C05 C13; do
#         python ./download-goes.py --bucket noaa-goes17 --year 2022 --month 8 --days $day $day --product ABI-L1b-RadC --channel $channel --bounds -125 46 -121 49 --dir /storage/cdalden/goes/washington
#         # python ../processing/ortho_batch.py /storage/cdalden/goes/washington/goes17/2022/08/$day/ABI-L1b-RadC/$channel/
#     done
# done

#!/bin/bash

# Loop through days and channels
for day in $(seq 1 28); do
    for channel in C02 C05 C13; do
        # Run the download script
        python ./download-goes.py --bucket noaa-goes16 \
                                  --year 2023 \
                                  --month 2 \
                                  --days $day $day \
                                  --product ABI-L1b-RadC \
                                  --channel $channel \
                                  --bounds -109 37 -104 41 \
                                  --dir /storage/cdalden/goes/colorado
    done
done