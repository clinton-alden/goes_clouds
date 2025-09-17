#!/bin/bash

# IMPORTANT - channels need to be 3 digits (ie C02, C05, or C13)
# Loop through days and channels
for day in $(seq 1 31); do
    for channel in C02 C05 C13; do
        # Run the download script
        python ./download-goes.py --bucket noaa-goes17 \
                                  --year 2022 \
                                  --month 8 \
                                  --days $day $day \
                                  --product ABI-L1b-RadC \
                                  --channel $channel \
                                  --bounds -125 45 -120 49 \
                                  --dir /storage/cdalden/goes/washington
    done
done