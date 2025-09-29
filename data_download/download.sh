#!/bin/bash

# IMPORTANT - channels need to be 3 digits (ie C02, C05, or C13)
# Loop through days and channels
for day in $(seq 1 30); do # don't need to add 1 to the end day since sh is inclusive
    for channel in C02 C05 C13; do
        # Run the download script
        python ./download-goes.py --bucket noaa-goes16 \
                                  --year 2022 \
                                  --month 6 \
                                  --days $day $day \
                                  --product ABI-L1b-RadC \
                                  --channel $channel \
                                  --bounds -109 37 -104 41 \
                                  --dir /storage/cdalden/goes/colorado
    done
done