#!/bin/bash

# IMPORTANT - channels need to be 3 digits (ie C02, C05, or C13)
# Loop through days and channels
# for day in $(seq 4 4); do # don't need to add 1 to the end day since sh is inclusive
#     for channel in C02 C05 C13; do
#         # Run the download script
#         python ./download-goes.py --bucket noaa-goes18 \
#                                   --year 2024 \
#                                   --month 12 \
#                                   --days $day $day \
#                                   --product ABI-L1-RadC \ 
#                                   -c $channel \
#                                   --bounds -125 45 -119 49 \
#                                   --dir /storage/cdalden/goes/washington
#     done
# done


# note - when downloading L2 products, feed a channel. It will be ignored anyways but the script still needs something here
# python ./download-goes.py -B noaa-goes16 -Y 2023 -M 12 -D 27 31 -p ABI-L2-ACMC -c C02 -b -109 37 -104 41 -d /storage/cdalden/goes/colorado/

for day in $(seq 4 4); do # don't need to add 1 to the end day since sh is inclusive
    for channel in C02 C05 C13; do
        # Run the download script
        python ./download-goes.py -B noaa-goes18 -Y 2024 -M 12 -D $day $day -p ABI-L1-RadC -c $channel -b -125 45 -119 49 -d /storage/cdalden/goes/washington
    done
done