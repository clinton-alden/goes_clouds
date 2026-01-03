#!/bin/bash

# python -u ../data_download/download-goes.py \
#   -B noaa-goes16 \
#   -Y 2022 \
#   -M 07 \
#   -D 1 30 \
#   -p ABI-L2-ACMC \
#   -c C02 \
#   -b -109 37 -104 41 \
#   -d /storage/cdalden/goes/colorado/ 2>&1 | tee download-test.log


# for channel in C02 C05 C13; do
#     python -u ../data_download/download-goes.py \
#         -B noaa-goes16 \
#         -Y 2022 \
#         -M 7 \
#         -D 1 30 \
#         -p ABI-L1-RadC \
#         -c $channel \
#         -b -109 37 -104 41 \
#         -d /storage/cdalden/goes/colorado/
# done



for day in $(seq 4 4); do # don't need to add 1 to the end day since sh is inclusive
    for channel in C02 C05 C13; do
        # Run the download script
        python ../data_download/download-goes.py -B noaa-goes18 -Y 2024 -M 12 -D $day $day -p ABI-L1-RadC -c $channel -b -125 45 -119 49 -d /storage/cdalden/goes/washington
    done
done