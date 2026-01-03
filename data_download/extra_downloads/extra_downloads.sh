#!/bin/bash

# IMPORTANT - channels need to be 3 digits (ie C02, C05, or C13)
# Loop through days and channels
# for day in $(seq 1 31); do # don't need to add 1 to the end day since sh is inclusive
#     for channel in C02 C05 C13; do
#         # Run the download script
#         python ../download-goes.py --bucket noaa-goes18 \
#                                   --year 2023 \
#                                   --month 7 \
#                                   --days $day $day \
#                                   --product ABI-L1b-RadC \
#                                   --channel $channel \
#                                   --bounds -118 32.5 -117 33.5 \
#                                   --dir /storage/cdalden/goes/scripps
#     done
# done


# note - when downloading L2 products, feed a channel. It will be ignored anyways but the script still needs something here
python ../download-goes.py -B noaa-goes18 -Y 2023 -M 5 -D 1 31 -p ABI-L2-ACMC -c C02 -b -118 32.5 -117 33.5 -d /storage/cdalden/goes/scripps/
python ../download-goes.py -B noaa-goes18 -Y 2023 -M 6 -D 1 30 -p ABI-L2-ACMC -c C02 -b -118 32.5 -117 33.5 -d /storage/cdalden/goes/scripps/
python ../download-goes.py -B noaa-goes18 -Y 2023 -M 7 -D 1 31 -p ABI-L2-ACMC -c C02 -b -118 32.5 -117 33.5 -d /storage/cdalden/goes/scripps/
python ../download-goes.py -B noaa-goes18 -Y 2023 -M 8 -D 1 31 -p ABI-L2-ACMC -c C02 -b -118 32.5 -117 33.5 -d /storage/cdalden/goes/scripps/
