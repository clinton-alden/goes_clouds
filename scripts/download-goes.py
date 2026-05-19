import argparse
import os
from glob import glob

import boto3
import numpy as np
from botocore import UNSIGNED
from botocore.config import Config
from goespy.Downloader import (
    ABI_Downloader,  # https://github.com/palexandremello/goes-py
)

from goes_ortho.clip import subsetNetCDF

##############################################################

# ---------------------------- COMMAND LINE ARGUMENTS ----------------------------#
# Set up argument and error handling
parser = argparse.ArgumentParser(
    description="Produces a timeseries of GOES ABI Radiance observations for a single location given a directory of GOES ABI files"
)
parser.add_argument(
    "-B",
    "--bucket",
    required=True,
    type=str,
    help="AWS S3 Bucket for GOES (e.g. noaa-goes16",
)
parser.add_argument(
    "-Y",
    "--year",
    required=True,
    type=int,
    help="Specify time range to search for GOES ABI imagery (year)",
)
parser.add_argument(
    "-M",
    "--month",
    required=True,
    type=int,
    help="Specify time range to search for GOES ABI imagery (month)",
)
parser.add_argument(
    "-D",
    "--days",
    required=True,
    type=int,
    nargs=2,
    help="Specify time range to search for GOES ABI imagery (start day, stop day)",
)
parser.add_argument(
    "-p",
    "--product",
    required=True,
    type=str,
    help="GOES ABI Product (e.g. ABI-L1b-RadC)",
)
parser.add_argument(
    "-c", "--channel", required=True, type=str, help="GOES ABI channel/band (e.g. C14)"
)
parser.add_argument(
    "-b",
    "--bounds",
    required=True,
    type=float,
    nargs=4,
    help="Bounds to crop GOES ABI image to (min_lat max_lat min_lon max_lon",
)
parser.add_argument(
    "-d", "--dir", required=True, help="Directory to save GOES ABI files (.nc)"
)
args = parser.parse_args()

# -----------------------SET ARGUMENTS TO VARIABLES----------------------------#
indir = args.dir


bucket = args.bucket  # AWS S3 Bucket for GOES
satellite = bucket[
    5:
]  # get the last part of the bucket name which contains satellite name (goes16 or goes17)
# Specify time range to search for GOES ABI imagery (year, month, days)
year = args.year
month = args.month
start_day = args.days[0]
stop_day = args.days[1]
days = np.linspace(
    start_day, stop_day, stop_day - start_day + 1, dtype=np.int16
)  # create list of days from start to stop date
hours = [
    "00",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23",
]  # all hours (can we just use linspace of ints here instead of list of strings?)

# Optional hour override via env var, e.g. GOES_HOURS=16-23 or GOES_HOURS=16,17,18
hours_env = os.environ.get("GOES_HOURS", "").strip()
if hours_env:
    parsed_hours = []
    if "-" in hours_env and "," not in hours_env:
        start_h, end_h = hours_env.split("-", 1)
        start_i = int(start_h)
        end_i = int(end_h)
        parsed_hours = [f"{h:02d}" for h in range(start_i, end_i + 1)]
    else:
        parsed_hours = [f"{int(h):02d}" for h in hours_env.split(",") if h.strip() != ""]

    if parsed_hours:
        hours = parsed_hours
        print(f"Using GOES_HOURS override: {hours}")

timesteps_per_hour_env = os.environ.get("GOES_TIMESTEPS_PER_HOUR", "").strip()
timesteps_per_hour = int(timesteps_per_hour_env) if timesteps_per_hour_env else None
# Specify GOES ABI product, channel, lat/lon bounds, directory path for storing files
product = args.product
channel = args.channel  # e.g. 'C14' is the 11.2 micron channel, "Longwave window"
bounds = args.bounds  # Define a bounding box to crop to: [Lat_min Lat_max Lon_min Lon_max] (e.g. 30, 50, -125, -105 for western CONUS)
storage_path = args.dir  # Local path where data will be stored (to do: check if files already exist in this directory)

# Show us the bounds we'll crop images to
print("\nFiles will be downloaded and then cropped to these bounds:")
print(
    "\t({w},{n}).\t.({e},{n})\n\n\n\n\t({w},{s}).\t.({e},{s})\n".format(
        n=bounds[1], w=bounds[2], e=bounds[3], s=bounds[0]
    )
)


def _download_limited_timesteps(
    storage_path, bucket, year, month, day, hour, product, channel, max_files
):
    """Download only the first N GOES files for an hour instead of the full hour."""
    import datetime as dt

    day = int(day)
    hour = f"{int(hour):02d}"
    julian_day = dt.date(int(year), int(month), day).timetuple().tm_yday
    prefix = f"{product}/{year}/{julian_day:03d}/{hour}/"
    out_dir = os.path.join(
        storage_path,
        bucket[5:],
        str(year),
        str(month),
        str(day),
        product,
        hour,
        channel,
    )
    os.makedirs(out_dir, exist_ok=True)

    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    keys = sorted(
        item["Key"]
        for item in response.get("Contents", [])
        if item["Key"].endswith(".nc") and f"M6{channel}" in os.path.basename(item["Key"])
    )[:max_files]
    if not keys:
        raise FileNotFoundError(f"No GOES files found at s3://{bucket}/{prefix}")

    print(f"Downloading {len(keys)} timestep(s) from s3://{bucket}/{prefix}")
    for key in keys:
        out_path = os.path.join(out_dir, os.path.basename(key))
        if os.path.exists(out_path):
            continue
        s3.download_file(bucket, key, out_path)


##############################################################
# For each S3 bucket, download the corresponding observations if we don't have them already
filepath = []  # store filepaths of the files we download
print("For each S3 bucket, download the corresponding observations")
i = 0
for d in range(len(days)):
    for h in range(len(hours)):
        filepath.append(
            "{}/{}/{}/{}/{}/{}/{}/{}/".format(
                storage_path,
                satellite,
                year,
                month,
                days[d],
                product,
                hours[h],
                channel,
            )
        )
        if timesteps_per_hour:
            _download_limited_timesteps(
                storage_path,
                bucket,
                year,
                month,
                days[d],
                hours[h],
                product,
                channel,
                timesteps_per_hour,
            )
        else:
            if not os.path.exists(filepath[i]):
                ABI = ABI_Downloader(
                    storage_path, bucket, year, month, days[d], hours[h], product, channel
                )

        # now try and crop these so they don't take up so much space - this is very inefficient but oh well it's what I have right now
        if os.path.exists(
            filepath[i]
        ):  # we have to make sure the path exists (meaning we downloaded something) before running the subsetNetCDF function
            print("\nSubsetting files in...{}".format(filepath[i]))
            for file in glob(filepath[i] + "*.nc"):
                if file.endswith("_ortho.nc"):
                    print(f"Skipping already orthorectified file: {file}")
                    continue
                subsetNetCDF(file, bounds)
        i += 1


print("Done")
