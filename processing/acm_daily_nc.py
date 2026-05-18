#!/usr/bin/env python3
"""Combine daily orthorectified ACM files into one NetCDF per day."""

import argparse
import calendar
import os
import re
from datetime import datetime

from netCDF4 import Dataset, date2num, num2date


STATIC_VARS = ("latitude", "longitude", "elevation", "abi_fixed_grid_x", "abi_fixed_grid_y")
OPTIONAL_STATIC_VARS = ("zone_labels",)
TIME_VAR = "ACM"
START_RE = re.compile(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write daily ACM NetCDF files from orthorectified 5-minute inputs."
    )
    parser.add_argument("--base-dir", required=True, help="Base GOES output directory.")
    parser.add_argument("--year", required=True, type=int)
    parser.add_argument("--month", required=True, type=int)
    parser.add_argument("--goes", required=True, help="GOES token, for example goes16.")
    parser.add_argument("--domain", required=True, help="Domain token for output names.")
    parser.add_argument(
        "--product",
        default="ABI-L2-ACMC",
        help="Product directory containing orthorectified files.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional override for daily NetCDF output directory.",
    )
    return parser.parse_args()


def day_input_dir(base_dir, goes, year, month, day, product):
    return os.path.join(base_dir, goes, str(year), str(month), str(day), product)


def find_ortho_files(input_dir):
    if not os.path.isdir(input_dir):
        return []

    files = []
    for root, _, names in os.walk(input_dir):
        for name in names:
            if name.endswith("_ortho.nc"):
                files.append(os.path.join(root, name))
    return sorted(files)


def read_source_time(nc_path):
    match = START_RE.search(os.path.basename(nc_path))
    if match is None:
        raise ValueError(f"Could not parse scan start time from filename: {nc_path}")

    year, julian_day, hour, minute, second = map(int, match.groups())
    return datetime.strptime(
        f"{year:04d} {julian_day:03d} {hour:02d} {minute:02d} {second:02d}",
        "%Y %j %H %M %S",
    )


def collect_valid_files(files):
    timed_files = []
    for path in files:
        try:
            timed_files.append((read_source_time(path), path))
        except Exception as exc:
            print(f"Skipping malformed ACM file: {path} ({exc})")
    return sorted(timed_files, key=lambda item: item[0])


def copy_attributes(src_var, dst_var, skip=None):
    skip = set(skip or ())
    for attr in src_var.ncattrs():
        if attr in skip:
            continue
        dst_var.setncattr(attr, src_var.getncattr(attr))


def init_output_dataset(out_path, first_file, year, month, day):
    with Dataset(first_file) as src:
        out = Dataset(out_path, "w", format="NETCDF4")

        lat_len = len(src.dimensions["latitude"])
        lon_len = len(src.dimensions["longitude"])
        out.createDimension("t", None)
        out.createDimension("latitude", lat_len)
        out.createDimension("longitude", lon_len)

        time_units = f"seconds since {year:04d}-{month:02d}-{day:02d} 00:00:00"
        time_var = out.createVariable("t", "f8", ("t",))
        time_var.units = time_units
        time_var.calendar = "proleptic_gregorian"

        lat_src = src.variables["latitude"]
        lat_dst = out.createVariable("latitude", lat_src.datatype, ("latitude",))
        lat_dst[:] = lat_src[:]
        copy_attributes(lat_src, lat_dst, skip={"_FillValue"})

        lon_src = src.variables["longitude"]
        lon_dst = out.createVariable("longitude", lon_src.datatype, ("longitude",))
        lon_dst[:] = lon_src[:]
        copy_attributes(lon_src, lon_dst, skip={"_FillValue"})

        for name in STATIC_VARS[2:] + OPTIONAL_STATIC_VARS:
            if name not in src.variables:
                continue
            src_var = src.variables[name]
            dims = src_var.dimensions
            fill_value = getattr(src_var, "_FillValue", None)
            dst_var = out.createVariable(
                name,
                src_var.datatype,
                dims,
                zlib=True,
                complevel=4,
                fill_value=fill_value,
            )
            dst_var[:] = src_var[:]
            copy_attributes(src_var, dst_var, skip={"_FillValue"})

        src_var = src.variables[TIME_VAR]
        fill_value = getattr(src_var, "_FillValue", None)
        dst_var = out.createVariable(
            TIME_VAR,
            src_var.datatype,
            ("t", "latitude", "longitude"),
            zlib=True,
            complevel=4,
            fill_value=fill_value,
        )
        copy_attributes(src_var, dst_var, skip={"_FillValue", "coordinates"})
        dst_var.coordinates = "t latitude longitude"

        for attr in src.ncattrs():
            if attr == "history":
                continue
            out.setncattr(attr, src.getncattr(attr))
        out.setncattr("source_product", "ABI-L2-ACMC orthorectified files")
        out.setncattr(
            "history",
            f"Combined daily from orthorectified 5-minute ACM files on {datetime.utcnow().isoformat()}Z",
        )

        return out


def output_path(base_dir, goes, domain, year, month, day, output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(base_dir, goes, "daily_nc")

    os.makedirs(output_dir, exist_ok=True)
    filename = f"{goes}_acm_{domain}_{year:04d}{month:02d}{day:02d}.nc"
    return os.path.join(output_dir, filename)


def append_one_timestep(out_ds, source_path, time_index):
    with Dataset(source_path) as src:
        if "t" not in src.variables or TIME_VAR not in src.variables:
            raise KeyError(f"Required variables missing in {source_path}")
        time_var = src.variables["t"]
        source_time = num2date(
            time_var[:].item(),
            time_var.units,
            getattr(time_var, "calendar", "standard"),
        )
        out_ds.variables["t"][time_index] = date2num(
            source_time,
            out_ds.variables["t"].units,
            out_ds.variables["t"].calendar,
        )
        out_ds.variables[TIME_VAR][time_index, :, :] = src.variables[TIME_VAR][:]


def write_day(base_dir, goes, domain, year, month, day, product, output_dir=None):
    input_dir = day_input_dir(base_dir, goes, year, month, day, product)
    timed_files = collect_valid_files(find_ortho_files(input_dir))
    date_str = f"{year:04d}-{month:02d}-{day:02d}"

    if not timed_files:
        print(f"{date_str}: no ortho ACM files found in {input_dir}")
        return False

    out_path = output_path(base_dir, goes, domain, year, month, day, output_dir)
    print(f"{date_str}: combining {len(timed_files)} ortho files")

    out_ds = init_output_dataset(out_path, timed_files[0][1], year, month, day)
    try:
        write_index = 0
        for _, source_path in timed_files:
            try:
                append_one_timestep(out_ds, source_path, write_index)
                write_index += 1
            except Exception as exc:
                print(f"Skipping malformed ACM file during write: {source_path} ({exc})")
    finally:
        out_ds.close()

    print(f"{date_str}: wrote {out_path}")
    return True


def main():
    args = parse_args()
    num_days = calendar.monthrange(args.year, args.month)[1]

    success_count = 0
    for day in range(1, num_days + 1):
        wrote = write_day(
            args.base_dir,
            args.goes,
            args.domain,
            args.year,
            args.month,
            day,
            args.product,
            args.output_dir,
        )
        success_count += int(wrote)

    print(
        f"Finished daily ACM NetCDF packaging for {args.year:04d}-{args.month:02d}: "
        f"{success_count}/{num_days} daily files written"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
