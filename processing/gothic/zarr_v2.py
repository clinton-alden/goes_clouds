import utils
import os
import sys
import calendar

# Gothic test version

# Check for command-line arguments
if len(sys.argv) not in (6, 8):
    print("Usage: python zarr_v2.py <directory> <year> <month> <goes_model> <domain> [start_day end_day]")
    sys.exit(1)

# Parse command-line arguments
in_dir = sys.argv[1]
year = int(sys.argv[2])
month = int(sys.argv[3])
goes_model = sys.argv[4]
domain = sys.argv[5]

# Default to the full month; optional bounds support isolated recovery retries.
num_days = calendar.monthrange(year, month)[1]
start_date = int(sys.argv[6]) if len(sys.argv) == 8 else 1
end_date = int(sys.argv[7]) if len(sys.argv) == 8 else num_days
if not (1 <= start_date <= end_date <= num_days):
    raise ValueError(f"Invalid day range {start_date}-{end_date} for {year}-{month:02d}")

channel_list = ['C02', 'C05', 'C13']

print(f"Processing {goes_model} for {year}-{month:02d} (days {start_date}-{end_date})")
print(f"Domain: {domain}")

utils.goes_nc_to_zarr(in_dir, channel_list, start_date, end_date, month, year,
                       domain, goes_model, surprise=True)

# Delete netcdf files after saving the zarrs
for day in range(start_date, end_date + 1):
    day_of_month = f'{day:02d}'
    nc_dir = in_dir + f'{goes_model}/{year}/{month}/{str(day)}/'
    # Recursively remove all .nc files in the directory
    for root, dirs, files in os.walk(nc_dir):
        print(f"Inspecting directory: {root}")
        print(f"Files found: {files}")
        for file in files:
            if file.endswith('.nc'):
                file_path = os.path.join(root, file)
                try:
                    os.remove(file_path)
                    print(f"Deleted: {file_path}")
                except Exception as e:
                    print(f"Failed to delete {file_path}: {e}")
                    raise
