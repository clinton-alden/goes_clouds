import utils
import os

in_dir = f'/storage/cdalden/goes/scripps/'
channel_list = ['C02', 'C05', 'C13']
start_date = 1
end_date = 25
month = 8
year = 2023
goes_model = 'goes18'


utils.goes_nc_to_zarr(in_dir, channel_list, start_date, end_date, month, year,
                       'scripps', goes_model, surprise=True)

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