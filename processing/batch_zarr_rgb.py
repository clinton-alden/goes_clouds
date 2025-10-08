import utils
import os
import shutil

in_dir = f'/storage/cdalden/goes/colorado/'
channel_list = ['C02', 'C05', 'C13']
start_date = 2
end_date = 3
month = 5
year = 2023
goes_model = 'goes16'
domain = 'colorado'

for day in range(start_date, end_date + 1):
    date = f'{year}{month:02d}{day:02d}'
    print(f'Starting zarr for {date}')
    try:
        utils.goes_nc_to_zarr(in_dir, channel_list, start_date, end_date, month, year,
                        'colorado', goes_model, surprise=True)
    except Exception as e:
        print(f"Error during zarr for {date}: {e}")
        continue

    # Delete netcdf files after saving the zarrs
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


    print(f'Starting RGB for {date}')
    try:
        utils.goes_rad_to_rgb(in_dir, date, goes_model, location=domain)
    except Exception as e:
        print(f"Error during RGB for {date}: {e}")
        continue
    print(f'Finished RGB for {date} 🟥🟩🟦')


    channel_list = ['C02', 'C05', 'C13']  # List of channels

    for channel in channel_list:
        # Construct the file path for the .zarr directory
        zarr_file = f'{in_dir}{channel}/{goes_model}_{channel}_{domain}_{date}.zarr'
        
        # Check if the path exists
        if os.path.exists(zarr_file):
            try:
                # Use shutil.rmtree to delete the directory
                shutil.rmtree(zarr_file)
                print(f'Deleted: {zarr_file}')
            except Exception as e:
                print(f'Failed to delete {zarr_file}: {e}')
        else:
            print(f'File not found: {zarr_file}')