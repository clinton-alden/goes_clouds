import utils
import os
import shutil

in_dir = f'/storage/cdalden/goes/washington/goes17/'
goes = 'goes17'
start_day = 26
end_day = 30
month = 9
year = 2022

for day in range(start_day, end_day + 1):
    date = f'{year}{month:02d}{day:02d}'
    utils.goes_rad_to_rgb(in_dir, date, goes, location='washington')
    print(f'Finished {date}')


    channel_list = ['C02', 'C05', 'C13']  # List of channels

    for channel in channel_list:
        # Construct the file path for the .zarr directory
        zarr_file = f'{in_dir}{channel}/{goes}_{channel}_washington_{date}.zarr'
        
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