import xarray as xr
import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# make a function to process ACM data to create daily CSV files with 3 binary cloud columns
# ACM binary mask, TSI binary mask, and RGB-derived binary mask
def concat_acm_rgb_files(year, month, day):
    base_path = f'/storage/cdalden/goes/colorado/goes16/{year}/{month}/{day}/ABI-L2-ACMC/'
    valid_hours = [str(hour) for hour in range(16, 24)]

    # Collect all .nc files in the valid subdirectories
    nc_files = [
        os.path.join(root, file)
        for root, dirs, files in os.walk(base_path)
        for file in files
        if os.path.basename(root) in valid_hours and file.endswith('.nc')
    ]

    # Sort the files to ensure proper time order
    nc_files.sort()

    # Use Dask for parallelized file reading and processing
    combined_ds = xr.open_mfdataset(
        nc_files,
        concat_dim='t',  # Concatenate along 't' instead of 'time'
        combine='nested',
        preprocess=lambda ds: ds.expand_dims(t=[ds['t'].values]).drop_vars('time', errors='ignore'),
        parallel=False
    )

    # Display the resulting dataset
    combined_ds = combined_ds['ACM']
    combined_ds = combined_ds.drop(['time', 'dem_px_angle_x', 'dem_px_angle_y'], errors='ignore')

    tsi_df = pd.read_csv('/storage/cdalden/goes/colorado/goes16/rgb_composite/rgb_tsi_east_river_202205_202305.csv')
    tsi_df['time'] = pd.to_datetime(tsi_df['time'])
    tsi_df = tsi_df[
        (tsi_df['time'].dt.date == pd.Timestamp(f'{year}-{month.zfill(2)}-{day.zfill(2)}').date()) &
        (tsi_df['time'].dt.hour >= 16)
    ]

    # Define masks for the DataFrame
    may_mask = ((tsi_df['red'] <= 0.23) & (tsi_df['blue'] >= 0.26)) | ((tsi_df['red'] > 0.23) & (tsi_df['blue'] >= 0.16))
    summer_mask = (tsi_df['red'] > 0.07) & (tsi_df['green'] > 0.19)
    winter_mask = (tsi_df['red'] > 0.39) & (tsi_df['blue'] > 0.13)

    # Apply the appropriate mask based on the month
    if month in ['1', '2']:
        tsi_df['rgb_cloud_binary'] = winter_mask.astype(int)
    elif month in ['4', '5']:
        tsi_df['rgb_cloud_binary'] = may_mask.astype(int)
    elif month in ['6', '7', '8']:
        tsi_df['rgb_cloud_binary'] = summer_mask.astype(int)
    else:
        tsi_df['rgb_cloud_binary'] = 0 # cloud (1) if all conds are met, not cloud (0) otherwise


    combined_ds_gothic = combined_ds.sel(
        latitude=slice(39.065, 38.904),
        longitude=slice(-107.08, -106.993))

    goes_acm_gothic_ds = xr.where((combined_ds_gothic == 3) | (combined_ds_gothic == 2), 1, 0)

    # Calculate the spatial average over latitude and longitude
    spatial_avg = goes_acm_gothic_ds.mean(dim=["latitude", "longitude"])

    # Convert spatial_avg to a Pandas Series or DataFrame to align with tsi_df
    spatial_avg_df = spatial_avg.to_dataframe().reset_index()

    # Merge spatial_avg with tsi_df to ensure alignment
    spatial_avg_df['t'] = pd.to_datetime(spatial_avg_df['t'])
    tsi_df = tsi_df.sort_values(by='time')
    spatial_avg_df = spatial_avg_df.sort_values(by='t')

    # Merge spatial_avg_df into tsi_df based on the nearest time
    tsi_df = pd.merge_asof(tsi_df, spatial_avg_df, left_on='time', right_on='t', direction='nearest')

    # Initialize the column with default values
    tsi_df['goes_acm_cloud'] = 4  # Default to 4 for values between 0.1 and 0.7

    # Assign 1 if spatial average > 0.7
    tsi_df.loc[tsi_df['ACM'] > 0.7, 'goes_acm_cloud'] = 1

    # Assign 0 if spatial average < 0.1
    tsi_df.loc[tsi_df['ACM'] < 0.1, 'goes_acm_cloud'] = 0

    tsi_df = tsi_df.drop(columns=['t','ACM'])

    out_path = '/storage/cdalden/goes/colorado/goes16/east_river_cloud_binary/'
    out_name = f'east_river_cloud_binary_{year}{month.zfill(2)}{day.zfill(2)}.csv'
    tsi_df.to_csv(out_path + out_name)


# loop through February 2023
for i in range(1,6):
    year = str(2022)
    month = str(5)
    day = str(i)
    
    concat_acm_rgb_files(year, month, day)
    print(f'Processing complete for {year}{month.zfill(2)}{day.zfill(2)} 🌟')