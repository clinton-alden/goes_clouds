####################################################################################
# Useful functions developed for analysis of GOES data
# updated September 2025
# Clinton Alden
####################################################################################

import xarray as xr
import pandas as pd

####################################################################################


def daily_cloud_frequency(year, month, state, goes):
    '''
    Function to calculate daily cloud frequency from GOES RGB composite data

    Parameters
    ----------
    year : str
        e.g., '2022'
    month : str
        e.g., '07' for July
    state : str
        e.g., 'washington'
    goes : str
        e.g., 'goes16' or 'goes17'
    '''

    if month == '07' or month == '08':
        start_day = 1
        end_day = 32
    elif month == '06' or month == '09':
        start_day = 1
        end_day = 31
    
    for day in range(start_day, end_day):    
        date = f'{year}{month}' + str(day).zfill(2)
        path = f'/storage/cdalden/goes/{state}/{goes}/rgb_composite/'
        file = f'{goes}_C02_C05_C13_rgb_{state}_{date}.nc'.format(date=date)

        ds = xr.open_dataset(path + file)

        # Make mask for cloud/no cloud
        clouds = xr.where((ds.blue >= 0.13) & (ds.green >= 0.15), 1, 0) # cloud (1) if all conds are met, not cloud (0) otherwise
        ds['clouds'] = clouds

        # Select the time range between 0000-0300 and 1400-2359
        cloud_frequency = ds.clouds.sel(
            t=((ds['t'].dt.hour >= 0) & (ds['t'].dt.hour < 3)) | (ds['t'].dt.hour >= 14)
        ).sum(dim='t')
        cloud_counts = xr.Dataset(
            {'cloud_frequency': (['latitude', 'longitude'], cloud_frequency.data)},
            coords={'latitude': ds.latitude, 'longitude': ds.longitude}
        )

        # Save the cloud frequency data to a NetCDF file
        out_path = f'/storage/cdalden/goes/{state}/{goes}/cloud_counts/'
        out_file = f'{goes}_cloud_frequency_{state}_{date}.nc'.format(date=date)
        cloud_counts.to_netcdf(out_path + out_file)
        print(f"Processed and saved cloud frequency for {date}")



def process_monthly_data(year, month, state, goes, save_file=False):
    '''
    Function to calculate daily cloud frequency from GOES RGB composite data

    Parameters
    ----------
    year : str
        e.g., '2022'
    month : str
        e.g., '07' for July
    state : str
        e.g., 'washington'
    goes : str
        e.g., 'goes16' or 'goes17'
    save_file : bool, optional
        False by default for quicker analysis.
        If True, saves the processed monthly data to a NetCDF file.


    Usage
    ----------
     months_dict = {'06': 'jun_cloud_freq', '07': 'jul_cloud_freq', 
                   '08': 'aug_cloud_freq', '09': 'sep_cloud_freq'}
        # Dictionary to store the datasets
        cloud_freq_datasets = {}

        # Loop through the months and process the data
        for month, var_name in months_dict.items():
            ds = process_monthly_data(year, month, state, goes, save_file=True)  # Process the data for the month
            cloud_freq_datasets[var_name] = ds  # Save the dataset with the variable name as the key

        # Access the datasets
        jun_cloud_freq = cloud_freq_datasets['jun_cloud_freq']


    '''
    ####################################################################################

    
    if month == '07' or month == '08':
        start_day = 1
        end_day = 32
    elif month == '06' or month == '09':
        start_day = 1
        end_day = 31

    # Define the date range and file path
    path = f'/storage/cdalden/goes/{state}/{goes}/cloud_counts/'
    file_name_template = f'{goes}_cloud_frequency_{state}_' + '{date}.nc'

    # Generate the list of file paths and corresponding times
    file_paths = []
    times = []
    for day in range(start_day, end_day):  # Days 1 to 30
        date = f'{year}{month}{str(day).zfill(2)}'  # Format the date as '202209DD'
        file_paths.append(path + file_name_template.format(date=date))
        times.append(pd.Timestamp(f'{year}-{month}-{str(day).zfill(2)}'))  # Create a timestamp for each day

    # Open all files and add the time dimension
    datasets = []
    for file, time in zip(file_paths, times):
        ds = xr.open_dataset(file)
        ds = ds.expand_dims(time=[time])  # Add the time dimension
        datasets.append(ds)

    # Combine all datasets along the time dimension
    combined_ds = xr.concat(datasets, dim='time')

    # Compute the sum across the time dimension
    monthly_sum = combined_ds['cloud_frequency'].sum(dim='time')

    # Add the monthly sum as a new variable to the dataset
    combined_ds['monthly_sum'] = monthly_sum

    # monthly frequency 
    combined_ds['monthly_frequency'] = combined_ds['monthly_sum'] / (156*30)  # Assuming 156 observations per day for 30 days

    out_name = f'{goes}_monthly_cloud_frequency_{state}_{year}{month}.nc'
    if save_file:
        combined_ds.to_netcdf(path + out_name)
        print(f"Saved monthly cloud frequency data for month {month}")
    
    return combined_ds