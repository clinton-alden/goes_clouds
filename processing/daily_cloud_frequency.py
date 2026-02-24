import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from pathlib import Path
import sys
import glob
import shutil

####################################################################################


def calculate_hourly_cloud_frequency(base_dir, year, month, domain, goes):
    '''
    Function to calculate hourly cloud frequency from GOES RGB composite data
    and delete RGB files after processing

    Parameters
    ----------
    base_dir : str
        Base directory path, e.g., '/storage/cdalden/goes/colorado/'
    year : str
        e.g., '2022'
    month : str
        e.g., '07' for July (can be '7' or '07')
    domain : str
        e.g., 'colorado', 'washington'
    goes : str
        e.g., 'goes16' or 'goes17'
    '''
    
    # Normalize inputs
    year = str(year)
    month = str(int(month)).zfill(2)
    
    # Calculate number of days in month
    if month in ['01', '03', '05', '07', '08', '10', '12']:
        num_days = 31
    elif month in ['04', '06', '09', '11']:
        num_days = 30
    elif month == '02':
        # Check for leap year
        num_days = 29 if (int(year) % 4 == 0 and (int(year) % 100 != 0 or int(year) % 400 == 0)) else 28
    
    rgb_composite_dir = Path(base_dir) / goes / 'rgb_composite'
    cloud_counts_dir = Path(base_dir) / goes / 'cloud_frequency'
    
    # Create output directory if it doesn't exist
    cloud_counts_dir.mkdir(parents=True, exist_ok=True)
    
    for day in range(1, num_days + 1):
        day_str = str(day).zfill(2)
        date = f'{year}{month}{day_str}'
        
        # Look for RGB composite file
        rgb_file_pattern = rgb_composite_dir / f'{goes}_C02_C05_C13_rgb_{domain}_{date}.nc'
        
        if not rgb_file_pattern.exists():
            print(f"RGB file not found: {rgb_file_pattern}")
            continue
        
        try:
            ds = xr.open_dataset(str(rgb_file_pattern))
            
            # Make mask for cloud/no cloud
            clouds = xr.where((ds.blue >= 0.13) & (ds.green >= 0.15), 1, 0)
            ds['clouds'] = clouds
            
            # Calculate hourly cloud frequency (group by hour and sum)
            hourly_freq = ds.clouds.groupby('t.hour').sum(dim='t')
            
            # Create output dataset
            cloud_counts = xr.Dataset(
                {'cloud_frequency': (['hour', 'latitude', 'longitude'], hourly_freq.data)},
                coords={
                    'hour': hourly_freq.hour,
                    'latitude': ds.latitude,
                    'longitude': ds.longitude
                }
            )
            
            # Save the hourly cloud frequency data to a NetCDF file
            out_file = cloud_counts_dir / f'{goes}_cloud_frequency_{domain}_{date}.nc'
            cloud_counts.to_netcdf(str(out_file))
            print(f"Processed and saved hourly cloud frequency for {date}")
            
        except Exception as e:
            print(f"ERROR processing {date}: {e}")
            continue
    
    # Delete RGB composite files after processing
    # print(f"Deleting RGB composite files in {rgb_composite_dir}")
    # try:
    #     # Remove the entire rgb_composite directory
    #     if rgb_composite_dir.exists():
    #         shutil.rmtree(rgb_composite_dir)
    #         print(f"Successfully deleted RGB composite directory")
    # except Exception as e:
    #     print(f"ERROR deleting RGB files: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python daily_cloud_frequency.py <base_dir> <year> <month> <domain> <goes>")
        print("Example: python daily_cloud_frequency.py /storage/cdalden/goes/colorado/ 2022 7 colorado goes16")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    year = sys.argv[2]
    month = sys.argv[3]
    domain = sys.argv[4]
    goes = sys.argv[5]
    
    calculate_hourly_cloud_frequency(base_dir, year, month, domain, goes)