#!/usr/bin/env python3
"""
Convert orthorectified ACM (Aerosol Cloud Mask) NetCDF files to daily Zarr format
Adapted from zarr_v2.py for L2-ACMC product
"""

import os
import sys
import glob
import gc
import shutil
import warnings
from datetime import datetime
from pathlib import Path

import xarray as xr

warnings.filterwarnings("ignore")

def acm_nc_to_zarr(in_dir, year, month, goes_model, domain):
    """
    Convert orthorectified ACM NetCDF files to daily Zarr files.
    
    Parameters
    ----------
    in_dir : str
        Base input directory containing GOES data
    year : int
        Year to process
    month : int
        Month to process (1-12)
    goes_model : str
        GOES satellite model (e.g., 'goes16', 'goes18')
    domain : str
        Domain identifier (e.g., 'gothic_acm', 'mammoth_acm')
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    
    # Determine number of days in month
    import calendar
    num_days = calendar.monthrange(year, month)[1]
    
    print(f"Converting ACM to zarr for {year}-{month:02d} ({num_days} days)")
    
    success_count = 0
    fail_count = 0
    
    for day in range(1, num_days + 1):
        day_str = f"{day:02d}"
        month_str = f"{month:02d}"
        
        # Path to orthorectified files for this day
        # Structure: goes_model/year/month/day/ABI-L2-ACMC/...
        nc_dir = os.path.join(in_dir, goes_model, str(year), str(month), day_str)
        
        if not os.path.exists(nc_dir):
            print(f"  {year}-{month_str}-{day_str}: Directory not found, skipping")
            continue
        
        # Find all orthorectified ACM files for this day
        nc_files = []
        for root, dirs, files in os.walk(os.path.join(nc_dir, 'ABI-L2-ACMC')):
            for file in files:
                if file.endswith('_ortho.nc'):
                    nc_files.append(os.path.join(root, file))
        
        if not nc_files:
            print(f"  {year}-{month_str}-{day_str}: No ortho files found")
            continue
        
        nc_files = sorted(nc_files)
        existing_nc_files = [f for f in nc_files if os.path.exists(f)]
        missing_count = len(nc_files) - len(existing_nc_files)
        
        if missing_count > 0:
            print(f"  {year}-{month_str}-{day_str}: WARNING {missing_count} files missing before open")
        
        if not existing_nc_files:
            print(f"  {year}-{month_str}-{day_str}: No readable files, skipping")
            continue
        
        out_name = f'{goes_model}_acm_{domain}_{year}{month_str}{day_str}.zarr'
        out_path = os.path.join(in_dir, goes_model, 'ACM', out_name)
        
        # Remove existing zarr
        if os.path.isdir(out_path):
            shutil.rmtree(out_path)
        
        try:
            wrote_any = False
            
            # Stream writes file-by-file to avoid keeping full day in memory
            for f in existing_nc_files:
                try:
                    with xr.open_dataset(f) as ds_single:
                        # Ensure time dimension exists
                        if 't' not in ds_single.dims:
                            ds_single = ds_single.expand_dims('t')
                        
                        # Drop unnecessary variables
                        drop_vars = [
                            'dem_px_angle_x', 'dem_px_angle_y',
                            'solar_zenith_angle', 'solar_azimuth_angle'
                        ]
                        ds_single = ds_single.drop_vars([v for v in drop_vars if v in ds_single.data_vars], errors='ignore')
                        
                        if not wrote_any:
                            ds_single.to_zarr(out_path, mode='w')
                            wrote_any = True
                            print(f"  {year}-{month_str}-{day_str}: Created {out_name}")
                        else:
                            ds_single.to_zarr(out_path, mode='a', append_dim='t')
                
                except FileNotFoundError:
                    print(f"    WARNING: File missing: {Path(f).name}")
                except Exception as e:
                    print(f"    WARNING: Failed to process {Path(f).name}: {e}")
                    fail_count += 1
            
            if wrote_any:
                success_count += 1
                gc.collect()
            else:
                print(f"  {year}-{month_str}-{day_str}: WARNING No datasets could be opened")
                fail_count += 1
        
        except Exception as e:
            print(f"  {year}-{month_str}-{day_str}: ERROR {e}")
            fail_count += 1
    
    print(f"ACM zarr conversion complete: {success_count} days succeeded, {fail_count} failed")
    return fail_count == 0


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <in_dir> <year> <month> <goes_model> <domain>")
        print(f"Example: {sys.argv[0]} /glade/derecho/scratch/cdalden/gothic_acm 2021 1 goes16 gothic_acm")
        sys.exit(1)
    
    in_dir = sys.argv[1]
    year = int(sys.argv[2])
    month = int(sys.argv[3])
    goes_model = sys.argv[4]
    domain = sys.argv[5]
    
    if not os.path.isdir(in_dir):
        print(f"ERROR: Input directory does not exist: {in_dir}")
        sys.exit(1)
    
    success = acm_nc_to_zarr(in_dir, year, month, goes_model, domain)
    sys.exit(0 if success else 1)
