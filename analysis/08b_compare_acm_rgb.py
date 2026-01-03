import xarray as xr
import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# # make a function to process ACM data to create daily CSV files with 3 binary cloud columns
# # ACM binary mask, TSI binary mask, and RGB-derived binary mask
# def concat_acm_rgb_files(year, month, day):
#     base_path = f'/storage/cdalden/goes/colorado/goes16/{year}/{month}/{day}/ABI-L2-ACMC/'
#     valid_hours = [str(hour) for hour in range(14, 24)]

#     # Collect all .nc files in the valid subdirectories
#     nc_files = [
#         os.path.join(root, file)
#         for root, dirs, files in os.walk(base_path)
#         for file in files
#         if os.path.basename(root) in valid_hours and file.endswith('.nc')
#     ]

#     # Sort the files to ensure proper time order
#     nc_files.sort()

#     # Use Dask for parallelized file reading and processing
#     combined_ds = xr.open_mfdataset(
#         nc_files,
#         concat_dim='t',  # Concatenate along 't' instead of 'time'
#         combine='nested',
#         preprocess=lambda ds: ds.expand_dims(t=[ds['t'].values]).drop_vars('time', errors='ignore'),
#         parallel=False
#     )

#     combined_ds = combined_ds.drop(['time', 'dem_px_angle_x', 'dem_px_angle_y'], errors='ignore')

#     path = '/storage/cdalden/goes/colorado/goes16/rgb_composite/'
#     file = f'combined_goes16_C02_C05_C13_rgb_colorado_{year}{month.zfill(2)}.nc'
#     rgb_ds = xr.open_dataset(path + file)
#     rgb_ds = rgb_ds.sel(t=slice(f'{year}-{month.zfill(2)}-{day.zfill(2)}'))

#     may_mask = ((rgb_ds['red'] <= 0.23) & (rgb_ds['blue'] >= 0.26)) | ((rgb_ds['red'] > 0.23) & (rgb_ds['blue'] >= 0.16))
#     summer_mask = (rgb_ds['red'] > 0.07) & (rgb_ds['green'] > 0.19)
#     winter_mask = (rgb_ds['red'] > 0.39) & (rgb_ds['blue'] > 0.13)

#     if month in ['1', '2']:
#         rgb_ds['rgb_cloud_frac'] = winter_mask.astype(int)
#     elif month in ['4', '5']:
#         rgb_ds['rgb_cloud_frac'] = may_mask.astype(int)
#     elif month in ['6', '7', '8']:
#         rgb_ds['rgb_cloud_frac'] = summer_mask.astype(int)

#     combined_ds_gothic = combined_ds.sel(
#         latitude=slice(39.065, 38.904),
#         longitude=slice(-107.08, -106.993))

#     goes_acm_gothic_ds = xr.where((combined_ds_gothic == 3) | (combined_ds_gothic == 2), 1, 0)
#     goes_acm_gothic_cf = goes_acm_gothic_ds.mean(dim=["latitude", "longitude"], skipna=True)

#     df_rgb = rgb_ds['rgb_cloud_frac'].mean(dim=["latitude", "longitude"], skipna=True).to_dataframe().reset_index()
#     df_goes = goes_acm_gothic_cf['ACM'].to_dataframe().reset_index()

#     df_opaque = pd.read_csv('/storage/cdalden/goes/surface_obs/sail_total_sky_imager/sail_tsi_cloud_frac.csv')

#     df_rgb = df_rgb.sort_values('t')
#     df_rgb = df_rgb[(df_rgb.t.dt.month == int(month)) & (df_rgb.t.dt.day == int(day))]
#     df_goes = df_goes.sort_values('t')
#     df_goes = df_goes[(df_goes.t.dt.month == int(month)) & (df_goes.t.dt.day == int(day))]
#     df_opaque['t'] = pd.to_datetime(df_opaque['t'])
#     df_opaque = df_opaque.sort_values('t')
#     df_opaque = df_opaque[(df_opaque.t.dt.month == int(month)) & (df_opaque.t.dt.day == int(day))]

#     # Merge on nearest time
#     df = pd.merge_asof(df_rgb, df_goes, on='t', direction='nearest', suffixes=('_rgb', '_goes'))
#     df = pd.merge_asof(df, df_opaque, on='t', direction='nearest')
#     df.set_index(pd.to_datetime(df.pop('t')), inplace=True)
#     df = df[df.index.hour >= 14]
#     df[df < 0] = np.nan # TSI sometimes has negative numbers so setting those to nan

#     out_path = '/storage/cdalden/goes/colorado/goes16/cloud_counts/east_river/'
#     out_name = f'east_river_cloud_fracs_{year}{month.zfill(2)}{day.zfill(2)}.csv'
#     df.to_csv(out_path + out_name, index=True)


# make a function to process ACM data to create daily CSV files with 3 binary cloud columns
# ACM binary mask, TSI binary mask, and RGB-derived binary mask
def concat_acm_rgb_files(year, month, day, domain, goes):
    base_path = f'/storage/cdalden/goes/{domain}/{goes}/{year}/{month}/{day}/ABI-L2-ACMC/'
    valid_hours = [str(hour) for hour in range(14, 24)]

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

    combined_ds = combined_ds.drop(['time', 'dem_px_angle_x', 'dem_px_angle_y'], errors='ignore')

    path = f'/storage/cdalden/goes/{domain}/{goes}/rgb_composite/'
    file = f'combined_{goes}_C02_C05_C13_rgb_{domain}_{year}{month.zfill(2)}.nc'
    rgb_ds = xr.open_dataset(path + file)
    rgb_ds = rgb_ds.sel(t=slice(f'{year}-{month.zfill(2)}-{day.zfill(2)}'))

    may_mask = ((rgb_ds['red'] <= 0.23) & (rgb_ds['blue'] >= 0.26)) | ((rgb_ds['red'] > 0.23) & (rgb_ds['blue'] >= 0.16))
    summer_mask = (rgb_ds['red'] > 0.07) & (rgb_ds['green'] > 0.19)
    winter_mask = (rgb_ds['red'] > 0.39) & (rgb_ds['blue'] > 0.13)
    scripps_mask = ((rgb_ds['red'] > 0) & (rgb_ds['green'] <= 0.13)) | (rgb_ds['green'] > 0.13)
    
    if domain == 'colorado':
        if month in ['1', '2']:
            rgb_ds['rgb_cloud_frac'] = winter_mask.astype(int)
        elif month in ['4', '5']:
            rgb_ds['rgb_cloud_frac'] = may_mask.astype(int)
        elif month in ['6', '7', '8']:
            rgb_ds['rgb_cloud_frac'] = summer_mask.astype(int)
        
        df_opaque = pd.read_csv('/storage/cdalden/goes/surface_obs/sail_total_sky_imager/sail_tsi_cloud_frac.csv')

        combined_ds_gothic = combined_ds.sel(
            latitude=slice(39.065, 38.904),
            longitude=slice(-107.08, -106.993))
    
    elif domain == 'scripps':
        df_opaque = pd.read_csv('/storage/cdalden/goes/surface_obs/scripps_total_sky_imager/scripps_tsi_cloud_frac.csv')

        rgb_ds['rgb_cloud_frac'] = scripps_mask.astype(int)

        combined_ds_gothic = combined_ds.sel(
            latitude=slice(32.97, 32.83),
            longitude=slice(-117.31, -117.23))
        

    goes_acm_gothic_ds = xr.where((combined_ds_gothic == 3) | (combined_ds_gothic == 2), 1, 0)
    goes_acm_gothic_cf = goes_acm_gothic_ds.mean(dim=["latitude", "longitude"], skipna=True)
    df_rgb = rgb_ds['rgb_cloud_frac'].mean(dim=["latitude", "longitude"], skipna=True).to_dataframe().reset_index()
    df_goes = goes_acm_gothic_cf['ACM'].to_dataframe().reset_index()



    df_rgb = df_rgb.sort_values('t')
    df_rgb = df_rgb[(df_rgb.t.dt.month == int(month)) & (df_rgb.t.dt.day == int(day))]
    df_goes = df_goes.sort_values('t')
    df_goes = df_goes[(df_goes.t.dt.month == int(month)) & (df_goes.t.dt.day == int(day))]
    df_opaque['t'] = pd.to_datetime(df_opaque['t'])
    df_opaque = df_opaque.sort_values('t')
    df_opaque = df_opaque[(df_opaque.t.dt.month == int(month)) & (df_opaque.t.dt.day == int(day))]

    # Merge on nearest time
    df = pd.merge_asof(df_rgb, df_goes, on='t', direction='nearest', suffixes=('_rgb', '_goes'))
    df = pd.merge_asof(df, df_opaque, on='t', direction='nearest')
    df.set_index(pd.to_datetime(df.pop('t')), inplace=True)
    df = df[df.index.hour >= 14]
    df[df < 0] = np.nan # TSI sometimes has negative numbers so setting those to nan

    if domain == 'colorado':
        out_path = f'/storage/cdalden/goes/{domain}/{goes}/cloud_counts/east_river/'
        out_name = f'east_river_cloud_fracs_{year}{month.zfill(2)}{day.zfill(2)}.csv'
    elif domain == 'scripps':
        out_path = f'/storage/cdalden/goes/{domain}/{goes}/cloud_counts/'
        out_name = f'scripps_cloud_fracs_{year}{month.zfill(2)}{day.zfill(2)}.csv'
    df.to_csv(out_path + out_name, index=True)




# loop through s?,21):
for i in range(21,31):
    year = str(2022)
    month = str(6)
    day = str(i)

    concat_acm_rgb_files(year, month, day, 'colorado', 'goes16')
    print(f'Processing complete for {year}{month.zfill(2)}{day.zfill(2)} 🌟')