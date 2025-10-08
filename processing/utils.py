# Infrastructure for processing GOES ABI data and creating the Day Cloud Phase RGB Composites
import xarray as xr
import numpy as np
import pandas as pd
import os
import gc
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import warnings
import imageio

warnings.filterwarnings("ignore")

#########################################################################################
# Processing to netcdfs to zarr files for each day
#########################################################################################

def goes_nc_to_zarr(in_dir, channels, startday, endday, month, year, 
                    location='washington', goes_model='goes17', surprise=False):
    """
    Convert multiple NetCDF files from a directory to a single Zarr file.
    Parameters
    ----------
    in_dir : str
        Input directory containing NetCDF files.
    channels : list
        List of channels to process (e.g., ['C13', 'C02', 'C05']).
    startday : int
        Starting day of the month for processing.
    endday : int
        Ending day of the month for processing.
    month : int
        Month of the data to process.
    year : int
        Year of the data to process.
    location : str
        Location identifier (e.g., 'washington').
    goes_model : str
        GOES satellite model (e.g., 'goes17').
    surprise : bool
        If True, print fun messages during processing.

    Example usage:
    >>> goes_nc_to_zarr('/storage/cdalden/goes/washington/', ['C13', 'C02', 'C05'], 
                        1, 2, 8, 2022, 'washington', 'goes17')
    """ 

    # if no date range is specified, run for full month
    if not startday:
        startday = 1
    if not endday:
        endday = 31

    # this is adds 0 to single digit days of the months
    # GOES sub dirs don't have it but it is needed for the rest of my processing
    for day in range(startday,endday+1):
        day_of_month = str(day)
        if day<10:
            out_day_of_month = '0' + day_of_month
        else:
            out_day_of_month = day_of_month
        # print(f'Starting {str(month)}/{day_of_month}')
        
        # Directory containing the NetCDF files
        nc_dir = in_dir + f'{goes_model}/{year}/{month}/{day_of_month}/'

        # Check if the directory exists
        if not os.path.exists(nc_dir):
            print(f'Directory {nc_dir} does not exist. Skipping...')
            continue

        # count number of files in nc_dir
        num_files = len([f for root, dirs, files in os.walk(nc_dir) for f in files if f.endswith('.nc')])
        expected_num_files = len(channels)*288  # 288 files per channel per day
        # if num_files != expected_num_files:
        #     print(f'WARNING: {num_files} files found in {nc_dir}, expected {expected_num_files}. Skipping...')
        #     continue
        if num_files != expected_num_files:
            print(f'WARNING: {num_files} files found in {nc_dir}, expected {expected_num_files}. Not skipping, just a warning')

        # loop through all needed channels
        for channel in channels:
            out_name = f'{str(goes_model)}_{channel}_{location}_{str(year)}0{str(month)}{out_day_of_month}.zarr'
            print(f'Processing {out_name}...')
            
            # Recursively list all NetCDF files in the directory and subdirectories
            nc_files = []
            for root, dirs, files in os.walk(nc_dir):
                for file in files:
                    # print(file)
                    if file.endswith('.nc') and channel in file:
                        nc_files.append(os.path.join(root, file))

            # Open multiple NetCDF files as a list of datasets
            datasets = [xr.open_dataset(f) for f in nc_files]


            # Concatenate datasets along the 't' coordinate
            combined_ds = xr.concat(datasets, dim='t')
            combined_ds = combined_ds.drop_vars(['dem_px_angle_x', 'dem_px_angle_y'])

            # Save the combined dataset to a Zarr file
            out_name = out_name
            combined_ds.to_zarr(in_dir+f'{goes_model}/{channel}/'+out_name)

            # Force garbage collection to free memory
            gc.collect()
            # if surprise:
                # print('Beep beep, here comes the garbage truck! 🚛')
            print('Finished')
            if surprise:
                print('🛰️')


#########################################################################################
# Processing to convert to RGB composite
#########################################################################################

# Calculate latitude and longitude from GOES ABI fixed grid projection data
# GOES ABI fixed grid projection is a map projection relative to the GOES satellite
# Units: latitude in °N (°S < 0), longitude in °E (°W < 0)
# See GOES-R Product User Guide (PUG) Volume 5 (L2 products) Section 4.2.8 for 
# details & example of calculations
# "file_id" is an ABI L1b or L2 .nc file opened using the netCDF4 library

# Clinton's note: I needed to change how attributes of the goes imager projection were accessed

def calculate_degrees(file_id):
    
    # Read in GOES ABI fixed grid projection variables and constants
    x_coordinate_1d = file_id.variables['x'][:]  # E/W scanning angle in radians
    y_coordinate_1d = file_id.variables['y'][:]  # N/S elevation angle in radians
    projection_info = file_id.variables['goes_imager_projection']
    lon_origin = projection_info.attrs['longitude_of_projection_origin']
    H = projection_info.attrs['perspective_point_height']+projection_info.attrs['semi_major_axis']
    r_eq = projection_info.attrs['semi_major_axis']
    r_pol = projection_info.attrs['semi_minor_axis']
    
    # Create 2D coordinate matrices from 1D coordinate vectors
    x_coordinate_2d, y_coordinate_2d = np.meshgrid(x_coordinate_1d, y_coordinate_1d)
    
    # Equations to calculate latitude and longitude
    lambda_0 = (lon_origin*np.pi)/180.0  
    a_var = np.power(np.sin(x_coordinate_2d),2.0) + (np.power(np.cos(x_coordinate_2d),2.0)*(np.power(np.cos(y_coordinate_2d),2.0)+(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(y_coordinate_2d),2.0))))
    b_var = -2.0*H*np.cos(x_coordinate_2d)*np.cos(y_coordinate_2d)
    c_var = (H**2.0)-(r_eq**2.0)
    r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)
    s_x = r_s*np.cos(x_coordinate_2d)*np.cos(y_coordinate_2d)
    s_y = - r_s*np.sin(x_coordinate_2d)
    s_z = r_s*np.cos(x_coordinate_2d)*np.sin(y_coordinate_2d)
    
    # Ignore numpy errors for sqrt of negative number; occurs for GOES-16 ABI CONUS sector data
    np.seterr(all='ignore')
    
    abi_lat = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
    abi_lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)
    
    return abi_lat, abi_lon

# From goes2go - annoying warnings with import so copying here manually
def goes_norm(value, lower_limit, upper_limit, clip=True, invert_norm=False):
    """
    Normalize values between an upper and lower limit between 0 and 1.

    # Brightness temperatures need to be inverted so 
    #   lowest vals are warm temps (blue and green when plotted)
    #   and highest vals are cold temps (red when plotted)

    Normalize between a lower and upper limit. In other words, it
    converts your number to a value in the range between 0 and 1.
    Follows `normalization formula
    <https://stats.stackexchange.com/a/70807/220885>`_

    This is the same concept as `contrast or histogram stretching
    <https://staff.fnwi.uva.nl/r.vandenboomgaard/IPCV20162017/LectureNotes/IP/PointOperators/ImageStretching.html>`_


    .. code:: python

        _normalizedValue = (OriginalValue-LowerLimit)/(UpperLimit-LowerLimit)

    Parameters
    ----------
    value :
        The original value. A single value, vector, or array.
    upper_limit :
        The upper limit.
    lower_limit :
        The lower limit.
    clip : bool
        - True: Clips values between 0 and 1 for RGB.
        - False: Retain the numbers that extends outside 0-1 range.
    """
    norm = (value - lower_limit) / (upper_limit - lower_limit)
    if clip:
        norm = np.clip(norm, 0, 1)
    if invert_norm:
        norm = 1-norm
    return norm

def radiance_to_brightness_temp(ds, band):
    """
    Convert radiance to brightness temperature using the inverse Planck formula.

    ds: xarray.Dataset with radiance data

    band: str, the band number (e.g., '13' for channel 13)

    """
    # Constants
    c1 = 1.191042e-5    # W/m^2/sr/m^4
    c2 = 1.4387752      # cm*K

    # Wavelength in meters for the specified band
    channel_dict = {'13':10.4e-4}
    wavelength = channel_dict[band] # cm
    wavenumber = 1 / wavelength     # cm^-1

    radiance = ds['Rad_C' + band]

    # Calculate brightness temperature
    ds['btemp_C' + band] = c2 * wavenumber / np.log((c1 * wavenumber**3) / radiance + 1)

    return ds

def goes_rad_to_rgb(path, date, goes, location):
    """
    Downscale GOES bands C02, C05, and C13 to the same grid and interpolate to match the target dataset.
    The function also calculates reflectivity and brightness temperature for the specified bands.

    Parameters:
    path (str): The path to the directory containing the GOES-16 data files.
    date (str): The date string in the format 'YYYYMMDD'.
    goes (str): The GOES satellite identifier (e.g., 'goes16').
    location (str): The location identifier (e.g., 'colorado', 'washington').

    Returns:
    xarray.Dataset: A dataset containing the downscaled reflectivity and brightness temperature.
    """

    # Load the GOES-16 ABI data
    C02_file = f'C02/{goes}_C02_{location}_' + date + '.zarr'
    ds_C02 = xr.open_dataset(path+C02_file)
    C05_file = f'C05/{goes}_C05_{location}_' + date + '.zarr'
    ds_C05 = xr.open_dataset(path+C05_file)
    C13_file = f'C13/{goes}_C13_{location}_' + date + '.zarr'
    ds_C13 = xr.open_dataset(path+C13_file)

    # convert to lat and lon from x and y coordinates
    # removed - Stevens code already does this
    # lat_C02, lon_C02 = calculate_degrees(ds_C02)
    # lat_C05, lon_C05 = calculate_degrees(ds_C05)
    # lat_C13, lon_C13 = calculate_degrees(ds_C13)
    # ds_C02 = ds_C02.assign_coords(y=("y", lat_C02[:,0]), x=("x", lon_C02[0,:]))
    # ds_C05 = ds_C05.assign_coords(y=("y", lat_C05[:,0]), x=("x", lon_C05[0,:]))
    # ds_C13 = ds_C13.assign_coords(y=("y", lat_C13[:,0]), x=("x", lon_C13[0,:]))

    # Some days having missing files, fill times manually
    # Extract the date information from ds_C02
    time_C02 = ds_C02['t']
    full_date = pd.to_datetime(time_C02.values[0])  # Get the first timestamp as reference
    year, month, day = full_date.year, full_date.month, full_date.day

    # Create a complete time range for the day at 5-minute increments
    start_time = pd.Timestamp(year=year, month=month, day=day, hour=0, minute=2, second=30)
    end_time = pd.Timestamp(year=year, month=month, day=day, hour=23, minute=59, second=59)
    complete_time_range = pd.date_range(start=start_time, end=end_time, freq='5T')

    # Interpolate missing timesteps along the 't' dimension
    ds_C02 = ds_C02.sortby('t')
    ds_C05 = ds_C05.sortby('t')
    ds_C13 = ds_C13.sortby('t')
    
    # Reindex all datasets to the complete time range, assigning to the nearest timesteps
    ds_C02 = ds_C02.reindex(t=complete_time_range, method='nearest')
    ds_C13 = ds_C13.reindex(t=complete_time_range, method='nearest')
    ds_C05 = ds_C05.reindex(t=complete_time_range, method='nearest')

    # ds_C02 = ds_C02.interpolate_na(dim='t', method='linear')
    # ds_C05 = ds_C05.interpolate_na(dim='t', method='linear')
    # ds_C13 = ds_C13.interpolate_na(dim='t', method='linear')

    # Ensure the coordinate ranges overlap
    target_lat = ds_C02['latitude']
    target_lon = ds_C02['longitude']

    # Align the coordinates of ds_C13 with ds_C02
    ds_C13_aligned = ds_C13.reindex_like(ds_C02, method='nearest')
    ds_C05_aligned = ds_C05.reindex_like(ds_C02, method='nearest')

    # Interpolate the source dataset to the target dataset's grid using a specified method
    ds_C13_resampled = ds_C13_aligned.interp(latitude=target_lat, longitude=target_lon, method='linear')
    ds_C05_resampled = ds_C05_aligned.interp(latitude=target_lat, longitude=target_lon, method='linear')
    # If there are still NaNs, try filling them with a method like nearest neighbor
    ds_C13_resampled_filled = ds_C13_resampled.fillna(ds_C13_aligned.interp(latitude=target_lat, longitude=target_lon, method='nearest'))
    ds_C05_resampled_filled = ds_C05_resampled.fillna(ds_C05_aligned.interp(latitude=target_lat, longitude=target_lon, method='nearest'))

    # Create a new dataset with the two Rad variables
    combined_ds = xr.Dataset(
        {'tb_C13': ds_C13_resampled_filled['tb'], 'ref_C02': ds_C02['ref'], 'ref_C05': ds_C05_resampled_filled['ref']},
        coords={'latitude': target_lat, 'longitude': target_lon, 't': ds_C02['t']})
    
    # following 2 blocks unnesessary, steven's code already does this
    # reflectivity = kappa factor * radiance
    # combined_ds['refl_C02'] = ds_C02['kappa0'] * combined_ds['Rad_C02']
    # combined_ds['refl_C05'] = ds_C05['kappa0'] * combined_ds['Rad_C05']

    # Convert radiance to brightness temperature
    # radiance_to_brightness_temp(combined_ds, '13')

    combined_ds['green'] = goes_norm(combined_ds['ref_C02'], 
                                     0, 0.78, clip=True)
    combined_ds['blue'] = goes_norm(combined_ds['ref_C05'], 
                                    0.01, 0.59, clip=True)
    # VERY IMPORTANT 
    # The RGB Day Cloud Phase Distinction GOES product 
    # anchors the warm temps to 0 and cold temps to 1
    # so that the coldest temps display red and warmest temps display blue (counterintuitve)
    # therefore: we have to invert the nomalized data (bright_temp=True)
    combined_ds['red'] = goes_norm(combined_ds['tb_C13'], 
                                   219.65, 280.65, clip=True, invert_norm=True)

    combined_ds = combined_ds.drop_vars(['ref_C02', 'ref_C05', 'tb_C13'])
    
    out_name = f'{goes}_C02_C05_C13_rgb_{location}_{date}.nc'
    combined_ds.to_netcdf(path + f'rgb_composite/{out_name}', 
                          mode='w', format='NETCDF4')


#########################################################################################
# One Stop Shop for processing
#########################################################################################

def goes_zarr_to_rgb(in_dir, out_dir, date, goes, gif=False):
    # ***********BROKEN***********
    """
    Convert GOES ABI NetCDF files to a single Zarr file and then create an RGB composite.

    Parameters:
    in_dir (str): Directory containing the NetCDF files.
    out_dir (str): Output directory for the Zarr file.
    date (str): Date in 'YYYYMMDD' format.
    goes (str): GOES satellite identifier (e.g., 'goes16').
    gif (bool): If True, will return the filename for GIF creation instead of displaying the plot.

    Returns:
    xarray.Dataset: Dataset containing the RGB composite.
    """
    
    # Convert NetCDF files to Zarr
    goes_nc_to_zarr(in_dir, out_dir, 'goes16_' + date + '.zarr')

    # Load the Zarr file
    zarr_file = out_dir + 'goes16_' + date + '.zarr'
    ds = goes_rad_to_rgb(zarr_file, date)

    # Save final ds as netcdf for 1 day with the 3 processed and normalized RGB channels
    out_name = goes + '_RGB_' + date + '.nc'

    # Generate time string
    time_of_day = '12:00:00'  # Example time of day
    time_str = generate_time(date, time_of_day)

    # Plot RGB image
    rgb_plot = plot_rgb_image(ds, date, time_of_day, gif=gif)

    return rgb_plot

#########################################################################################
# Making GIF loops
#########################################################################################

# Function to generate the time string based on date and time of day
def generate_time(date, time_of_day):
    # Convert date string to datetime object
    date_obj = datetime.strptime(date, '%Y%m%d')
    # Format the date part
    date_str = date_obj.strftime('%Y-%m-%d')
    # Combine date and time of day
    time_str = date_str + 'T' + time_of_day
    return time_str

def plot_rgb_image(ds, date, time_of_day, gif=False):
    """""
    Plot RGB image from xarray dataset and save as PNG.
    Parameters:
    ds (xarray.Dataset): The dataset containing the RGB channels.
    date (str): The date in 'YYYYMMDD' format.
    time_of_day (str): The time of day in 'HH:MM:SS' format.
    gif (bool): Whether to save the plot as a GIF.
    """

    # Extract the values as NumPy arrays
    red = ds['red'].values
    green = ds['green'].values
    blue = ds['blue'].values

    # Find the minimum shape among the arrays
    min_shape = np.min([red.shape, green.shape, blue.shape], axis=0)

    # Resize the arrays to the minimum shape
    red_resized = red[:min_shape[0], :min_shape[1]]
    green_resized = green[:min_shape[0], :min_shape[1]]
    blue_resized = blue[:min_shape[0], :min_shape[1]]

    # Ensure the arrays have the same dimensions
    assert red_resized.shape == green_resized.shape == blue_resized.shape, "Arrays must have the same shape"

    # Stack the arrays along the last dimension to create an RGB image
    rgb_image = np.stack([red_resized, green_resized, blue_resized], axis=-1)

    # Extract longitude and latitude values
    lon = ds['x'].values
    lat = ds['y'].values

    # Plot the RGB image using matplotlib's imshow
    rgb_plot = plt.imshow(rgb_image, extent=[lon.min(), lon.max(), lat.min(), lat.max()])
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('GOES Day Cloud Phase RGB Composite - ' + date + ' ' + time_of_day + ' UTC')
    plt.axis('on')  # Show the axis

    # Save the plot as a PNG file
    filename = f'./plots/goes_RGB_{date}_{time_of_day}.png'
    plt.savefig(filename)

    
    if gif:
        plt.close()
        return filename
    else:
        return rgb_plot
        plt.show()


def create_gif_from_pngs(png_dir, output_gif):
    # Get a list of all PNG files in the directory
    png_files = sorted([f for f in os.listdir(png_dir) if f.endswith('.png')])

    # Read each image and append it to the images list
    images = []
    for png_file in png_files:
        image_path = os.path.join(png_dir, png_file)
        images.append(imageio.imread(image_path))

    # Create and save the GIF
    imageio.mimsave(output_gif, images, duration=0.1)  # Adjust the duration as needed

def cloud_mask(ds):
    """
    create a "cloud" mask where the green channel is greater than 0.5 and blue channel 
    is less than 0.2 set these pixels to red=0, blue=1, green=1

    can be modified for future processing and for different cloud types

    """

    green = ds['green'].values
    blue = ds['blue'].values
    red = ds['red'].values
    mask = (green > 0.4) & (blue > 0.4) & (red > 0.4)
    ds_masked = ds.copy()
    ds_masked['red'] = xr.where(mask, 0, ds['red'])
    ds_masked['green'] = xr.where(mask, 1, ds['green'])
    ds_masked['blue'] = xr.where(mask, 1, ds['blue'])

    return ds_masked

def make_gif(ds, date, start_time, end_time, mask=False):
    """
    Create a GIF from GOES data for a specific date and time range.

    Parameters:
    ds(xarray.Dataset): The dataset containing the GOES data.
    date (str): The date in 'YYYYMMDD' format.
    start_time (str): The start time in 'HHMM' format.
    end_time (str): The end time in 'HHMM' format.
    mask (bool): Whether to apply a cloud mask.
    """
    # input_file = f'/storage/cdalden/goes/goes16/RGB_composite/goes16_C02_C05_C13_RGB_colorado_{date}.nc'
    # ds = xr.open_dataset(input_file)

    if mask:
        print('Applying cloud mask...')
        ds = cloud_mask(ds)


    start_time = datetime.strptime(f"{date}T{start_time}00", '%Y%m%dT%H%M%S')
    end_time = datetime.strptime(f"{date}T{end_time}00", '%Y%m%dT%H%M%S')

    ds = ds.sortby('t')

    # List to store the filenames of the generated plots
    filenames = []

    # Loop through every 10-minute chunk
    current_time = start_time
    while current_time <= end_time:
        time_str = current_time.strftime('%Y-%m-%dT%H:%M:%S')
        ds_i = ds.sel(t=time_str, method='nearest')
        hour_of_day = current_time.strftime('%H:%M')
        filename = plot_rgb_image(ds_i, date, hour_of_day, gif=True)
        filenames.append(filename)
        current_time += timedelta(minutes=10)
        print(f'Generated RGB image for {time_str}')

    # Create the GIF
    start_time_out = start_time.strftime('%H%M')
    end_time_out =  end_time.strftime('%H%M')
    if mask:
        output_gif = f'./gifs/masked_goes_RGB_{date}_{start_time_out}_{end_time_out}.gif'
    else:
        output_gif = f'./gifs/goes_RGB_{date}_{start_time_out}_{end_time_out}.gif'
    with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Clean up the temporary files
    for filename in filenames:
        os.remove(filename)