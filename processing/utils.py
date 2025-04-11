# Infrastructure for processing GOES ABI data and creating the Day Cloud Phase RGB Composites
import xarray as xr
import numpy as np
import os

def goes_nc_to_zarr(in_dir, out_dir, out_name): 
    """
    Convert multiple NetCDF files to a single Zarr file.
    
    Parameters:
    in_dir (str): Directory containing the NetCDF files.
    out_dir (str): Output directory for the Zarr file.
    out_name (str): Name of the output Zarr file.
    
    Returns:
    None
    """
    
 
    # Directory containing the NetCDF files
    
    # Recursively list all NetCDF files in the directory and subdirectories
    nc_files = []
    for root, dirs, files in os.walk(in_dir):
        for file in files:
            # print(file)
            if file.endswith('.nc'):
                nc_files.append(os.path.join(root, file))

    # Open multiple NetCDF files as a list of datasets
    datasets = [xr.open_dataset(f) for f in nc_files]

    # Concatenate datasets along the 't' coordinate
    combined_ds = xr.concat(datasets, dim='t')

    # Save the combined dataset to a Zarr file
    out_name = out_name
    combined_ds.to_zarr(out_dir+out_name)

    return print("Zarr file saved to " + out_dir + out_name)


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
def goes_norm(value, lower_limit, upper_limit, clip=True):
    """
    Normalize values between an upper and lower limit between 0 and 1.

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

def goes_rad_to_rgb(path, date):
    """
    Downscale GOES-16 bands C02, C05, and C13 to the same grid and interpolate to match the target dataset.
    The function also calculates reflectivity and brightness temperature for the specified bands.

    Parameters:
    path (str): The path to the directory containing the GOES-16 data files.
    date (str): The date string in the format 'YYYYMMDD'.

    Returns:
    xarray.Dataset: A dataset containing the downscaled reflectivity and brightness temperature.
    """

    # Load the GOES-16 ABI data
    file = 'goes16_C02_colorado_' + date + '.zarr'
    ds_C02 = xr.open_dataset(path+'channel02/'+file)
    file = 'goes16_C05_colorado_' + date + '.zarr'
    ds_C05 = xr.open_dataset(path+'channel05/'+file)
    file = 'goes16_C13_colorado_' + date + '.zarr'
    ds_C13 = xr.open_dataset(path+'channel13/'+file)

    # convert to lat and lon from x and y coordinates
    lat_C02, lon_C02 = calculate_degrees(ds_C02)
    lat_C05, lon_C05 = calculate_degrees(ds_C05)
    lat_C13, lon_C13 = calculate_degrees(ds_C13)
    ds_C02 = ds_C02.assign_coords(y=("y", lat_C02[:,0]), x=("x", lon_C02[0,:]))
    ds_C05 = ds_C05.assign_coords(y=("y", lat_C05[:,0]), x=("x", lon_C05[0,:]))
    ds_C13 = ds_C13.assign_coords(y=("y", lat_C13[:,0]), x=("x", lon_C13[0,:]))

    # Extract time coordinates from ds_C02
    time_C02 = ds_C02['t']
    ds_C13 = ds_C13.assign_coords(t=time_C02)
    ds_C05 = ds_C05.assign_coords(t=time_C02)

    # Ensure the coordinate ranges overlap
    target_lat = ds_C02['y']
    target_lon = ds_C02['x']

    # Align the coordinates of ds_C13 with ds_C02
    ds_C13_aligned = ds_C13.reindex_like(ds_C02, method='nearest')
    ds_C05_aligned = ds_C05.reindex_like(ds_C02, method='nearest')

    # Interpolate the source dataset to the target dataset's grid using a specified method
    ds_C13_resampled = ds_C13_aligned.interp(y=target_lat, x=target_lon, method='linear')
    ds_C05_resampled = ds_C05_aligned.interp(y=target_lat, x=target_lon, method='linear')
    # If there are still NaNs, try filling them with a method like nearest neighbor
    ds_C13_resampled_filled = ds_C13_resampled.fillna(ds_C13_aligned.interp(y=target_lat, x=target_lon, method='nearest'))
    ds_C05_resampled_filled = ds_C05_resampled.fillna(ds_C05_aligned.interp(y=target_lat, x=target_lon, method='nearest'))

    # Create a new dataset with the two Rad variables
    combined_ds = xr.Dataset(
        {'Rad_C13': ds_C13_resampled_filled['Rad'], 'Rad_C02': ds_C02['Rad'], 'Rad_C05': ds_C05_resampled_filled['Rad']},
        coords={'y': target_lat, 'x': target_lon, 't': ds_C02['t']})
    
    # reflectivity = kappa factor * radiance
    combined_ds['refl_C02'] = ds_C02['kappa0'] * combined_ds['Rad_C02']
    combined_ds['refl_C05'] = ds_C05['kappa0'] * combined_ds['Rad_C05']

    # Convert radiance to brightness temperature
    radiance_to_brightness_temp(combined_ds, '13')

    combined_ds['green'] = goes_norm(combined_ds['refl_C02'], 0, 0.78, clip=True)
    combined_ds['blue'] = goes_norm(combined_ds['refl_C05'], 0.01, 0.59, clip=True)
    combined_ds['red'] = goes_norm(combined_ds['btemp_C13'], 219.65, 280.65, clip=True)

    combined_ds = combined_ds.drop_vars(['Rad_C02', 'Rad_C05', 'Rad_C13'])

    return combined_ds


#########################################################################################
# Making GIF loops
#########################################################################################

