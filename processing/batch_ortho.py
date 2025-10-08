import os
import sys
import orthorectify_modded

domain = input('What domain are you ortho-ing?    ')
# domain = 'scripps'
### VERY IMPORTANT ###
# CHANGE THE BOUNDS
def process_files(root_dir):
    # Loop through all subdirectories and files
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.nc'):  # Check if the file is a NetCDF file
                netcdf_path = os.path.join(subdir, file)
                print('working on ' + str(subdir))

                # Define all args needed for the ortho function
                goes_image_path = netcdf_path
                data_vars = ["Rad"]
                new_goes_filename = netcdf_path.replace('.nc', '_ortho.nc')
                if domain =='washington':
                    bounds = (-125, 45, -120, 49)
                elif domain == 'colorado':
                    bounds = (-109, 37, -104, 41)
                elif domain == 'scripps':
                    bounds = (-118, 32.5, -117, 33.5)
                api_key = "41d14aae7e761c0de3e8f99aa4fd24d9"

                if 'ortho' in netcdf_path:
                    print(f"File {goes_image_path} already ortho'd, skipping.")
                else:
                    orthorectify_modded.ortho(
                        goes_image_path,
                        data_vars,
                        bounds,
                        api_key,
                        new_goes_filename,
                        dem_filepath=None, # make this None to download domain first time, otherwise use 'temp_SRTMGL3_DEM.tif'
                        demtype="SRTMGL3",
                        keep_dem=True,
                    )
                    # Optionally delete the original NetCDF file
                    os.remove(netcdf_path)

if __name__ == "__main__":
    # Check if the root directory is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <root_directory>")
        sys.exit(1)

    # Get the root directory from the command-line argument
    root_dir = sys.argv[1]

    # Process the files
    process_files(root_dir)