import os
import sys
import orthorectify_modded

# domain = input('What domain are you ortho-ing?    ')
# domain = 'scripps'

def require_opentopography_api_key():
    api_key = os.environ.get("OPENTOPOGRAPHY_API_KEY")
    if not api_key:
        raise RuntimeError(
            "OPENTOPOGRAPHY_API_KEY is required for DEM downloads. "
            "Create your own OpenTopography API key and export it before running ortho."
        )
    return api_key

### VERY IMPORTANT ###
# CHANGE THE BOUNDS
def process_files(root_dir, domain):
    print('DOMAIN: ' + str(domain))
    api_key = require_opentopography_api_key()
    # Loop through all subdirectories and files
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.nc'):  # Check if the file is a NetCDF file
                netcdf_path = os.path.join(subdir, file)
                if not os.path.exists(netcdf_path):
                    print(f"File disappeared before processing, skipping: {netcdf_path}")
                    continue
                print('working on ' + str(subdir))

                # Define all args needed for the ortho function
                goes_image_path = netcdf_path
                data_vars = ["Rad"]
                # data_vars = ['ACM']
                new_goes_filename = netcdf_path.replace('.nc', '_ortho.nc')
                if domain =='washington':
                    bounds = (-125, 45, -120, 49)
                elif domain == 'colorado':
                    bounds = (-109, 37, -104, 41)
                elif domain == 'colorado_west23':
                    bounds = (-109, 37, -104, 41)
                elif domain == 'sierra':
                    bounds = (-121, 35, -118, 40)
                elif domain == 'scripps':
                    bounds = (-118, 32.5, -117, 33.5)

                if goes_image_path.endswith('_ortho.nc'):
                    print(f"File {goes_image_path} already ortho'd, skipping.")
                    continue

                if os.path.exists(new_goes_filename):
                    print(f"Output already exists, skipping: {new_goes_filename}")
                    continue

                try:
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
                except FileNotFoundError:
                    print(f"Missing during ortho open, skipping: {goes_image_path}")
                    continue
                except PermissionError:
                    print(f"Permission denied during ortho write, skipping: {new_goes_filename}")
                    continue

                # Delete the original file if it still exists. Some ortho paths
                # may already remove or move it.
                try:
                    os.remove(netcdf_path)
                except FileNotFoundError:
                    print(f"Source file already removed, skipping delete: {netcdf_path}")

if __name__ == "__main__":
    # Check if the root directory is provided as a command-line argument
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <root_directory> <domain>")
        sys.exit(1)

    # Get the root directory from the command-line argument
    root_dir = sys.argv[1]
    domain = sys.argv[2]

    # Process the files
    process_files(root_dir, domain)
