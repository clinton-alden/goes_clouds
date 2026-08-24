import os
import subprocess
import sys
import orthorectify_v3 as orthorectify_modded

# domain = input('What domain are you ortho-ing?    ')
domain = 'colorado'

### VERY IMPORTANT ###
# CHANGE THE BOUNDS


def require_opentopography_api_key():
    api_key = os.environ.get("OPENTOPOGRAPHY_API_KEY")
    if not api_key:
        raise RuntimeError(
            "OPENTOPOGRAPHY_API_KEY is required for DEM downloads. "
            "Create your own OpenTopography API key and export it before running ortho."
        )
    return api_key


def get_bounds(domain):
    if domain == 'washington':
        return (-125, 45, -120, 49)
    if domain == 'colorado':
        return (-109, 37, -104, 41)
    if domain == 'scripps':
        return (-118, 32.5, -117, 33.5)
    raise ValueError(f"Unknown domain '{domain}'. Expected washington, colorado, or scripps.")


def process_one_file(netcdf_path, domain):
    if not os.path.exists(netcdf_path):
        print(f"File disappeared before processing, skipping: {netcdf_path}")
        return 0

    print('working on ' + str(os.path.dirname(netcdf_path)))

    goes_image_path = netcdf_path
    data_vars_env = os.environ.get("GOES_DATA_VARS", "").strip()
    if data_vars_env:
        data_vars = [part.strip() for part in data_vars_env.split(",") if part.strip()]
    else:
        data_vars = ['ACM']
    new_goes_filename = netcdf_path.replace('.nc', '_ortho.nc')
    bounds = get_bounds(domain)
    api_key = require_opentopography_api_key()

    if goes_image_path.endswith('_ortho.nc'):
        print(f"File {goes_image_path} already ortho'd, skipping.")
        return 0

    if os.path.exists(new_goes_filename):
        print(f"Output already exists, skipping: {new_goes_filename}")
        return 0

    try:
        orthorectify_modded.ortho(
            goes_image_path,
            data_vars,
            bounds,
            api_key,
            new_goes_filename,
            dem_filepath=None,  # None downloads DEM to scratch path in v3.
            demtype="SRTMGL3",
            keep_dem=True,
        )
    except FileNotFoundError:
        print(f"Missing during ortho open, skipping: {goes_image_path}")
        return 0
    except PermissionError:
        print(f"Permission denied during ortho write, skipping: {new_goes_filename}")
        return 0
    except RuntimeError as e:
        print(f"Runtime error during ortho, skipping: {new_goes_filename} ({e})")
        return 0

    try:
        os.remove(netcdf_path)
    except FileNotFoundError:
        print(f"Source file already removed, skipping delete: {netcdf_path}")

    return 0


def process_files(root_dir, domain):
    print('DOMAIN: ' + str(domain))
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if not file.endswith('.nc') or file.endswith('_ortho.nc'):
                continue

            netcdf_path = os.path.join(subdir, file)
            cmd = [
                sys.executable,
                os.path.abspath(__file__),
                '--single-file',
                netcdf_path,
                domain,
            ]
            result = subprocess.run(cmd)
            if result.returncode != 0:
                print(
                    f"Child process failed for {netcdf_path} with return code {result.returncode}; continuing."
                )

if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == '--single-file':
        sys.exit(process_one_file(sys.argv[2], sys.argv[3]))

    if len(sys.argv) != 3:
        print("Usage: python script_name.py <root_directory> <domain>")
        print("   or: python script_name.py --single-file <netcdf_path> <domain>")
        sys.exit(1)

    root_dir = sys.argv[1]
    domain = sys.argv[2]
    process_files(root_dir, domain)
