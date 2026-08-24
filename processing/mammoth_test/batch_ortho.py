import os
import sys
import orthorectify_modded
from multiprocessing import Pool, cpu_count

# Mammoth Mountain test version - 50km x 50km domain

def require_opentopography_api_key():
    api_key = os.environ.get("OPENTOPOGRAPHY_API_KEY")
    if not api_key:
        raise RuntimeError(
            "OPENTOPOGRAPHY_API_KEY is required for DEM downloads. "
            "Create your own OpenTopography API key and export it before running ortho."
        )
    return api_key


def process_single_file(args):
    """Process a single NetCDF file - worker function for multiprocessing"""
    netcdf_path, domain, api_key = args
    
    if not os.path.exists(netcdf_path):
        return f"File disappeared before processing, skipping: {netcdf_path}"
    
    goes_image_path = netcdf_path
    data_vars = ["Rad"]
    new_goes_filename = netcdf_path.replace('.nc', '_ortho.nc')
    
    if domain == 'mammoth':
        bounds = (-119.32, 37.41, -118.75, 37.86)
    elif domain =='washington':
        bounds = (-125, 45, -120, 49)
    elif domain == 'colorado':
        bounds = (-109, 37, -104, 41)
    elif domain == 'scripps':
        bounds = (-118, 32.5, -117, 33.5)
    else:
        raise ValueError(f"Unknown domain: {domain}")

    if goes_image_path.endswith('_ortho.nc'):
        return f"File {goes_image_path} already ortho'd, skipping."

    if os.path.exists(new_goes_filename):
        return f"Output already exists, skipping: {new_goes_filename}"

    try:
        orthorectify_modded.ortho(
            goes_image_path,
            data_vars,
            bounds,
            api_key,
            new_goes_filename,
            dem_filepath='temp_SRTMGL3_DEM.tif',
            demtype="SRTMGL3",
            keep_dem=True,
        )
        try:
            os.remove(netcdf_path)
        except FileNotFoundError:
            pass
        return f"SUCCESS: {netcdf_path}"
    except FileNotFoundError:
        return f"Missing during ortho open, skipping: {goes_image_path}"
    except PermissionError:
        return f"Permission denied during ortho write, skipping: {new_goes_filename}"
    except Exception as e:
        return f"ERROR processing {netcdf_path}: {e}"

def process_files(root_dir, domain, max_workers=None):
    print('DOMAIN: ' + str(domain))
    api_key = require_opentopography_api_key()
    
    # Collect all .nc files first
    nc_files = []
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.nc') and not file.endswith('_ortho.nc'):
                netcdf_path = os.path.join(subdir, file)
                nc_files.append((netcdf_path, domain, api_key))
    
    print(f"Found {len(nc_files)} files to process")
    
    if not nc_files:
        print("No files to process")
        return
    
    # Download DEM once with first file before parallel processing
    if nc_files:
        first_file = nc_files[0]
        print(f"Downloading DEM with first file: {first_file[0]}")
        result = process_single_file(first_file)
        print(result)
        nc_files = nc_files[1:]
    
    # Process remaining files in parallel if any left
    if nc_files:
        if max_workers is None:
            # Default to all visible CPUs unless caller constrains via env.
            max_workers = cpu_count()
        
        print(f"Processing {len(nc_files)} remaining files with {max_workers} workers")
        with Pool(processes=max_workers) as pool:
            results = pool.map(process_single_file, nc_files)
        
        for result in results:
            print(result)

if __name__ == "__main__":
    # Check if the root directory is provided as a command-line argument
    if len(sys.argv) != 3:
        print("Usage: python batch_ortho.py <root_directory> <domain>")
        sys.exit(1)

    # Get the root directory from the command-line argument
    root_dir = sys.argv[1]
    domain = sys.argv[2]

    # Optional worker override (e.g., ORTHO_MAX_WORKERS=8)
    max_workers_env = os.environ.get("ORTHO_MAX_WORKERS")
    max_workers = int(max_workers_env) if max_workers_env else None

    # Process the files
    process_files(root_dir, domain, max_workers=max_workers)
