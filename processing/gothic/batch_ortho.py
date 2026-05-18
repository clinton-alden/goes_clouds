import os
import sys
import orthorectify_modded
from multiprocessing import Pool, cpu_count

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# Gothic test version - 15km x 15km domain


def _resolve_data_vars():
    """Read comma-separated GOES_DATA_VARS from environment, fallback to Rad."""
    env_value = os.environ.get("GOES_DATA_VARS", "Rad")
    resolved = [v.strip() for v in env_value.split(",") if v.strip()]
    return resolved if resolved else ["Rad"]


def _resolve_worker_count(requested_workers, data_vars):
    """Choose a safe worker count for the current data variables.

    ACM orthorectification can be memory-heavy; cap workers unless explicitly forced.
    """
    visible_cpus = cpu_count()
    workers = requested_workers if requested_workers is not None else visible_cpus
    workers = max(1, min(workers, visible_cpus))

    is_acm = any(v.upper() == "ACM" for v in data_vars)
    force_workers = os.environ.get("FORCE_ORTHO_WORKERS", "0") == "1"
    if is_acm and not force_workers:
        acm_cap = int(os.environ.get("ORTHO_ACM_MAX_WORKERS", "4"))
        workers = max(1, min(workers, acm_cap))

    return workers


def _domain_bounds(domain):
    if domain == "gothic":
        return (-107.08, 38.89, -106.90, 39.03)
    if domain == "mammoth":
        return (-119.32, 37.41, -118.75, 37.86)
    if domain == "cues":
        return (-119.32, 37.41, -118.75, 37.86)
    if domain == "washington":
        return (-125, 45, -120, 49)
    if domain == "colorado":
        return (-109, 37, -104, 41)
    if domain == "scripps":
        return (-118, 32.5, -117, 33.5)
    raise ValueError(f"Unknown domain: {domain}")


def process_single_file(args):
    """Process a single NetCDF file - worker function for multiprocessing"""
    netcdf_path, domain, api_key, data_vars, dem_filepath = args
    
    if not os.path.exists(netcdf_path):
        return f"File disappeared before processing, skipping: {netcdf_path}"
    
    goes_image_path = netcdf_path
    new_goes_filename = netcdf_path.replace('.nc', '_ortho.nc')
    
    bounds = _domain_bounds(domain)

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
            dem_filepath=dem_filepath,
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

def process_files(root_dir, domain, data_vars=None, max_workers=None):
    print('DOMAIN: ' + str(domain))
    api_key = "41d14aae7e761c0de3e8f99aa4fd24d9"
    if data_vars is None:
        data_vars = ["Rad"]
    
    # Collect all .nc files first
    nc_files = []
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.nc') and not file.endswith('_ortho.nc'):
                netcdf_path = os.path.join(subdir, file)
                nc_files.append((netcdf_path, domain, api_key, data_vars))
    
    print(f"Found {len(nc_files)} files to process")
    
    if not nc_files:
        print("No files to process")
        return
    
    # Download DEM once with first file before parallel processing
    if nc_files:
        dem_filepath = os.path.join(root_dir, "temp_SRTMGL3_DEM.tif")
        if not os.path.exists(dem_filepath):
            print(f"Preparing DEM cache: {dem_filepath}")
            orthorectify_modded.get_dem(
                demtype="SRTMGL3",
                bounds=_domain_bounds(domain),
                api_key=api_key,
                out_fn=dem_filepath,
                proj="+proj=lonlat +ellps=GRS80",
            )

        first_file = nc_files[0]
        first_file = (first_file[0], domain, api_key, data_vars, dem_filepath)
        print(f"Processing first file: {first_file[0]}")
        result = process_single_file(first_file)
        print(result)
        nc_files = [
            (p[0], domain, api_key, data_vars, dem_filepath) for p in nc_files[1:]
        ]
    
    # Process remaining files in parallel if any left
    if nc_files:
        max_workers = _resolve_worker_count(max_workers, data_vars)
        print(f"Processing {len(nc_files)} remaining files with {max_workers} workers")
        with Pool(processes=max_workers, maxtasksperchild=8) as pool:
            for result in pool.imap_unordered(process_single_file, nc_files, chunksize=1):
                print(result)

        # Keep a clear stage marker for downstream log checks.
        print("Orthorectification stage complete")

if __name__ == "__main__":
    # Check if the root directory is provided as a command-line argument
    if len(sys.argv) != 3:
        print("Usage: python batch_ortho.py <root_directory> <domain>")
        sys.exit(1)

    # Get the root directory from the command-line argument
    root_dir = sys.argv[1]
    domain = sys.argv[2]
    data_vars = _resolve_data_vars()

    # Optional worker override (e.g., ORTHO_MAX_WORKERS=8)
    max_workers_env = os.environ.get("ORTHO_MAX_WORKERS")
    max_workers = int(max_workers_env) if max_workers_env else None

    # Process the files
    process_files(root_dir, domain, data_vars=data_vars, max_workers=max_workers)
