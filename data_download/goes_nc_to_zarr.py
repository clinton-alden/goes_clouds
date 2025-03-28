# %% [markdown]
# ## Combine GOES netCDFs into one .zarr file
# 
# This could also be done using the `goes-ortho.build_zarr()` function but that is currently not working.
# 
# Clinton Alden
# 27 March 2025

# %%
import xarray as xr
import os

# %%
# Directory containing the NetCDF files
nc_dir = '/storage/cdalden/goes/goes16/2023/6/'

# Recursively list all NetCDF files in the directory and subdirectories
nc_files = []
for root, dirs, files in os.walk(nc_dir):
    for file in files:
        print(file)
        if file.endswith('.nc'):
            nc_files.append(os.path.join(root, file))

# Open multiple NetCDF files as a list of datasets
datasets = [xr.open_dataset(f) for f in nc_files]

# Concatenate datasets along the 't' coordinate
combined_ds = xr.concat(datasets, dim='t')

# Save the combined dataset to a Zarr file
combined_ds.to_zarr('/storage/cdalden/goes/goes16/2023/june_CO.zarr')

combined_ds

# %%



