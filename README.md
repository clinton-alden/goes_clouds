# Using GOES to study orographic clouds and valley scale winds in the East River of Colorado

### Steps to install goes_ortho
1. create new env
2. I could not install goes_ortho from Steven Pestana's github at first. The following 3 installs allowed me to successfully install goes_ortho via pip.
`conda install -c conda-forge gdal`
`conda install -c conda-forge gcc_linux-64`
`conda install -c conda-forge gxx_linux-64`
3. `pip install goes_ortho`
4. install goespy by cloning the github repo (pip installing uses a weird/broken version that does not download files correctly) - https://github.com/palexandremello/goes-py/tree/master


### Current processing workflow
1. Create environment using above steps
2. Download specified channels and dates for GOES data using `./data_download/download.sh`
    - usage: `nohup ./download.sh`
3. Orthorectify raw GOES .nc files using `./processing/ortho_batch.py`
    - usage: `python ./ortho_batch.py /path/to/GOES/files/`
4. Combine .nc files to daily .zarr for each day by channel using `utils.goes_nc_to_zarr`
    - usage: see `./processing/00_goes_rad_corrections.ipynb`
    - ***WANT TO MAKE .sh FILE*** 
5. Create RGB composite file for each day using `utils.goes_rad_to_rgb`
    - usage: see `./processing/00_goes_rad_corrections.ipynb`
    - ***WANT TO MAKE .sh FILE*** 
