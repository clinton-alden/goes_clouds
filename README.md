# Using GOES to study orographic clouds and valley scale winds in the East River of Colorado

### Steps to install goes_ortho
1. create new env
2. I could not install goes_ortho from Steven Pestana's github at first. The following 3 installs allowed me to successfully install goes_ortho via pip.
`conda install -c conda-forge gdal`
`conda install -c conda-forge gcc_linux-64`
`conda install -c conda-forge gxx_linux-64`
3. `pip install goes_ortho`


### Current to-do list
1. [X] Need to convert from x and y to lat lon
2. [X] Resample to match resolutions
3. [ ] Correct radiance to reflectivity for bands 2 and 5
4. [ ] Correct radiance to brightness temperature for band 13
5. [ ] Adjust RGB normalization based on the NOAA guidance from Matt Jochum
6. [ ] Ortho rectify images using goes-ortho
7. [ ] Combine 3 RGB bands into one zarr file 