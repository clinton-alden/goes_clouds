# Quick Reference Guide - GOES Processing Code

**Location**: `/glade/u/home/cdalden/goes_work/`

---

## 📑 Start Here

1. **[GITHUB_FILES_INDEX.md](GITHUB_FILES_INDEX.md)** ← Complete index of all files
2. **[GITHUB_REPOSITORY_FILES.md](GITHUB_REPOSITORY_FILES.md)** ← Detailed documentation

---

## 🔧 Main Functions (in utils.py)

### Convert NetCDF to Zarr
```python
import sys
sys.path.append('/glade/u/home/cdalden/goes_work')
from GITHUB_SOURCES_utils import goes_nc_to_zarr

goes_nc_to_zarr(
    in_dir='/path/to/data/',
    channels=['C02', 'C05', 'C13'],
    startday=1,
    endday=31,
    month=6,
    year=2022,
    location='colorado',
    goes_model='goes16',
    surprise=True
)
```

### Create RGB Composite
```python
from GITHUB_SOURCES_utils import goes_rad_to_rgb

goes_rad_to_rgb(
    path='/path/to/zarr/files/',
    date='20220601',
    goes='goes16',
    location='colorado'
)
# Outputs: goes16_C02_C05_C13_rgb_colorado_20220601.nc
```

### Detect Clouds
```python
from GITHUB_SOURCES_utils import cloud_mask
import xarray as xr

ds = xr.open_dataset('goes16_C02_C05_C13_rgb_colorado_20220601.nc')
ds_with_clouds = cloud_mask(ds)
# Pixels with green > 0.4 AND blue > 0.4 AND red > 0.4 are marked as clouds
```

### Create Animated GIF
```python
from GITHUB_SOURCES_utils import make_gif

make_gif(
    ds=ds,
    date='20220601',
    start_time='0600',
    end_time='1800',
    mask=False
)
# Outputs:./gifs/goes_RGB_20220601_0600_1800.gif
```

---

## 🚀 Quick Scripts

### Process a full month (Zarr)
```bash
cd /glade/u/home/cdalden/goes_work
python GITHUB_SOURCES_zarr_v2.py /storage/data/ 2022 6 goes16 colorado
```

### Process a full month (RGB)
```bash
python GITHUB_SOURCES_rgb_v2.py /storage/data/ 2022 6 colorado goes16
```

### Calculate cloud frequency
```bash
python GITHUB_SOURCES_daily_cloud_frequency.py /storage/data/ 2022 06 colorado goes16
```

---

## 📊 Data Workflow

```
Step 1: Download GOES NetCDF
  └─> 288 × 3 channels per day = 864 files/day
  └─> ~200 MB per file = ~170 GB/day

Step 2: Convert to Zarr (faster access)
  └─> Concatenate 288 5-min timesteps
  └─> Save as 3 Zarr files/day (~2 GB per set)
  └─> Delete NetCDF to save space

Step 3: Generate RGB Composite
  └─> Load C02 (visible red), C05 (NIR), C13 (thermal IR)
  └─> Normalize to 0-1 range
  └─> Combine into single RGB NetCDF (~200 MB)
  └─> Delete Zarr to save space

Step 4: Analyze Cloud Frequency
  └─> Detect clouds: blue >= 0.13 AND green >= 0.15
  └─> Calculate hourly frequency
  └─> Save cloud frequency dataset
```

---

## 🎨 RGB Color Mapping

```
Red Channel (C13):     Thermal IR 10.3 μm - Brightness Temperature
  Range: 219.65 - 280.65 K (inverted)
  Cold clouds (high altitude) → Bright Red
  Warm surface (low altitude) → Dark Blue

Green Channel (C02):   Visible Red 0.64 μm - Reflectivity
  Range: 0 - 0.78
  Bright clouds/snow → Bright Green
  Dark surface/ocean → Dark

Blue Channel (C05):    Near-IR 1.38 μm - Reflectivity
  Range: 0.01 - 0.59
  Cirrus clouds → Bright Blue
  Water/surface → Dark
```

---

## 📦 Dependencies

```python
import xarray as xr       # NetCDF/Zarr I/O
import numpy as np        # Array operations
import pandas as pd       # Time series
import matplotlib.pyplot  # Plotting
import imageio            # GIF creation
import os, gc, shutil     # System operations
```

---

## ⚙️ Key Configuration Parameters

**Channels used for RGB:**
- C02: 0.64 μm (visible red)
- C05: 1.38 μm (near-infrared)
- C13: 10.3 μm (thermal IR)

**Time resolution:** 5 minutes (288 timesteps/day)

**Domains available:**
- colorado, washington, scripps, (others)

**GOES satellites:**
- goes16 (GOES-East)
- goes17 (GOES-West)
- goes18

**Cloud detection threshold:**
- Cloud: blue >= 0.13 AND green >= 0.15

---

## 📍 File Organization

```
/glade/u/home/cdalden/goes_work/
│
├─ GITHUB_FILES_INDEX.md                 ← You are here
├─ GITHUB_REPOSITORY_FILES.md            ← Full documentation
│
├─ GITHUB_SOURCES_utils.py               ← Main utilities (import these)
├─ GITHUB_SOURCES_rgb_v2.py              ← Run as script
├─ GITHUB_SOURCES_zarr_v2.py             ← Run as script
├─ GITHUB_SOURCES_daily_cloud_frequency.py ← Run as script
├─ GITHUB_SOURCES_batch_zarr.py          ← Legacy (hardcoded)
└─ GITHUB_SOURCES_batch_rgb.py           ← Legacy (hardcoded)
```

---

## 🔍 Example: Full Processing Pipeline

```python
#!/usr/bin/env python3

# Set up paths
import sys
sys.path.append('/glade/u/home/cdalden/goes_work')

from GITHUB_SOURCES_utils import (
    goes_nc_to_zarr, goes_rad_to_rgb, 
    cloud_mask, make_gif
)
import xarray as xr

# Configuration
BASE_DIR = '/storage/cdalden/goes/'
YEAR, MONTH, DAY = 2022, 6, 15
DOMAIN, GOES = 'colorado', 'goes16'
DATE = f'{YEAR}{MONTH:02d}{DAY:02d}'

# Step 1: Convert NetCDF to Zarr
print("Converting NetCDF to Zarr...")
goes_nc_to_zarr(
    BASE_DIR, ['C02', 'C05', 'C13'],
    DAY, DAY, MONTH, YEAR, DOMAIN, GOES
)

# Step 2: Create RGB composite
print("Creating RGB composite...")
goes_rad_to_rgb(BASE_DIR, DATE, GOES, DOMAIN)

# Step 3: Analyze clouds and create GIF
print("Analyzing clouds and creating GIF...")
rgb_file = f'{BASE_DIR}{GOES}/rgb_composite/{GOES}_C02_C05_C13_rgb_{DOMAIN}_{DATE}.nc'
ds = xr.open_dataset(rgb_file)

ds_clouds = cloud_mask(ds)
make_gif(ds_clouds, DATE, '0600', '1800', mask=True)

print("✅ Complete!")
```

---

## 🐛 Troubleshooting

### "Directory does not exist"
- Check input path format (should end with `/`)
- Verify GOES subdirectory structure exists

### "Expected X files, found Y"
- Some days may have missing data
- Script continues with warning (not fatal)

### Out of memory
- Zarr conversion processes one day at a time
- RGB generation processes one timestep at a time
- Usually completes in minutes with 100GB+ files

### Missing output files
- Check directory exists: `mkdir -p ./rgb_composite`
- Verify input Zarr files were created
- Check network/storage connectivity

---

## 📞 Reference

**GitHub**: https://github.com/clinton-alden/goes_clouds  
**Author**: Clinton Alden  
**Last Updated**: February 24, 2025

**Relevant Papers/Docs:**
- GOES-R Product User Guide (PUG) Volume 5 (L2 products)
- GOES ABI Fixed Grid Projection specifications

---

## 🎯 Next Steps

1. ✅ Review `GITHUB_REPOSITORY_FILES.md` for full API documentation
2. ✅ Try the example pipeline above
3. ✅ Adapt scripts for your specific domain/date range
4. ✅ Check GitHub repository for latest updates

Good luck with your GOES analysis! 🛰️
