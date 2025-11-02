# This just runs the code developed in ./analysis/01_plot_RGB.ipynb in a .py file
import xarray as xr
from PIL import Image
from analysis_utils import make_gif

# suppress warnings
import warnings
warnings.filterwarnings("ignore")

for i in range(1,32): 
    day = str(i).zfill(2)
    date = '202305' + day
    start_time = '1400'
    end_time = '2355'
    state = 'colorado'

    path = '/storage/cdalden/goes/colorado/goes16/rgb_composite/'
    file = f'goes16_C02_C05_C13_rgb_{state}_{date}.nc'
    ds = xr.open_dataset(path+file)

    make_gif(ds, date, start_time, end_time, subset=state, mask=False)
    ds.close()