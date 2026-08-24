#!/usr/bin/env python3
"""Extract inclusive-footprint ACM fractions for three Colorado sites."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS, Transformer

SITES = {
    "gothic": (38.96, -106.99),
    "table_mountain": (40.12498, -105.23680),
    "senator_beck": (37.90, -107.72),
}

def crs_for(ds):
    a=ds.goes_imager_projection.attrs
    return CRS.from_proj4(f"+proj=geos +h={a['perspective_point_height']} +lon_0={a['longitude_of_projection_origin']} +sweep={a['sweep_angle_axis']} +a={a['semi_major_axis']} +b={a['semi_minor_axis']} +units=m +no_defs")

def footprint_selection(ds, lat0, lon0):
    lat_min=lat0-5/111.32; lat_max=lat0
    half=5/(111.32*np.cos(np.deg2rad(lat0))); lon_min=lon0-half; lon_max=lon0+half
    crs=crs_for(ds); h=float(ds.goes_imager_projection.attrs['perspective_point_height'])
    fwd=Transformer.from_crs(4326,crs,always_xy=True); inv=Transformer.from_crs(crs,4326,always_xy=True)
    x=np.asarray(ds.x,float); y=np.asarray(ds.y,float); dx=abs(float(np.median(np.diff(x)))); dy=abs(float(np.median(np.diff(y))))
    cx,cy=fwd.transform([lon_min,lon_min,lon_max,lon_max],[lat_min,lat_max,lat_min,lat_max])
    ix=np.where((x>=min(cx)/h-2*dx)&(x<=max(cx)/h+2*dx))[0]; iy=np.where((y>=min(cy)/h-2*dy)&(y<=max(cy)/h+2*dy))[0]
    xx,yy=np.meshgrid(x[ix],y[iy]); los=[]; las=[]
    for sx,sy in ((-.5,-.5),(-.5,.5),(.5,-.5),(.5,.5)):
        lo,la=inv.transform((xx+sx*dx)*h,(yy+sy*dy)*h); los.append(lo); las.append(la)
    mask=(np.max(las,axis=0)>=lat_min)&(np.min(las,axis=0)<=lat_max)&(np.max(los,axis=0)>=lon_min)&(np.min(los,axis=0)<=lon_max)
    return iy,ix,mask

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input-dir',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--date-key',required=True); a=p.parse_args()
    files=sorted(a.input_dir.rglob('OR_ABI-L2-ACMC-*.nc')); rows={s:[] for s in SITES}; cache={}
    for path in files:
        try:
            with xr.open_dataset(path,mask_and_scale=True) as ds:
                key=(len(ds.x),len(ds.y),float(ds.x[0]),float(ds.y[0]),float(ds.goes_imager_projection.attrs['longitude_of_projection_origin']))
                if key not in cache: cache[key]={s:footprint_selection(ds,*coord) for s,coord in SITES.items()}
                for site,(iy,ix,mask) in cache[key].items():
                    bcm=np.asarray(ds.BCM.isel(y=iy,x=ix),float); valid=np.isfinite(bcm)&mask
                    rows[site].append({'time':pd.Timestamp(ds.t.values),'cloud_pixels':int(((bcm==1)&valid).sum()),'valid_pixels':int(valid.sum()),'grid_pixels':int(mask.sum())})
        except (OSError,ValueError,KeyError) as e: print('skip',path,e,flush=True)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    for site,data in rows.items():
        if not data: continue
        d=pd.DataFrame(data).sort_values('time').drop_duplicates('time'); d['hour']=d.time.dt.floor('h')
        h=d.groupby('hour',as_index=False).agg(acm_cloud_pixels=('cloud_pixels','sum'),acm_valid_pixels=('valid_pixels','sum'),acm_n_scans=('time','size'),acm_grid_pixels=('grid_pixels','median'))
        h['acm_cloud_fraction']=h.acm_cloud_pixels/h.acm_valid_pixels
        h.to_csv(a.output_dir/f'{site}_acm_hourly_{a.date_key}.csv',index=False)
    print(f'{a.date_key}: {len(files)} scans -> {a.output_dir}',flush=True)
if __name__=='__main__': main()
