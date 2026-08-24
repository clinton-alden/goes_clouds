#!/usr/bin/env python3
"""Apply Gothic rules pixelwise and aggregate hourly for a Colorado site/month."""
from __future__ import annotations
import argparse,re
from pathlib import Path
import numpy as np,pandas as pd,xarray as xr
from pyproj import CRS,Transformer

ROOT=Path(__file__).resolve().parent
RULES=ROOT/'output_11d_gothic/gothic_rgb_tempbin_decision_tree_rules.csv'
COND=re.compile(r'(red|green|blue)\s*(<=|>)\s*([-+0-9.eE]+)')
CONFIG={
 'table_mountain':(40.12498,-105.23680,Path('/glade/derecho/scratch/cdalden/colorado/goes16/rgb_composite'),ROOT/'output_13b_boulder_domain_cloud_binary/boulder_rgb_binary_all_available.csv','time'),
 'senator_beck':(37.90,-107.72,Path('/glade/derecho/scratch/cdalden/senator_beck/goes16/rgb_composite'),ROOT/'output_13c_senator_beck_domain_cloud_binary/senator_beck_rgb_binary_all_available.csv','time'),
}
def rule_mask(v,rule):
 out=np.ones_like(v['red'],bool)
 for f,op,t in COND.findall(rule): out &= v[f]<=float(t) if op=='<=' else v[f]>float(t)
 return out
def native_mask(ds,lat0,lon0):
 h=35786023.; crs=CRS.from_proj4('+proj=geos +h=35786023 +lon_0=-75 +sweep=x +a=6378137 +b=6356752.31414 +units=m +no_defs')
 xx,yy=np.meshgrid(np.asarray(ds.x,float)*h,np.asarray(ds.y,float)*h); lon,lat=Transformer.from_crs(crs,4326,always_xy=True).transform(xx,yy)
 la=lat0-5/111.32; half=5/(111.32*np.cos(np.deg2rad(lat0)))
 return (lat>=la)&(lat<=lat0)&(lon>=lon0-half)&(lon<=lon0+half)
def latlon_mask(ds,lat0,lon0):
 lat=np.asarray(ds.latitude,float);lon=np.asarray(ds.longitude,float); dlat=abs(np.median(np.diff(lat)));dlon=abs(np.median(np.diff(lon)))
 la=lat0-5/111.32;half=5/(111.32*np.cos(np.deg2rad(lat0)))
 iy=np.where((lat+dlat/2>=la)&(lat-dlat/2<=lat0))[0];ix=np.where((lon+dlon/2>=lon0-half)&(lon-dlon/2<=lon0+half))[0]
 return iy,ix
def main():
 p=argparse.ArgumentParser();p.add_argument('--site',choices=CONFIG,required=True);p.add_argument('--year',type=int,required=True);p.add_argument('--month',type=int,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
 lat0,lon0,rgbdir,tempfile,_=CONFIG[a.site]; files=sorted(rgbdir.glob(f'*_{a.year}{a.month:02d}??.nc')); temps=pd.read_csv(tempfile,usecols=['time','era5_temp_c'],parse_dates=['time']).dropna().drop_duplicates('time').sort_values('time');temps['time']=pd.to_datetime(temps.time).astype('datetime64[ns]')
 rules=pd.read_csv(RULES);rules=rules[rules.status.eq('trained')]; cloudy={str(k):g[g.prediction.eq(1)].rule.astype(str).tolist() for k,g in rules.groupby('temp_bin')}; rows=[]
 for path in files:
  with xr.open_dataset(path) as ds:
   internal=pd.DatetimeIndex(pd.to_datetime(ds.t.values));date=pd.to_datetime(path.stem.rsplit('_',1)[-1]);times=pd.DatetimeIndex(date+(internal-internal.normalize())).astype('datetime64[ns]')
   td=pd.DataFrame({'time':times});td=pd.merge_asof(td,temps,on='time',direction='nearest',tolerance=pd.Timedelta('61min'))
   bins=np.full(len(times),None,object)
   for label,g in rules.groupby('temp_bin'): bins[(td.era5_temp_c>=float(g.temp_left_c.iloc[0]))&(td.era5_temp_c<float(g.temp_right_c.iloc[0]))]=str(label)
   if 'latitude' in ds.coords:
    iy,ix=latlon_mask(ds,lat0,lon0); arrays={b:np.asarray(ds[b].isel(latitude=iy,longitude=ix),np.float32).reshape(len(times),-1) for b in ('red','green','blue')}
   else:
    mask=native_mask(ds,lat0,lon0); arrays={b:np.asarray(ds[b],np.float32)[:,mask].reshape(len(times),-1) for b in ('red','green','blue')}
   valid=np.isfinite(arrays['red'])&np.isfinite(arrays['green'])&np.isfinite(arrays['blue']); cloud=np.zeros_like(valid,bool)
   for label in pd.unique(bins):
    if label is None:continue
    ind=np.where(bins==label)[0];v={b:arrays[b][ind] for b in arrays}
    for rule in cloudy.get(label,[]):cloud[ind]|=rule_mask(v,rule)
   cloud&=valid; rows.append(pd.DataFrame({'time':times,'cloud':cloud.sum(axis=1),'valid':valid.sum(axis=1)}))
 if not rows:return
 d=pd.concat(rows).sort_values('time').drop_duplicates('time');d['hour']=d.time.dt.floor('h');h=d.groupby('hour',as_index=False).agg(rgb_cloud_pixels=('cloud','sum'),rgb_valid_pixels=('valid','sum'),rgb_n_scans=('time','size'));h['rgb_cloud_fraction']=h.rgb_cloud_pixels/h.rgb_valid_pixels
 a.output_dir.mkdir(parents=True,exist_ok=True);out=a.output_dir/f'{a.site}_rgb_hourly_{a.year}{a.month:02d}.csv';h.to_csv(out,index=False);print(out,len(h),flush=True)
if __name__=='__main__':main()
