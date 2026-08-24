#!/usr/bin/env python3
"""Assemble Colorado hourly RGB/ACM/SW comparisons after batch arrays finish."""
from pathlib import Path
import sys
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent; OUT=ROOT/'output_20_colorado_sites_hourly'; BINS=np.linspace(-1,1,41)

def load_acm(site):
 ps=sorted((OUT/'acm_daily').glob(f'{site}_acm_hourly_????????.csv')); return pd.concat([pd.read_csv(p,parse_dates=['hour']) for p in ps]).sort_values('hour').drop_duplicates('hour')
def load_rgb(site):
 if site=='gothic':
  d=pd.read_csv(ROOT/'output_08c_11d_rules/gothic_rgb_11d_pixel_cloud_fraction.csv',parse_dates=['t']);d['hour']=d.t.dt.floor('h')
  return d.groupby('hour',as_index=False).agg(rgb_cloud_fraction=('rgb_11d_cloud_frac','mean'),rgb_n_scans=('t','size'))
 ps=sorted((OUT/'rgb_monthly').glob(f'{site}_rgb_hourly_*.csv'));return pd.concat([pd.read_csv(p,parse_dates=['hour']) for p in ps]).sort_values('hour').drop_duplicates('hour')
def load_sw(site):
 if site=='table_mountain': p=ROOT/'output_13b_boulder_domain_cloud_binary/boulder_rgb_binary_all_available.csv';d=pd.read_csv(p,usecols=['time','k_t','cos_sza'],parse_dates=['time']);d['hour']=d.time.dt.floor('h')
 elif site=='senator_beck': p=ROOT/'output_13c_senator_beck_domain_cloud_binary/senator_beck_rgb_binary_all_available.csv';d=pd.read_csv(p,usecols=['time','k_t','cos_sza'],parse_dates=['time']);d['hour']=d.time.dt.floor('h')
 else:
  sys.path.insert(0,str(ROOT));import train_gothic_sw_rgb_tempbin_trees as g
  raw=g.read_sw_daily_files();s=raw.set_index('time_local').sw_obs.sort_index();times=pd.date_range(s.index.min().floor('D'),s.index.max().ceil('D'),freq='5min');obs=s.reindex(s.index.union(times)).sort_index().interpolate('time').reindex(times);d=pd.DataFrame({'time_local':times,'sw_obs':obs.values}).dropna();d=d.set_index('time_local').between_time(g.LOCAL_START,g.LOCAL_END).reset_index();clear=g.build_metsim_clear_sky(d.time_local,g.GOTHIC_LAT,g.GOTHIC_LON,g.GOTHIC_ELEV_M,g.METSIM_LAPSE_RATE,g.TIME_STEP_MIN,g.METSIM_TBASE,g.LOCAL_TZ);d=d.merge(clear,on='time_local');d['k_t']=d.sw_obs/d.sw_clear;d['cos_sza']=1.;d['time']=d.time_local.dt.tz_localize(g.LOCAL_TZ,ambiguous='NaT',nonexistent='shift_forward').dt.tz_convert('UTC').dt.tz_localize(None);d['hour']=d.time.dt.floor('h')
 d=d[np.isfinite(d.k_t)&(d.cos_sza>=.25)].copy();d['cloud']=(d.k_t<.65).astype(float);h=d.groupby('hour',as_index=False).agg(sw_cloud_fraction=('cloud','mean'),sw_n_valid=('cloud','size'))
 required=1 if site=='senator_beck' else 12;h.loc[h.sw_n_valid!=required,'sw_cloud_fraction']=np.nan;return h
def figures(site,pair):
 pair['season']=pair.hour.dt.month.map({12:'DJF',1:'DJF',2:'DJF',3:'MAM',4:'MAM',5:'MAM',6:'JJA',7:'JJA',8:'JJA',9:'SON',10:'SON',11:'SON'});pair['utc_hour']=pair.hour.dt.hour
 for col,label,color,token in [('rgb_residual','Gothic RGB','tab:blue','rgb'),('acm_residual','GOES ACM','tab:orange','acm')]:
  fig,axs=plt.subplots(2,2,figsize=(13,9),sharex=True,sharey=True)
  for ax,s in zip(axs.ravel(),['DJF','MAM','JJA','SON']):
   v=pair.loc[pair.season==s,col];ax.hist(v,BINS,weights=np.ones(len(v))/len(v),color=color,edgecolor='white');ax.axvline(0,color='k',ls='--');ax.set_title(f'{s} n={len(v):,} bias={v.mean():+.3f}')
  fig.tight_layout();fig.savefig(OUT/f'{site}_{token}_residual_by_season.png',dpi=180);plt.close(fig)
  hours=sorted(pair.utc_hour.unique());fig,axs=plt.subplots(3,4,figsize=(16,11),sharex=True,sharey=True)
  for ax,h in zip(axs.ravel(),hours):
   v=pair.loc[pair.utc_hour==h,col];ax.hist(v,BINS,weights=np.ones(len(v))/len(v),color=color,edgecolor='white');ax.axvline(0,color='k',ls='--');ax.set_title(f'{h:02d}Z n={len(v):,} bias={v.mean():+.3f}')
  for ax in axs.ravel()[len(hours):]:ax.set_visible(False)
  fig.tight_layout();fig.savefig(OUT/f'{site}_{token}_residual_by_utc_hour.png',dpi=180);plt.close(fig)
def main():
 summaries=[]
 for site in ('gothic','table_mountain','senator_beck'):
  rgb,acm,sw=load_rgb(site),load_acm(site),load_sw(site);p=rgb.merge(acm,on='hour').merge(sw,on='hour');p=p.dropna(subset=['rgb_cloud_fraction','acm_cloud_fraction','sw_cloud_fraction']);p=p[(p.rgb_n_scans==12)&(p.acm_n_scans==12)]
  p['rgb_residual']=p.rgb_cloud_fraction-p.sw_cloud_fraction;p['acm_residual']=p.acm_cloud_fraction-p.sw_cloud_fraction;p.to_csv(OUT/f'{site}_rgb_acm_sw_hourly_pairs.csv',index=False)
  for model,col in [('Gothic RGB','rgb_residual'),('GOES ACM','acm_residual')]:
   r=p[col];summaries.append({'site':site,'model':model,'n':len(p),'bias':r.mean(),'mae':r.abs().mean(),'rmse':np.sqrt(np.mean(r*r))})
  figures(site,p)
 pd.DataFrame(summaries).to_csv(OUT/'colorado_sites_rgb_vs_acm_metrics.csv',index=False);print(pd.DataFrame(summaries).to_string(index=False))
if __name__=='__main__':main()
