import numpy as np, pandas as pd
from scipy.io import netcdf_file
from scipy.interpolate import RegularGridInterpolator
D='downloads/'
def load(v):
    nc=netcdf_file(D+f'agg_met_{v}_1979_CurrentYear_CONUS.nc','r',mmap=False); k=[k for k in nc.variables if len(nc.variables[k].shape)==3][0]; a=nc.variables[k]
    x=a.data.astype('float32'); x[a.data==32767]=np.nan; x=x*a.scale_factor+a.add_offset
    lat=nc.variables['lat'][:].copy(); lon=nc.variables['lon'][:].copy(); t=pd.Timestamp('1900-01-01')+pd.to_timedelta(nc.variables['day'][:],unit='D')
    return x,lat,lon,t
pr,lat,lon,t=load('pr'); srad,_,_,_=load('srad'); tn,_,_,_=load('tmmn'); tx,_,_,_=load('tmmx'); tavg=(tn+tx)/2
def roll(x,w,how):
    c=np.nancumsum(np.nan_to_num(x),0); n=np.cumsum(np.isfinite(x),0)
    s=c.copy(); s[w:]=c[w:]-c[:-w]; nn=n.copy(); nn[w:]=n[w:]-n[:-w]
    out=np.where(nn>=int(0.8*w),(s if how=='sum' else s/np.maximum(nn,1)),np.nan); out[:w-1]=np.nan; return out.astype('float32')
F={'P30':roll(pr,30,'sum'),'P90':roll(pr,90,'sum'),'P180':roll(pr,180,'sum'),'P365':roll(pr,365,'sum'),'R30':roll(srad,30,'mean'),'T90':roll(tavg,90,'mean')}
np.savez_compressed('gridmet_features_2012_2025.npz',lat=lat,lon=lon,dates=np.array([str(d.date()) for d in t]),**F)
# consistency check at stations vs training features
a=pd.read_pickle('all_features.pkl'); st=pd.read_csv('sitevars/stations_coast.csv').set_index('station')
la=a.station.map(st.latitude).values; lo=a.station.map(st.longitude).values
ti=pd.Index(t).get_indexer(a.date)
print('tavg training range',a.tavg.describe()[['min','mean','max']].round(1).to_dict())
lat_s=lat[::-1] if lat[0]>lat[-1] else lat
for k,arr in F.items():
    A=arr[:, ::-1,:] if lat[0]>lat[-1] else arr
    vals=np.full(len(a),np.nan)
    for tt in np.unique(ti):
        m=ti==tt
        if tt<0: continue
        f=RegularGridInterpolator((lat_s,lon),A[tt],bounds_error=False,fill_value=np.nan)
        vals[m]=f(np.column_stack([la[m],lo[m]]))
    ok=np.isfinite(vals)&a[k].notna().values
    r=np.corrcoef(vals[ok],a[k].values[ok])[0,1]; rel=np.nanmedian(np.abs(vals[ok]-a[k].values[ok])/np.maximum(np.abs(a[k].values[ok]),1e-6))
    print(k,'n',ok.sum(),'r',round(r,4),'median rel diff',round(rel,4),'grid NaN at stations with training value',int((~np.isfinite(vals)&a[k].notna().values).sum()),'training NaN but grid value',int((np.isfinite(vals)&a[k].isna().values).sum()))
