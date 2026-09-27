# EDEN quarterly NetCDFs: /Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/data/EDEN_GRIDS_Q (mounted here as SERVER/...)
import numpy as np, pandas as pd, os
from scipy.io import netcdf_file
H=os.path.expanduser; DL=H('SERVER/data/EDEN_GRIDS_Q'); OUT=H('downloads/eden_monthly'); os.makedirs(OUT,exist_ok=True)
cache={}
def load(q):
    if q in cache: return cache[q]
    f=f'{DL}/{q}.nc'
    if not os.path.exists(f): return None
    nc=netcdf_file(f,'r',mmap=False); s=nc.variables['stage'][:].astype('float32'); tm=nc.variables['time'][:]; u=nc.variables['time'].units.decode(); nc.close()
    d=pd.Timestamp(u.split('since ')[1][:10])+pd.to_timedelta(tm,unit='D')
    if len(cache)>=3: cache.pop(next(iter(cache)))
    cache[q]=(s,d); return cache[q]
def qname(t): return f'{t.year}_q{(t.month-1)//3+1}'
for end in pd.date_range('2013-01-01','2025-12-01',freq='MS')+pd.Timedelta(days=14):
    o=f'{OUT}/{end.date()}.npz'
    if os.path.exists(o): continue
    res={}
    for w in [30,90]:
        start=end-pd.Timedelta(days=w-1); qs=sorted({qname(t) for t in pd.date_range(start,end)})
        acc=np.zeros((405,287),'float64'); n=np.zeros((405,287),'float64')
        for q in qs:
            L=load(q)
            if L is None: continue
            s,d=L; m=(d>=start)&(d<=end)
            if m.any(): v=s[m]; ok=np.isfinite(v); acc+=np.where(ok,v,0).sum(0); n+=ok.sum(0)
        res[f'm{w}']=np.where(n>0,acc/np.maximum(n,1),np.nan).astype('float32')
    np.savez_compressed(o,**res)
print('done', len(os.listdir(OUT)))
