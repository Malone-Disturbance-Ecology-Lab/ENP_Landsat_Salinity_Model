# EDEN quarterly NetCDFs: /Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/data/EDEN_GRIDS_Q (mounted here as SERVER/...)
import numpy as np, pandas as pd, os
from scipy.io import netcdf_file
H=os.path.expanduser; DL=H('SERVER/data/EDEN_GRIDS_Q')
def load(q):
    nc=netcdf_file(f'{DL}/{q}.nc','r',mmap=False); s=nc.variables['stage'][:].astype('float32'); tm=nc.variables['time'][:]; u=nc.variables['time'].units.decode()
    t0=pd.Timestamp(u.split('since ')[1][:10]); return s,t0+pd.to_timedelta(tm,unit='D')
cache={}
def stack(qs):
    S=[];D=[]
    for q in qs:
        s,d=load(q); S.append(s); D.append(d)
    return np.concatenate(S),pd.DatetimeIndex(np.concatenate([x.values for x in D]))
out={}
for tag,end,qs in [('2025-05-15',pd.Timestamp('2025-05-15'),['2025_q1','2025_q2']),('2025-10-15',pd.Timestamp('2025-10-15'),['2025_q3','2025_q4'])]:
    s,d=stack(qs)
    for w in [30,90]:
        m=(d>end-pd.Timedelta(days=w))&(d<=end)
        out[f'{tag}_m{w}']=np.nanmean(s[m],0).astype('float32'); print(tag,w,m.sum())
np.savez_compressed('downloads/eden_mean_grids_2025.npz',**out)
print('saved')
