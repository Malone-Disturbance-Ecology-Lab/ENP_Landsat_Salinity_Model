# EDEN quarterly NetCDFs: /Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/data/EDEN_GRIDS_Q (mounted here as SERVER/...)
import numpy as np, pandas as pd, os, sys, glob
from scipy.io import netcdf_file
from scipy import ndimage as ndi
from pyproj import Transformer
H=os.path.expanduser
DL=H('SERVER/data/EDEN_GRIDS_Q'); OUT=H('downloads')
st=pd.read_csv(H('SERVER/DBHydro_lonlat.csv')).drop_duplicates('station')
t=Transformer.from_crs(4326,26917,always_xy=True)
X,Y=t.transform(st.longitude.values,st.latitude.values)
def cells(valid):
    # nearest valid cell for each station (EDEN grid 400 m, x0=463400, y0=2790200)
    dist,(iy,ix)=ndi.distance_transform_edt(~valid,return_indices=True)
    xi=np.clip(np.round((X-463400)/400).astype(int),0,286); yi=np.clip(np.round((Y-2790200)/400).astype(int),0,404)
    # distance from station to nearest valid cell centre (m), accounting for stations outside the extent
    ny=iy[yi,xi]; nx=ix[yi,xi]
    d=np.hypot(463400+nx*400-X, 2790200+ny*400-Y)
    return ny,nx,d
meta=None
for f in sorted(glob.glob(DL+'/20[12][0-9]_q[1-4].nc')):
    q=os.path.basename(f)[:-3]; o=f'{OUT}/{q}.csv'
    if os.path.exists(o): continue
    try:
        nc=netcdf_file(f,'r',mmap=False); s=nc.variables['stage'][:].astype('float32'); tm=nc.variables['time'][:]; units=nc.variables['time'].units.decode()
    except Exception as e:
        print('skip',q,e); continue
    t0=pd.Timestamp(units.split('since ')[1][:10]); dates=t0+pd.to_timedelta(tm,unit='D')
    valid=np.isfinite(s).any(0)
    if meta is None:
        ny,nx,d=cells(valid); meta=pd.DataFrame({'station':st.station,'lon':st.longitude,'lat':st.latitude,'cell_y':ny,'cell_x':nx,'dist_to_cell_m':d.round(0)})
        meta.to_csv(f'{OUT}/station_cells.csv',index=False)
    # 5 km neighbourhood mean (valid cells within 12.5 cells radius)
    r=12; yy,xx=np.mgrid[-r:r+1,-r:r+1]; k=((yy**2+xx**2)<=(12.5**2)).astype('float32')
    rows=[]
    for di in range(s.shape[0]):
        a=s[di]; v=np.isfinite(a); num=ndi.convolve(np.where(v,a,0),k,mode='constant'); den=ndi.convolve(v.astype('float32'),k,mode='constant')
        m5=np.where(den>0,num/np.maximum(den,1),np.nan)
        rows.append(pd.DataFrame({'station':meta.station,'date':dates[di].strftime('%Y-%m-%d'),'stage_cell_cm':a[meta.cell_y,meta.cell_x],'stage_5km_cm':m5[meta.cell_y,meta.cell_x]}))
    pd.concat(rows).round(2).to_csv(o,index=False); print('done',q,s.shape)
