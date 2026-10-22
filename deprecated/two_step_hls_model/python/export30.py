"""30 m monthly mean salinity (2013-2025) on the native HLS/Landsat grid.
monthly30(m) = bilinear(250 m monthly mean of daily maps) + mean over the month's days of the 30 m detail C,
where C is linearly interpolated in time between mid-month full-model runs (map30/C_int16.npy).
Writes map30/monthly30_int16.npy (156,3839,4171; psu*10; -1 nodata) and one GeoTIFF per month."""
import numpy as np, pandas as pd, rasterio, os, time
t0=time.time(); H,W=3839,4171; NM=156
OUT='/mnt/user-data/outputs/salinity30m'; os.makedirs(f'{OUT}/monthly',exist_ok=True)
wfr=rasterio.open('/mnt/user-data/uploads/ENP_Landsat_Salinity_Model/site_variables/ENP_WaterFrequency_L30.tif'); inside=wfr.read(1)>-1
C=np.load('map30/C_int16.npy',mmap_mode='r'); MED=np.load('daily_maps/median_int16.npy',mmap_mode='r')
edates=pd.date_range('2013-01-01','2025-12-31'); months=pd.date_range('2013-01-01',periods=NM,freq='MS'); mid=months+pd.Timedelta(days=14)
midn=mid.values.astype('datetime64[D]').astype(float); dn=edates.values.astype('datetime64[D]').astype(float)
j=np.clip(np.searchsorted(midn,dn)-1,0,NM-2); f=np.clip((dn-midn[j])/(midn[j+1]-midn[j]),0,1)
rr=np.arange(H); cc=np.arange(W)
fi=(rr-4)/8.0; iy=np.clip(np.floor(fi).astype(int),0,478); fy=np.clip(fi-iy,0,1)[:,None]
fj=(cc-4)/8.0; ix=np.clip(np.floor(fj).astype(int),0,519); fx=np.clip(fj-ix,0,1)[None,:]
def up(g):   # NaN-aware bilinear 250 m (480x521) -> 30 m (3839x4171)
    a=g[iy][:,ix]; b=g[iy][:,ix+1]; c=g[iy+1][:,ix]; d=g[iy+1][:,ix+1]
    w=[(1-fy)*(1-fx),(1-fy)*fx,fy*(1-fx),fy*fx]; v=[a,b,c,d]
    num=sum(np.where(np.isfinite(x),x*y,0) for x,y in zip(v,w)); den=sum(np.where(np.isfinite(x),y,0) for x,y in zip(v,w))
    return np.where(den>0,num/np.maximum(den,1e-9),np.nan).astype('float32')
MO=np.lib.format.open_memmap('map30/monthly30_int16.npy','w+','int16',(NM,H,W))
prof=wfr.profile.copy(); prof.update(count=1,dtype='int16',nodata=-1,compress='deflate',predictor=2,tiled=True,blockxsize=256,blockysize=256,driver='GTiff')
for m,M0 in enumerate(months):
    days=np.flatnonzero((edates>=M0)&(edates<M0+pd.offsets.MonthBegin(1)))
    g=np.asarray(MED[days[0]:days[-1]+1],dtype='float32')/10; g[g<0]=np.nan; g=np.nanmean(g,0)
    # time weights of the mid-month detail layers for this month's days
    wts={}
    for d in days:
        wts[j[d]]=wts.get(j[d],0)+(1-f[d]); wts[j[d]+1]=wts.get(j[d]+1,0)+f[d]
    cm=np.zeros((H,W),'float32'); ws=np.zeros((H,W),'float32')
    for k,wk in wts.items():
        if wk<=0: continue
        c=np.asarray(C[k],dtype='float32'); ok=c>-32768
        cm+=np.where(ok,c/10*wk,0); ws+=np.where(ok,wk,0)
    cm=np.where(ws>0,cm/np.maximum(ws,1e-9),0)
    s=np.clip(up(g)+cm,0,None); out=np.where(inside&np.isfinite(s),np.round(s*10),-1).astype('int16')
    MO[m]=out
    with rasterio.open(f'{OUT}/monthly/ENP_salinity30m_monthly_{M0.year}_{M0.month:02d}.tif','w',**prof) as o:
        o.write(out,1); o.set_band_description(1,f'{M0.year}-{M0.month:02d}'); o.update_tags(units='monthly mean salinity, psu x 10 (divide by 10); -1 = no data',grid='HLS/Landsat 30 m (EPSG:4326, 0.000270923 deg)')
    if m%12==0: print(M0.date(),round(time.time()-t0),'s',flush=True)
MO.flush(); print('done',round(time.time()-t0))
