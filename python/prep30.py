"""Inputs for 30 m mapping on the native HLS/Landsat grid (3839 x 4171, EPSG:4326, 0.000271 deg).
1) static30.npy  : (n_static, n_inside) float32 static features at every 30 m pixel inside the domain
2) spec_grid.npy : (32, 156, 240, 521) float16 vegetation-spectral composites on the sampling grid
                   (rows every 16 px, 8-px blocks), bilinearly resampled to 30 m in run30.py"""
import numpy as np, pandas as pd, rasterio, time, json
from scipy import ndimage as ndi
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from pyproj import Transformer
t0=time.time()
names=['SRS','TaylorSlough','C111','BigCypress']; SRC=['src_dist_SRS_km','src_dist_TaylorSlough_km','src_dist_C111_km','src_dist_BigCypressGulf_km']
U='/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/site_variables/'
wfr=rasterio.open(U+'ENP_WaterFrequency_L30.tif'); wf=wfr.read(1); T=wfr.transform; H,W=wf.shape
inside=(wf>-1).ravel(); idx=np.flatnonzero(inside); NP=len(idx); np.save('px30_index.npy',idx.astype(np.int64))
print('30 m grid',H,W,'inside pixels',NP,flush=True)
st={}
st['distkm']=(rasterio.open('sitevars/ENP_DistCoast_signed.tif').read(1).ravel()[idx]/1000).astype('float32')
Cn=np.load('connectivity_rasters.npz')
for k in ['chan_dist_strict_km','chan_dist_soft_km','chan_density_500m','chan_density_1km','inland_water_1km']: st[k]=Cn[k].ravel()[idx].astype('float32')
st['chan_sinuosity']=np.clip(st['chan_dist_soft_km']/np.maximum(np.clip(st['distkm'],0,None),0.1),None,10); st['chan_dist_strict_km']=np.clip(st['chan_dist_strict_km'],None,100)
Sm=np.load('offshore_smoothed_v3.npz')
for n in names: st['w_'+n]=Sm['w_'+n].ravel()[idx].astype('float32'); st[f'basin_cost_{n}_km']=np.clip(Sm['cost_'+n].ravel()[idx],None,150).astype('float32')
for k in SRC: st[k]=np.clip(Sm[k].ravel()[idx],None,150).astype('float32')
np.save('basin30.npy',Sm['basin'].ravel()[idx].astype('int8'))
print('connectivity/basins',round(time.time()-t0),flush=True)
V=rasterio.open('sitevars/ENP_VegMap_30m.tif').read(1); cd=pd.read_csv('sitevars/vegmap_codes.csv').set_index('id')
cls_of=np.array(['']+list(cd.cls.reindex(range(1,cd.index.max()+1)).fillna('')),dtype=object); Vc=cls_of[V]
def frac(mask,rad):
    r_=int(round(rad/28.6)); yy,xx=np.mgrid[-r_:r_+1,-r_:r_+1]; k=((yy**2+xx**2)<=r_**2).astype('float32'); k/=k.sum(); return fftconvolve(mask.astype('float32'),k,mode='same')
den={1000:frac(V>0,1000).ravel()[idx],250:frac(V>0,250).ravel()[idx]}
mapped=den[1000]
for nm,cl,rad in [('p_Open water_1000m','Open water',1000),('p_Mangrove_1000m','Mangrove',1000),('p_Freshwater marsh/prairie_1000m','Freshwater marsh/prairie',1000),('p_Swamp scrub/shrubland_1000m','Swamp scrub/shrubland',1000),('p_Salt marsh_1000m','Salt marsh',1000),('p_Open water_250m','Open water',250),('p_Mangrove_250m','Mangrove',250)]:
    f=frac(Vc==cl,rad).ravel()[idx]; st[nm]=np.where(den[rad]>0,f/np.maximum(den[rad],1e-6),0).astype('float32')
unm=mapped<0.5; st['p_Open water_1000m'][unm]=1; st['p_Open water_250m'][unm]=1
np.save('vegclass30.npy',V.ravel()[idx].astype('int16')); del V,Vc
print('veg fractions',round(time.time()-t0),flush=True)
EM=np.load('downloads/eden_monthly_mean_grids_2013_2025.npz'); anyv=np.zeros((405,287),bool)
for k in EM.files[:24]: anyv|=np.isfinite(EM[k])
vy,vx=np.nonzero(anyv); tree=cKDTree(np.column_stack([463400+vx*400,2790200+vy*400]))
rr,cc=np.unravel_index(idx,(H,W)); lon,lat=rasterio.transform.xy(T,rr,cc); lon=np.asarray(lon,'float64'); lat=np.asarray(lat,'float64')
np.save('lonlat30.npy',np.stack([lon,lat]).astype('float64'))
X,Y=Transformer.from_crs(4326,26917,always_xy=True).transform(lon,lat); np.save('utm30.npy',np.stack([X,Y]).astype('float64'))
dd,_=tree.query(np.column_stack([X,Y]),workers=2); st['dist_to_eden_cell_km']=(dd/1000).astype('float32')
keys=sorted(st); np.save('static30.npy',np.stack([st[k] for k in keys])); json.dump(keys,open('static30_keys.json','w'))
print('static saved',len(keys),round(time.time()-t0),flush=True)
# ---------------- vegetation spectral composites on the sampling grid (same processing as the 250 m maps)
IDX=['NDVI','EVI2','NDMI','SWIR1']; II=[0,1,2,5]; vrows=list(range(4,3839,16)); VR=np.load('veg_rows_L30.npz')
M=np.stack([VR[f'row_{r:04d}'][:,II,:].astype('float32') for r in vrows],2)
for q,(lo_,hi_) in enumerate([(-1,1),(-1,1.5),(-1,1),(0,1)]): M[:,q][(M[:,q]<lo_)|(M[:,q]>hi_)]=np.nan
med=np.nanmedian(M,0,keepdims=True); mad=np.nanmedian(np.abs(M-med),0,keepdims=True); M[np.abs(M-med)>4*1.4826*np.maximum(mad,0.01)]=np.nan; del med,mad
NM=156; dy=16*abs(T.e)*111.0; dx=8*T.a*111.0*np.cos(np.deg2rad(25.3))
def disc(rkm):
    ry=int(np.ceil(rkm/dy)); rx=int(np.ceil(rkm/dx)); yy,xx=np.mgrid[-ry:ry+1,-rx:rx+1]; return (((yy*dy)**2+(xx*dx)**2)<=rkm**2+1e-9).astype('float32')
SPK=[]; SPG=np.lib.format.open_memmap('spec_grid.npy','w+','float16',(32,NM)+M.shape[2:])
q_=0
for rn,rk in {'1km':1.0,'3km':3.0}.items():
    k=disc(rk)
    for q,nm in enumerate(IDX):
        g=M[:,q]; v=np.isfinite(g)
        num=np.stack([ndi.convolve(np.where(v[m],g[m],0),k,mode='constant') for m in range(NM)]); dn=np.stack([ndi.convolve(v[m].astype('float32'),k,mode='constant') for m in range(NM)])
        nb=np.where(dn>0,num/np.maximum(dn,1e-6),np.nan).astype('float32')
        flat=pd.DataFrame(nb.reshape(NM,-1))
        comp={'c12':flat.rolling(12,center=True,min_periods=3).median().values,'c3':flat.rolling(3,center=True,min_periods=1).median().values}
        lt=np.nanmedian(nb.reshape(NM,-1),0); comp['lt']=np.broadcast_to(lt,(NM,lt.size)); comp['d12']=comp['c12']-lt[None]
        for kk,a in comp.items():
            SPG[q_]=np.asarray(a,dtype='float32').reshape(nb.shape).astype('float16'); SPK.append(f'{nm}_{kk}_{rn}'); q_+=1
        print('spectra',rn,nm,round(time.time()-t0),flush=True)
SPG.flush(); json.dump(SPK,open('spec_grid_keys.json','w')); print('done',round(time.time()-t0))
