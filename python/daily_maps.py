"""Daily 250 m salinity maps 2013-01-01..2025-12-31 from the final blend
   pred = 0.7*(B_hat + anomaly_hat) + 0.3*v3,  B_hat = smoothed spectral stage-1 mean (monthly -> daily interp).
Resumable: daily medians are written to daily_maps/median_int16.npy (psu*10, -1 = no data) with a progress file."""
import pandas as pd, numpy as np, pickle, rasterio, os, sys, time, warnings; warnings.filterwarnings('ignore')
from scipy import ndimage as ndi
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from pyproj import Transformer
from fillers import laplace_solver
FM=pickle.load(open('final_models.pkl','rb')); m1,m2,m3=FM['m1'],FM['m2'],FM['m3']
F1=list(m1.feature_names_in_); F2=list(m2.feature_names_in_); F3=list(m3.feature_names_in_)
OUT='daily_maps'; os.makedirs(OUT,exist_ok=True)
names=['SRS','TaylorSlough','C111','BigCypress']; SRC=['src_dist_SRS_km','src_dist_TaylorSlough_km','src_dist_C111_km','src_dist_BigCypressGulf_km']
U='SERVER/site_variables/'
wfr=rasterio.open(U+'ENP_WaterFrequency_L30.tif'); wf=wfr.read(1); T=wfr.transform; K=8
RR,CC=np.meshgrid(np.arange(K//2,wf.shape[0],K),np.arange(K//2,wf.shape[1],K),indexing='ij'); SH=RR.shape
inside=(wf[RR,CC]>-1).ravel(); NP=int(inside.sum()); print('map grid',SH,'pixels',NP,flush=True)
lon,lat=rasterio.transform.xy(T,RR.ravel()[inside],CC.ravel()[inside]); lon=np.array(lon); lat=np.array(lat)
# ---------------- static pixel features (same definitions as v3)
st={}
st['distkm']=rasterio.open('sitevars/ENP_DistCoast_signed.tif').read(1)[RR,CC].ravel()[inside]/1000
Cn=np.load('connectivity_rasters.npz')
for k in ['chan_dist_strict_km','chan_dist_soft_km','chan_density_500m','chan_density_1km','inland_water_1km']: st[k]=Cn[k][RR,CC].ravel()[inside]
st['chan_sinuosity']=np.clip(st['chan_dist_soft_km']/np.maximum(np.clip(st['distkm'],0,None),0.1),None,10); st['chan_dist_strict_km']=np.clip(st['chan_dist_strict_km'],None,100)
Sm=np.load('offshore_smoothed_v3.npz')
for n in names: st['w_'+n]=Sm['w_'+n][RR,CC].ravel()[inside]; st[f'basin_cost_{n}_km']=np.clip(Sm['cost_'+n][RR,CC].ravel()[inside],None,150)
for k in SRC: st[k]=np.clip(Sm[k][RR,CC].ravel()[inside],None,150)
basin=np.array(names)[Sm['basin'][RR,CC].ravel()[inside]-1]
V=rasterio.open('sitevars/ENP_VegMap_30m.tif').read(1); cd=pd.read_csv('sitevars/vegmap_codes.csv').set_index('id')
cls_of=np.array(['']+list(cd.cls.reindex(range(1,cd.index.max()+1)).fillna('')),dtype=object); Vc=cls_of[V]
def frac(mask,rad):
    r_=int(round(rad/28.6)); yy,xx=np.mgrid[-r_:r_+1,-r_:r_+1]; k=((yy**2+xx**2)<=r_**2).astype('float32'); k/=k.sum(); return fftconvolve(mask.astype('float32'),k,mode='same')
mapped=frac(V>0,1000)[RR,CC].ravel()[inside]
for nm,cl,rad in [('p_Open water_1000m','Open water',1000),('p_Mangrove_1000m','Mangrove',1000),('p_Freshwater marsh/prairie_1000m','Freshwater marsh/prairie',1000),('p_Swamp scrub/shrubland_1000m','Swamp scrub/shrubland',1000),('p_Salt marsh_1000m','Salt marsh',1000),('p_Open water_250m','Open water',250),('p_Mangrove_250m','Mangrove',250)]:
    f=frac(Vc==cl,rad)[RR,CC].ravel()[inside]; den=frac(V>0,rad)[RR,CC].ravel()[inside]; st[nm]=np.where(den>0,f/np.maximum(den,1e-6),0)
unm=mapped<0.5; st['p_Open water_1000m'][unm]=1; st['p_Open water_250m'][unm]=1
del V,Vc
EM=np.load('downloads/eden_monthly_mean_grids_2013_2025.npz'); anyv=np.zeros((405,287),bool)
for k in EM.files[:24]: anyv|=np.isfinite(EM[k])
vy,vx=np.nonzero(anyv); tree=cKDTree(np.column_stack([463400+vx*400,2790200+vy*400]))
tU=Transformer.from_crs(4326,26917,always_xy=True); X,Y=tU.transform(lon,lat); dd,_=tree.query(np.column_stack([X,Y])); st['dist_to_eden_cell_km']=dd/1000
print('static done',flush=True)
BSF=f'{OUT}/stage1_monthly_mean_smoothed.npy'; RESUME=os.path.exists(BSF)
NM=156; months=pd.date_range('2013-01-01',periods=NM,freq='MS'); mid=months+pd.Timedelta(days=14)
if not RESUME:
    # ---------------- vegetation spectral features on the veg grid -> map grid
    IDX=['NDVI','EVI2','NDMI','SWIR1']; II=[0,1,2,5]; vrows=list(range(4,3839,16)); VR=np.load('veg_rows_L30.npz')
    M=np.stack([VR[f'row_{r:04d}'][:,II,:].astype('float32') for r in vrows],2)
    for q,(lo_,hi_) in enumerate([(-1,1),(-1,1.5),(-1,1),(0,1)]): M[:,q][(M[:,q]<lo_)|(M[:,q]>hi_)]=np.nan
    med=np.nanmedian(M,0,keepdims=True); mad=np.nanmedian(np.abs(M-med),0,keepdims=True); M[np.abs(M-med)>4*1.4826*np.maximum(mad,0.01)]=np.nan; del med,mad
    NM=156; dy=16*abs(T.e)*111.0; dx=8*T.a*111.0*np.cos(np.deg2rad(25.3))
    def disc(rkm):
        ry=int(np.ceil(rkm/dy)); rx=int(np.ceil(rkm/dx)); yy,xx=np.mgrid[-ry:ry+1,-rx:rx+1]; return (((yy*dy)**2+(xx*dx)**2)<=rkm**2+1e-9).astype('float32')
    # map-pixel -> veg-grid (row fraction, col)
    vi=(RR.ravel()[inside]-4)/16.0; i0=np.clip(np.floor(vi).astype(int),0,len(vrows)-1); i1=np.clip(i0+1,0,len(vrows)-1); fi=np.clip(vi-i0,0,1).astype('float32')
    jc=np.clip(CC.ravel()[inside]//8,0,M.shape[3]-1)
    SPECPX={}   # name -> (156, NP) float16
    for rn,rk in {'1km':1.0,'3km':3.0}.items():
        k=disc(rk)
        for q,nm in enumerate(IDX):
            g=M[:,q]; v=np.isfinite(g)
            num=np.stack([ndi.convolve(np.where(v[m],g[m],0),k,mode='constant') for m in range(NM)]); den=np.stack([ndi.convolve(v[m].astype('float32'),k,mode='constant') for m in range(NM)])
            nb=np.where(den>0,num/np.maximum(den,1e-6),np.nan).astype('float32')        # (156, nr, nc)
            flat=pd.DataFrame(nb.reshape(NM,-1))
            comp={'c12':flat.rolling(12,center=True,min_periods=3).median().values,'c3':flat.rolling(3,center=True,min_periods=1).median().values}
            lt=np.nanmedian(nb.reshape(NM,-1),0); comp['lt']=np.broadcast_to(lt,(NM,lt.size)); comp['d12']=comp['c12']-lt[None]
            for kk,a in comp.items():
                a=np.asarray(a,dtype='float32').reshape(nb.shape)
                v0=a[:,i0,jc]; v1=a[:,i1,jc]; val=np.where(np.isfinite(v0)&np.isfinite(v1),v0*(1-fi)+v1*fi,np.where(np.isfinite(v0),v0,v1))
                SPECPX[f'{nm}_{kk}_{rn}']=val.astype('float16')
            print('spectra',rn,nm,flush=True)
    del M
# ---------------- gridMET (padded + harmonic fill, as v3) with pixel bilinear weights
G=np.load('gridmet_features_2012_2025.npz'); glat=G['lat']; glon=G['lon']; gd=pd.to_datetime(G['dates'])
flip=glat[0]>glat[-1]; glat0=glat[::-1] if flip else glat; step=np.diff(glat0).mean(); npad=7; glat_s=np.r_[glat0[0]-step*np.arange(npad,0,-1),glat0]
CLIMV=['P30','P90','P180','P365','R30','T90']; GA={}
for k in CLIMV:
    A_=G[k][:,::-1,:] if flip else G[k]; A_=np.concatenate([np.full((A_.shape[0],npad,A_.shape[2]),np.nan,'float32'),A_],1)
    kn=np.isfinite(A_[-1]); GA[k]=laplace_solver(kn)(np.where(kn[None],A_,0)).astype('float32'); GA[k][~np.isfinite(A_).any(axis=(1,2))]=np.nan
def cen(x,w=365,minp=60):
    v=np.isfinite(x); c=np.concatenate([np.zeros((1,)+x.shape[1:]),np.cumsum(np.where(v,x,0),0,dtype='float64')]); n=np.concatenate([np.zeros((1,)+x.shape[1:]),np.cumsum(v,0)])
    i=np.arange(len(x)); hi=np.minimum(i+w//2+1,len(x)); lo=np.maximum(i-w//2,0); return np.where(n[hi]-n[lo]>=minp,(c[hi]-c[lo])/np.maximum(n[hi]-n[lo],1),np.nan).astype('float32')
for k in ['P30','P90','P180']: GA[k+'_dep']=GA[k]-cen(GA[k])
iy=np.clip(np.searchsorted(glat_s,lat)-1,0,len(glat_s)-2); ix=np.clip(np.searchsorted(glon,lon)-1,0,len(glon)-2)
fy=np.clip((lat-glat_s[iy])/(glat_s[iy+1]-glat_s[iy]),0,1); fx=np.clip((lon-glon[ix])/(glon[ix+1]-glon[ix]),0,1)
def gsamp(A2): return (1-fy)*(1-fx)*A2[iy,ix]+(1-fy)*fx*A2[iy,ix+1]+fy*(1-fx)*A2[iy+1,ix]+fy*fx*A2[iy+1,ix+1]
gidx=pd.Index(gd)
# ---------------- EDEN daily fields (cropped padded grid), pixel bilinear weights
meta=np.load('eden_daily/crop_meta.npy'); ey0,ex0,eh,ew=meta; eh,ew=int(eh),int(ew)
EF={k:np.load(f'eden_daily/{k}.npy',mmap_mode='r') for k in ['cell','cell30','Rcell','5km30','5km90','R5']}
ry=(Y-ey0)/400; rx=(X-ex0)/400; ey=np.clip(np.floor(ry).astype(int),0,eh-2); ex=np.clip(np.floor(rx).astype(int),0,ew-2); efy=np.clip(ry-ey,0,1); efx=np.clip(rx-ex,0,1)
def esamp(A2):
    A2=np.asarray(A2,dtype='float32'); return (1-efy)*(1-efx)*A2[ey,ex]+(1-efy)*efx*A2[ey,ex+1]+efy*(1-efx)*A2[ey+1,ex]+efy*efx*A2[ey+1,ex+1]
edates=pd.date_range('2013-01-01','2025-12-31')
# ---------------- flows + anchors (global daily)
AF=pd.read_pickle('all_features.pkl'); AF['date']=pd.to_datetime(AF.date)
FL=['SRS30','SRS90','SRS180','TSC11130','TSC11190','TSC111180','S18C30','S18C90','S19790']
fl=AF.groupby('date')[FL].first().reindex(edates).interpolate(limit=5)
for c in FL: fl[c+'_dep']=fl[c]-fl[c].rolling(365,center=True,min_periods=180).mean()
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
piv=AF.pivot_table(index='date',columns='station',values='salinity').reindex(edates)
anc=pd.DataFrame(index=edates)
for s in ANCH:
    a30=piv[s].rolling(30,min_periods=10).mean(); anc[f'anc_{s}_30']=a30; anc[f'anc_{s}_3']=piv[s].rolling(3,min_periods=1).mean()
    anc[f'anc_{s}_dep']=a30-piv[s].rolling(365,center=True,min_periods=120).mean()
if RESUME:
    Bs=np.load(BSF).astype('float32'); print('resumed: loaded stage-1 means',flush=True)
else:
    # ---------------- stage 1: monthly B_hat per pixel, smoothed, then daily interpolation
    months=pd.date_range('2013-01-01',periods=NM,freq='MS'); mid=months+pd.Timedelta(days=14)
    Bm=np.full((NM,NP),np.nan,'float32')
    for m,Dt in enumerate(mid):
        cols={}
        cols.update(st)
        for k in SPECPX: cols[k]=SPECPX[k][m].astype('float32')
        ti=gidx.get_loc(Dt); cols['P365']=gsamp(GA['P365'][ti]); cols['P180']=gsamp(GA['P180'][ti])
        cols['stage_5km_cm90']=esamp(EF['5km90'][edates.get_loc(Dt)])
        mf=fl.loc[months[m]:months[m]+pd.offsets.MonthEnd(0)]; cols['SRS180']=np.full(NP,mf.SRS180.mean()); cols['TSC111180']=np.full(NP,mf.TSC111180.mean())
        Bm[m]=m1.predict(np.column_stack([cols[c] for c in F1]))
    Bs=pd.DataFrame(Bm).rolling(12,center=True,min_periods=1).mean().values.astype('float32')
    np.save(f'{OUT}/stage1_monthly_mean_smoothed.npy',Bs.astype('float16'))
    print('stage1 monthly done; B range',np.nanpercentile(Bs,[1,50,99]).round(1),flush=True)
    del SPECPX
midn=mid.values.astype('datetime64[D]').astype(float)
# ---------------- daily loop
Q=FM['Q']; bins=FM['bins']; lab=FM['lab']; qlev=np.array([Q.get(l,FM['Qall']) for l in lab],dtype='float32')
MED=np.lib.format.open_memmap(f'{OUT}/median_int16.npy','r+' if os.path.exists(f'{OUT}/median_int16.npy') else 'w+','int16',(len(edates),)+SH)
prog=f'{OUT}/progress.txt'; start=int(open(prog).read()) if os.path.exists(prog) else 0
# monthly accumulators for mean / 90% bounds (recomputed from stored medians at the end)
w=np.array([0.7,0.3],dtype='float32'); base2={c:st[c] for c in F2 if c in st}; base3={c:st[c] for c in F3 if c in st}
t0=time.time()
for di in range(start,len(edates)):
    Dt=edates[di]; ti=gidx.get_loc(Dt); cols={}
    for k in CLIMV+['P30_dep','P90_dep','P180_dep']: cols[k]=gsamp(GA[k][ti])
    doy=Dt.dayofyear; cols['doy_cos']=np.full(NP,np.cos(2*np.pi*doy/365.25)); cols['doy_sin']=np.full(NP,np.sin(2*np.pi*doy/365.25))
    e={k:esamp(EF[k][di]) for k in EF}
    cols['stage_cell_cm']=e['cell']; cols['stage_cell_cm30']=e['cell30']; cols['stage_5km_cm30']=e['5km30']; cols['stage_5km_cm90']=e['5km90']
    cols['stage_cell_cm_dep']=e['cell']-e['Rcell']; cols['stage_cell_cm30_dep']=e['cell30']-e['Rcell']; cols['stage_5km_cm30_dep']=e['5km30']-e['R5']; cols['stage_5km_cm90_dep']=e['5km90']-e['R5']
    for c in fl.columns: cols[c]=np.full(NP,fl.at[Dt,c])
    for c in anc.columns: cols[c]=np.full(NP,anc.at[Dt,c])
    cols['src_flow90_w']=(st['w_SRS']+st['w_BigCypress'])*cols['SRS90']+(st['w_TaylorSlough']+st['w_C111'])*cols['TSC11190']
    x=float(Dt.to_datetime64().astype('datetime64[D]').astype(float)); j=np.clip(np.searchsorted(midn,x)-1,0,NM-2); f=np.clip((x-midn[j])/(midn[j+1]-midn[j]),0,1)
    B=Bs[j]*(1-f)+Bs[j+1]*f
    a=m2.predict(np.column_stack([cols[c] if c in cols else base2[c] for c in F2]).astype('float32'))
    v3=m3.predict(np.column_stack([cols[c] if c in cols else base3[c] for c in F3]).astype('float32'))
    p=np.clip(0.7*(B+a)+0.3*v3,0,None)
    out=np.full(SH[0]*SH[1],-1,'int16'); out[inside]=np.round(p*10).astype('int16'); MED[di]=out.reshape(SH)
    if di%30==0:
        MED.flush(); open(prog,'w').write(str(di+1)); print(Dt.date(),'median psu',round(float(np.median(p)),1),'elapsed',round(time.time()-t0),'s',flush=True)
MED.flush(); open(prog,'w').write(str(len(edates)))
np.save(f'{OUT}/interval_q_by_level.npy',qlev); np.save(f'{OUT}/basin_px.npy',basin)
print('done')
