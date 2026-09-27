import pandas as pd, numpy as np, rasterio, pickle, os, warnings; warnings.filterwarnings('ignore')
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
from pyproj import Transformer
exec(open('conn_cv.py').read().split("res=[]; zt={}")[0])
from sklearn.ensemble import HistGradientBoostingRegressor
names=['SRS','TaylorSlough','C111','BigCypress']
CHAN=['chan_dist_soft_km','chan_dist_strict_km','chan_sinuosity','chan_density_500m','chan_density_1km','inland_water_1km']
W=['w_'+n for n in names]; BC=[f'basin_cost_{n}_km' for n in names]
SRC=['src_dist_SRS_km','src_dist_TaylorSlough_km','src_dist_C111_km','src_dist_BigCypressGulf_km']
CLIMV=['P30','P90','P180','P365','R30','T90']
FEATS=F0+CHAN+SRC+W+BC+['src_flow90_w']+AF

from fillers import laplace_solver, nearest_inside
print('stage feats in FEATS:',[f for f in FEATS if 'stage' in f or 'eden' in f],flush=True)
# ---- gridMET: pad south to 24.78N, harmonic fill of ocean / outside cells
G=np.load('gridmet_features_2012_2025.npz'); glat=G['lat']; glon=G['lon']; gd=pd.to_datetime(G['dates'])
flip=glat[0]>glat[-1]; glat0=glat[::-1] if flip else glat; step=np.diff(glat0).mean(); npad=7
glat_s=np.r_[glat0[0]-step*np.arange(npad,0,-1),glat0]
GA={}
for k in CLIMV:
    A=G[k][:,::-1,:] if flip else G[k]; A=np.concatenate([np.full((A.shape[0],npad,A.shape[2]),np.nan,'float32'),A],1)
    kn=np.isfinite(A[-1]); GA[k]=laplace_solver(kn)(np.where(kn[None],A,0)).astype('float32'); GA[k][~np.isfinite(A).any(axis=(1,2))]=np.nan
print('gridMET filled',{k:float(np.isnan(GA[k][-1]).mean()) for k in CLIMV},flush=True)
def gfield(k,ti): return GA[k][ti]
stco=pd.read_csv('sitevars/stations_coast.csv').set_index('station')
def bilin_w(la,lo):
    iy=int(np.clip(np.searchsorted(glat_s,la)-1,0,len(glat_s)-2)); ix=int(np.clip(np.searchsorted(glon,lo)-1,0,len(glon)-2))
    fy=(la-glat_s[iy])/(glat_s[iy+1]-glat_s[iy]); fx=(lo-glon[ix])/(glon[ix+1]-glon[ix])
    inb=(glat_s[0]<=la<=glat_s[-1]) and (glon[0]<=lo<=glon[-1])
    return iy,ix,fy,fx,inb
def station_clim(df):
    df=df.copy(); ti=pd.Index(gd).get_indexer(df.date)
    for s in df.station.unique():
        m=(df.station==s).values; iy,ix,fy,fx,inb=bilin_w(stco.loc[s,'latitude'],stco.loc[s,'longitude'])
        if not inb: continue
        tt=ti[m]; okt=tt>=0
        for k in CLIMV:
            A=GA[k]; v=(1-fy)*(1-fx)*A[:,iy,ix]+(1-fy)*fx*A[:,iy,ix+1]+fy*(1-fx)*A[:,iy+1,ix]+fy*fx*A[:,iy+1,ix+1]
            col=df.loc[m,k].values.astype(float); new=np.full(int(m.sum()),np.nan); new[okt]=v[tt[okt]]
            df.loc[m,k]=np.where(np.isfinite(new),new,col)
    return df
def add_basin(df,SB):
    df=df.copy()
    for c in W+BC+SRC: df[c]=df.station.map(SB[c]).astype(float)
    for c in BC+SRC: df[c]=df[c].clip(upper=150)
    df['basin']=df.station.map(SB.basin)
    df['src_flow90_w']=(df.w_SRS+df.w_BigCypress)*df.SRS90+(df.w_TaylorSlough+df.w_C111)*df.TSC11190
    return df

SB3=pd.read_csv('station_basins_v3.csv',index_col=0)
# ---- EDEN: harmonic extrapolation beyond the EDEN domain (replaces nearest-cell fill)
EM=np.load('downloads/eden_monthly_mean_grids_2013_2025.npz')
anyv=np.zeros((405,287),bool)
for k in EM.files[:24]: anyv|=np.isfinite(EM[k])
r_=12; yy,xx=np.mgrid[-r_:r_+1,-r_:r_+1]; kk=((yy**2+xx**2)<=12.5**2).astype('float32')
reg5=ndi.convolve(anyv.astype('float32'),kk,mode='constant')>0
PS,PW,PE,PN=140,70,10,10; y0p=2790200-PS*400; x0p=463400-PW*400
def pad(g): return np.pad(g,((PS,PN),(PW,PE)),constant_values=np.nan)
fill_c=laplace_solver(pad(anyv)); fill_5=laplace_solver(pad(reg5))
mkeys=sorted(set(k[:10] for k in EM.files)); EXT={}; OLD={}
vy,vx=np.nonzero(anyv); tree=cKDTree(np.column_stack([463400+vx*400,2790200+vy*400]))
def m5(g):
    v=np.isfinite(g); num=ndi.convolve(np.where(v,g,0),kk,mode='constant'); den=ndi.convolve(v.astype('float32'),kk,mode='constant'); return np.where(den>0,num/np.maximum(den,1),np.nan)
def nearest_valid(g):
    v=np.isfinite(g)
    if v.all() or not v.any(): return g
    _,(iy,ix)=ndi.distance_transform_edt(~v,return_indices=True); return g[iy,ix]
for d_ in mkeys:
    for w_ in ('m30','m90'):
        g=EM[f'{d_}_{w_}']; gi=nearest_inside(g,anyv); g5=m5(gi)
        EXT[(d_,w_,'cell')]=fill_c(np.nan_to_num(pad(gi))).astype('float32'); EXT[(d_,w_,'5km')]=fill_5(np.nan_to_num(pad(g5))).astype('float32')
        OLD[(d_,w_,'cell')]=nearest_valid(g).astype('float32'); OLD[(d_,w_,'5km')]=m5(g).astype('float32')
print('EDEN extrapolated',len(mkeys),flush=True)
from scipy.ndimage import map_coordinates
tU=Transformer.from_crs(4326,26917,always_xy=True)
def samp(ext,X,Y): return map_coordinates(ext,[(np.asarray(Y)-y0p)/400,(np.asarray(X)-x0p)/400],order=1,mode='nearest')
def station_stage(df):
    df=df.copy(); mdates=pd.to_datetime(mkeys); t_m=mdates.values.astype('datetime64[D]').astype(float)
    for s in df.station.unique():
        m=(df.station==s).values
        if s not in stco.index: continue
        X,Y=tU.transform(stco.loc[s,'longitude'],stco.loc[s,'latitude']); _,ii=tree.query([X,Y]); ny_,nx_=vy[ii],vx[ii]
        td=df.loc[m,'date'].values.astype('datetime64[D]').astype(float)
        for w_ in ('m30','m90'):
            for typ,col in (('cell',f'stage_cell_cm{w_[1:]}'),('5km',f'stage_5km_cm{w_[1:]}')):
                if col not in df: continue
                off=np.array([float(samp(EXT[(d_,w_,typ)],[X],[Y])[0])-OLD[(d_,w_,typ)][ny_,nx_] for d_ in mkeys])
                ok=np.isfinite(off)
                if ok.any(): df.loc[m,col]=df.loc[m,col].values+np.interp(td,t_m[ok],off[ok])
    return df
trn_c=station_stage(station_clim(trn)); ten_c=station_stage(station_clim(ten)); print('station climate + stage done',flush=True)
A=add_basin(trn_c,SB3); B_=add_basin(ten_c,SB3)
pickle.dump(dict(A=A,B=B_,FEATS=FEATS),open('v3_tables.pkl','wb')); print('saved',A.shape,B_.shape,flush=True)
