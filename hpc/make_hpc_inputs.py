"""Small inputs for the HPC daily 30 m run: filled gridMET fields (+ departures) and daily global drivers (flows, anchors).
Run from the pipeline working directory (needs gridmet_features_2012_2025.npz and all_features.pkl)."""
import numpy as np, pandas as pd, sys
sys.path.insert(0,'.'); from fillers import laplace_solver
OUT=sys.argv[1] if len(sys.argv)>1 else 'hpc_inputs'
import os; os.makedirs(OUT,exist_ok=True)
G=np.load('gridmet_features_2012_2025.npz'); glat=G['lat']; glon=G['lon']
flip=glat[0]>glat[-1]; glat0=glat[::-1] if flip else glat; step=np.diff(glat0).mean(); npad=7; glat_s=np.r_[glat0[0]-step*np.arange(npad,0,-1),glat0]
GA={}
for k in ['P30','P90','P180','P365','R30','T90']:
    A_=G[k][:,::-1,:] if flip else G[k]; A_=np.concatenate([np.full((A_.shape[0],npad,A_.shape[2]),np.nan,'float32'),A_],1)
    kn=np.isfinite(A_[-1]); GA[k]=laplace_solver(kn)(np.where(kn[None],A_,0)).astype('float32'); GA[k][~np.isfinite(A_).any(axis=(1,2))]=np.nan
def cen(x,w=365,minp=60):
    v=np.isfinite(x); c=np.concatenate([np.zeros((1,)+x.shape[1:]),np.cumsum(np.where(v,x,0),0,dtype='float64')]); n=np.concatenate([np.zeros((1,)+x.shape[1:]),np.cumsum(v,0)])
    i=np.arange(len(x)); hi=np.minimum(i+w//2+1,len(x)); lo=np.maximum(i-w//2,0); return np.where(n[hi]-n[lo]>=minp,(c[hi]-c[lo])/np.maximum(n[hi]-n[lo],1),np.nan).astype('float32')
for k in ['P30','P90','P180']: GA[k+'_dep']=GA[k]-cen(GA[k])
np.savez_compressed(f'{OUT}/gridmet_filled.npz',glat_s=glat_s,glon=glon,dates=G['dates'],**GA)
AF=pd.read_pickle('all_features.pkl'); AF['date']=pd.to_datetime(AF.date); ed=pd.date_range('2013-01-01','2025-12-31')
FL=['SRS30','SRS90','SRS180','TSC11130','TSC11190','TSC111180','S18C30','S18C90','S19790']
fl=AF.groupby('date')[FL].first().reindex(ed).interpolate(limit=5)
for c in FL: fl[c+'_dep']=fl[c]-fl[c].rolling(365,center=True,min_periods=180).mean()
piv=AF.pivot_table(index='date',columns='station',values='salinity').reindex(ed)
for s in ['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']:
    a30=piv[s].rolling(30,min_periods=10).mean(); fl[f'anc_{s}_30']=a30; fl[f'anc_{s}_3']=piv[s].rolling(3,min_periods=1).mean(); fl[f'anc_{s}_dep']=a30-piv[s].rolling(365,center=True,min_periods=120).mean()
fl.index.name='date'; fl.to_csv(f'{OUT}/daily_global_drivers.csv'); print('ok',fl.shape)
