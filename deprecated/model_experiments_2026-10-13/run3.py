import sys, numpy as np, pandas as pd, time; sys.path.insert(0,'.')
from harness import *
from sklearn.ensemble import HistGradientBoostingRegressor as HGB
TA=pd.read_parquet('TA_index.parquet'); T=pd.read_parquet('T.parquet'); F=pd.read_parquet('fresh_station_features.parquet')
ref=T.drop_duplicates('date').set_index('date')
for c in ['msl_1m','msl_3m','msl_12m']: F[c]=F.date.map(ref[c])
F['year_frac']=F.date.dt.year+F.date.dt.dayofyear/365.25; F['block']=TA[TA.is_new].block.values
A=pd.concat([T,F[T.columns.intersection(F.columns)]],ignore_index=True); assert (A.station.values==TA.station.values).all()
orig=~TA.is_new.values; MS=['msl_1m','msl_3m','msl_12m']
def lin_basis(D):
    coast=np.exp(-np.clip(D.distkm.values,0,None)/5.0)          # 1 offshore/at coast, ~0 > 15 km inland
    s=D.msl_12m.values*100                                        # cm
    return np.c_[s, s*coast]
def run(feats, mode, tag, prm=None):
    obs=A.salinity.notna().values; X=A[feats].astype('float32'); y=A.salinity.values; p=np.full(len(A),np.nan); t0=time.time(); coefs=[]
    for b in range(10):
        trn=obs&(A.block.values!=b); h=A.block.values==b
        if mode=='lin':
            # salinity response to sea level (cm), common + coastal part, after removing station means (within-station OLS)
            Bt=lin_basis(A[trn]); st=A.station.values[trn]
            yd=y[trn]-pd.Series(y[trn]).groupby(st).transform('mean').values; Bd=Bt-pd.DataFrame(Bt).groupby(st).transform('mean').values
            beta=np.linalg.lstsq(Bd,yd,rcond=None)[0]; coefs.append(beta)
            off_t=Bt@beta; m=HGB(**(prm or PRM)).fit(X[trn],y[trn]-off_t); p[h]=m.predict(X[h])+lin_basis(A[h])@beta
        else:
            m=HGB(**(prm or PRM)).fit(X[trn],y[trn]); p[h]=m.predict(X[h])
    print(tag,'CV',round(time.time()-t0),'s', 'beta (psu per cm: all, coastal extra) mean', np.round(np.mean(coefs,0),3) if coefs else '', flush=True)
    p=np.clip(p,0,None); m=metrics(T,p[orig],tag); nn=~orig; m['newfresh_bias']=(p[nn]-y[nn]).mean(); return m,p
R=[]; P={}
for feats,mode,tag in [(F0,'lin','E8 E4 + linear sea-level term'),(F0+MS,'lin','E9 E6 + linear sea-level term')]:
    m,p=run(feats,mode,tag); R.append(m); P[tag[:2]]=p; print(pd.DataFrame(R).round(3).T.to_string(),flush=True)
pd.DataFrame(R).to_csv('results3.csv',index=False); np.savez_compressed('preds3.npz',**P)
