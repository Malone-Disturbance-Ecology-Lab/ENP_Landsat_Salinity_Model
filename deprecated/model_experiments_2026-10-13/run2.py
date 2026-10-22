import sys, numpy as np, pandas as pd; sys.path.insert(0,'.')
from harness import *
from sklearn.cluster import KMeans
T=pd.read_parquet('T.parquet'); F=pd.read_parquet('fresh_station_features.parquet')
# sea-level + year on fresh rows (same definitions as harness)
mo=F.date.dt.to_period('M').dt.to_timestamp(); ref=T.drop_duplicates('date').set_index('date')
for c in ['msl_1m','msl_3m','msl_12m']: F[c]=F.date.map(ref[c])
F['year_frac']=F.date.dt.year+F.date.dt.dayofyear/365.25
# blocks: same k-means on original station coordinates; new stations go to the nearest centroid
S=pd.read_csv(f'{U}/site_variables/station_site_variables.csv').set_index('station'); st=sorted(set(T.station))
km=KMeans(10,n_init=10,random_state=0).fit(S.loc[st,['longitude','latitude']].values); assert (km.predict(S.loc[st,['longitude','latitude']].values)==T.drop_duplicates('station').set_index('station').loc[st,'block'].values).all()
NS=pd.read_csv('cond_station_summary.csv',index_col=0); F['block']=km.predict(NS.loc[F.station,['lon','lat']].values)
print(F.groupby('station').block.first().to_string())
TA=pd.concat([T,F[T.columns.intersection(F.columns)]],ignore_index=True); TA['is_new']=TA.station.isin(F.station.unique())
orig=~TA.is_new.values
def M(p,tag):
    m=metrics(TA[orig].reset_index(drop=True),p[orig],tag)
    nn=TA.is_new.values&np.isfinite(p); m['newfresh_bias']=(p[nn]-TA.salinity.values[nn]).mean(); m['newfresh_pred_mean']=p[nn].mean(); return m
MS=['msl_1m','msl_3m','msl_12m']
exps=[('E1 train 2000-25',F0,'raw',None),
      ('E4 E1 + 13 fresh sites',F0,'raw','fresh'),
      ('E5 E4 + sqrt target',F0,'sqrt','fresh'),
      ('E6 E4 + sea level',F0+MS,'raw','fresh'),
      ('E7 E4 + year (diagnostic)',F0+['year_frac'],'raw','fresh')]
R=[]; P={}
for tag,f,tg,data in exps:
    sub=TA if data=='fresh' else TA[orig]
    p=np.full(len(TA),np.nan); q=run_cv(sub.reset_index(drop=True),f,(2000,2025),tg,tag=tag); p[sub.index.values]=q; P[tag[:2]]=p
    R.append(M(p,tag)); print(pd.DataFrame(R).round(3).T.to_string(),flush=True)
pd.DataFrame(R).to_csv('results2.csv',index=False); np.savez_compressed('preds2.npz',**P); TA[['station','date','salinity','block','is_new']].to_parquet('TA_index.parquet')
