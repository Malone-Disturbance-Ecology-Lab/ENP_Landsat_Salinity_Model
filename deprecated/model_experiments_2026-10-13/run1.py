import sys, numpy as np, pandas as pd, pickle; sys.path.insert(0,'.')
from harness import *
T=load(); T.to_parquet('T.parquet'); print(T.shape, T.salinity.notna().sum(), flush=True)
print(T[['msl_1m','msl_12m']].describe().round(3))
exps=[('E0 released setup, train 2013-25',F0,(2013,2025),'raw'),
      ('E1 train 2000-25',F0,(2000,2025),'raw'),
      ('E2 E1 + sea level',F0+['msl_1m','msl_3m','msl_12m'],(2000,2025),'raw'),
      ('E3 E2 + log target',F0+['msl_1m','msl_3m','msl_12m'],(2000,2025),'log')]
R=[]; P={}
for tag,f,yrs,tg in exps:
    p=run_cv(T,f,yrs,tg,tag=tag); P[tag]=p; R.append(metrics(T,p,tag)); print(pd.DataFrame(R).round(3).T.to_string(),flush=True)
pd.DataFrame(R).to_csv('results1.csv',index=False); np.savez_compressed('preds1.npz',**{k[:2]:v for k,v in P.items()})
