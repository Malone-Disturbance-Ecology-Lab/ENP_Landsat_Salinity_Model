import sys, numpy as np, pandas as pd; sys.path.insert(0,'.')
from harness import *
TA=pd.read_parquet('TA_index.parquet'); T=pd.read_parquet('T.parquet'); F=pd.read_parquet('fresh_station_features.parquet')
ref=T.drop_duplicates('date').set_index('date')
for c in ['msl_1m','msl_3m','msl_12m']: F[c]=F.date.map(ref[c])
F['year_frac']=0; F['block']=TA[TA.is_new].block.values
A=pd.concat([T,F[T.columns.intersection(F.columns)]],ignore_index=True); orig=~TA.is_new.values; MS=['msl_1m','msl_3m','msl_12m']
R=[]
for seed in [1,2]:
    prm=dict(PRM); prm['random_state']=seed
    for tag,f,data in [('E1',F0,'orig'),('E4',F0,'all'),('E6',F0+MS,'all')]:
        sub=A if data=='all' else A[orig]
        p=np.full(len(A),np.nan); p[sub.index.values]=run_cv(sub.reset_index(drop=True),f,(2000,2025),'raw',prm=prm,tag=f'{tag} seed {seed}')
        m=metrics(T,p[orig],f'{tag} s{seed}'); nn=~orig; m['newfresh_bias']=np.nanmean(p[nn]-A.salinity.values[nn]); R.append(m)
        pd.DataFrame(R).to_csv('results4.csv',index=False)
print(pd.DataFrame(R)[['exp','R2_2000-12','bias_2000-12','R2_2013-25','bias_2013-25','fresh_bias','lt10_bias','station_mean_R2','trend_sign_agree','trend_rho','change_mod','newfresh_bias']].round(3).to_string())
