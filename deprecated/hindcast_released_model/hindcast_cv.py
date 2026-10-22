"""Hindcast 2000-2013 with the one-step (no-spectra) model and test it against the pre-2013 gauges.
A) the released one-step model (m3, trained on all 84 model stations 2013-2025)  -> held out in time only
B) spatial-block CV refits of the same model (10 k-means blocks)                  -> held out in space (and in time before 2013)
C) reproduction check: CV on 2013-2024 + 2025 out-of-time, to compare with Table 2 (one-step R2 0.786 / 0.701)."""
import pandas as pd, numpy as np, pickle, warnings, time; warnings.filterwarnings('ignore')
from sklearn.ensemble import HistGradientBoostingRegressor as HGB
from sklearn.cluster import KMeans
U='/mnt/user-data/uploads/MaloneLab--Research--ENP_Landsat_Salinity_Model'
T=pd.read_parquet(f'{U}/dbhydro_hindcast/features/station_features_2000_2025_v2.parquet')
FM=pickle.load(open(f'{U}/best_two_step/salinity30m/hpc_inputs/final_models.pkl','rb')); m3=FM['m3']; F=list(m3.feature_names_in_)
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
T=T[(T.date>='2000-04-01')&~T.station.isin(ANCH)].reset_index(drop=True)
X=T[F].astype('float32')
T['p_final']=np.clip(m3.predict(X),0,None)
S=pd.read_csv('station_site_variables.csv').set_index('station'); st=sorted(set(T.station))
km=KMeans(10,n_init=10,random_state=0).fit(S.loc[st,['longitude','latitude']].values); blk=dict(zip(st,km.labels_)); T['block']=T.station.map(blk)
print('block sizes',np.bincount(km.labels_))
prm=m3.get_params(); obs=T.salinity.notna()
tr_all=obs&(T.src=='training')
T['p_cv']=np.nan; T['p_cv24']=np.nan; t0=time.time()
for b in range(10):
    trn=tr_all&(T.block!=b); h=T.block==b
    m=HGB(**prm).fit(X[trn],T.salinity[trn]); T.loc[h,'p_cv']=np.clip(m.predict(X[h]),0,None)
    trn24=trn&(T.date<'2025-01-01'); m=HGB(**prm).fit(X[trn24],T.salinity[trn24]); T.loc[h,'p_cv24']=np.clip(m.predict(X[h]),0,None)
    print('block',b,round(time.time()-t0),'s',flush=True)
m=HGB(**prm).fit(X[tr_all&(T.date<'2025-01-01')],T.salinity[tr_all&(T.date<'2025-01-01')]); T['p_all24']=np.clip(m.predict(X),0,None)
T.to_parquet('hindcast_station_predictions.parquet',index=False); print('saved')
