import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import pandas as pd, numpy as np, pickle, warnings; warnings.filterwarnings('ignore')
U=f'{ROOT}'
T=pd.read_parquet(f'{U}/dbhydro_hindcast/features/station_features_2000_2025.parquet')
FM=pickle.load(open(f'{U}/best_two_step/salinity30m/hpc_inputs/final_models.pkl','rb')); m3=FM['m3']; F=list(m3.feature_names_in_)
print('missing features:',[f for f in F if f not in T])
DRV=pd.read_csv(f'{U}/dbhydro_hindcast/features/daily_drivers_1998_2025.csv',parse_dates=['date']).set_index('date')
G=pd.read_csv(f'{U}/best_two_step/salinity30m/hpc_inputs/daily_global_drivers.csv',parse_dates=['date']).set_index('date')
for c in [c for c in G.columns if c.startswith('anc_') and not c.endswith('_dep')]:
    j=pd.concat([DRV[c].rename('new'),G[c].rename('old')],axis=1).dropna(); d=(j.new-j.old).abs()
    print(c,len(j),'days >0.05 psu:',int((d>0.05).sum()),'max',round(d.max(),3))
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
T['p3']=np.clip(m3.predict(T[F].astype('float32')),0,None)
o=T[T.salinity.notna()&~T.station.isin(ANCH)].copy(); o['period']=np.where(o.date>='2013-04-01','2013-25 (training)','2000-13 (hindcast)')
o=o[o.date>='2000-04-01']
def met(g):
    y,p=g.salinity,g.p3; return pd.Series(dict(n=len(g),stations=g.station.nunique(),R2=1-((y-p)**2).sum()/((y-y.mean())**2).sum(),RMSE=np.sqrt(((y-p)**2).mean()),bias=(p-y).mean()))
print(o.groupby('period').apply(met).round(3))
# reference: in-sample fit reported for the training years with the original features is not stored; compare with held-out blend
