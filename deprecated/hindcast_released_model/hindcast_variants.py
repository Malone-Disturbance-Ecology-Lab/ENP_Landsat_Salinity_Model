import pandas as pd, numpy as np, pickle, warnings, time; warnings.filterwarnings('ignore')
from sklearn.ensemble import HistGradientBoostingRegressor as HGB
U='/mnt/user-data/uploads/MaloneLab--Research--ENP_Landsat_Salinity_Model'
T=pd.read_parquet('hindcast_station_predictions.parquet')
FM=pickle.load(open(f'{U}/best_two_step/salinity30m/hpc_inputs/final_models.pkl','rb')); m3=FM['m3']; F=list(m3.feature_names_in_); prm=m3.get_params()
NOFLOW=[f for f in F if f not in ('SRS90','TSC11190','src_flow90_w')]
tr=T.salinity.notna()&(T.src=='training'); T['p_cv_noflow']=np.nan
for b in range(10):
    trn=tr&(T.block!=b); h=T.block==b
    m=HGB(**prm).fit(T.loc[trn,NOFLOW].astype('float32'),T.salinity[trn]); T.loc[h,'p_cv_noflow']=np.clip(m.predict(T.loc[h,NOFLOW].astype('float32')),0,None)
T.to_parquet('hindcast_station_predictions.parquet',index=False)
o=T[T.salinity.notna()].copy(); o['per']=np.where(o.date<'2013-04-01','2000-13','2013-25')
B=pd.read_csv(f'{U}/best_two_step/timeseries/station_timeseries_metrics.csv').set_index('station').basin; o['basin']=o.station.map(B)
def m(g,c): y,p=g.salinity,g[c]; return pd.Series(dict(R2=1-((y-p)**2).sum()/((y-y.mean())**2).sum(),RMSE=np.sqrt(((y-p)**2).mean()),bias=(p-y).mean()))
for c in ['p_cv','p_cv_noflow']:
    print(c); print(o.groupby('per').apply(lambda g:m(g,c)).round(3).to_string()); print(o.groupby(['basin','per']).apply(lambda g:m(g,c)).round(2).unstack().to_string())
# interannual: annual means at sonde stations with >= 300 obs-days in the year
o['y']=o.date.dt.year; A=o.groupby(['station','y']).agg(n=('salinity','size'),obs=('salinity','mean'),cv=('p_cv','mean'),fin=('p_final','mean')).reset_index(); A=A[A.n>=300]
for per,g in [('2000-12',A[A.y<=2012]),('2014-25',A[A.y>=2014])]:
    g=g.copy(); 
    for c in ['obs','cv','fin']: g[c+'_a']=g[c]-g.groupby('station')[c].transform('mean')
    print(per,'station-years',len(g),'stations',g.station.nunique(),'r(annual anomaly) CV',round(np.corrcoef(g.obs_a,g.cv_a)[0,1],3),'released',round(np.corrcoef(g.obs_a,g.fin_a)[0,1],3))
yr=A.groupby('y')[['obs','cv','fin']].mean().round(1); print(yr.to_string())
