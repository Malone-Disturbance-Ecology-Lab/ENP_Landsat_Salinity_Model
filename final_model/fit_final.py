"""Final one-step model (no HLS): released predictors + sea level (msl_1m, msl_3m, msl_12m), trained 2000-04..2025-12 at all
non-anchor stations plus 13 fresh DBHydro conductivity sites (equal weights). Block CV -> conformal levels; final fit."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import sys, json, pickle, numpy as np, pandas as pd; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from harness import *
E=f'{WORK}/exp'; MS=['msl_1m','msl_3m','msl_12m']; FE=F0+MS
TA=pd.read_parquet(f'{E}/TA_index.parquet'); T=pd.read_parquet(f'{E}/T.parquet'); F=pd.read_parquet(f'{E}/fresh_station_features.parquet')
ref=T.drop_duplicates('date').set_index('date')
for c in MS: F[c]=F.date.map(ref[c])
F['block']=TA[TA.is_new].block.values
A=pd.concat([T,F[T.columns.intersection(F.columns)]],ignore_index=True); new=TA.is_new.values; orig=~new
p=run_cv(A,FE,(2000,2025),'raw',prm=PRM,tag='final CV')
m=metrics(T,p[orig],'final (v2 + sea level)'); m['newfresh_bias']=np.nanmean(p[new]-A.salinity.values[new])
pd.DataFrame([m]).to_csv('cv_metrics_final.csv',index=False); print(pd.Series(m).to_string())
A[['date','station','salinity','block']].assign(pred_cv=p,fresh_site=new).to_parquet('cv_predictions_final.parquet')
o=A.salinity.notna().values; e=np.abs(p[o]-A.salinity.values[o]); bins=[-1,5,15,30,1e9]; lev=np.digitize(p[o],bins[1:-1])
q=[float(np.quantile(e[lev==k],0.9)) for k in range(4)]
json.dump(dict(bins=bins,q=q,source='2000-2025 block-CV residuals of the final one-step model (stations + 13 fresh sites, sea level; held out in space)'),open('q_levels.json','w'))
yr=A.date.dt.year.values[o]; cov=e<=np.array(q)[lev]; print('q',np.round(q,2))
for nm,s in [('2000-12',yr<2013),('2013-25',yr>=2013),('fresh',new[o])]: print('coverage',nm,round(cov[s].mean(),3))
mdl=HGB(**PRM).fit(A.loc[o,FE].astype('float32'),A.salinity.values[o]); print('iters',mdl.n_iter_)
pickle.dump(dict(m3=mdl,F=FE,PRM=PRM,q_levels=dict(bins=bins,q=q),
    sea_level=dict(stations=[8724580,8723970,8723214],definition='mean of station monthly MSL minus its 2000-2025 mean (m); msl_1m month, msl_3m/msl_12m trailing 3/12-month means'),
    train=dict(dates='2000-04-01..2025-12-31',n_rows=int(o.sum()),n_stations=int(A.station[o].nunique()),fresh_sites=sorted(F.station.unique()))),open('final_model_onestep.pkl','wb'))
ref[MS].to_csv('msl_by_date.csv')
import treelite, tl2cgen
tl2cgen.export_lib(treelite.sklearn.import_model(mdl), toolchain='gcc', libpath='m3_final.so', params={'parallel_comp':8}, nthread=2)
X=T.sample(20000,random_state=0)[FE].astype('float32').values
print('max |compiled - sklearn|', float(np.max(np.abs(tl2cgen.Predictor('m3_final.so',nthread=2).predict(tl2cgen.DMatrix(X)).ravel()-mdl.predict(X)))))
