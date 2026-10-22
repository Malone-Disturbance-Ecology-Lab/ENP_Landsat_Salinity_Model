"""Held-out (10-block spatial CV) validation of the final one-step model, 2000-2025: overall, by station, by region, by salinity class, intervals."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import json, numpy as np, pandas as pd
U=f'{ROOT}'
C=pd.read_parquet(f'{WORK}/final/cv_predictions_final.parquet'); C=C[C.salinity.notna()].copy()
S=pd.read_csv(f'{U}/site_variables/station_site_variables.csv').set_index('station')
F=pd.read_csv(f'{WORK}/exp/fresh13_coords.csv').set_index('station').rename(columns={'lat':'latitude','lon':'longitude'})
LL=pd.concat([S[['latitude','longitude']],F[['latitude','longitude']]])
IDX=np.load(f'{WORK}/map/in/px30_index.npy'); REG=np.load(f'{U}/best_two_step/events_masked/in/region30.npy'); W=4171
D30=0.00027092271166751034; LEFT,TOP=-81.51874748177192,25.889645249658873
REGN=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']
def region(st):
    la,lo=LL.loc[st,['latitude','longitude']]; p=int(np.floor((TOP-la)/D30))*W+int(np.floor((lo-LEFT)/D30)); i=np.searchsorted(IDX,p)
    return REGN[REG[i]] if i<len(IDX) and IDX[i]==p and 0<=REG[i]<6 else 'outside'
st=sorted(C.station.unique()); RG={s:region(s) for s in st}; C['region']=C.station.map(RG)
Q=json.load(open(f'{WORK}/final/q_levels.json')); bins=np.array(Q['bins']); q=np.array(Q['q'])
C['hw']=q[np.clip(np.digitize(C.pred_cv,bins[1:-1]),0,3)]; C['cov']=(C.salinity-C.pred_cv).abs()<=C.hw
def kge(o,p):
    r=np.corrcoef(o,p)[0,1]; return 1-np.sqrt((r-1)**2+(p.std()/o.std()-1)**2+(p.mean()/o.mean()-1)**2)
def pooled(g):
    e=g.pred_cv-g.salinity; return dict(n=len(g),stations=g.station.nunique(),R2=1-(e**2).sum()/((g.salinity-g.salinity.mean())**2).sum(),RMSE=np.sqrt((e**2).mean()),bias=e.mean(),coverage=g['cov'].mean())
rows=[]
for s,g in C.groupby('station'):
    if len(g)<30 or g.salinity.std()==0: continue
    e=g.pred_cv-g.salinity; rows.append(dict(station=s,region=RG[s],fresh=g.fresh_site.iloc[0],days=len(g),mean_obs=g.salinity.mean(),r=np.corrcoef(g.salinity,g.pred_cv)[0,1],KGE=kge(g.salinity.values,g.pred_cv.values),RMSE=np.sqrt((e**2).mean()),bias=e.mean(),coverage=g['cov'].mean()))
ST=pd.DataFrame(rows); ST.round(3).to_csv('station_skill_final.csv',index=False)
M=C[~C.fresh_site]; Sm=ST[~ST.fresh]
out={'pooled_model_stations':pooled(M),'pooled_fresh':pooled(C[C.fresh_site]),
     'station_median':dict(n=len(Sm),r=Sm.r.median(),KGE=Sm.KGE.median(),RMSE=Sm.RMSE.median(),abs_bias=Sm.bias.abs().median(),pct_r07=100*(Sm.r>0.7).mean(),pct_kge05=100*(Sm.KGE>0.5).mean(),pct_kge0=100*(Sm.KGE>0).mean(),pct_cov_lt70=int((Sm.coverage<0.7).sum()))}
print(json.dumps(out,indent=1,default=float))
R=[]
for k,g in M.groupby('region'):
    s=Sm[Sm.region==k]; R.append(dict(region=k,stations=len(s),days=len(g),R2=pooled(g)['R2'],median_r=s.r.median(),median_KGE=s.KGE.median(),median_RMSE=s.RMSE.median(),median_bias=s.bias.median()))
R=pd.DataFrame(R); R.round(3).to_csv('skill_by_region_final.csv',index=False); print(R.round(2).to_string())
Sm=Sm.assign(cls=pd.cut(Sm.mean_obs,[-1,2,10,25,99],labels=['< 2','2–10','10–25','> 25']))
print(Sm.groupby('cls').agg(stations=('r','size'),r=('r','median'),RMSE=('RMSE','median'),bias=('bias','median')).round(2))
M=M.assign(lev=pd.cut(M.pred_cv,bins,labels=['0–5','5–15','15–30','>30']))
cv=M.groupby('lev')['cov'].mean(); print('coverage by bin',cv.round(3).to_dict(),'overall',round(M['cov'].mean(),3),'median width',np.median(2*M.hw))
cv.to_csv('interval_coverage_final.csv')
