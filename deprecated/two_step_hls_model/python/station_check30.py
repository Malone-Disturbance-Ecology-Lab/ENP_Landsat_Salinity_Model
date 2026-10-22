"""Station checks for the 30 m product: daily skill (30 m vs 250 m), annual trend agreement, event-day agreement."""
import numpy as np, pandas as pd
from scipy.stats import norm, spearmanr
D='change30'
S=pd.read_csv('/mnt/user-data/outputs/daily_salinity_maps/station_daily_model_series_2013_2025.csv.gz',parse_dates=['date'])
M=pd.read_csv(f'{D}/station_daily_model_series_30m.csv.gz',parse_dates=['date'])
S=S.merge(M,on=['station','date'],how='left').dropna(subset=['gauge_salinity','model_salinity','model_salinity_30m']); S['year']=S.date.dt.year
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
def sk(o,m): return dict(n=len(o),r2=1-((o-m)**2).sum()/((o-o.mean())**2).sum(),rmse=np.sqrt(((o-m)**2).mean()),bias=(m-o).mean(),r=np.corrcoef(o,m)[0,1])
rows=[]
for lab,s in [('all',S),('non-anchor',S[~S.station.isin(ANCH)])]:
    for col,res in [('model_salinity','250 m'),('model_salinity_30m','30 m')]: rows.append(dict(set=lab,resolution=res,stations=s.station.nunique(),**sk(s.gauge_salinity,s[col])))
SK=pd.DataFrame(rows).round(3); SK.to_csv(f'{D}/station_daily_skill_30m_vs_250m.csv',index=False); print(SK.to_string())
def sen(y,x):
    i,j=np.triu_indices(len(x),1); sl=np.median((y[j]-y[i])/(x[j]-x[i])); s=np.sign(y[j]-y[i]).sum(); n=len(x); v=n*(n-1)*(2*n+5)/18
    z=(s-np.sign(s))/np.sqrt(v) if s!=0 else 0; return sl,2*norm.sf(abs(z))
rows=[]
for st,g in S.groupby('station'):
    a=g.groupby('year').agg(obs=('gauge_salinity','mean'),mod=('model_salinity_30m','mean'),n=('gauge_salinity','size')); a=a[a.n>=6]
    if len(a)<10: continue
    so,po=sen(a.obs.values,a.index.values.astype(float)); sm,pm=sen(a['mod'].values,a.index.values.astype(float))
    rows.append(dict(station=st,anchor=st in ANCH,years=len(a),obs_mean=a.obs.mean(),obs_slope=so,obs_p=po,model_slope=sm,model_p=pm))
R=pd.DataFrame(rows); R.round(4).to_csv(f'{D}/station_trend_check_30m.csv',index=False); na=R[~R.anchor]
print('trend: stations',len(R),'spearman all',round(spearmanr(R.obs_slope,R.model_slope)[0],2),'non-anchor',round(spearmanr(na.obs_slope,na.model_slope)[0],2),
      'sign agree',round((np.sign(R.obs_slope)==np.sign(R.model_slope)).mean(),2),'obs sig up/down',int(((R.obs_p<.05)&(R.obs_slope>0)).sum()),int(((R.obs_p<.05)&(R.obs_slope<0)).sum()),
      'model sig up/down',int(((R.model_p<.05)&(R.model_slope>0)).sum()),int(((R.model_p<.05)&(R.model_slope<0)).sum()),'median slopes',round(R.obs_slope.median(),3),round(R.model_slope.median(),3))
# events: >=7 consecutive (calendar) days above the series' own 95th percentile, on gauge days only
def evdays(s):
    thr=s.quantile(.95); ex=(s>thr); d=s.index.to_series().diff().dt.days.ne(1); grp=(d|(ex!=ex.shift())).cumsum()
    out=pd.Series(False,index=s.index)
    for _,gg in ex.groupby(grp):
        if gg.iloc[0] and len(gg)>=7: out[gg.index]=True
    return out
rows=[]
for st,g in S.groupby('station'):
    g=g.set_index('date').sort_index()
    if len(g)<365: continue
    eo=evdays(g.gauge_salinity); em=evdays(g.model_salinity_30m); n=g.groupby('year').size()
    for y in n.index[n>=180]: rows.append(dict(station=st,anchor=st in ANCH,year=y,mod=int(em[g.year==y].sum()),obs=int(eo[g.year==y].sum())))
E=pd.DataFrame(rows); E.to_csv(f'{D}/event_days_gauge_vs_model_30m.csv',index=False)
print('events: stations',E.station.nunique(),'station-years',len(E)); print(E.groupby('year')[['obs','mod']].mean().round(1).T.to_string())
print('corr of annual means',round(np.corrcoef(*E.groupby('year')[['obs','mod']].mean().values.T)[0,1],2))
